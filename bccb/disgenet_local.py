#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2022
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: Dénes Türei (turei.denes@gmail.com)
#           Nicolàs Palacio
#           Sebastian Lobentanzer
#           Erva Ulusoy
#           Olga Ivanova
#           Ahmet Rifaioglu
#           Melih Darcan
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

import collections
import csv
import json
from getpass import getpass

from typing import List, Union, NamedTuple, Dict, Tuple, Optional, Literal

import pypath.share.curl as curl
import pypath.share.session as session
import pypath.resources.urls as urls
import pypath.utils.mapping as mapping

_logger = session.Logger(name="disgenet_input")
_log = _logger._log


class DisgenetApi:
    
    _name = "DisGeNET"
    _api_url = urls.urls["disgenet"]["api_url"]
    _authenticated: bool = False
    _api_key: str = None
    e_mail: str = None
    password: str = None

    def authenticate(self) -> bool:
        """
        Starts an authorization process in DisGeNET API.
        Returns a boolean which is success of authentication.
        """

        if self._authenticated and self._api_key != None:
            return True

        print(f"Authorizing in {self._name} API...")
        # DisGeNET migrated from email/password login (returning a session
        # token) to a static per-user API key generated on their website.
        # There is no /auth/ endpoint to call anymore, so we just prompt
        # for the key and store it directly instead of making a POST
        # request and parsing a token from the response.
        api_key: str = getpass("API Key: ")

        self._api_key = api_key
        self._authenticated = True

        return self._authenticated and (self._api_key != None)

    def _if_authenticated(f):
        """
        Simple wrapper to get rid of the burden of
        checking authentication status and authenticating
        if already haven't.
        """

        def wrapper(self, *args, **kwargs):
            if self.authenticate():
                return f(self, *args, **kwargs)

            else:
                _log("DisGeNET failure in authorization, check your credentials.")

        return wrapper

    def _delete_cache(f):
        """
        A necessary wrapper as databases may be updated
        after the initial download of a particular data.
        """

        def wrapper(*args, **kwargs):
            with curl.cache_delete_on():
                return f(*args, **kwargs)

        return wrapper

    # get_ddas_that_share_genes() and get_ddas_that_share_variants() were
    # removed: the new DisGeNET DDA API has a single /dda endpoint that
    # always returns both gene-sharing and variant-sharing metrics
    # together in one query, so the separate genes/variants split no
    # longer applies.

    
    # _get_vdas(), get_vdas_by_variants(), get_vdas_by_genes(),
    # get_vdas_by_diseases(), and get_vdas_by_source() were removed. These
    # five functions existed to select a "by" routing mode (gene, disease,
    # variant, source) for the old path-based /vda/{by}/... endpoint,
    # where each identifier type required a different URL. The new API
    # has a single fixed endpoint (/vda/summary) where variant, gene,
    # disease, source, etc. are all just optional query parameters on the
    # same call, so there's no more separate URL/routing per identifier
    # type. Replaced by the single get_vda_summary() function below, which
    # accepts all identifier types and filters at once.

    def get_vda_summary(
    self,
    variant: Union[str, List[str]] = None,
    gene_ncbi_id: Union[str, List[str]] = None,
    gene_symbol: Union[str, List[str]] = None,
    disease: Union[str, List[str]] = None,
    source: Union[str, List[str]] = None,
    min_score: float = None,
    max_score: float = None,
    min_ei: float = None,
    max_ei: float = None,
    min_dsi: float = None,
    max_dsi: float = None,
    min_dpi: float = None,
    max_dpi: float = None,
    dis_class_list: Union[str, List[str]] = None,
    page_number: int = None,
) -> NamedTuple(
    "VariantDiseaseAssociation",
    [
        ("variantid", str),
        ("gene_symbol", Tuple[str]),
        ("variant_dsi", float),
        ("variant_dpi", float),
        ("variant_consequence_type", str),
        ("diseaseid", str),
        ("disease_name", str),
        ("disease_classes_do", Tuple[str]),
        ("disease_classes_hpo", Tuple[str]),
        ("disease_classes_msh", Tuple[str]),
        ("disease_classes_umls_st", Tuple[str]),
        ("disease_type", str),
        ("score", float),
        ("ei", float),
        ("year_initial", int),
        ("year_final", int),
        ("source", str),
    ],
):
        """
        Returns Variant-Disease Associations (aggregated summary), narrowed
        to the fields the old get_vdas_by_variants/genes/diseases/source()
        functions returned.

        disease_type is now populated from the API's "diseaseType" field.
        It was previously hardcoded to None based on an incorrect reading
        of the response schema; confirmed via real test data (e.g.
        disease_type='phenotype' for a real query) that this field is
        actually present in /vda/summary's response.

        disease_class_name is not included: same situation as in
        get_dda/get_gda_summary - the old API had a separate class-name
        field, the new API only returns the combined "name (code)" string
        in disease_class.

        @variant: Union[str, List[str]]
            Variant (dbSNP Identifier) or list of variants, up to 100.
        @gene_ncbi_id / @gene_symbol: Union[str, List[str]]
            Gene identifier(s); the old "gene" param accepted either NCBI
            ID or HGNC symbol, the new API splits these into two params.
        @disease: Union[str, List[str]]
            Disease id(s) with vocabulary prefix, e.g. "UMLS_C0005745".
        @source: Union[str, List[str]]
            Source of the VDA.
        @min_score/@max_score, @min_ei/@max_ei, @min_dsi/@max_dsi,
        @min_dpi/@max_dpi: float
            Score ranges, in [0,1].
        @dis_class_list: Union[str, List[str]]
            MeSH Disease Classes - corresponds to the old "disease_class"
            param.
        @page_number: int
            Page number - corresponds to the old "limit" param (100 per
            page; TRIAL accounts are capped at the top-30 results and do
            not support pagination).
        """

        url = f"{self._api_url}/vda/summary"
        get_params = dict()

        if variant != None:
            get_params["variant"] = self._list_to_str(variant, "Variant ID", limit=100)

        if gene_ncbi_id != None:
            get_params["gene_ncbi_id"] = self._list_to_str(gene_ncbi_id, "Gene NCBI ID", limit=100)

        if gene_symbol != None:
            get_params["gene_symbol"] = self._list_to_str(gene_symbol, "Gene Symbol", limit=100)

        if disease != None:
            get_params["disease"] = self._list_to_str(disease, "Disease ID", limit=100)

        if source != None:
            get_params["source"] = source

        if min_score != None:
            get_params["min_score"] = str(min_score)

        if max_score != None:
            get_params["max_score"] = str(max_score)

        if min_ei != None:
            get_params["min_ei"] = str(min_ei)

        if max_ei != None:
            get_params["max_ei"] = str(max_ei)

        if min_dsi != None:
            get_params["min_dsi"] = str(min_dsi)

        if max_dsi != None:
            get_params["max_dsi"] = str(max_dsi)

        if min_dpi != None:
            get_params["min_dpi"] = str(min_dpi)

        if max_dpi != None:
            get_params["max_dpi"] = str(max_dpi)

        if dis_class_list != None:
            get_params["dis_class_list"] = dis_class_list

        if page_number != None:
            get_params["page_number"] = str(page_number)

        result = self._retrieve_data(url, get_params)

        if result == None:
            return None

        VariantDiseaseAssociation = collections.namedtuple(
            "VariantDiseaseAssociation",
            [
                "variantid",
                "gene_symbol",
                "variant_dsi",
                "variant_dpi",
                "variant_consequence_type",
                "diseaseid",
                "disease_name",
                "disease_classes_do",
                "disease_classes_hpo",
                "disease_classes_msh",
                "disease_classes_umls_st",
                "disease_type",
                "score",
                "ei",
                "year_initial",
                "year_final",
                "source",
            ],
        )

        for index, entry in enumerate(result):
            result[index] = VariantDiseaseAssociation(
                self._get_string(entry.get("variantStrID")),
                tuple(entry["geneSymbol_keyword"]) if entry.get("geneSymbol_keyword") else None,
                self._get_float(entry.get("variantDSI")),
                self._get_float(entry.get("variantDPI")),
                self._get_string(entry.get("mostSevereConsequences")),
                self._get_string(entry.get("diseaseUMLSCUI")),
                self._get_string(entry.get("diseaseName")),
                tuple(entry["diseaseClasses_DO"]) if entry.get("diseaseClasses_DO") else None,
                tuple(entry["diseaseClasses_HPO"]) if entry.get("diseaseClasses_HPO") else None, # Likely None for most records, but keeping the field since it may be populated depending on the variant/gene - didn't want to silently drop this classification data.
                tuple(entry["diseaseClasses_MSH"]) if entry.get("diseaseClasses_MSH") else None,
                tuple(entry["diseaseClasses_UMLS_ST"]) if entry.get("diseaseClasses_UMLS_ST") else None,
                self._get_string(entry.get("diseaseType")),
                self._get_float(entry.get("score")),
                self._get_float(entry.get("ei")),
                self._get_int(entry.get("yearInitial")),
                self._get_int(entry.get("yearFinal")),
                self._get_string(entry.get("source")),
            )

        return result

    
    # get_gdas_by_genes(), get_gdas_by_diseases(), get_gdas_by_uniprots(),
    # get_gdas_by_source(), and _get_gdas() were removed. These functions
    # existed to select a "by" routing mode (gene, disease, uniprot,
    # source) for the old path-based /gda/{by}/... endpoint, where each
    # identifier type required a different URL. The new API has a single
    # fixed endpoint (/gda/summary) where gene, disease, uniprot, source,
    # etc. are all just optional query parameters on the same call, so
    # there's no more separate URL/routing per identifier type. Replaced
    # by the single get_gda_summary() function below, which accepts all
    # identifier types and filters at once. This mirrors the same removal
    # already done for the VDA summary functions (get_vdas_by_*/_get_vdas).
    def get_gda_summary(
        self,
        gene_ncbi_id: Union[str, List[str]] = None,
        gene_ensembl_id: Union[str, List[str]] = None,
        gene_symbol: Union[str, List[str]] = None,
        uniprot_id: Union[str, List[str]] = None,
        disease: Union[str, List[str]] = None,
        source: Union[str, List[str]] = None,
        min_score: float = None,
        max_score: float = None,
        min_ei: float = None,
        max_ei: float = None,
        min_dsi: float = None,
        max_dsi: float = None,
        min_dpi: float = None,
        max_dpi: float = None,
        min_pli: float = None,
        max_pli: float = None,
        type: str = None,
        dis_class_list: Union[str, List[str]] = None,
        page_number: int = None,
    ) -> NamedTuple(
        "GeneDiseaseAssociation",
        [
            ("geneid", int),
            ("gene_symbol", str),
            ("uniprotid", Tuple[str]),
            ("gene_dsi", float),
            ("gene_dpi", float),
            ("gene_pli", float),
            ("protein_class", Tuple[str]),
            ("protein_class_name", Tuple[str]),
            ("diseaseid", str),
            ("disease_name", str),
            ("disease_classes_do", Tuple[str]),
            ("disease_classes_hpo", Tuple[str]),
            ("disease_classes_msh", Tuple[str]),
            ("disease_classes_umls_st", Tuple[str]),
            ("disease_type", str),
            ("score", float),
            ("ei", float),
            ("el", str),
            ("year_initial", int),
            ("year_final", int),
            ("source", str),
        ],
    ):
        """
        Returns Gene-Disease Associations (aggregated summary), narrowed
        to the fields the old get_gdas_by_genes/diseases/uniprots/source()
        functions returned.

        @gene_ncbi_id / @gene_ensembl_id / @gene_symbol: Union[str, List[str]]
            Gene identifier(s) in the respective vocabulary, up to 100.
        @uniprot_id: Union[str, List[str]]
            Uniprot accession(s), up to 100.
        @disease: Union[str, List[str]]
            Disease id(s) with vocabulary prefix, e.g. "UMLS_C0005745".
        @source: Union[str, List[str]]
            Source of the GDA.
        @min_score/@max_score, @min_ei/@max_ei, @min_dsi/@max_dsi,
        @min_dpi/@max_dpi, @min_pli/@max_pli: float
            Score ranges, in [0,1].
        @type: str
            DisGeNET Disease Type ("disease", "phenotype", "group").
        @dis_class_list: Union[str, List[str]]
            MeSH Disease Classes - corresponds to the old "disease_class"
            param.
        @page_number: int
            Page number - corresponds to the old "limit" param (100 per
            page; TRIAL accounts are capped at the top-30 results and do
            not support pagination).
        """

        url = f"{self._api_url}/gda/summary"
        get_params = dict()

        if gene_ncbi_id != None:
            get_params["gene_ncbi_id"] = self._list_to_str(gene_ncbi_id, "Gene NCBI ID", limit=100)

        if gene_ensembl_id != None:
            get_params["gene_ensembl_id"] = self._list_to_str(gene_ensembl_id, "Gene Ensembl ID", limit=100)

        if gene_symbol != None:
            get_params["gene_symbol"] = self._list_to_str(gene_symbol, "Gene Symbol", limit=100)

        if uniprot_id != None:
            get_params["uniprot_id"] = self._list_to_str(uniprot_id, "Uniprot ID", limit=100)

        if disease != None:
            get_params["disease"] = self._list_to_str(disease, "Disease ID", limit=100)

        if source != None:
            get_params["source"] = source

        if min_score != None:
            get_params["min_score"] = str(min_score)

        if max_score != None:
            get_params["max_score"] = str(max_score)

        if min_ei != None:
            get_params["min_ei"] = str(min_ei)

        if max_ei != None:
            get_params["max_ei"] = str(max_ei)

        if min_dsi != None:
            get_params["min_dsi"] = str(min_dsi)

        if max_dsi != None:
            get_params["max_dsi"] = str(max_dsi)

        if min_dpi != None:
            get_params["min_dpi"] = str(min_dpi)

        if max_dpi != None:
            get_params["max_dpi"] = str(max_dpi)

        if min_pli != None:
            get_params["min_pli"] = str(min_pli)

        if max_pli != None:
            get_params["max_pli"] = str(max_pli)

        if type != None:
            get_params["type"] = type

        if dis_class_list != None:
            get_params["dis_class_list"] = dis_class_list

        if page_number != None:
            get_params["page_number"] = str(page_number)

        result = self._retrieve_data(url, get_params)

        if result == None:
            return None

        GeneDiseaseAssociation = collections.namedtuple(
            "GeneDiseaseAssociation",
            [
                "geneid",
                "gene_symbol",
                "uniprotid",
                "gene_dsi",
                "gene_dpi",
                "gene_pli",
                "protein_class",
                "protein_class_name",
                "diseaseid",
                "disease_name",
                "disease_classes_do",
                "disease_classes_hpo",
                "disease_classes_msh",
                "disease_classes_umls_st",
                "disease_type",
                "score",
                "ei",
                "el",
                "year_initial",
                "year_final",
                "source",
            ],
        )

        for index, entry in enumerate(result):
            result[index] = GeneDiseaseAssociation(
                self._get_int(entry.get("geneNcbiID")),
                self._get_string(entry.get("symbolOfGene")),
                tuple(entry["geneProteinStrIDs"]) if entry.get("geneProteinStrIDs") else None,
                self._get_float(entry.get("geneDSI")),
                self._get_float(entry.get("geneDPI")),
                self._get_float(entry.get("genepLI")),
                tuple(entry["geneProteinClassIDs"]) if entry.get("geneProteinClassIDs") else None,
                tuple(entry["geneProteinClassNames"]) if entry.get("geneProteinClassNames") else None,
                self._get_string(entry.get("diseaseUMLSCUI")),
                self._get_string(entry.get("diseaseName")),
                tuple(entry["diseaseClasses_DO"]) if entry.get("diseaseClasses_DO") else None,
                tuple(entry["diseaseClasses_HPO"]) if entry.get("diseaseClasses_HPO") else None,
                tuple(entry["diseaseClasses_MSH"]) if entry.get("diseaseClasses_MSH") else None,
                tuple(entry["diseaseClasses_UMLS_ST"]) if entry.get("diseaseClasses_UMLS_ST") else None,
                self._get_string(entry.get("diseaseType")),
                self._get_float(entry.get("score")),
                self._get_float(entry.get("ei")),
                self._get_string(entry.get("el")),
                self._get_int(entry.get("yearInitial")),
                self._get_int(entry.get("yearFinal")),
                self._get_string(entry.get("source")),
            )

        return result

    
        
    # get_vda_evidences_by_variant() and get_vda_evidences_by_disease() were
    # removed. They existed only to select a "by" routing mode for the old
    # path-based /vda/evidences/{by}/... endpoint. The new API has a single
    # fixed endpoint (/vda/evidence). Replaced by get_vda_evidence() below,
    # which keeps only the params the old functions/​_get_evidences() used
    # (variant/gene/disease, source, min/max_year, min/max_score,
    # limit/offset) rather than the full param set the new endpoint
    # supports. Like the old code, this does not convert results into a
    # namedtuple - it returns the raw list of record dicts as-is.
    def get_vda_evidence(
        self,
        variant: Union[str, List[str]] = None,
        gene_ncbi_id: Union[str, List[str]] = None,
        gene_symbol: Union[str, List[str]] = None,
        disease: Union[str, List[str]] = None,
        source: Union[str, List[str]] = None,
        min_pmYear: str = None,
        max_pmYear: str = None,
        min_score: float = None,
        max_score: float = None,
        page_number: int = None,
    ):
        """
        Returns evidences that support Variant-Disease Associations.
        Params narrowed to what the old get_vda_evidences_by_variant()/
        get_vda_evidences_by_disease() (via _get_evidences()) used.

        @variant: Union[str, List[str]]
            Variant (dbSNP Identifier) or list of variants, up to 100.
        @gene_ncbi_id / @gene_symbol: Union[str, List[str]]
            Gene identifier(s); the old "gene" param accepted either.
        @disease: Union[str, List[str]]
            Disease id(s) with vocabulary prefix, e.g. "UMLS_C0005745".
        @source: Union[str, List[str]]
            Source of the VDA.
        @min_pmYear/@max_pmYear: str
            Publication year range - corresponds to the old
            "min_year"/"max_year" params. The new API has both
            min_timestamp/max_timestamp (evidence timestamp) and
            min_pmYear/max_pmYear (publication year); mapped to the latter
            as the closer match, but this is an assumption, not confirmed.
        @min_score/@max_score: float
            Variant-disease normalized score range, in [0,1].
        @page_number: int
            Page number - corresponds to the old "limit"/"offset" params.
            The old code auto-looped through cursor-based pages
            (data["next"]) when get_all=True; the new API paginates via
            page_number instead of a cursor, and TRIAL accounts don't
            support pagination at all, so the auto-loop was removed. This
            returns a single page; call again with an incremented
            page_number for more.
        """

        url = f"{self._api_url}/vda/evidence"
        get_params = dict()

        if variant != None:
            get_params["variant"] = self._list_to_str(variant, "Variant ID", limit=100)

        if gene_ncbi_id != None:
            get_params["gene_ncbi_id"] = self._list_to_str(gene_ncbi_id, "Gene NCBI ID", limit=100)

        if gene_symbol != None:
            get_params["gene_symbol"] = self._list_to_str(gene_symbol, "Gene Symbol", limit=100)

        if disease != None:
            get_params["disease"] = self._list_to_str(disease, "Disease ID", limit=100)

        if source != None:
            get_params["source"] = source

        if min_pmYear != None:
            get_params["min_pmYear"] = min_pmYear

        if max_pmYear != None:
            get_params["max_pmYear"] = max_pmYear

        if min_score != None:
            get_params["min_score"] = str(min_score)

        if max_score != None:
            get_params["max_score"] = str(max_score)

        if page_number != None:
            get_params["page_number"] = str(page_number)

        return self._retrieve_data(url, get_params)


    # get_gda_evidences_by_gene() and get_gda_evidences_by_disease() were
    # removed, for the same reason as the VDA evidence functions above.
    # Replaced by get_gda_evidence() below, same narrowed-param approach,
    # same raw-passthrough (no namedtuple) behavior as the old code.
    def get_gda_evidence(
        self,
        gene_ncbi_id: Union[str, List[str]] = None,
        gene_symbol: Union[str, List[str]] = None,
        disease: Union[str, List[str]] = None,
        source: Union[str, List[str]] = None,
        min_pmYear: str = None,
        max_pmYear: str = None,
        min_score: float = None,
        max_score: float = None,
        page_number: int = None,
    ):
        """
        Returns evidences that support Gene-Disease Associations.
        Params narrowed to what the old get_gda_evidences_by_gene()/
        get_gda_evidences_by_disease() (via _get_evidences()) used.

        @gene_ncbi_id / @gene_symbol: Union[str, List[str]]
            Gene identifier(s); the old "gene" param accepted either.
        @disease: Union[str, List[str]]
            Disease id(s) with vocabulary prefix, e.g. "UMLS_C0005745".
        @source: Union[str, List[str]]
            Source of the GDA.
        @min_pmYear/@max_pmYear: str
            Publication year range - corresponds to the old
            "min_year"/"max_year" params (same assumption noted in
            get_vda_evidence above).
        @min_score/@max_score: float
            Gene-disease normalized score range, in [0,1].
        @page_number: int
            Page number - corresponds to the old "limit"/"offset" params
            (same auto-loop removal noted in get_vda_evidence above).
        """

        url = f"{self._api_url}/gda/evidence"
        get_params = dict()

        if gene_ncbi_id != None:
            get_params["gene_ncbi_id"] = self._list_to_str(gene_ncbi_id, "Gene NCBI ID", limit=100)

        if gene_symbol != None:
            get_params["gene_symbol"] = self._list_to_str(gene_symbol, "Gene Symbol", limit=100)

        if disease != None:
            get_params["disease"] = self._list_to_str(disease, "Disease ID", limit=100)

        if source != None:
            get_params["source"] = source

        if min_pmYear != None:
            get_params["min_pmYear"] = min_pmYear

        if max_pmYear != None:
            get_params["max_pmYear"] = max_pmYear

        if min_score != None:
            get_params["min_score"] = str(min_score)

        if max_score != None:
            get_params["max_score"] = str(max_score)

        if page_number != None:
            get_params["page_number"] = str(page_number)

        return self._retrieve_data(url, get_params)
    
    
    # get_gda_evidences_by_gene() and get_gda_evidences_by_disease() (and
    # the "gda" branch of the shared _get_evidences()) were removed. They
    # existed only to select a "by" routing mode (gene vs disease) for the
    # old path-based /gda/evidences/{by}/... endpoint. The new API has a
    # single fixed endpoint (/gda/evidence) where gene, disease, chemical,
    # etc. are all just optional query parameters on the same call, so
    # there's no more separate URL/routing per identifier type. Replaced
    # by the single get_gda_evidence() function below, which accepts all
    # identifier types at once. This also mirrors the same removal already
    # done for get_vda_evidences_by_variant()/by_disease().

    def get_dda(
        self,
        disease_1: Union[str, List[str]],
        disease_2: Union[str, List[str]],
        source: Union[str, List[str]] = None,
        page_number: int = None,
    ) -> NamedTuple(
        "DiseaseDiseaseAssociation",
        [
            ("disease1_name", str),
            ("disease2_name", str),
            ("disease1_classes_do", Tuple[str]),
            ("disease1_classes_hpo", Tuple[str]),
            ("disease1_classes_msh", Tuple[str]),
            ("disease1_classes_umls_st", Tuple[str]),
            ("disease2_classes_do", Tuple[str]),
            ("disease2_classes_hpo", Tuple[str]),
            ("disease2_classes_msh", Tuple[str]),
            ("disease2_classes_umls_st", Tuple[str]),
            ("jaccard_genes", float),
            ("pvalue_jaccard_genes", float),
            ("jaccard_variants", float),
            ("pvalue_jaccard_variants", float),
            ("source", str),
            ("ngenes1", int),
            ("ngenes2", int),
            ("nvariants1", int),
            ("nvariants2", int),
            ("shared_genes", int),
            ("shared_variants", int),
            ("diseaseid1", str),
            ("diseaseid2", str),
        ],
    ):
        """
        Returns Disease-Disease Associations between disease_1 and disease_2.
        Fields narrowed to the union of what the old
        get_ddas_that_share_genes()/get_ddas_that_share_variants() returned;
        the new API returns both gene-sharing and variant-sharing metrics
        together in a single query, so both jaccard_genes and
        jaccard_variants are always populated (the old code only returned
        one or the other depending on which function was called). The old
        duplicate fields disease1_ngenes/disease2_ngenes/disease1_nvariants/
        disease2_nvariants (which held the same values as ngenes1/ngenes2/
        nvariants1/nvariants2) were consolidated into the single generic
        fields, since the new API only exposes one copy of each.

        @disease_1 / @disease_2: Union[str, List[str]]
            Disease id(s) with vocabulary prefix, e.g. "UMLS_C0005745".
        @source: Union[str, List[str]]
            Source of the DDA.
        @page_number: int
            Page number - corresponds to the old "limit" param (100 per
            page; TRIAL accounts are capped at the top-10 results and do
            not support pagination).
        """

        disease_1 = self._list_to_str(disease_1, "Disease 1 ID", limit=100)
        disease_2 = self._list_to_str(disease_2, "Disease 2 ID", limit=100)

        url = f"{self._api_url}/dda"
        get_params = dict()

        get_params["disease_1"] = disease_1
        get_params["disease_2"] = disease_2

        if source != None:
            get_params["source"] = source

        if page_number != None:
            get_params["page_number"] = str(page_number)

        result = self._retrieve_data(url, get_params)

        if result == None:
            return None

        DiseaseDiseaseAssociation = collections.namedtuple(
            "DiseaseDiseaseAssociation",
            [
                "disease1_name",
                "disease2_name",
                "disease1_classes_do",
                "disease1_classes_hpo",
                "disease1_classes_msh",
                "disease1_classes_umls_st",
                "disease2_classes_do",
                "disease2_classes_hpo",
                "disease2_classes_msh",
                "disease2_classes_umls_st",
                "jaccard_genes",
                "pvalue_jaccard_genes",
                "jaccard_variants",
                "pvalue_jaccard_variants",
                "source",
                "ngenes1",
                "ngenes2",
                "nvariants1",
                "nvariants2",
                "shared_genes",
                "shared_variants",
                "diseaseid1",
                "diseaseid2",
            ],
        )

        for index, entry in enumerate(result):
            result[index] = DiseaseDiseaseAssociation(
                self._get_string(entry.get("disease1_Name")),
                self._get_string(entry.get("disease2_Name")),
                tuple(entry["disease1_Classes_DO"]) if entry.get("disease1_Classes_DO") else None,
                tuple(entry["disease1_Classes_HPO"]) if entry.get("disease1_Classes_HPO") else None,
                tuple(entry["disease1_Classes_MSH"]) if entry.get("disease1_Classes_MSH") else None,
                tuple(entry["disease1_Classes_UMLS_ST"]) if entry.get("disease1_Classes_UMLS_ST") else None,
                tuple(entry["disease2_Classes_DO"]) if entry.get("disease2_Classes_DO") else None,
                tuple(entry["disease2_Classes_HPO"]) if entry.get("disease2_Classes_HPO") else None,
                tuple(entry["disease2_Classes_MSH"]) if entry.get("disease2_Classes_MSH") else None,
                tuple(entry["disease2_Classes_UMLS_ST"]) if entry.get("disease2_Classes_UMLS_ST") else None,
            
                self._get_float(entry.get("jaccard_genes")),
                self._get_float(entry.get("pvalue_jaccard_genes")),
                self._get_float(entry.get("jaccard_variants")),
                self._get_float(entry.get("pvalue_jaccard_variants")),
                self._get_string(entry.get("source")),
                self._get_int(entry.get("ngenes_diseaseID_1")),
                self._get_int(entry.get("ngenes_diseaseID_2")),
                self._get_int(entry.get("nvariants_diseaseID_1")),
                self._get_int(entry.get("nvariants_diseaseID_2")),
                self._get_int(entry.get("shared_genes")),
                self._get_int(entry.get("shared_variants")),
                self._get_string(entry.get("disease1_UMLSCUI")),
                self._get_string(entry.get("disease2_UMLSCUI")),
            )

        return result
    # _get_evidences() was removed. It tried to share one generic function
    # between GDA and VDA evidence retrieval using "of"/"by" flags for
    # path-based routing (e.g. /vda/evidences/variant/...). The new API
    # has separate, fixed endpoints (/gda/evidence, /vda/evidence) with
    # much richer, endpoint-specific parameters (e.g. chromcoord, hgvsc,
    # hgvsp, min_polyphen for VDA; uniprot_id, nct_phase for GDA), so a
    # single shared function is no longer a good fit. 
   
    

    def _list_to_str(
        self, list_obj: List[str], name: str, limit: Optional[int] = None
    ) -> List[str]:
        """
        Joins the list object and returns a string

        @list_obj : List[str]
            List object to be processed
        @name: str
            Name of items, like Gene ID or something
        @limit: Optional[int]
            Maximum number of list items
        """

        if isinstance(list_obj, list):
            if limit != None and len(list_obj) > limit:
                _log(
                    f"DisGeNET maximum length of {name}'s are {limit}, "
                    f"first {limit} {name}'s will be used."
                )

                return ",".join(list_obj[:limit])

            return ",".join(list_obj)

        return list_obj

    # Re-enabling this decorator: it was previously commented out, which let
    # _retrieve_data run with no valid api_key and silently send
    # "Authorization: Bearer None" to the server. Now that authenticate()
    # is fixed for the new API, this guard is safe to use again.
    @_if_authenticated
    @_delete_cache
    def _retrieve_data(
        self, url: str, get_params: Union[List[str], Dict[str, str]]
    ) -> List[Dict[str, str]]:
        """
        Retrieves the data with given request body

        @url : str
            Query url
        @get_params: Union[List[str], Dict[str, str]]
            GET request parameters
        """

        headers = ["accept: */*", f"Authorization: Bearer {self._api_key}"]

        # Any param whose value is a Python list (e.g. source,
        # dda_relation, dis_class_list) is expanded into repeated
        # "key=value" entries, one per list item. This replaces the old
        # code, which only special-cased a single param named
        # "disease_class" (which no longer exists in the new API; it's
        # now "dis_class_list") and would otherwise have serialized any
        # other list param as a literal Python list string, e.g.
        # "source=['CURATED', 'CLINVAR']", which the DisGeNET API does
        # not parse correctly.
        get_params_list = []
        for key, value in get_params.items():
            if isinstance(value, list):
                get_params_list.extend([f"{key}={item}" for item in value])
            else:
                get_params_list.append(f"{key}={value}")

        get_params = get_params_list

        c = curl.Curl(url=url, get=get_params, req_headers=headers)

        if c.status == 0 or c.status == 200:
            result = c.result
            result = json.loads(result)

            # The new API wraps the actual records in a "payload" field
            # alongside status/paging/warnings metadata. The old code
            # returned the raw parsed JSON as-is, which worked because
            # the old API's response body WAS the list of records
            # directly. Every function below expects a plain list of
            # record dicts, so we extract "payload" here.
            return result.get("payload")

        _log(f"DisGeNET: an error occurred with the code {c.status}")
        _log(f"DisGeNET response body: {repr(c.result)}")
        _log(f"DisGeNET curl object: {vars(c)}")

    def _get_int(self, str_obj) -> int:
        """
        Returns an int if the object is not None

        @str_obj : str
            String to be processed
        """

        if str_obj != None and not isinstance(str_obj, int):
            return int(str_obj)

        return str_obj

    def _get_float(self, str_obj) -> float:
        """
        Returns a float if the object is not None

        @str_obj : str
            String to be processed
        """

        if str_obj != None and not isinstance(str_obj, float):
            return float(str_obj)

        return str_obj

    def _get_string(self, obj) -> str:
        """
        Returns a string if the object is not None

        @obj : object
            Object to be processed
        """

        if obj != None and not isinstance(obj, str):
            return str(obj).strip()

        return obj

    def _get_tuple(self, str_obj: str, delim: str) -> Tuple[str]:
        """
        Returns a splitted tuple with given delimiter

        @str_obj : str
            String object to be processed
        @delim : str
            Char to split the str_obj
        """

        if str_obj != None and not isinstance(str_obj, tuple):
            return tuple([item.strip() for item in str_obj.split(delim)])

        return str_obj


@DisgenetApi._delete_cache
def variant_gene_mappings(
    api: "DisgenetApi",
    gene_ncbi_ids: List[str],
    batch_size: int = 10,
) -> Dict[str, "VariantGeneMapping"]:
    """
    Builds the same {snpId: [VariantGeneMapping(geneId, geneSymbol,
    sourceIds), ...]} structure the old bulk-download version produced,
    but queries the new DisGeNET API instead (no bulk mapping file
    exists anymore).

    Unlike the old version (which took no arguments and downloaded the
    full DisGeNET variant-gene mapping file), this now requires a list
    of known NCBI Gene IDs to query against - e.g. from
    uniprot_adapter.py's xref_geneid field (same source used for
    disgenet_annotations()). This queries /entity/variant using
    gene_ncbi_id, which returns each matching variant's own
    variantToGenes field - a list of {geneNcbiID, symbolOfGene,
    geneEnsemblID, sources} already grouped per variant, so no manual
    grouping across rows is needed (unlike the old bulk file, where the
    same snpId/geneId pair could appear on multiple rows with different
    sourceId values that had to be accumulated by hand).

    @api: DisgenetApi
        An already-authenticated DisgenetApi instance.
    @gene_ncbi_ids: List[str]
        Known NCBI Gene (Entrez) IDs to query, e.g. ["672", "675", ...].
    @batch_size: int
        Number of gene IDs to send per request. Kept small (10) by
        default to stay safe under TRIAL account limits while testing;
        raise this once running under a full academic account.
    """

    VariantGeneMapping = collections.namedtuple(
        "VariantGeneMapping",
        [
            "geneId",
            "geneSymbol",
            "sourceIds",
        ],
    )

    mapping = dict()

    for i in range(0, len(gene_ncbi_ids), batch_size):
        batch = gene_ncbi_ids[i : i + batch_size]
        page_number = 0

        while True:
            url = f"{api._api_url}/entity/variant"
            get_params = {
                "gene_ncbi_id": batch,
                "page_number": str(page_number),
            }

            result = api._retrieve_data(url, get_params)

            if not result:
                break

            for entry in result:
                snp_id = entry.get("strID")
                variant_to_genes = entry.get("variantToGenes")

                if snp_id == None or not variant_to_genes:
                    continue

                if snp_id not in mapping:
                    mapping[snp_id] = []

                for vtg in variant_to_genes:
                    gene_id = vtg.get("geneNcbiID")
                    gene_symbol = vtg.get("symbolOfGene")
                    sources = vtg.get("sources")

                    mapping[snp_id].append(
                        VariantGeneMapping(
                            api._get_string(gene_id) if gene_id != None else None,
                            gene_symbol,
                            tuple(sources) if sources else None,
                        )
                    )

            if len(result) < 100:
                break

            page_number += 1

    return mapping


@DisgenetApi._delete_cache
def disease_id_mappings(
    api: "DisgenetApi",
    mondo_ids: List[str],
    batch_size: int = 10,
) -> Dict[str, "DiseaseIdMapping"]:
    """
    Builds the same {diseaseId: DiseaseIdMapping(name, vocabularies)}
    structure the old bulk-download version produced, but queries the
    new DisGeNET API instead (no bulk download endpoint exists anymore).

    Unlike the old version (which took no arguments and downloaded the
    full DisGeNET disease list), this now requires a list of known
    MONDO disease IDs to query against - e.g. from disease_adapter.py's
    MONDO ontology download. Each ID must be in the API's expected
    format, e.g. "MONDO_0007254".

    This is a standalone function (not a DisgenetApi method) so it
    doesn't change how the class itself is called elsewhere; it takes
    an already-authenticated DisgenetApi instance as a parameter
    instead, and calls _retrieve_data directly rather than going
    through get_gda_summary(), because get_gda_summary()'s narrowed
    output (matching the old API's fields) doesn't include
    diseaseVocabularies, which is exactly the field this function needs.

    Note on data shape: /gda/summary returns one row per matching gene
    for a queried disease, and every row repeats the same disease-level
    diseaseVocabularies value. We only need that value once per
    disease, so once a disease_id has been recorded we skip it on
    subsequent rows - this does not lose any information, since this
    function was never collecting gene data in the first place (same
    as the old bulk-file version, which only ever returned name +
    vocabularies per disease, no gene info).

    @api: DisgenetApi
        An already-authenticated DisgenetApi instance.
    @mondo_ids: List[str]
        Known MONDO disease IDs to query, e.g. ["MONDO_0007254", ...].
    @batch_size: int
        Number of disease IDs to send per request. Kept small (10) by
        default to stay safe under TRIAL account limits while testing;
        raise this once running under a full academic account.
    """

    Vocabulary = collections.namedtuple(
        "Vocabulary",
        [
            "vocabulary",
            "code",
            "vocabularyName",
        ],
    )

    DiseaseIdMapping = collections.namedtuple(
        "DiseaseIdMapping",
        [
            "name",
            "vocabularies",
        ],
    )

    mapping = dict()

    for i in range(0, len(mondo_ids), batch_size):
        batch = mondo_ids[i : i + batch_size]
        page_number = 0

        while True:
            url = f"{api._api_url}/gda/summary"
            get_params = {
                "disease": batch,
                "page_number": str(page_number),
            }

            result = api._retrieve_data(url, get_params)

            if not result:
                break

            for entry in result:
                disease_id = entry.get("diseaseUMLSCUI")
                raw_vocabs = entry.get("diseaseVocabularies")

                if disease_id == None or not raw_vocabs:
                    continue

                # A disease query returns one row per matching gene;
                # every row repeats the same disease-level vocabulary
                # info, so we only process the first occurrence.
                if disease_id in mapping:
                    continue

                name = entry.get("diseaseName")
                vocab_list = []

                for raw in raw_vocabs:
                    parts = raw.split("_", 1)

                    if len(parts) != 2:
                        continue

                    vocabulary, code = parts
                    vocab_list.append(
                        Vocabulary(
                            vocabulary,
                            code,
                            # The new API only returns short prefixes
                            # (e.g. "MESH", "MONDO"), not the full
                            # vocabulary name the old bulk file had
                            # (e.g. "Medical Subject Headings"). No
                            # equivalent field exists in the new API
                            # or the web interface (checked disease
                            # detail page - it also only shows short
                            # labels like "MeSH Disease Class", not
                            # the full vocabulary name), so this is
                            # left as None rather than guessing.
                            None,
                        )
                    )

                mapping[disease_id] = DiseaseIdMapping(name, tuple(vocab_list))

            if len(result) < 100:
                break

            page_number += 1

    return mapping


@DisgenetApi._delete_cache
def disgenet_annotations(
    api: "DisgenetApi",
    gene_ncbi_ids: List[str],
    dataset: str = "curated",
    batch_size: int = 10,
) -> Dict[str, set]:
    """
    Builds the same {uniprot_id: {DisGeNetAnnotation(...), ...}} structure
    the old bulk-download version produced, but queries the new DisGeNET
    API instead (no bulk annotation files exist anymore - confirmed via
    AIDA that curated_gene_disease_associations is no longer available).

    Unlike the old version (which took no gene list and downloaded the
    full DisGeNET file, then mapped gene symbols to UniProt via pypath's
    own mapping.map_name()), this now requires a list of known NCBI Gene
    IDs to query against - e.g. from uniprot_adapter.py's xref_geneid
    field. The gene-symbol-to-UniProt translation step is no longer
    needed: the new API already returns UniProt IDs directly in
    geneProteinStrIDs, so we read those instead of calling pypath's
    mapping module.

    This calls _retrieve_data directly rather than get_gda_summary(),
    because get_gda_summary()'s narrowed output (matching the old
    _get_gdas() fields) doesn't include numPMIDs or
    numDBSNPsupportingAssociation, which map to this function's
    nof_pmids/nof_snps fields and aren't part of the old GDA shape.

    @api: DisgenetApi
        An already-authenticated DisgenetApi instance.
    @gene_ncbi_ids: List[str]
        Known NCBI Gene (Entrez) IDs to query, e.g. ["672", "675", ...].
    @dataset: str
        Only "curated" and "all" are supported (mapped to
        source=["CURATED"] and no source filter, respectively). The old
        "literature" and "befree" datasets have no clear 1:1 equivalent
        in the new API's source list (BEFREE no longer exists; the
        closest match, TEXTMINING_HUMAN/TEXTMINING_MODELS, isn't the
        same grouping) - flagging this rather than guessing a mapping.
    @batch_size: int
        Number of gene IDs to send per request. Kept small (10) by
        default to stay safe under TRIAL account limits while testing;
        raise this once running under a full academic account.
    """

    dataset_source_map = {
        "curated": ["CURATED"],
        "all": None,
    }

    if dataset not in dataset_source_map:
        _log(
            f"DisGeNET dataset='{dataset}' is not supported. Only 'curated' and "
            "'all' are mapped to the new API's source filter; 'literature' "
            "and 'befree' have no clear equivalent in the new source list."
        )
        return None

    source = dataset_source_map[dataset]

    DisGeNetAnnotation = collections.namedtuple(
        "DisGeNetAnnotation",
        [
            "disease",
            "type",
            "score",
            "dsi",
            "dpi",
            "nof_pmids",
            "nof_snps",
            "source",
        ],
    )

    data = collections.defaultdict(set)

    for i in range(0, len(gene_ncbi_ids), batch_size):
        batch = gene_ncbi_ids[i : i + batch_size]
        page_number = 0

        while True:
            url = f"{api._api_url}/gda/summary"
            get_params = {
                "gene_ncbi_id": batch,
                "page_number": str(page_number),
            }

            if source != None:
                get_params["source"] = source

            result = api._retrieve_data(url, get_params)

            if not result:
                break

            for entry in result:
                uniprot_ids = entry.get("geneProteinStrIDs")

                if not uniprot_ids:
                    continue

                disease = entry.get("diseaseName")
                disease_type = entry.get("diseaseType")
                score = entry.get("score")
                dsi = entry.get("geneDSI")
                dpi = entry.get("geneDPI")
                nof_pmids = entry.get("numPMIDs")
                nof_snps = entry.get("numDBSNPsupportingAssociation")

                # The new /gda/summary response does not include a
                # per-record source field (unlike the old API's
                # GeneDiseaseAssociation, which had one). When a single
                # source filter was requested (dataset="curated"), that
                # filter is the only source these aggregated results
                # could have come from, so we record it as such. For
                # dataset="all" (no filter), the true per-record source
                # is unknown, so this is left as None.
                record_source = tuple(source) if source != None else None

                annotation = DisGeNetAnnotation(
                    disease,
                    disease_type,
                    api._get_float(score) if score != None else None,
                    api._get_float(dsi) if dsi != None else None,
                    api._get_float(dpi) if dpi != None else None,
                    api._get_int(nof_pmids) if nof_pmids != None else None,
                    api._get_int(nof_snps) if nof_snps != None else None,
                    record_source,
                )

                for uniprot in uniprot_ids:
                    data[uniprot].add(annotation)

            if len(result) < 100:
                break

            page_number += 1

    return dict(data)