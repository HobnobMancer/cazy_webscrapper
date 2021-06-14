#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
# Author:
# Emma E. M. Hobbs
#
# Contact
# eemh1@st-andrews.ac.uk
#
# Emma E. M. Hobbs,
# Biomolecular Sciences Building,
# University of St Andrews,
# North Haugh Campus,
# St Andrews,
# KY16 9ST
# Scotland,
# UK
#
# The MIT License
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""Submodule for retieving sequences for proteins in dict (JSON file)"""


import logging
import json
import math
import sys

from tqdm import tqdm

from scraper.expand import get_accession_chunks
from scraper.expand.get_genbank_sequences.ncbi import query_entrez
from scraper.utilities import parse_configuration


def sequences_for_proteins_from_dict(date_today, args):
    """Coordinate retrievel of protein sequences for proteins in a CAZy dict (JSON file).
    
    :param date_today: str, date script was invoked, used for naming files
    :param args: cmd-line args parser
    
    Return nothing.
    """
    logger = logging.getLogger(__name__)

    # retrieve dict of CAZy family classifications of proteins
    cazy_dict = get_cazy_dict(args)

    # retrieve configuration data, as a dict of CAZy classes and families to retrieve seqs for
    config_dict = parse_configuration.parse_configuration_for_cazy_dict(args)

    # retrieve sequences for proteins in the cazy_dict matching user specified critera
    protein_list = get_qualifying_proteins(cazy_dict, config_dict)

    try:
        protein_list.remove("NA")  # incase 'NA' was added when scraping CAZy
    except ValueError:
        pass
    
    logger.warning(f"Retrieving sequences for {len(protein_list)} proteins")

    # break up protein_list into multiple, smaller lists for batch querying Entrez
    # batches of greater than 200 can be rejected by Entrez during busy periods
    # args.epost=size of chunks

    accessions_lists_for_individual_queries = []

    for accession_list in tqdm(
        get_accession_chunks(protein_list, args.epost),
        desc="Batch retrieving sequences from NCBI",
        total=(math.ceil(len(protein_list) / args.epost)),
    ):
        try:
            query_entrez.get_sequences_for_dict(accession_list, args)

        except RuntimeError as err:  # typically Some IDs have invalid value and were omitted.
            logger.error(
                "RuntimeError raised for accession list. Will query accessions individualy after\n"
                "RuntimeError:\n"
                f"{err}"
            )
            accessions_lists_for_individual_queries.append(accession_list)

    if len(accessions_lists_for_individual_queries) != 0:
        for accession_list in tqdm(
            accessions_lists_for_individual_queries,
            desc="Performing individual queries to parse GenBank accessions without records",
        ):
            for accession in tqdm(accession_list, desc="Retrieving individual sequences"):
                try:
                    query_entrez.get_sequences_for_dict([accession], args)
                except RuntimeError as err:
                    logger.error(
                        f"Querying NCBI for {accession} raised the following RuntimeError:\n"
                        f"{err}"
                    )

    return


def get_cazy_dict(args):
    """Retrieve dict of CAZy family annotations of proteins."""
    logger = logging.getLogger(__name__)

    try:
        with open(args.dict, "r") as fh:
            cazy_dict = json.load(fh)

    except FileNotFoundError:
        logger.error(
            "Did not find the local CAZy dict (JSON) file.\n"
            "Check the path is correct.\n"
            "Terminating programme"
        )
        sys.exit(1)

    return cazy_dict


def get_qualifying_proteins(cazy_dict, config_data):
    """Identify proteins to retrieve sequences for, those that meet at least one config criteria.

    :param cazy_dict: dict of proteins catalogued in CAZy
        Keyed by protein GenBank accession, valued by list of CAZy family annotations
    :param config_data: dict of two sets: CAZy classes and families to retrieve seqs for

    Return list of proteins (1 protein = 1 GenBank acession) to retrieve seqs for.
    """
    logger = logging.getLogger(__name__)

    logger.info("Retrieving GenBank accessions of proteins matching provided criteria")

    proteins = set()  # proteins to retrieve sequences for, set prevents duplicates

    if (len(list(config_data["classes"])) == 0) and (len(list(config_data["families"])) != 0):
        # check only in family configuration data
        for protein in tqdm(cazy_dict, desc="Identifying CAZymes matching config data"):
            for fam in cazy_dict[protein]:
                if fam in config_data["families"]:
                    proteins.add(protein)
                    continue
    
    elif (len(list(config_data["classes"])) != 0) and (len(list(config_data["families"])) == 0):
        # check only class configuration data
        for protein in tqdm(cazy_dict, desc="Identifying CAZymes matching config data"):
            for cazy_class in config_data["classes"]:
                for fam in cazy_dict[protein]:
                    if fam.startswith(cazy_class):
                        proteins.add(protein)
                        continue
    
    else:
        # check both configuration data
        for protein in tqdm(cazy_dict, desc="Identifying CAZymes matching config data"):
            for fam in cazy_dict[protein]:
                if fam in config_data["families"]:
                    proteins.add(protein)
                    continue
            for cazy_class in config_data["classes"]:
                for fam in cazy_dict[protein]:
                    if fam.startswith(cazy_class):
                        proteins.add(protein)
                        continue
    
    return list(proteins)
