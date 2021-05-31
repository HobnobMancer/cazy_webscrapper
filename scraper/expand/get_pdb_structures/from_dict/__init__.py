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
"""Submodule for retieving structure files for proteins in dict (JSON file)"""


import logging
import math

from tqdm import tqdm

from scraper.expand import get_genbank_sequences, get_cazy_dict, get_qualifying_proteins
from scraper.expand.get_genbank_sequences.ncbi import query_entrez
from scraper.utilities import file_io, parse_configuration


def structures_for_proteins_from_dict(date_today, args):
    """Coordinate retrievel of protein sequences for proteins in a CAZy dict (JSON file).
    
    :param date_today: str, date script was invoked, used for naming files
    :param args: cmd-line args parser
    
    Return nothing.
    """
    logger = logging.getLogger(__name__)

    # create output directory
    if args.outdir is None:
        # save structure files to the cwd
        outdir = os.getcwd()
    else:
        outdir = args.outdir
        file_io.make_output_directory(outdir, args.force, args.nodelete)

    # retrieve configuration data, as a dict of CAZy classes and families to retrieve seqs for
    parse_configuration_path = parse_configuration.__file__
    config_dict = parse_configuration.parse_configuration_for_cazy_dict(
        parse_configuration_path,
        args,
    )

    # retrieve dict of CAZy family classifications of proteins
    cazy_dict = get_cazy_dict(args)

    # retrieve sequences for proteins in the cazy_dict matching user specified critera
    protein_list = get_qualifying_proteins(cazy_dict, config_dict)

    try:
        protein_list.remove("NA")  # incase 'NA' was added when scraping CAZy
    except ValueError:
        pass
    
    logger.warning(f"Retrieving structures for {len(protein_list)} proteins")

    # break up protein_list into multiple, smaller lists to reduce burden on PDB

    accessions_lists_for_individual_queries = []

    for accession_list in tqdm(
        get_genbank_sequences.get_accession_chunks(protein_list, args.epost),
        desc="Batch retrieving structures from PDB",
        total=(math.ceil(len(protein_list) / args.epost)),
    ):
        try:
            query_entrez.get_sequences_for_dict(accession_list, date_today, args)
        except RuntimeError as err:  # typically Some IDs have invalid value and were omitted.
            logger.warning(
                "RuntimeError raised for accession list. Will query accessions individualy after"
            )
            with open("legihton_error.txt", "a") as fh:
                fh.write(f"{err}\n{str(accession_list)}")
            accessions_lists_for_individual_queries.append(accession_list)

    if len(accessions_lists_for_individual_queries) != 0:
        for accession_list in tqdm(
            accessions_lists_for_individual_queries,
            desc="Performing individual queries to parse GenBank accessions without records",
        ):
            for accession in tqdm(accession_list, desc="Retrieving individual sequences"):
                try:
                    query_entrez.get_sequences_for_dict([accession], date_today, args)
                except RuntimeError as err:
                    logger.warning(
                        f"Querying NCBI for {accession} raised the following RuntimeError:\n"
                        f"{err}"
                    )
    return
