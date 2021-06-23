#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
# Author:
# Emma E. M. Hobbs

# Contact
# eemh1@st-andrews.ac.uk

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
#
# Bio.PDB reference:
# Hamelryck, T., Manderick, B. (2003) PDB parser and structure class 
# implemented in Python. Bioinformatics 19: 2308–2310
"""Retrieve protein sequences from the local CAZyme database and write to FASTA files"""


import logging
import sys

import pandas as pd

from datetime import datetime
from pathlib import Path
from typing import List, Optional

from tqdm import tqdm

from scraper.expand.get_genbank_sequences import from_sql_db
from scraper.expand.query_sql_db import query_extract_seq
from scraper.sql.sql_orm import get_db_session
from scraper.utilities import config_logger, file_io, parse_configuration
from scraper.utilities.parsers import build_extract_sequences_parser


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Set up programme and initate run."""
    start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # used in terminating message
    start_time = pd.to_datetime(start_time)

    # parse cmd-line arguments
    if argv is None:
        parser = build_extract_sequences_parser()
        args = parser.parse_args()
    else:
        args = build_extract_sequences_parser(argv).parse_args()

    if logger is None:
        logger = logging.getLogger(__package__)
        config_logger(args)

    if (args.fasta is None) and (args.blastdb is None):
        logger.error(
            "Please use --fasta and/or --blastdb, otherwise the mode of operation is unknown\n"
            "Terminating program"
        )
        sys.exit(1)

    if args.fasta is not None:
        if (args.fasta == Path('.')) or (str(args.fasta) == '.'):
                logger.info(f"Writing out FASTA file {args.fasta} to the current working directory")

        # check if passed a path to a file or dir
        elif len((args.fasta).suffix) == 0:  # passed path to a dir
            logger.info(f"Building directory for FASTA files: {args.fasta}")
            logger.info(f"Writing out to dir if already exists: {args.force}")
            logger.info(f"Nucking dir for the FASTA file: {args.nodelete}")
            file_io.make_output_directory(args.fasta, args.force, args.nodelete)
        
        else:  # passed path to file
            # check if need to create path to dir
            logger.info("Provided path to a FASTA file")
            logger.info("Compiling path to directory for the FASTA file")
            dir_path = Path("/".join((args.fasta).parts[:-1]))

            logger.info(f"Building directory for the FASTA file: {dir_path}")
            logger.info(f"Writing out to FASTA dir if already exists: {args.force}")
            logger.info(f"Nucking dir for the FASTA file: {args.nodelete}")
        
            file_io.make_output_directory(dir_path, args.force, args.nodelete)

    if args.blastdb is not None:  # build directory to store FASTA file for BLAST db
        logger.info("Enabled creating BLAST database containing retrieved protein sequences")

        if len((args.blastdb).parts) != 1:  # == 1 when database will be built in the cwd
            logger.info("Compiling path to directory for BLAST db")
            dir_path = Path("/".join((args.blastdb).parts[:-1]))

            logger.info(f"Building directory for BLAST db: {dir_path}")
            logger.info(f"Writing out to BLAST db dir if already exists: {args.force}")
            logger.info(f"Nucking dir for the BLAST db: {args.nodelete}")
            file_io.make_output_directory(dir_path, args.force, args.nodelete)

    # get db session
    session = get_db_session(args)

    # get configuration
    (
        config_dict,
        taxonomy_filters,
        kingdoms,
        ec_filters,
    ) = parse_configuration.parse_configuration_for_cazy_database(args)

    genbank_records = get_genbank_records(
        args,
        session,
        config_dict,
        taxonomy_filters,
        kingdoms,
        ec_filters,
    )

    for record in tqdm(genbank_records, desc="Writing out seqs to FASTA"):
        file_io.write_out_extract_seq_to_fasta(
            record[0].sequence, 
            record[0].genbank_accession,
            args,
        )

    # genbank_accessions = [query_result.genbank_accession for query_result in filtered_query_results]

    if (args.blastdb is not None):  # build a local BLAST database
        for record in tqdm(genbank_records, desc="Writing out seqs to FASTA for BLAST db"):
            file_io.write_extracted_fasta_for_db(
                record[0].sequence,
                record[0].genbank_accession,
                args,
            )

        # write out fasta files
        file_io.build_blast_db(args)
    
    end_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # used in terminating message
    end_time = pd.to_datetime(start_time)
    end_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    end_time = pd.to_datetime(end_time)
    total_time = end_time - start_time

    logger.info(
        "Finished retrieving protein sequences.\n"
        "Terminating program.\n"
        f"Scrape initated at {start_time}\n"
        f"Scrape finished at {end_time}\n"
        f"Total run time: {total_time}"
    )

    print(
        "=====================cazy_webscraper-expand-get_genank_sequences=====================\n"
        "Finished retrieving protein sequences.\n"
        "Terminating program.\n"
        f"Scrape initated at {start_time}\n"
        f"Scrape finished at {end_time}\n"
        f"Total run time: {total_time}\n"
    )


def get_genbank_records(
    args,
    session,
    config_dict,
    taxonomy_filters,
    kingdoms,
    ec_filters,
):
    """Retrieve the GenBank records from a local CAZyme database for all proteins from which 
    a protein sequence is to be extracted.

    :param args: cmd-line argument parser
    :param session: open SQLite db session
    :param config_dict: dict, defines CAZy classes and families to get sequences for
    :param taxonomy_filters: set of genera, species and strains to retrieve sequences for
    :param kingdoms: set of taxonomy Kingdoms to retrieve sequences for
    :param ec_filters: set of EC numbers annotations CAZymes must have at least one to retrieve
        a sequence

    Return a list of GenBank accessions, containing no duplicate GenBank accessions
    """
    logger = logging.getLogger(__name__)

    if config_dict is not None:  # there are specific CAZy classes/families to retrieve sequences for

        if args.primary:
            logger.warning(
                "Retrieving sequences for PRIMARY GenBank accessions that:\n"
                "belong to specific CAZy classes/families AND\n"
                "have a sequence in the db"
            )
            (
                genbank_query_class,
                genbank_query_family,
            ) = query_extract_seq.get_prim_gnbk_acc_from_clss_fams_with_seq(
                session,
                config_dict,
            )

        else:
            logger.warning(
                "Retrieving sequences for ALL GenBank accessions that:\n"
                "belong to specific CAZy classes/families AND\n"
                "have a sequence in the db"
            )
            (
                genbank_query_class,
                genbank_query_family,
            ) = query_extract_seq.get_all_gnbk_acc_from_clss_fams_with_seq(
                session,
                config_dict,
            )

        query_results = genbank_query_class + genbank_query_family
        logger.info(
            f"Retrieved {len(query_results)} records from the specified CAZy classes and families\n"
            "prior to applying any taxonomy and EC number filters"
        )

    else:

        if args.primary:
            logger.warning(
                "Retrieving sequences for all PRIMARY GenBank accessions that\n"
                "have a sequence in the db"
            )
            query_results = query_extract_seq.get_prim_genbank_accessions_with_seq(session)

        else:
            logger.warning(
                "Retrieving sequences for ALL GenBank accessions that\n"
                "have a sequence in the db"
            )
            query_results = query_extract_seq.get_genbank_accessions_with_seq(session)
            
        logger.info(
            f"Retrieved {len(query_results)} records from the local CAZyme database\n"
            "Prior to apply any taxonomy and EC number filters"
        )

    if args.accessions is not None:
        logger.info(
            "Retrieving records from the local CAZyme database for accessions "
            "provided using args.accessions"
        )
        accessions_list = (args.accessions).split(",")
        user_accessions = query_extract_seq.get_user_accessions_with_seq(accessions_list, session)

        query_results += user_accessions
    
    if args.accessions_path is not None:
        logger.info(
            f"Retrieveing accessions listed in\n{args.accessions_path}"
        )
        accessions_list = file_io.get_accessions_from_file(args)

        logger.info(
            "Retrieving records from the local CAZyme database for accessions "
            "provided using args.accessions_path"
        )
        user_accessions = query_extract_seq.get_user_accessions_with_seq(
            accessions_list,
            session,
        )

        query_results += user_accessions

    # check if any records were retrived from the querying of the local CAZyme database
    if len(query_results) == 0:
        logger.warning(
            "Retrieved no records from the local CAZyme database mathcing provided criteria\n"
            "and with sequences"
        )
        return

    # apply taxonomic and EC number filters
    logger.info(
        "Applying any provided taxonomic and EC number filters to records"
        "retrieved from the local CAZyme database"
    )

    filtered_query_results = from_sql_db.parse_genbank_query(
        query_results,
        taxonomy_filters,
        kingdoms,
        ec_filters,
        session,
    )

    return list(set(filtered_query_results))  # prevent quering the same accession multiple times


if __name__ == "__main__":
    main()
