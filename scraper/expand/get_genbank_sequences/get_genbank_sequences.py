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

# The MIT License
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""Retrieve proteins sequences from GenBank and populate the local database and write to FASTA"""


import logging
import os
import sys

import pandas as pd

from datetime import datetime
from pathlib import Path
from typing import List, Optional

from Bio import Entrez

from scraper.expand.get_genbank_sequences import (
    from_dict,
    from_sql_db,
)
from scraper.utilities import config_logger, file_io
from scraper.utilities.parsers import build_genbank_sequences_parser


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Set up programme and initate run."""
    start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # used in terminating message
    start_time = pd.to_datetime(start_time)
    date_today = datetime.now().strftime("%Y/%m/%d")  # used as seq_update_date in the db

    # parse cmd-line arguments
    if argv is None:
        parser = build_genbank_sequences_parser()
        args = parser.parse_args()
    else:
        args = build_genbank_sequences_parser(argv).parse_args()

    if logger is None:
        logger = logging.getLogger(__package__)
        config_logger(args)
    
    Entrez.email = args.email

    if args.fasta is not None:
        logger.info("Enabled writing out FASTA files")

        # check if passed a path to a file or dir
        if len((args.fasta).suffix) == 0:  # passed path to a dir
            logger.info(f"Building directory for FASTA files: {args.fasta}")
            logger.info(f"Writing out to dir if already exists: {args.force}")
            logger.info(f"Nucking dir for the FASTA file: {args.nodelete}")
            file_io.make_output_directory(args.fasta, args.force, args.nodelete)
        
        else:  # passed path to file
            # check if need to create path to dir
            logger.info("Provided path to a FASTA file")
            logger.info("Compiling path to directory for the FASTA file")
            dir_path = Path("/".join((args.fasta).parts[:-1]))

            if dir_path == Path('.'):
                logger.info(f"Writing out FASTA file {args.fasta} to the current working directory")
            
            else:
                logger.info(f"Building directory for the FASTA file: {dir_path}")
                logger.info(f"Writing out to FASTA dir if already exists: {args.force}")
                logger.info(f"Nucking dir for the FASTA file: {args.nodelete}")
            
                file_io.make_output_directory(dir_path, args.force, args.nodelete)

    if args.fasta_only is not None:
        logger.info(
            "Enabled writing out FASTA files ONLY\n"
            "Seqs will be retrieved irrespective of the local CAZyme db seq storage status\n"
            "Annd seqs will be written to FASTA files ONLY, "
            "they will NOT be written to the local db"
        )

        # check if passed a path to a file or dir
        if len((args.fasta_only).suffix) == 0:  # passed path to a dir
            logger.info(f"Building directory for FASTA files: {args.fasta_only}")
            logger.info(f"Writing out to dir if already exists: {args.force}")
            logger.info(f"Nucking dir for the FASTA file: {args.nodelete}")
            file_io.make_output_directory(args.fasta_only, args.force, args.nodelete)
        
        else:  # passed path to file
            # check if need to create path to dir
            logger.info("Provided path to a FASTA file")
            logger.info("Compiling path to directory for the FASTA file")
            dir_path = Path("/".join((args.fasta_only).parts[:-1]))

            if dir_path == Path('.'):
                logger.info(
                    f"Writing out FASTA file {args.fasta_only} to the current working directory"
                )
            
            else:
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

    if str(args.database).endswith(".json"):
        logger.info("CAZy dictionary (JSON file) provided")
        # move to script that handles retrieving sequences for proteins in dict (JSON file)
        err = from_dict.sequences_for_proteins_from_dict(date_today, args)

    else:
        logger.info("CAZy database provided")
        # move to script that handles retrieving sequences for proteins in a SQL database
        err = from_sql_db.sequences_for_proteins_from_db(date_today, args)

    if (args.blastdb is not None) and (err is None):  # build a local BLAST database
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


if __name__ == "__main__":
    main()
