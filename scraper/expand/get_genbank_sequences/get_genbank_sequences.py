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
from scraper.expand.get_genbank_sequences.ncbi import blast_db
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
        logger = logging.getLogger(__name__)
        config_logger(args)
    
    Entrez.email = args.email

    if args.fasta is not None:
        logger.info("Enabled writing out FASTA files")
        if len((args.fasta).parts) != 1:
            logger.info("Compiling path to directory for FASTA files")
            dir_path = Path("/".join((args.fasta).parts))
            
            logger.info(f"Building directory for FASTA files: {dir_path}")
            file_io.make_output_directory(dir_path, args.force, args.nodelete)

    if args.blastdb is not None:  # build directory to store FASTA file for BLAST db
        logger.info("Enabled creating BLAST database containing retrieved protein sequences")
        if len((args.blastdb).parts) != 1:
            logger.info("Compile path to directory for BLAST db")
            dir_path = Path("/".join((args.blastdb).parts))

            logger.info(f"Building directory for BLAST db: {dir_path}")
            file_io.make_output_directory(dir_path, args.force, args.nodelete)

    if str(args.database).endswith(".json"):
        logger.info("CAZy dictionary (JSON file) provided")
        # move to script that handles retrieving sequences for proteins in dict (JSON file)
        from_dict.sequences_for_proteins_from_dict(date_today, args)

    else:
        logger.info("CAZy database provided")
        # move to script that handles retrieving sequences for proteins in a SQL database
        from_sql_db.sequences_for_proteins_from_db(date_today, args)

    if args.blastdb is not None:  # build a local BLAST database
        blast_db.build_blast_db(args)
    
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
