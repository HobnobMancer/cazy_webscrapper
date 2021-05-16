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
import sys

import pandas as pd

from datetime import datetime
from typing import List, Optional

from scraper.expand.get_genbank_sequences import (
    from_dict,
    from_sql_db,
)
from scraper.utilities import config_logger
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

    if (args.database is None) and (args.dict is None):
        logger.warning(
            "No CAZy data provided.\n"
            "Please provide a path to a local dict (JSON file) or SQL db.\n"
            "Terminating program."
        )
    
    elif (args.database is not None) and (args.dict is not None):
        logger.warning(
            "CAZy dict (JSON file) AND SQL db provided.\n"
            "Please provide ONLY ONE to path to a local dict (JSON file) OR SQL db.\n"
            "Terminating program."
        )
    
    if args.database is not None:
        # check path to the local db is valid
        if os.path.isfile(args.database) is False:
            logger.error(
                "Could not find local CAZy database.\n"
                "Check path is correct.\nTerminating programme."
            )
            sys.exit(1)
        
        # move to script that handles retrieving sequences for proteins in a SQL database
        from_sql_db.sequences_for_proteins_from_db(args)
    
    elif args.dict is not None:
        # check the path to the local CAZy dict is valid
        if os.path.isfile(args.dict) is False:
            logger.error(
                "Could not find local CAZy dict (JSON file).\n"
                "Check path is correct.\nTerminating programme."
            )
            sys.exit(1)
        
        # move to script that handles retrieving sequences for proteins in dict (JSON file)
        from_dict.sequences_for_proteins_from_dict(args)
    
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
