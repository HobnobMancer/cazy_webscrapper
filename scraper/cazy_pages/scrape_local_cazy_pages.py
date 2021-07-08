#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
# (c) James Hutton Institute 2020-2021
#
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
"""
Web scraper to scrape CAZy website and retrieve all protein data.

:cmd_args --config: path to configruration file
:cmd_args --classes: specify CAZy classes to scrape
:cmd_args --database: provide path to a local SQL database to add additional data to
:cmd_args --families: specify CAZy families to retrieve CAZymes from
:cmd_args --force: force overwriting content in exisiting output directory
:cmd_args --genera: specify Genera to retrieve CAZymes from
:cmd_args --kingdoms: specify taxonomy Kingdoms to scrape proteins from
:cmd_args --log: path to log file, enables writing out log messages to a log file
:cmd_args --nodelete: if true does not delete content in pre-existing output directory
:cmd_args --output: path to output directory
:cmd_args --retries: specify the number of times to try scraping a page if connection fails
:cmd_args --subfamilies: enable retrieval of subfamilies from CAZy
:cmd_args --species: specify species to retrieve CAZymes from
:cmd_args --strains: specify specific strains of species to retrieve CAZymes from
:cmd_args --timeout: specify the maximum time (in seconds) before determining connection timed out
:cmd_args --verbose: change logger level from warning to info, verbose logging
"""


import logging
import sys

import pandas as pd

from datetime import datetime
from typing import List, Optional

from scraper.crawler.cazy_html_pages import parse_local_pages
from scraper.sql import sql_orm, sql_interface
from scraper.utilities import (
    config_logger,
    file_io,
    parse_configuration,
    termcolour,
)
from scraper.utilities.parsers import parser_scrape_local_pages


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Set up parser, logger and coordinate overall scrapping of CAZy."""
    # Program preparation
    time_stamp = datetime.now().strftime("%Y-%m-%d--%H-%M-%S")  # used in naming files
    start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # used in terminating message
    start_time = pd.to_datetime(start_time)

    if argv is None:
        parser = parser_scrape_local_pages.build_parser()
        args = parser.parse_args()
    else:
        parser = parser_scrape_local_pages.build_parser(argv)
        args = parser.parse_args()

    if logger is None:
        logger = logging.getLogger(__package__)
        config_logger(args)

    # check output setup
    logger.info("Preparing output directory")
    logger.info(
        f"Output dir:{args.output}\n"
        f"Force writing to exiting output dir: {args.force}\n"
        f"Delete content of exiting output dir: {args.nodelete}\n"
    )
    file_io.make_output_directory(args.output, args.force, args.nodelete)

    cazy_home = "http://www.cazy.org"

    # retrieve configuration data
    logger.info("Parsing configuration")
    (
        excluded_classes,
        config_dict,
        cazy_dict,
        taxonomy_filters_dict,
        kingdoms,
        ec_filters,
    ) = parse_configuration.parse_configuration(args)

    # convert taxonomy_filters to a set for quicker identification of species to scrape
    taxonomy_filters = get_filter_set(taxonomy_filters_dict)
    
    # ADD IN GETTING SESSION !!!!!!!!!!!!!!!!!

    # log scraping of CAZy in local db
    logger.info("Add log of scrape to the local CAZyme database")
    sql_interface.log_scrape_in_db(
        "CAZy scrape",
        time_stamp,
        config_dict,
        taxonomy_filters_dict,
        kingdoms,
        ec_filters,
        session,
        args,
    )

    logger.info("Scraping data from local CAZy HTML pages")
    parse_local_pages.parse_local_pages(
        args,
        cazy_home,
        start_time,
        time_stamp,
        session,
        taxonomy_filters,
        ec_filters,
    )

    end_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    end_time = pd.to_datetime(end_time)
    total_time = end_time - start_time

    logger.info(
        "Finished scraping local CAZy pages\n"
        f"Scrape initated at {start_time}\n"
        f"Scrape finished at {end_time}\n"
        f"Total run time: {total_time}"
    )
    
    end_message = [
        termcolour(
            "=====================cazy_webscraper=====================",
            "green",
        ),
        termcolour(
            "=====================scrape local CAZy page library=====================",
            "yellow",
        ),
        termcolour(
            f"Scrape initated at {start_time}",
            "cyan",
        ),
        termcolour(
            f"Scrape finished at {end_time}",
            "cyan",
        ),
        termcolour(
            f"Total run time: {total_time}",
            "cyan",
        ),
        termcolour(
            "For citation information use: cazy_webscraper -C",
            "white",
        ),
    ]
    sys.stderr.write("\n".join(end_message) + "\n")


def get_filter_set(taxonomy_filters_dict):
    """Create a set of all taxonomy filters from a dictionary.

    :param taxonomy_filers: dict of taxonomy filters

    Return a set.
    """
    taxonomy_filters = []

    for key in taxonomy_filters_dict:
        try:
            if len(taxonomy_filters_dict[key]) != 0:
                taxonomy_filters += taxonomy_filters_dict[key]
        except TypeError:
            pass

    if len(taxonomy_filters) == 0:
        taxonomy_filters = None

    else:
        taxonomy_filters = set(taxonomy_filters)

    return taxonomy_filters


if __name__ == "__main__":
    main()
