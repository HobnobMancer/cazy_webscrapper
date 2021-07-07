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
"""Submodule for building cmd-line parsers"""


import argparse

from pathlib import Path
from typing import List, Optional


def build_parser(argv: Optional[List] = None):
    """Return ArgumentParser parser for script."""
    # Create parser object
    parser = argparse.ArgumentParser(
        prog="get_cazy_pages",
        description="Scrapes the CAZy database and writes HTML pages to local HTML files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "output",
        type=Path,
        help="Path to output directory to store HTML files",
    )

    # Add optional arguments to parser

    # Add option to specify path to configuration file
    parser.add_argument(
        "-c",
        "--config",
        type=Path,
        metavar="config file",
        default=None,
        help="Path to configuration file. Default: None, scrapes entire database",
    )

    # Add option to define complete classes to scrape
    parser.add_argument(
        "--classes",
        type=str,
        default=None,
        help="Classes from which all families are to be scraped. Separate classes by ','"
    )

    # Add option to use own CAZy class synoymn dict
    parser.add_argument(
        "--cazy_synonyms",
        type=Path,
        default=None,
        help="Path to JSON file containing CAZy class synoymn names",
    )

    # Add option to scrape all of CAZy
    parser.add_argument(
        "--complete_download",
        dest="complete_download",
        action="store_true",
        default=False,
        help=(
            "Download all of CAZy. Note this is very demanding on the CAZy server\n"
            "Please only what is required if possible"
        ),
    )

    # Add option to force file over writting
    parser.add_argument(
        "-f",
        "--force",
        dest="force",
        action="store_true",
        default=False,
        help="Force writing in output directory that already exists",
    )

    # Add option to specify families to scrape
    parser.add_argument(
        "--families",
        type=str,
        default=None,
        help=(
            "Families to retrieve HTML pages for.\n"
            "Separate families by commas 'GH1,GH2' (case sensitive)"
        ),
    )

    # Add option to restrict the scrape to specific kingdoms
    parser.add_argument(
        "--kingdoms",
        type=str,
        default=None,
        help=(
            "Kingdoms to retrieve pages for each CAZy family. Separate by a single comma.\n"
            "Options= archaea, bacteria, eukaryota, viruses, unclassified (not case sensitive)"
        ),
    )

    # Add log file name option
    # If not given, no log file will be written out
    parser.add_argument(
        "-l",
        "--log",
        type=Path,
        metavar="log file name",
        default=None,
        help="Defines log file name and/or path",
    )

    # Add option to prevent over writing of existing files
    # and cause addition of files to output directory
    parser.add_argument(
        "-n",
        "--nodelete",
        dest="nodelete",
        action="store_true",
        default=False,
        help="enable/disable deletion of files already present in output directory",
    )

    # Add option to enable number of times to retry scraping
    parser.add_argument(
        "-r",
        "--retries",
        type=int,
        default=10,
        help="Number of times to retry scraping a family or class page if error encountered",
    )

    # Add option to enable retrieval of subfamilies
    parser.add_argument(
        "-s",
        "--subfamilies",
        dest="subfamilies",
        action="store_true",
        default=False,
        help="Enable retrieval of pages for CAZy subfamilies from CAZy",
    )

    # Add option to define time out limit for trying to connect to CAZy
    parser.add_argument(
        "-t",
        "--timeout",
        type=int,
        default=45,
        help="Connection timeout limit (seconds)"
    )

    # Add option for more detail (verbose) logging
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        default=False,
        help="Set logger level to 'INFO'",
    )

    if argv is None:
        # parse command-line
        return parser
    else:
        # return namespace
        return parser.parse_args(argv)
