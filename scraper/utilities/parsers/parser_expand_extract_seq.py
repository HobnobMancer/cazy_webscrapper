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
import sys

from pathlib import Path
from typing import List, Optional


def build_extract_sequences_parser(argv: Optional[List] = None):
    """Return ArgumentParser parser for the script 'expand.genbank_sequences.py'."""
    # Create parser object
    parser = argparse.ArgumentParser(
        prog="genbank_sequences.py",
        description=(
            "Extract protein sequences from a local CAZyme database write to FASTA and/or\n"
            "Build a BLAST database"
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Add positional/required arguments
    parser.add_argument(
        "database",
        type=Path,
        metavar="local CAZy database",
        help="Path to local CAZy database",
    )

    # Add optional arguments to parser

    # add option for user to specify their own genbank accessions
    parser.add_argument(
        "--accessions",
        type=str,
        default=None,
        help=(
            "List of GenBank protein accessions to retrieve protein sequences for\n"
            "Separate the accessions with a single comma\n"
            "For example 'WP00001.1,PYN15678.2"
        ),
    )

    parser.add_argument(
        "--accessions_path",
        type=Path,
        default=None,
        help=(
            "Path to a plain text file containing a list of GenBank protein accessions "
            "to extract protein sequences for\n"
            "Write one unique Genbank accession per line of the text file\n"
        ),
    )

    # Add option for building a BLAST database of retrieved protein sequences
    parser.add_argument(
        "-b",
        "--blastdb",
        type=Path,
        default=None,
        help=(
            "Enable creating a BLAST db of extracted protein sequences.\n"
            "Pass a path to a FILE (without file extension) for a BLAST database.\n"
            "This path will be used to create the FASTA file for the database and the database."
        ),
    )

    # Add option to specify path to configuration file
    parser.add_argument(
        "-c",
        "--config",
        type=Path,
        metavar="config file",
        default=None,
        help="Path to configuration file. Default: None, scrapes entire database",
    )

    # Add option to use own CAZy class synoymn dict
    parser.add_argument(
        "--cazy_synonyms",
        type=Path,
        default=None,
        help="Path to JSON file containing CAZy class synoymn names",
    )

    # Add option to define classes to retrieve protein sequences for
    parser.add_argument(
        "--classes",
        type=str,
        default=None,
        help="CAZy Classes to extract protein sequences from. Separate classes by ','"
    )

    parser.add_argument(
        "--ec",
        type=str,
        default=None,
        help="Defines EC numbers to restrict the extraction of protein sequences to"
    )

    # Add option to force file over writting
    parser.add_argument(
        "-f",
        "--force",
        dest="force",
        action="store_true",
        default=False,
        help="Force file over writting",
    )

    # Add option to specify families to retrieve protein sequences for
    parser.add_argument(
        "--families",
        type=str,
        default=None,
        help=(
            "CAZy familes and subfamilies to extract protein sequences from.\n"
            "Separate families by commas 'GH1,GH2,GH5_1'"
        ),
    )

    # Add option to enable writing sequences to FASTA file or files, or not at all
    parser.add_argument(
        "--fasta",
        type=Path,
        default=None,
        help=(
            "Enable writing out retrieved sequences to FASTA file(s).\n"
            "Pass a path to a dir to write out one FASTA file for each retrieved seq\n"
            "OR a path to a FASTA file to write all retrieved seqs to"
        ),
    )

    # Add option to restrict the scrape to specific kingdoms
    parser.add_argument(
        "--kingdoms",
        type=str,
        default=None,
        help=(
            "Kingdoms to scrape. Separate by a single comma.\n"
            "Options= archaea, bacteria, eukaryota, viruses, unclassified (not case sensitive)"
        ),
    )

    # Add option to restrict scrape to specific genera
    parser.add_argument(
        "--genera",
        type=str,
        default=None,
        help="Genera to restrict the extraction of protein sequences to"
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
        help="Enable not deleteing content in already existing output directory",
    )

    # enable retrieving protein sequences for only primary GenBank accessions
    parser.add_argument(
        "-p",
        "--primary",
        dest="primary",
        action="store_true",
        default=False,
        help="Enable retrieveing protein sequences for only primary GenBank accessions",
    )

    # Add option to restrict the scrape to specific species. This will scrape CAZymes from
    # all strains belonging to each listed species
    parser.add_argument(
        "--species",
        type=str,
        default=None,
        help="Species (written as Genus Species) to restrict the extraction of protein sequences to"
    )

    # Add option to restrict scraping to specific strains of organisms
    parser.add_argument(
        "--strains",
        type=str,
        default=None,
        help=(
            "Specific strains of organisms to restrict the extraction of protein sequences to"
            "(written as Genus Species Strain)"
        ),
    )

    # Add option for more detail (verbose) logging
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        default=False,
        help="Set logger level to 'INFO' (default level is warning)",
    )

    if argv is None:
        # parse command-line
        return parser
    else:
        # return namespace
        return parser.parse_args(argv)
