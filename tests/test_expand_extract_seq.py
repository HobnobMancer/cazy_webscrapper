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
"""Tests the script scraper.expand.extract_db_sequences.extract_db_sequences.py

These test are intened to be run from the root of the repository using:
pytest -v
"""


import pytest

from argparse import Namespace, ArgumentParser
from pathlib import Path

from scraper import utilities
from scraper.expand.extract_db_sequences import extract_db_sequences
from scraper.expand.query_sql_db import query_extract_seq
from scraper.expand.get_genbank_sequences import from_sql_db
from scraper.sql.sql_orm import Genbank
from scraper.utilities import file_io, parsers, parse_configuration


# test main()


def test_main_both_none(monkeypatch):
    """Test main FASTA and BLASTDB are None"""

    db_path = Path("cazy_database.db")
    email = "dummy.email@domain"
    fasta_file = Path("fasta_dir/fasta_file.fasta")
    fasta_only = Path("fasta_only_dir/dir1/fasta_only.fasta")

    args_namespace = {"args": Namespace(
        database=db_path,
        fasta=None,
        blastdb=None,
        verbose=False,
        log=None,
        force=False,
        nodelete=False,
    )}

    def mock_building_parser(*args, **kwargs):
        parser_args = ArgumentParser(
            prog="extract_db_sequences",
            usage=None,
            description="Retrieve protein sequences from local database",
            conflict_handler="error",
            add_help=True,
        )
        return parser_args

    def mock_parser(*args, **kwargs):
        return args_namespace["args"]

    def mock_config_logger(*args, **kwargs):
        return

    def mock_making_output_dir(*args, **kwargs):
        return
    
    def mock_from_sql(*args, **kwargs):
        return

    def mock_build_db(*args, **kwargs):
        return

    monkeypatch.setattr(parsers, "build_extract_sequences_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(utilities, "config_logger", mock_config_logger)

    with pytest.raises(SystemExit):
        extract_db_sequences.main()


def test_main_output_files(monkeypatch):
    """Test main when writing proteins to one fasta file"""

    db_path = Path("cazy_database.db")
    email = "dummy.email@domain"
    fasta_file = Path("fasta_dir/fasta_file.fasta")

    args_namespace = {"args": Namespace(
        database=db_path,
        fasta=fasta_file,
        verbose=False,
        log=None,
        blastdb=Path("blastdb_dir/dir1/dir2"),
        force=False,
        nodelete=False,
    )}

    def mock_building_parser(*args, **kwargs):
        parser_args = ArgumentParser(
            prog="extract_db_sequences",
            usage=None,
            description="Retrieve protein sequences from local database",
            conflict_handler="error",
            add_help=True,
        )
        return parser_args

    def mock_parser(*args, **kwargs):
        return args_namespace["args"]

    def mock_config_logger(*args, **kwargs):
        return

    def mock_making_output_dir(*args, **kwargs):
        return

    def mock_build_db(*args, **kwargs):
        return
    
    def mock_configuration(*args, **kwargs):
        return {}, set(), set(), set()
    
    def mock_accession(*args, **kwargs):
        genbank = Genbank(genbank_accession="accession")
        return [[genbank], [genbank]]

    monkeypatch.setattr(parsers, "build_extract_sequences_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(utilities, "config_logger", mock_config_logger)
    monkeypatch.setattr(file_io, "make_output_directory", mock_making_output_dir)
    monkeypatch.setattr(extract_db_sequences, "get_db_session", mock_build_db)
    monkeypatch.setattr(parse_configuration, "parse_configuration_for_cazy_database", mock_configuration)
    monkeypatch.setattr(extract_db_sequences, "get_genbank_records", mock_accession)
    monkeypatch.setattr(file_io, "write_out_extract_seq_to_fasta", mock_making_output_dir)
    monkeypatch.setattr(file_io, "write_extracted_fasta_for_db", mock_making_output_dir)
    monkeypatch.setattr(file_io, "build_blast_db", mock_build_db)

    extract_db_sequences.main()


# test get_genbank_records()


def test_get_acc_config_prim(monkeypatch):
    """Test get_genbank_records(), config_dict provided, primary is True"""
    config_dict = {}

    args_dict = {"args": Namespace(primary=True, accessions="acc,acc", accessions_path="acc")}

    def mock_db_query(*args, **kwargs):
        genbank = Genbank(genbank_accession="accession")
        return [[genbank], [genbank]], [[genbank], [genbank]]
    
    def mock_acc_file(*args, **kwargs):
        return [1,2,3,4]
    
    monkeypatch.setattr(query_extract_seq, "get_prim_gnbk_acc_from_clss_fams_with_seq", mock_db_query)
    monkeypatch.setattr(query_extract_seq, "get_user_accessions_with_seq", mock_db_query)
    monkeypatch.setattr(file_io, "get_accessions_from_file", mock_acc_file)
    monkeypatch.setattr(from_sql_db, "parse_genbank_query", mock_acc_file)
