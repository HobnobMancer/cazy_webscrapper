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
"""Tests the script scraper.expand.get_gebank_sequences.from_sql_db.__init__.py submodule

These test are intened to be run from the root of the repository using:
pytest -v
"""


import pytest

from argparse import Namespace, ArgumentParser
from pathlib import Path

from scraper.expand.get_genbank_sequences import from_sql_db, ncbi
from scraper.sql import sql_orm
from scraper.utilities import parse_configuration


@pytest.fixture
def args_fasta_only():
    args = {"args": Namespace(
        fasta_only=Path("fasta_only_path"),
        database=Path("db.db"),
        epost=4,
    )}
    return args


@pytest.fixture
def args_not_fasta_only():
    args = {"args": Namespace(
        fasta_only=None,
        database=Path("db.db"),
        epost=4,
    )}
    return args


def test_sequences_for_proteins_from_db_none(db_session, args_fasta_only, monkeypatch):
    """Test sequences_for_proteins_from _db() when None is returned from get_accessions"""

    def mock_get_session(*args, **kwargs):
        return db_session
    
    def mock_parse_config(*args, **kwargs):
        return {}, set(), set(), set()

    def mock_get_accessions(*args, **kwargs):
        return
    
    monkeypatch.setattr(sql_orm, "get_db_session", mock_get_session)
    monkeypatch.setattr(parse_configuration, "parse_configuration_for_cazy_database", mock_parse_config)
    monkeypatch.setattr(from_sql_db, "get_genbank_accessions", mock_get_accessions)

    from_sql_db.sequences_for_proteins_from_db("data_today", args_fasta_only["args"])


def test_sequences_for_proteins_from_db_empty_list(db_session, args_fasta_only, monkeypatch):
    """Test sequences_for_proteins_from _db() when an empty list is returned from get_accessions"""

    def mock_get_session(*args, **kwargs):
        return db_session
    
    def mock_parse_config(*args, **kwargs):
        return {}, set(), set(), set()

    def mock_get_accessions(*args, **kwargs):
        return []
    
    monkeypatch.setattr(sql_orm, "get_db_session", mock_get_session)
    monkeypatch.setattr(parse_configuration, "parse_configuration_for_cazy_database", mock_parse_config)
    monkeypatch.setattr(from_sql_db, "get_genbank_accessions", mock_get_accessions)

    from_sql_db.sequences_for_proteins_from_db("data_today", args_fasta_only["args"])


def test_sequences_for_proteins_from_db_acc_retrieved_fasta_only(db_session, args_fasta_only, monkeypatch):
    """Test sequences_for_proteins_from _db() when an empty list is returned from get_accessions"""

    def mock_get_session(*args, **kwargs):
        return db_session
    
    def mock_parse_config(*args, **kwargs):
        return {}, set(), set(), set()

    def mock_get_accessions(*args, **kwargs):
        return ["GCA_01234567", "WP123456789"]  # does not contain NA
    
    def mock_runtime_error(*args, **kwargs):
        raise RuntimeError
    
    monkeypatch.setattr(sql_orm, "get_db_session", mock_get_session)
    monkeypatch.setattr(parse_configuration, "parse_configuration_for_cazy_database", mock_parse_config)
    monkeypatch.setattr(from_sql_db, "get_genbank_accessions", mock_get_accessions)
    monkeypatch.setattr(ncbi, "get_sequences", mock_runtime_error)

    from_sql_db.sequences_for_proteins_from_db("data_today", args_fasta_only["args"])


def test_sequences_for_proteins_from_db_acc_retrieved(db_session, args_not_fasta_only, monkeypatch):
    """Test sequences_for_proteins_from _db() when an empty list is returned from get_accessions"""

    def mock_get_session(*args, **kwargs):
        return db_session
    
    def mock_parse_config(*args, **kwargs):
        return {}, set(), set(), set()

    def mock_get_accessions(*args, **kwargs):
        return ["GCA_01234567", "WP123456789"]  # does not contain NA
    
    def mock_runtime_error(*args, **kwargs):
        raise RuntimeError
    
    monkeypatch.setattr(sql_orm, "get_db_session", mock_get_session)
    monkeypatch.setattr(parse_configuration, "parse_configuration_for_cazy_database", mock_parse_config)
    monkeypatch.setattr(from_sql_db, "get_genbank_accessions", mock_get_accessions)
    monkeypatch.setattr(ncbi, "get_sequences_add_to_db", mock_runtime_error)

    from_sql_db.sequences_for_proteins_from_db("data_today", args_not_fasta_only["args"])
