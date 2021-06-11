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
"""Tests the script scraper.expand.get_genbank_sequences.py

These test are intened to be run from the root of the repository using:
pytest -v
"""


import os
import pytest

from argparse import Namespace, ArgumentParser

from scraper import utilities
from scraper.expand.get_genbank_sequences import get_genbank_sequences
from scraper.sql import sql_orm
from scraper.utilities import file_io, parse_configuration, parsers


@pytest.fixture
def args_both_none():
    argsdict = {
        "args": Namespace(
            email="dummyemail@domain.sec",
            database=None,
            dict=None,
            fasta=os.getcwd(),
            blast_db="Not None",
            verbose=False,
            log=None,
        )
    }
    return argsdict


@pytest.fixture
def args_both_not_none():
    argsdict = {
        "args": Namespace(
            email="dummyemail@domain.sec",
            database="db",
            dict="dict",
            fasta=os.getcwd(),
            blast_db="Not None",
            verbose=False,
            log=None,
        )
    }
    return argsdict


@pytest.fixture
def args_false_db():
    argsdict = {
        "args": Namespace(
            email="dummyemail@domain.sec",
            database="db",
            dict=None,
            fasta=os.getcwd(),
            blast_db="Not None",
            verbose=False,
            log=None,
        )
    }
    return argsdict



# test main()


def test_gbk_seq_main_both_none(args_both_none, monkeypatch):
    """Test main() when args database and dict are None."""

    def mock_building_parser(*args, **kwargs):
        parser_args = ArgumentParser(
            prog="get_genbank_sequences.py",
            usage=None,
            description="Retrieve protein sequnces from GenBank",
            conflict_handler="error",
            add_help=True,
        )
        return parser_args
    
    def mock_parser(*args, **kwargs):
        return args_both_none['args']

    def mock_config_logger(*args, **kwargs):
        return

    monkeypatch.setattr(parsers, "build_genbank_sequences_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(utilities, "config_logger", mock_config_logger)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        get_genbank_sequences.main()
    assert pytest_wrapped_e.type == SystemExit


def test_gbk_seq_main_both_non_none(args_both_not_none, monkeypatch):
    """Test main() when args database and dict are both not None."""

    def mock_building_parser(*args, **kwargs):
        parser_args = ArgumentParser(
            prog="get_genbank_sequences.py",
            usage=None,
            description="Retrieve protein sequnces from GenBank",
            conflict_handler="error",
            add_help=True,
        )
        return parser_args
    
    def mock_parser(*args, **kwargs):
        return args_both_not_none['args']

    def mock_config_logger(*args, **kwargs):
        return

    monkeypatch.setattr(parsers, "build_genbank_sequences_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(utilities, "config_logger", mock_config_logger)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        get_genbank_sequences.main()
    assert pytest_wrapped_e.type == SystemExit


def test_gbk_seq_main_db_false(args_false_db, monkeypatch):
    """Test main() when args database is not none but the path does not exist."""

    def mock_building_parser(*args, **kwargs):
        parser_args = ArgumentParser(
            prog="get_genbank_sequences.py",
            usage=None,
            description="Retrieve protein sequnces from GenBank",
            conflict_handler="error",
            add_help=True,
        )
        return parser_args
    
    def mock_parser(*args, **kwargs):
        return args_false_db['args']

    def mock_config_logger(*args, **kwargs):
        return

    monkeypatch.setattr(parsers, "build_genbank_sequences_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(utilities, "config_logger", mock_config_logger)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        get_genbank_sequences.main()
    assert pytest_wrapped_e.type == SystemExit
