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
"""Tests the script scraper.expand.get_pdb_structures.py

These test are intened to be run from the root of the repository using:
pytest -v
"""


import pytest

from argparse import Namespace, ArgumentParser

from scraper import utilities
from scraper.expand.get_pdb_structures import get_pdb_structures
from scraper.sql import sql_orm
from scraper.utilities import file_io, parse_configuration, parsers


@pytest.fixture
def output_dir(test_dir):
    return test_dir / "test_outputs/test_expand_pdb"


@pytest.fixture
def args_parser(db_path, output_dir):
    args = {"args": Namespace(
        database=db_path,
        pdb="pdb,xml",
        contig=None,
        classes=None,
        force=True,
        families=None,
        kingdoms=None,
        genera=None,
        log=None,
        nodelete=False,
        outdir=output_dir,
        species=None,
        strains=None,
        verbose=None,
    )}
    return args


def args_no_db():
    args = {"args": Namespace(
        database="Fake_path"
    )}
    return args


def test_main_args(args_parser, monkeypatch):
    """Test main() when args is None."""

    def mock_building_parser(*args, **kwargs):
        parser_args = ArgumentParser(
            prog="get_pdb_structures",
            usage=None,
            description="Retrieve structures for CAZymes",
            conflict_handler="error",
            add_help=True,
        )
        return parser_args

    def mock_parser(*args, **kwargs):
        return args_parser["args"]

    def mock_config_logger(*args, **kwargs):
        return

    def mock_get_sess(*args, **kwargs):
        return

    def mock_making_output_dir(*args, **kwargs):
        return
    
    def mock_pdb_acc(*args, **kwargs):
        return [1,2,3,4,5,6]
    
    def mock_download(*args, **kwargs):
        return

    monkeypatch.setattr(parsers, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(utilities, "config_logger", mock_config_logger)
    monkeypatch.setattr(get_pdb_structures, "get_database_session", mock_get_sess)
    monkeypatch.setattr(file_io, "make_output_directory", mock_making_output_dir)
    monkeypatch.setattr(get_pdb_structures, "get_pdb_accessions", mock_pdb_acc)
    monkeypatch.setattr(get_pdb_structures, "download_pdb_structures", mock_download)
    
    get_pdb_structures.main()


def test_main_argv(args_parser, monkeypatch):
    """Test main() when args is not None."""

    def mock_building_parser(*args, **kwargs):
        parser_args = ArgumentParser(
            prog="get_pdb_structures",
            usage=None,
            description="Retrieve structures for CAZymes",
            conflict_handler="error",
            add_help=True,
        )
        return parser_args

    def mock_parser(*args, **kwargs):
        return args_parser["args"]

    def mock_config_logger(*args, **kwargs):
        return

    def mock_get_sess(*args, **kwargs):
        return

    def mock_making_output_dir(*args, **kwargs):
        return
    
    def mock_pdb_acc(*args, **kwargs):
        return [1,2,3,4,5,6]
    
    def mock_download(*args, **kwargs):
        return

    monkeypatch.setattr(parsers, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(utilities, "config_logger", mock_config_logger)
    monkeypatch.setattr(get_pdb_structures, "get_database_session", mock_get_sess)
    monkeypatch.setattr(file_io, "make_output_directory", mock_making_output_dir)
    monkeypatch.setattr(get_pdb_structures, "get_pdb_accessions", mock_pdb_acc)
    monkeypatch.setattr(get_pdb_structures, "download_pdb_structures", mock_download)
    
    get_pdb_structures.main(argv=args_parser["args"])


def test_get_sesh_false_path(args_no_db):
    """Test get_database_session when the path does not exist"""
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        get_pdb_structures.get_database_session(args_no_db["args"])
    assert pytest_wrapped_e.type == SystemExit


def test_get_sesh_error(args_parser, monkeypatch):
    """Test get_database_session when an error arises when opening the session."""
    
    def mock_raise_error(*args, **kawrgs):
        raise ValueError

    monkeypatch.setattr(sql_orm, "get_db_session", mock_raise_error)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        get_pdb_structures.get_database_session(args_parser["args"])
    assert pytest_wrapped_e.type == SystemExit


def test_get_sesh_success(args_parser, monkeypatch):
    """Test get_database_session when successfully retrieves session."""
    
    def mock_get_sesh(*args, **kawrgs):
        raise

    monkeypatch.setattr(sql_orm, "get_db_session", mock_get_sesh)

    get_pdb_structures.get_database_session(args_parser["args"])
