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


@pytest.fixture
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
        acc1 = {'acc': Namespace(pdb_accession="accession1")}
        return [acc1['acc']]
    
    def mock_download(*args, **kwargs):
        return

    monkeypatch.setattr(parsers, "build_pdb_structures_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(utilities, "config_logger", mock_config_logger)
    monkeypatch.setattr(get_pdb_structures, "get_database_session", mock_get_sess)
    monkeypatch.setattr(file_io, "make_output_directory", mock_making_output_dir)
    monkeypatch.setattr(get_pdb_structures, "get_pdb_accessions", mock_pdb_acc)
    monkeypatch.setattr(get_pdb_structures, "download_pdb_structures", mock_download)
    
    get_pdb_structures.main()


@pytest.mark.skip(reason="Issue mocking parsing args, test to be fixed")
def test_main_argv(db_path, args_parser, monkeypatch):
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
        parser = Namespace(
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
        )
        return parser

    def mock_config_logger(*args, **kwargs):
        return

    def mock_get_sess(*args, **kwargs):
        return

    def mock_making_output_dir(*args, **kwargs):
        return
    
    def mock_pdb_acc(*args, **kwargs):
        acc1 = {'acc': Namespace(pdb_accession="accession1")}
        return [acc1['acc']]
    
    def mock_download(*args, **kwargs):
        return

    monkeypatch.setattr(parsers, "build_pdb_structures_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(utilities, "config_logger", mock_config_logger)
    monkeypatch.setattr(get_pdb_structures, "get_database_session", mock_get_sess)
    monkeypatch.setattr(file_io, "make_output_directory", mock_making_output_dir)
    monkeypatch.setattr(get_pdb_structures, "get_pdb_accessions", mock_pdb_acc)
    monkeypatch.setattr(get_pdb_structures, "download_pdb_structures", mock_download)
    
    get_pdb_structures.main(["args"])


# test get_database_session()


def test_get_sesh_false_path(args_no_db):
    """Test get_database_session when the path does not exist"""
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        get_pdb_structures.get_database_session(args_no_db["args"])
    assert pytest_wrapped_e.type == SystemExit


def test_get_sesh_error(args_parser, monkeypatch):
    """Test get_database_session when an error arises when opening the session."""
    
    def mock_raise_error(*args, **kawrgs):
        return ValueError

    monkeypatch.setattr(sql_orm, "get_db_session", mock_raise_error)

    get_pdb_structures.get_database_session(args_parser["args"])


def test_get_sesh_success(args_parser, monkeypatch):
    """Test get_database_session when successfully retrieves session."""
    
    def mock_get_sesh(*args, **kawrgs):
        raise

    monkeypatch.setattr(sql_orm, "get_db_session", mock_get_sesh)

    get_pdb_structures.get_database_session(args_parser["args"])


# test get_pdb_accessions()


def test_no_config(args_parser, monkeypatch):
    """Test get_pdb_accessions when no configuration is given."""

    def mock_parse_config(*args, **kwargs):
        return None, None, None, None
    
    def mock_db_query(*args, **kwargs):
        [['item','acc'], ['item1','acc'], ['item','acc']]
    
    monkeypatch.setattr(
        parse_configuration,
        "parse_configuration_for_cazy_database",
        mock_parse_config,
    )
    monkeypatch.setattr(get_pdb_structures, "get_all_pdb_accessions", mock_db_query)

    get_pdb_structures.get_pdb_accessions(args_parser["args"], "session")


def test_get_accessions_config_tax_only(args_parser, monkeypatch):
    """Test get_pdb_accessions when configuration is given, only apply tax filter."""

    def mock_parse_config(*args, **kwargs):
        return {"classes": ["GH"], "families": ["PL28"]}, {"Aspergillus"}, None, None
    
    def mock_db_query(*args, **kwargs):
        genus = Namespace(genus='Aspergillus', species='Fumigatus')
        kingdom = Namespace(kingoms='Bacteria')
        ecs = Namespace(ec_number='1.2.3.4')
        [
            ['NA','acc'],
            ['item1', genus, kingdom, ecs],
        ]
    
    monkeypatch.setattr(
        parse_configuration,
        "parse_configuration_for_cazy_database",
        mock_parse_config,
    )
    monkeypatch.setattr(get_pdb_structures, "get_pdb_acc_from_clss_fams", mock_db_query)

    get_pdb_structures.get_pdb_accessions(args_parser["args"], "session")


def test_get_accessions_config_kngdm_only(args_parser, monkeypatch):
    """Test get_pdb_accessions when configuration is given, only apply kingdom filter."""

    def mock_parse_config(*args, **kwargs):
        return {"classes": ["GH"], "families": ["PL28"]}, None, {'Bacteria'}, None
    
    def mock_db_query(*args, **kwargs):
        genus = Namespace(genus='Aspergillus', species='Fumigatus')
        kingdom = Namespace(kingoms='Bacteria')
        ecs = Namespace(ec_number='1.2.3.4')
        [
            ['NA','acc'],
            ['item1', genus, kingdom, ecs],
        ]
    
    monkeypatch.setattr(
        parse_configuration,
        "parse_configuration_for_cazy_database",
        mock_parse_config,
    )
    monkeypatch.setattr(get_pdb_structures, "get_pdb_acc_from_clss_fams", mock_db_query)

    get_pdb_structures.get_pdb_accessions(args=args_parser["args"], session="session")


def test_get_accessions_config_ec_only(args_parser, monkeypatch):
    """Test get_pdb_accessions when configuration is given, only apply EC number filter."""

    def mock_parse_config(*args, **kwargs):
        return {"classes": ["GH"], "families": ["PL28"]}, None, None, {'1.2.3.4'}
    
    def mock_db_query(*args, **kwargs):
        genus = Namespace(genus='Aspergillus', species='Fumigatus')
        kingdom = Namespace(kingoms='Bacteria')
        ecs = Namespace(ec_number='1.2.3.4')
        [
            ['NA','acc'],
            ['item1', genus, kingdom, ecs],
        ]
    
    monkeypatch.setattr(
        parse_configuration,
        "parse_configuration_for_cazy_database",
        mock_parse_config,
    )
    monkeypatch.setattr(get_pdb_structures, "get_pdb_acc_from_clss_fams", mock_db_query)

    get_pdb_structures.get_pdb_accessions(args=args_parser["args"], session="session")


def test_get_accession_config_tax_kngdm(args_parser, monkeypatch):
    """Test get_pdb_accessions when configuration is given, Taxa and Kingdom filters."""

    def mock_parse_config(*args, **kwargs):
        return {"classes": ["GH"], "families": ["PL28"]}, {'Aspergillus'}, {'Bacteria'}, None
    
    def mock_db_query(*args, **kwargs):
        genus = Namespace(genus='Aspergillus', species='Fumigatus')
        kingdom = Namespace(kingoms='Bacteria')
        ecs = Namespace(ec_number='1.2.3.4')
        [
            ['NA','acc'],
            ['item1', genus, kingdom, ecs],
        ]
    
    monkeypatch.setattr(
        parse_configuration,
        "parse_configuration_for_cazy_database",
        mock_parse_config,
    )
    monkeypatch.setattr(get_pdb_structures, "get_pdb_acc_from_clss_fams", mock_db_query)

    get_pdb_structures.get_pdb_accessions(args=args_parser["args"], session="session")


def test_get_accession_config_tax_ec(args_parser, monkeypatch):
    """Test get_pdb_accessions when configuration is given, Taxa and EC number filters."""

    def mock_parse_config(*args, **kwargs):
        return {"classes": ["GH"], "families": ["PL28"]}, {'Aspergillus'}, None, {'1.2.3.4'}
    
    def mock_db_query(*args, **kwargs):
        genus = Namespace(genus='Aspergillus', species='Fumigatus')
        kingdom = Namespace(kingoms='Bacteria')
        ecs = Namespace(ec_number='1.2.3.4')
        [
            ['NA','acc'],
            ['item1', genus, kingdom, ecs],
        ]
    
    monkeypatch.setattr(
        parse_configuration,
        "parse_configuration_for_cazy_database",
        mock_parse_config,
    )
    monkeypatch.setattr(get_pdb_structures, "get_pdb_acc_from_clss_fams", mock_db_query)

    get_pdb_structures.get_pdb_accessions(args=args_parser["args"], session="session")
