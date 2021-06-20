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

from scraper.expand.get_genbank_sequences import from_sql_db
from scraper.expand.get_genbank_sequences.from_sql_db  import query_sql_db, ncbi
from scraper.sql import sql_orm
from scraper.sql.sql_orm import Genbank, Cazyme, Taxonomy, Kingdom
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


@pytest.fixture
def genbanks_to_parse():
    genbank = Genbank(genbank_accession="WP123456")
    caz_gbk = "item"
    cazyme = Cazyme(cazyme_name="cazyme_name")
    taxonomy = Taxonomy(genus="Aspergillus", species="Fumagatus strain123")
    kingdom = Kingdom(kingdom="Bacteria")

    query_result = [genbank, caz_gbk, cazyme, taxonomy, kingdom]
    query_results = [query_result, query_result, query_result]
    
    return query_results


# test sequences_for_proteins_from_db()


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


# test get_genbank_accessions() with config dict


def test_get_acc_update_primary_none(db_session, monkeypatch):
    """Test get_genbank_accessions, when a config is given and args.update  and args.primary are True"""

    def mock_query_db(*args, **kwargs):
        return [], []

    monkeypatch.setattr(query_sql_db, "get_prim_gnbk_acc_from_clss_fams", mock_query_db)

    args_mocker = {'args': Namespace(update=True, primary=True, fasta_only=False)}

    from_sql_db.get_genbank_accessions(
        args_mocker["args"],
        db_session,
        "data_today",
        {"classes": [], "GH": ["GH1"]},
        set(),
        set(),
        set()
    )


def test_get_acc_update_primary(db_session, monkeypatch):
    """Test get_genbank_accessions, when a config is given and args.update  and args.primary are True"""

    def mock_query_db(*args, **kwargs):
        return ["acc1"], ["acc2", "acc3"]
    
    def mock_retrieve_accessions(*args, **kwargs):
        return ["acc1", "acc2", "acc3"]
    
    monkeypatch.setattr(query_sql_db, "get_prim_gnbk_acc_from_clss_fams", mock_query_db)
    monkeypatch.setattr(from_sql_db, "parse_genbank_query", mock_retrieve_accessions)
    monkeypatch.setattr(from_sql_db, "check_if_to_update", mock_retrieve_accessions)

    args_mocker = {'args': Namespace(update=True, primary=True, fasta_only=False)}

    from_sql_db.get_genbank_accessions(
        args_mocker["args"],
        db_session,
        "data_today",
        {"classes": [], "GH": ["GH1"]},
        set(),
        set(),
        set()
    )


def test_get_acc_update(db_session, monkeypatch):
    """Test get_genbank_accessions, when a config is given and args.update is True"""

    def mock_query_db(*args, **kwargs):
        return ["acc1"], ["acc2", "acc3"]
    
    def mock_retrieve_accessions(*args, **kwargs):
        return ["acc1", "acc2", "acc3"]
    
    monkeypatch.setattr(query_sql_db, "get_all_gnbk_acc_from_clss_fams", mock_query_db)
    monkeypatch.setattr(from_sql_db, "parse_genbank_query", mock_retrieve_accessions)
    monkeypatch.setattr(from_sql_db, "check_if_to_update", mock_retrieve_accessions)

    args_mocker = {'args': Namespace(update=True, primary=False, fasta_only=False)}

    from_sql_db.get_genbank_accessions(
        args_mocker["args"],
        db_session,
        "data_today",
        {"classes": [], "GH": ["GH1"]},
        set(),
        set(),
        set()
    )


def test_get_acc_primary(db_session, monkeypatch):
    """Test get_genbank_accessions, when a config is given and args.primary is True"""

    genbank_1 = Genbank(genbank_accession="GBK123456")
    genbank_2 = Genbank(genbank_accession="GBK123456")

    def mock_query_db(*args, **kwargs):
        return [genbank_1], [genbank_2]
    
    def mock_retrieve_accessions(*args, **kwargs):
        return [genbank_1, genbank_2]
    
    monkeypatch.setattr(query_sql_db, "get_prim_gnbk_acc_from_clss_fams_no_seq", mock_query_db)
    monkeypatch.setattr(from_sql_db, "parse_genbank_query", mock_retrieve_accessions)
    monkeypatch.setattr(from_sql_db, "check_if_to_update", mock_retrieve_accessions)

    args_mocker = {'args': Namespace(update=False, primary=True, fasta_only=False)}

    from_sql_db.get_genbank_accessions(
        args_mocker["args"],
        db_session,
        "data_today",
        {"classes": [], "GH": ["GH1"]},
        set(),
        set(),
        set()
    )


def test_get_acc_all(db_session, monkeypatch):
    """Test get_genbank_accessions, when a config is given and args.primary and args.update are False"""

    genbank_1 = Genbank(genbank_accession="GBK123456")
    genbank_2 = Genbank(genbank_accession="GBK123456")

    def mock_query_db(*args, **kwargs):
        return [genbank_1], [genbank_2]
    
    def mock_retrieve_accessions(*args, **kwargs):
        return [genbank_1, genbank_2]
    
    monkeypatch.setattr(query_sql_db, "get_all_gnbk_acc_from_clss_fams_no_seq", mock_query_db)
    monkeypatch.setattr(from_sql_db, "parse_genbank_query", mock_retrieve_accessions)
    monkeypatch.setattr(from_sql_db, "check_if_to_update", mock_retrieve_accessions)

    args_mocker = {'args': Namespace(update=False, primary=False, fasta_only=False)}

    from_sql_db.get_genbank_accessions(
        args_mocker["args"],
        db_session,
        "data_today",
        {"classes": [], "GH": ["GH1"]},
        set(),
        set(),
        set()
    )


# test get_genbank_accessions() with no config dict


def test_get_acc_update_primary_none_no_config(db_session, monkeypatch):
    """Test get_genbank_accessions, when a config is None is given and args.update  and args.primary are True"""

    def mock_query_db(*args, **kwargs):
        return []

    monkeypatch.setattr(query_sql_db, "get_prim_genbank_acc_for_update", mock_query_db)

    args_mocker = {'args': Namespace(update=True, primary=True, fasta_only=False)}

    from_sql_db.get_genbank_accessions(
        args_mocker["args"],
        db_session,
        "data_today",
        None,
        set(),
        set(),
        set()
    )


def test_get_acc_update_primary_no_config(db_session, monkeypatch):
    """Test get_genbank_accessions, when a config is None and args.update  and args.primary are True"""

    def mock_query_db(*args, **kwargs):
        return ["acc1", "acc2", "acc3"]
    
    def mock_retrieve_accessions(*args, **kwargs):
        return ["acc1", "acc2", "acc3"]
    
    monkeypatch.setattr(query_sql_db, "get_prim_genbank_acc_for_update", mock_query_db)
    monkeypatch.setattr(from_sql_db, "parse_genbank_query", mock_retrieve_accessions)
    monkeypatch.setattr(from_sql_db, "check_if_to_update", mock_retrieve_accessions)

    args_mocker = {'args': Namespace(update=True, primary=True, fasta_only=False)}

    from_sql_db.get_genbank_accessions(
        args_mocker["args"],
        db_session,
        "data_today",
        None,
        set(),
        set(),
        set()
    )


def test_get_acc_update_no_config(db_session, monkeypatch):
    """Test get_genbank_accessions, when a config is None and args.update is True"""

    def mock_query_db(*args, **kwargs):
        return ["acc1", "acc2", "acc3"]
    
    def mock_retrieve_accessions(*args, **kwargs):
        return ["acc1", "acc2", "acc3"]
    
    monkeypatch.setattr(query_sql_db, "get_all_genbank_acc_for_update", mock_query_db)
    monkeypatch.setattr(from_sql_db, "parse_genbank_query", mock_retrieve_accessions)
    monkeypatch.setattr(from_sql_db, "check_if_to_update", mock_retrieve_accessions)

    args_mocker = {'args': Namespace(update=True, primary=False, fasta_only=False)}

    from_sql_db.get_genbank_accessions(
        args_mocker["args"],
        db_session,
        "data_today",
        None,
        set(),
        set(),
        set()
    )


def test_get_acc_primary_no_config(db_session, monkeypatch):
    """Test get_genbank_accessions, when a config is given and args.primary is True"""

    genbank_1 = Genbank(genbank_accession="GBK123456")
    genbank_2 = Genbank(genbank_accession="GBK123456")

    def mock_query_db(*args, **kwargs):
        return [genbank_1, genbank_2]
    
    def mock_retrieve_accessions(*args, **kwargs):
        return [genbank_1, genbank_2]
    
    monkeypatch.setattr(query_sql_db, "get_prim_genbank_accessions_with_no_seq", mock_query_db)
    monkeypatch.setattr(from_sql_db, "parse_genbank_query", mock_retrieve_accessions)
    monkeypatch.setattr(from_sql_db, "check_if_to_update", mock_retrieve_accessions)

    args_mocker = {'args': Namespace(update=False, primary=True, fasta_only=False)}

    from_sql_db.get_genbank_accessions(
        args_mocker["args"],
        db_session,
        "data_today",
        None,
        set(),
        set(),
        set()
    )


def test_get_acc_all_no_config(db_session, monkeypatch):
    """Test get_genbank_accessions, when a config is given and args.primary and args.update are False"""

    genbank_1 = Genbank(genbank_accession="GBK123456")
    genbank_2 = Genbank(genbank_accession="GBK123456")

    def mock_query_db(*args, **kwargs):
        return [genbank_1, genbank_2]
    
    def mock_retrieve_accessions(*args, **kwargs):
        return [genbank_1, genbank_2]
    
    monkeypatch.setattr(query_sql_db, "get_genbank_accessions_with_no_seq", mock_query_db)
    monkeypatch.setattr(from_sql_db, "parse_genbank_query", mock_retrieve_accessions)
    monkeypatch.setattr(from_sql_db, "check_if_to_update", mock_retrieve_accessions)

    args_mocker = {'args': Namespace(update=False, primary=False, fasta_only=False)}

    from_sql_db.get_genbank_accessions(
        args_mocker["args"],
        db_session,
        "data_today",
        None,
        set(),
        set(),
        set()
    )


# test parse_genbank_query


def test_parse_gbk_no_filters(genbanks_to_parse):
    """Test parse_genbank() when no filters are applied"""
    taxonomy_filters = None
    kingdoms = None
    ec_filters = None
    session = None

    from_sql_db.parse_genbank_query(
        genbanks_to_parse,
        taxonomy_filters,
        kingdoms,
        ec_filters,
        session,
    )


def test_parse_gbk_only_tax_filter(genbanks_to_parse):
    """Test parse_genbank() when only a taxonomy filter is applied are applied"""
    taxonomy_filters = {"Aspergillus", "Fumagatus"}
    kingdoms = None
    ec_filters = None
    session = None

    from_sql_db.parse_genbank_query(
        genbanks_to_parse,
        taxonomy_filters,
        kingdoms,
        ec_filters,
        session,
    )



def test_parse_gbk_only_kingdom_filter(genbanks_to_parse):
    """Test parse_genbank() when only a kingdom filter is applied are applied"""
    taxonomy_filters = None
    kingdoms = {"Bacteria"}
    ec_filters = None
    session = None

    from_sql_db.parse_genbank_query(
        genbanks_to_parse,
        taxonomy_filters,
        kingdoms,
        ec_filters,
        session,
    )


def test_parse_gbk_all_filters(genbanks_to_parse, monkeypatch):
    """Test parse_genbank() when all filters is applied are applied"""
    taxonomy_filters = {"Aspergillus"}
    kingdoms = {"Bacteria"}
    ec_filters = {"2.4.1.-"}
    session = None


    def mock_get_ec_numbers(*args, **kwargs):
        return ["2.4.1.-"]
    
    monkeypatch.setattr(query_sql_db, "query_ec_number", mock_get_ec_numbers)

    from_sql_db.parse_genbank_query(
        genbanks_to_parse,
        taxonomy_filters,
        kingdoms,
        ec_filters,
        session,
    )

