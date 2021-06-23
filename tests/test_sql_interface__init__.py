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
"""Tests the submodule that interacts with the SQL database.

These test are intened to be run from the root of the repository using:
pytest -v
"""


import pytest

from argparse import Namespace
from datetime import datetime, time

from sqlalchemy.orm.exc import ObjectDeletedError

from scraper.sql import sql_interface
from scraper.sql.sql_interface import add_cazyme_data


def test_custom_error():
    """Test custom error."""
    with pytest.raises(sql_interface.SqlInterfaceException) as pytest_wrapped_err:
        raise sql_interface.SqlInterfaceException("message")
    assert pytest_wrapped_err.type == sql_interface.SqlInterfaceException


# tests for adding log to the local database


def test_add_db_log_no_config(db_session):
    config_dict = None
    taxonomy_filters = {"genera": None, "species": None, "strains": None}
    kingdoms = ["Archaea", "Bacteria", "Eukaryota", "Viruses", "Unclassified"]
    ec_filters = ['EC1.2.3.4', 'EC3.4.5.6']
    args = {
        "args": Namespace(
            classes=None,
            families=None,
            genera=None,
            species=None,
            strains=None,
            kingdoms=None,
            ec=None,
            streamline=None,
        )
    }

    sql_interface.log_scrape_in_db(
        "CAZy scrape",
        "YYYY-MM-DD--HH-MM-SS",
        config_dict,
        taxonomy_filters,
        kingdoms,
        ec_filters,
        db_session,
        args["args"],
    )


def test_add_db_log_with_config(db_session):
    config_dict = {}
    taxonomy_filters = {
        "genera": ["Caldivirga", "Cuniculiplasma"],
        "species": ["Pyrococcus furiosus"],
        "strains": ["Saccharolobus solfataricus POZ149", "Saccharolobus solfataricus SULB"]
    }
    ec_filters = ['EC1.2.3.4', 'EC3.4.5.6']
    kingdoms = None
    args = {
        "args": Namespace(
            classes="GH,PL",
            families="AA1,AA2",
            genera="Trichoderma",
            species="Aspergillus Niger",
            strains="Acidianus ambivalens LEI 10",
            kingdoms="Archaea,Bacteria",
            ec="EC1.2.3.4,EC5.6.4.7",
            streamline="genbank",
        )
    }

    sql_interface.log_scrape_in_db(
        "CAZy scrape",
        "YYYY-MM-DD--HH-MM-SS",
        config_dict,
        taxonomy_filters,
        kingdoms,
        ec_filters,
        db_session,
        args["args"],
    )


def test_add_db_log_with_full_config(db_session):
    config_dict = {
        'classes': ['GH','PL','CE'],
        'GH': ['GH1', 'GH2'], 
        'PL': ['PL1_1'],
    }
    taxonomy_filters = {
        "genera": ["Caldivirga", "Cuniculiplasma"],
        "species": ["Pyrococcus furiosus"],
        "strains": ["Saccharolobus solfataricus POZ149", "Saccharolobus solfataricus SULB"]
    }
    ec_filters = ['EC1.2.3.4', 'EC3.4.5.6']
    kingdoms = ["Archaea"]
    args = {
        "args": Namespace(
            classes="GH,PL",
            families="AA1,AA2",
            genera="Trichoderma",
            species="Aspergillus Niger",
            strains="Acidianus ambivalens LEI 10",
            kingdoms="Archaea,Bacteria",
            ec="EC1.2.3.4,EC5.6.4.7",
            streamline="genbank",
        )
    }

    sql_interface.log_scrape_in_db(
        "CAZy scrape",
        "YYYY-MM-DD--HH-MM-SS",
        config_dict,
        taxonomy_filters,
        kingdoms,
        ec_filters,
        db_session,
        args["args"],
    )


def test_add_db_log_no_ec(db_session):
    config_dict = None
    taxonomy_filters = {"genera": None, "species": None, "strains": None}
    kingdoms = ["Archaea", "Bacteria", "Eukaryota", "Viruses", "Unclassified"]
    ec_filters = None
    args = {
        "args": Namespace(
            classes=None,
            families=None,
            genera=None,
            species=None,
            strains=None,
            kingdoms=None,
            ec=None,
            streamline=None,
        )
    }

    sql_interface.log_scrape_in_db(
        "CAZy scrape",
        "YYYY-MM-DD--HH-MM-SS",
        config_dict,
        taxonomy_filters,
        kingdoms,
        ec_filters,
        db_session,
        args["args"],
    )


def test_add_db_log_with_gbk_seq(db_session):
    config_dict = {
        'classes': ['GH','PL','CE'],
        'GH': ['GH1', 'GH2'], 
        'PL': ['PL1_1'],
    }
    taxonomy_filters = {
        "genera": ["Caldivirga", "Cuniculiplasma"],
        "species": ["Pyrococcus furiosus"],
        "strains": ["Saccharolobus solfataricus POZ149", "Saccharolobus solfataricus SULB"]
    }
    ec_filters = ['EC1.2.3.4', 'EC3.4.5.6']
    kingdoms = ["Archaea"]
    args = {
        "args": Namespace(
            classes="GH,PL",
            families="AA1,AA2",
            genera="Trichoderma",
            species="Aspergillus Niger",
            strains="Acidianus ambivalens LEI 10",
            kingdoms="Archaea,Bacteria",
            ec=None,
            update="overwrite",
        )
    }

    sql_interface.log_scrape_in_db(
        "GenBank sequences",
        "YYYY-MM-DD--HH-MM-SS",
        config_dict,
        taxonomy_filters,
        kingdoms,
        ec_filters,
        db_session,
        args["args"],
    )


# tests for add_proteins_to_db


def test_adding_a_new_protein(db_session):
    """Test adding a new protein to the local database."""
    time_stamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    args = {'args': Namespace(streamline=None)}

    try:
        sql_interface.add_protein_to_db(
            "test_cazyme",
            "test_fam",
            "test_genus test_species",
            "kingdom",
            time_stamp,
            db_session,
            args['args'],
            ["EC1.2.3.4"],
            ["Genbank1", "Genbank2"],
            ["Uni1P"],
            ["UniNP1"],
            ["PDB1", "PDB2"],
        )
    except ObjectDeletedError as e:
        pass


def test_add_data_to_an_existing_record_in_db(db_session):
    """Test adding data to an existing CAZyme in the local database."""
    args = {'args': Namespace(streamline=None)}
    sql_interface.add_protein_to_db(
        "JCM19301_1832",
        "PL28",
        "pallidilutea JCM 19301",
        "Bacteria",
        "GAL67220.1",
        db_session,
        args['args'],
        ec_numbers=["EC4.2.2.-"],
        gbk_nonprimary=["Genbank1", "Genbank2"],
        uni_nonprimary=["Uni1", "Uni2"],
    )


def test_genbank_no_cazymes(db_session, monkeypatch):
    """test adding a protein to a database, GenBank accession is found with no linked CAZymes."""

    def mock_adding_a_new_protein(*args, **kwargs):
        return

    monkeypatch.setattr(sql_interface, "add_new_protein_to_db", mock_adding_a_new_protein)

    existing_genbank_with_no_cazyme = "test_genbank_no_cazyme"
    args = {'args': Namespace(streamline=None)}
    sql_interface.add_protein_to_db(
        "test_cazyme_name",
        "cazy_family",
        "source_genus organism",
        "kingdom",
        existing_genbank_with_no_cazyme,
        db_session,
        args['args'],
        ec_numbers=["EC4.2.2.-"],
    )


def test_genbank_no_cazymes_raise_error(db_session, monkeypatch):
    """test adding a protein to a database, GenBank accession is found with no linked CAZymes."""

    def mock_adding_a_new_protein(*args, **kwargs):
        return "error_message"

    monkeypatch.setattr(sql_interface, "add_new_protein_to_db", mock_adding_a_new_protein)

    existing_genbank_with_no_cazyme = "test_genbank_no_cazyme"
    args = {'args': Namespace(streamline=None)}
    sql_interface.add_protein_to_db(
        "test_cazyme_name",
        "cazy_family",
        "source_genus organism",
        "kingdom",
        existing_genbank_with_no_cazyme,
        db_session,
        args['args'],
        ec_numbers=["EC4.2.2.-"],
        gbk_nonprimary=["Genbank1", "Genbank2"],
        uni_nonprimary=["Uni1", "Uni2"],
    )


def test_one_genbank_multiple_cazymes(db_session, monkeypatch):
    """test adding protein to db when genbank is found linked to muliple CAZymes."""

    def mock_add_protein_to_db(*args, **kwargs):
        return

    monkeypatch.setattr(sql_interface, "add_data_to_protein_record", mock_add_protein_to_db)

    genbank_with_multiple_cazymes = "one_genbank_multi_cazymes"
    args = {'args': Namespace(streamline=None)}
    sql_interface.add_protein_to_db(
        "test_cazyme_name",
        "cazy_family",
        "source_genus organism",
        "kingdom",
        genbank_with_multiple_cazymes,
        db_session,
        args['args'],
    )


def test_multiple_genbanks_multiple_cazymes(db_session, monkeypatch):
    """test adding protein to db when finding multiple identical CAZymes and GenBank accesisons."""

    def mock_add_protein_to_db(*args, **kwargs):
        return

    monkeypatch.setattr(sql_interface, "add_data_to_protein_record", mock_add_protein_to_db)

    identical_genbank_accession = "identical_accession"

    args = {'args': Namespace(streamline=None)}

    sql_interface.add_protein_to_db(
        "test_cazyme_name",
        "cazy_family",
        "source_genus organism",
        "kingdom",
        identical_genbank_accession,
        db_session,
        args['args'],
    )


def test_multiple_genbanks_multiple_cazymes_streamline(db_session, monkeypatch):
    """test adding protein to db when finding multiple identical CAZymes and GenBank accesisons.
    Streamline scraping is enabled."""

    def mock_add_protein_to_db(*args, **kwargs):
        return

    monkeypatch.setattr(sql_interface, "parse_unique_genbank_conflict", mock_add_protein_to_db)
    monkeypatch.setattr(sql_interface, "add_cazy_family", mock_add_protein_to_db)
    monkeypatch.setattr(add_cazyme_data, "add_cazy_family", mock_add_protein_to_db)

    identical_genbank_accession = "identical_accession"

    args = {'args': Namespace(streamline='uniprot,ec')}

    sql_interface.add_protein_to_db(
        "test_cazyme_name",
        "cazy_family",
        "source_genus organism",
        "kingdom",
        identical_genbank_accession,
        db_session,
        args['args'],
    )


def test_multiple_genbanks_one_cazyme(db_session, monkeypatch):
    """test adding protien to db when identical GenBank accessions with one CAZyme link."""

    args = {'args': Namespace(streamline=None)}

    def mock_add_protein_to_db(*args, **kwargs):
        return

    monkeypatch.setattr(sql_interface, "add_data_to_protein_record", mock_add_protein_to_db)

    multiple_accession = "multiple_accession"

    sql_interface.add_protein_to_db(
        "test_cazyme_name",
        "cazy_family",
        "source_genus organism",
        "kingdom",
        multiple_accession,
        db_session,
        args['args'],
    )


def test_multiple_genbanks_no_cazymes(db_session, monkeypatch):
    """test adding protien to db when identical GenBank accessions are linked to no CAZymes."""

    args = {'args': Namespace(streamline=None)}

    identical_genbanks_no_cazymes = "multi_accession_no_caz"

    def mock_add_protein_to_db(*args, **kwargs):
        return

    monkeypatch.setattr(sql_interface, "add_new_protein_to_db", mock_add_protein_to_db)

    sql_interface.add_protein_to_db(
        "test_cazyme_name",
        "cazy_family",
        "source_genus organism",
        "kingdom",
        identical_genbanks_no_cazymes,
        db_session,
        args['args'],
    )


def test_add_protein_error_message(db_session, monkeypatch):
    """Test add_protein_to_db() when an error message is returned."""

    def mock_adding_a_new_protein(*args, **kwargs):
        return "error message"

    monkeypatch.setattr(sql_interface, "add_new_protein_to_db", mock_adding_a_new_protein)
    monkeypatch.setattr(sql_interface, "parse_unique_genbank_conflict", mock_adding_a_new_protein)

    existing_genbank_with_no_cazyme = "test_genbank_no_cazyme"
    args = {'args': Namespace(streamline=None)}
    sql_interface.add_protein_to_db(
        "test_cazyme_name",
        "cazy_family",
        "source_genus organism",
        "kingdom",
        existing_genbank_with_no_cazyme,
        db_session,
        args['args'],
        ec_numbers=["EC4.2.2.-"],
    )


# test streamline_addition()


def test_streamline_addition(db_session, monkeypatch):
    """Test streamline addition"""
    time_stamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    cazyme_name = 'ulvan lyase (BN863_22190)'
    primary_genbank = 'WP_038530530.1'
    gbk_nonprimary = ['CDF79931.1', f'nonprimary_{time_stamp}']
    ec_numbers = [f'EC{time_stamp}']
    uni_primary = [f'primaryUNI{time_stamp}']
    uni_nonprimary = [f'nonprimaryUNI{time_stamp}']
    pdb_accessions = [f'PDB{time_stamp}']
    family = f"FAM{time_stamp}"
    source_organism = "Formosa agariphila KMM 3901"
    kingdom = "Bacteria"

    args = {"args": Namespace(streamline="")}

    def mock_none(*args, **kwargs):
        return
    
    monkeypatch.setattr(sql_interface, "add_nonprimary_gbk_accessions", mock_none)
    monkeypatch.setattr(sql_interface, "add_cazy_family", mock_none)
    monkeypatch.setattr(sql_interface, "add_ec_numbers", mock_none)
    monkeypatch.setattr(sql_interface, "add_uniprot_accessions", mock_none)
    monkeypatch.setattr(sql_interface, "add_pdb_accessions", mock_none)

    sql_interface.streamline_addition(
        cazyme_name,
        family,
        source_organism,
        kingdom,
        primary_genbank,
        db_session,
        args['args'],
        ec_numbers,
        gbk_nonprimary,
        uni_primary,
        uni_nonprimary,
        pdb_accessions,
    )


def test_streamline_addition_no_gbk(db_session, monkeypatch):
    """Test streamline addition"""
    time_stamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    cazyme_name = f'cazyme_{time_stamp}'
    primary_genbank = 'previous_primary'
    gbk_nonprimary = []
    ec_numbers = [f'EC{time_stamp}']
    uni_primary = [f'primaryUNI{time_stamp}']
    uni_nonprimary = [f'nonprimaryUNI{time_stamp}']
    pdb_accessions = [f'PDB{time_stamp}']
    family = f"FAM{time_stamp}"
    source_organism = "Formosa agariphila KMM 3901"
    kingdom = "Bacteria"

    args = {"args": Namespace(streamline="")}

    def mock_none(*args, **kwargs):
        return
    
    def mock_error_message(*args, **kwargs):
        return "error message"
    
    monkeypatch.setattr(sql_interface, "parse_unique_genbank_conflict", mock_error_message)
    monkeypatch.setattr(sql_interface, "add_nonprimary_gbk_accessions", mock_none)
    monkeypatch.setattr(sql_interface, "add_cazy_family", mock_none)
    monkeypatch.setattr(sql_interface, "add_ec_numbers", mock_none)
    monkeypatch.setattr(sql_interface, "add_uniprot_accessions", mock_none)
    monkeypatch.setattr(sql_interface, "add_pdb_accessions", mock_none)

    sql_interface.streamline_addition(
        cazyme_name,
        family,
        source_organism,
        kingdom,
        primary_genbank,
        db_session,
        args['args'],
        ec_numbers,
        gbk_nonprimary,
        uni_primary,
        uni_nonprimary,
        pdb_accessions,
    )


# test parse_unique_genbank_conflict()


def test_parse_unique_genbank_conflict(db_session, monkeypatch):
    """Test parse_unique_genbank_conflict() when no primary genbank accessions returned."""
    time_stamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    cazyme_name = "cazyme_test_name"
    family = "sub_famfamily"
    source_organism = "test source organism strain"
    kingdom = "Bacteria"
    primary_genbank = f"unique_primary_genbank_{time_stamp}"

    def mock_none(*args, **kwargs):
        return
    
    monkeypatch.setattr(sql_interface, "add_new_protein_to_db", mock_none)

    sql_interface.parse_unique_genbank_conflict(
        cazyme_name,
        family,
        source_organism,
        kingdom,
        primary_genbank,
        db_session,
    )


# test add deleted cazy family()


def test_add_existing_deleted_fam(db_session):
    """Test add_deteleted_cazy_family() when the fam is already present."""
    fam = "deletedFAM"

    sql_interface.add_deleted_cazy_family(fam, db_session)


def test_add_new_deleted_subfam(db_session):
    """Test add_deteleted_cazy_family() when the fam is already present."""

    time_stamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    fam = f"GH4_{time_stamp}"

    sql_interface.add_deleted_cazy_family(fam, db_session)
