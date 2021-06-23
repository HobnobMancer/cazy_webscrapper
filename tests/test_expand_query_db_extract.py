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
"""Tests the script scraper.expand.get_gebank_sequences.from_sql_db.query_sql.py submodule

These test are intened to be run from the root of the repository using:
pytest -v
"""

import pytest

from scraper.expand.query_sql_db import query_extract_seq


def test_expnd_extract_seq_prim_gbk_class_fams(db_session):
    """Test get_prim_gnbk_acc_from_clss_fams()"""

    config_data = {
        "classes": ['PL'],
        "GH": ["GH5", "GH3_1"],
        "PL": None,
    }

    query_extract_seq.get_prim_gnbk_acc_from_clss_fams_with_seq(db_session, config_data)


def test_expnd_extract_seq_all_gbk_class_fams(db_session):
    """Test get_all_gnbk_acc_from_clss_fams()"""

    config_data = {
        "classes": ['PL'],
        "GH": ["GH5", "GH3_1"],
        "PL": None,
    }

    query_extract_seq.get_all_gnbk_acc_from_clss_fams_with_seq(db_session, config_data)


def test_expnd_extract_seq_prim_gbk(db_session):

    query_extract_seq.get_prim_genbank_accessions_with_seq(db_session)


def test_expnd_extract_seq_all_gbk_update(db_session):

    query_extract_seq.get_genbank_accessions_with_seq(db_session)


def test_user_accession_empty_extracted(db_session):
    """Test get_user_accessions_with_seq() when no records are retireved."""

    gbks = ["fake", "test"]
    assert [] == query_extract_seq.get_user_accessions_with_seq(gbks, db_session)


def test_user_accessions_susccess_extracted(db_session):
    """Test get_user_accessions_with_seq()."""

    gbks = ["WP_019387681.1", "CDF79931.1"]
    query_extract_seq.get_user_accessions_with_seq(gbks, db_session)
