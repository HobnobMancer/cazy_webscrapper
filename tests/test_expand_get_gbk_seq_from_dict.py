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
"""Tests the script scraper.expand.get_gebank_sequences.from_dict.__init__.py

These test are intened to be run from the root of the repository using:
pytest -v
"""

import pytest

from argparse import Namespace, ArgumentParser
from pathlib import Path

from scraper.utilities import parse_configuration
from scraper.expand.get_genbank_sequences import from_dict, ncbi


@pytest.fixture
def cazy_dict():
    cazy_dict = {
        "Protein1": ['GH1', 'CBM1'],
        "Protein2": ['PL1', 'PL2'],
        "Protien3": ['CE5', 'AA7'],
        "Protein4": ['GT5', 'CBM9'],
    }
    return cazy_dict


# test sequences_for_proteins_from_dict()


def test_seq_for_proteins_no_proteins(monkeypatch):
    """Test sequences_for_proteins_from_dict() when no proteins are retrieved"""

    def mock_get_cazy_dict(*args, **kwargs):
        return {}

    def mock_parse_config(*args, **kwargs):
        return {}

    def mock_get_qualifying_proteins(*args, **kwargs):
        return ["NA"]
    
    monkeypatch.setattr(from_dict, "get_cazy_dict", mock_get_cazy_dict)
    monkeypatch.setattr(parse_configuration, "parse_configuration_for_cazy_dict", mock_parse_config)
    monkeypatch.setattr(from_dict, "get_qualifying_proteins", mock_get_qualifying_proteins)

    assert "EXIT" == from_dict.sequences_for_proteins_from_dict("date_today", "args")


def test_seq_for_proteins(monkeypatch):
    """Test sequences_for_proteins_from_dict()"""

    def mock_get_cazy_dict(*args, **kwargs):
        return {}

    def mock_parse_config(*args, **kwargs):
        return {}

    def mock_get_qualifying_proteins(*args, **kwargs):
        return ["Protein1", "Protein2", "Protein3"]
    
    def mock_runtime_error(*args, **kwargs):
        raise RuntimeError

    args_namespace = {"args": Namespace(
        epost=2,
    )}
    
    monkeypatch.setattr(from_dict, "get_cazy_dict", mock_get_cazy_dict)
    monkeypatch.setattr(parse_configuration, "parse_configuration_for_cazy_dict", mock_parse_config)
    monkeypatch.setattr(from_dict, "get_qualifying_proteins", mock_get_qualifying_proteins)
    monkeypatch.setattr(ncbi, "get_sequences_for_dict", mock_runtime_error)
    monkeypatch.setattr(ncbi, "get_sequences_for_dict", mock_runtime_error)

    from_dict.sequences_for_proteins_from_dict("date_today", args_namespace["args"])


# test get_cazy_dict()


def test_get_cazy_dict_expand(test_input_dir):
    """Test get_cazy_dict() when successful"""
    dict_path = test_input_dir / "cazy_dictionary.json"

    args_namespace = {"args": Namespace(database=dict_path)}

    from_dict.get_cazy_dict(args_namespace["args"])


def test_get_cazy_dict_expand_error(test_input_dir):
    """Test get_cazy_dict() when FileNotFound error is raised"""
    dict_path = test_input_dir / "cazy_dictaaaaaaionary.json"

    args_namespace = {"args": Namespace(database=dict_path)}

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        from_dict.get_cazy_dict(args_namespace["args"])
    assert pytest_wrapped_e.type == SystemExit


# test get_qualifying_proteins()


def test_get_qual_proteins(cazy_dict):
    """Test get_qualifying_proteins() when no fams or classes specified"""
    config_dict = {"classes": [], "families": []}
    from_dict.get_qualifying_proteins(cazy_dict, config_dict)


def test_get_qual_proteins_fams(cazy_dict):
    """Test get_qualifying_proteins() when only families are specified"""
    config_dict = {"classes": [], "families": ['GH1', 'PL1']}
    from_dict.get_qualifying_proteins(cazy_dict, config_dict)


def test_get_qual_proteins_classes(cazy_dict):
    """Test get_qualifying_proteins() when only classes are specified"""
    config_dict = {"classes": ['GH', 'CBM'], "families": []}
    from_dict.get_qualifying_proteins(cazy_dict, config_dict)


def test_get_qual_proteins_classes_fams(cazy_dict):
    """Test get_qualifying_proteins() when classes and fams are specified"""
    config_dict = {"classes": ['PL'], "families": ['GH1']}
    from_dict.get_qualifying_proteins(cazy_dict, config_dict)
