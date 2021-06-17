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
from scraper.expand.get_genbank_sequences import from_dict


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
    