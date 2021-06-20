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
"""Tests the script scraper.expand.get_gebank_sequences.ncbi submodule

These test are intened to be run from the root of the repository using:
pytest -v
"""

import pytest

from argparse import Namespace
from pathlib import Path

from scraper.utilities import file_io
from scraper.expand.get_genbank_sequences import ncbi


@pytest.fixture
def fasta_path(test_dir):
    fasta_path = test_dir / "test_inputs"
    fasta_path = fasta_path / "test_inputs_expand"
    fasta_path = fasta_path / "expand_ncbi_fasta.fasta"
    return fasta_path


@pytest.fixture()
def gbk_accessions():
    return ["WP_001307453.1", "PNY18054.1", "accession_not_retrievied"]


@pytest.fixture
def args_namespace():
    args_dict = {"args": Namespace(fasta="path", fasta_only="path", blastdb="path")}
    return args_dict


def test_get_dict_ncbi(gbk_accessions, args_namespace, fasta_path, monkeypatch):
    """Test get_sequences_for_dict()"""

    with open(fasta_path, "r") as fh:
        retrieved_record = fh

        def mock_epost(*args, **kwargs):
            return {'QueryKey': '1', 'WebEnv': 'MCID_60cf5b5ef457b676172c93c4'}

        def mock_entrez(*args, **kwargst):
            return retrieved_record

        def mock_writing_fastas(*args, **kwargs):
            return

        monkeypatch.setattr(ncbi, "bulk_query_ncbi", mock_epost)
        monkeypatch.setattr(ncbi, "entrez_retry", mock_entrez)
        monkeypatch.setattr(file_io, "write_out_fasta", mock_writing_fastas)
        monkeypatch.setattr(file_io, "write_out_fasta_only", mock_writing_fastas)
        monkeypatch.setattr(file_io, "write_fasta_for_db", mock_writing_fastas)

        ncbi.get_sequences_for_dict(gbk_accessions, args_namespace["args"])


