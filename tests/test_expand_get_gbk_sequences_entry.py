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
"""Tests the script scraper.expand.get_gebank_sequences.get_gebank_sequences.py
which is the entry point for retrieving protein sequences from GenBank (NCBI)

These test are intened to be run from the root of the repository using:
pytest -v
"""


import pytest

from argparse import Namespace, ArgumentParser
from pathlib import Path

from scraper import utilities
from scraper.expand.get_genbank_sequences import get_genbank_sequences, from_sql_db, from_dict
from scraper.utilities import file_io, parsers


def test_main_db(monkeypatch):
    """Test main when a path to a CAZy database is provided,
    and writing seqs to single FASTA files"""

    db_path = Path("cazy_database.db")
    email = "dummy.email@domain"
    fasta_file = Path("fasta_dir/fasta_file.fasta")
    fasta_only = Path("fasta_only_dir/dir1/fasta_only.fasta")

    args_namespace = {"args": Namespace(
        email=email,
        database=db_path,
        fasta=fasta_file,
        fasta_only=fasta_only,
        verbose=False,
        log=None,
        blastdb=Path("blastdb_dir/dir1/dir2"),
        force=False,
        nodelete=False,
        update="uptate_only",
    )}

    def mock_building_parser(*args, **kwargs):
        parser_args = ArgumentParser(
            prog="get_genbank_sequences",
            usage=None,
            description="Retrieve protein sequences for CAZymes",
            conflict_handler="error",
            add_help=True,
        )
        return parser_args

    def mock_parser(*args, **kwargs):
        return args_namespace["args"]

    def mock_config_logger(*args, **kwargs):
        return

    def mock_making_output_dir(*args, **kwargs):
        return
    
    def mock_from_sql(*args, **kwargs):
        return

    def mock_build_db(*args, **kwargs):
        return

    monkeypatch.setattr(parsers, "build_pdb_structures_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(utilities, "config_logger", mock_config_logger)
    monkeypatch.setattr(from_sql_db, "sequences_for_proteins_from_db", mock_from_sql)
    monkeypatch.setattr(file_io, "make_output_directory", mock_making_output_dir)
    monkeypatch.setattr(file_io, "build_blast_db", mock_build_db)

    get_genbank_sequences.main()


def test_main_fastas_into_cwd(monkeypatch):
    """Test main when a path to a CAZy database is provided,
    and writing seqs to single FASTA file in the cwd"""

    db_path = Path("cazy_database.db")
    email = "dummy.email@domain"
    fasta_file = Path("fasta_file.fasta")
    fasta_only = Path("fasta_only.fasta")

    args_namespace = {"args": Namespace(
        email=email,
        database=db_path,
        fasta=fasta_file,
        fasta_only=fasta_only,
        verbose=False,
        log=None,
        blastdb=Path("blastdb_dir/dir1/dir2"),
        force=False,
        nodelete=False,
        update="overwrite",
    )}

    def mock_building_parser(*args, **kwargs):
        parser_args = ArgumentParser(
            prog="get_genbank_sequences",
            usage=None,
            description="Retrieve protein sequences for CAZymes",
            conflict_handler="error",
            add_help=True,
        )
        return parser_args

    def mock_parser(*args, **kwargs):
        return args_namespace["args"]

    def mock_config_logger(*args, **kwargs):
        return

    def mock_making_output_dir(*args, **kwargs):
        return
    
    def mock_from_sql(*args, **kwargs):
        return

    def mock_build_db(*args, **kwargs):
        return

    monkeypatch.setattr(parsers, "build_pdb_structures_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(utilities, "config_logger", mock_config_logger)
    monkeypatch.setattr(from_sql_db, "sequences_for_proteins_from_db", mock_from_sql)
    monkeypatch.setattr(file_io, "make_output_directory", mock_making_output_dir)
    monkeypatch.setattr(file_io, "build_blast_db", mock_build_db)

    get_genbank_sequences.main()


def test_main_dict(monkeypatch):
    """Test main when a path to a CAZy dict is provided,
    and writing seqs to single a FASTA per retrieved seq"""

    db_path = Path("cazy_database.json")
    email = "dummy.email@domain"
    fasta_dir = Path("fasta_dir/dir1/dir2")
    fasta_only = Path("fasta_only/dir1/dir2")

    args_namespace = {"args": Namespace(
        email=email,
        database=db_path,
        fasta=fasta_dir,
        fasta_only=fasta_only,
        verbose=False,
        log=None,
        blastdb=Path("blastdb_dir/dir1/dir2"),
        force=False,
        nodelete=False,
        update=None,
    )}

    def mock_building_parser(*args, **kwargs):
        parser_args = ArgumentParser(
            prog="get_genbank_sequences",
            usage=None,
            description="Retrieve protein sequences for CAZymes",
            conflict_handler="error",
            add_help=True,
        )
        return parser_args

    def mock_parser(*args, **kwargs):
        return args_namespace["args"]

    def mock_config_logger(*args, **kwargs):
        return

    def mock_making_output_dir(*args, **kwargs):
        return
    
    def mock_from_dict(*args, **kwargs):
        return "error message"

    def mock_build_db(*args, **kwargs):
        return

    monkeypatch.setattr(parsers, "build_pdb_structures_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(utilities, "config_logger", mock_config_logger)
    monkeypatch.setattr(from_dict, "sequences_for_proteins_from_dict", mock_from_dict)
    monkeypatch.setattr(file_io, "make_output_directory", mock_making_output_dir)
    monkeypatch.setattr(file_io, "build_blast_db", mock_build_db)

    get_genbank_sequences.main()
