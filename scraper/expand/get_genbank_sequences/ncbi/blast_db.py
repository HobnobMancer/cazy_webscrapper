#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
# (c) James Hutton Institute 2020-2021
#
# Author:
# Emma E. M. Hobbs

# Contact
# eemh1@st-andrews.ac.uk

# Emma E. M. Hobbs,
# Biomolecular Sciences Building,
# University of St Andrews,
# North Haugh Campus,
# St Andrews,
# KY16 9ST
# Scotland,
# UK

# The MIT License
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""Submodule for creating a local BLAST database"""


import logging

from Bio import SeqIO
from Bio.Blast.Applications import NcbimakeblastdbCommandline


def build_blast_db(args):
    """Build BLAST database of sequences retrieved from GenBank.

    :param args: cmd-line arguments parser

    Return nothing.
    """
    logger = logging.getLogger(__name__)

    blast_database_dir_name = str(args.blastdb).split("/")[-1]

    db_fasta_name = args.blastdb
    db_fasta_name = db_fasta_name / f"blast_db_{blast_database_dir_name}.fasta"

    logger.warning(f"Building local BLAST database containing proteins in:\n{db_fasta_name}")

    # build the command
    if args.blastdb_name is None:
        cmd_makedb = NcbimakeblastdbCommandline(
            cmd='makeblastdb',
            dbtype='prot',
            input_file=db_fasta_name,
        )
    
    else:
        cmd_makedb = NcbimakeblastdbCommandline(
            cmd='makeblastdb',
            dbtype='prot',
            input_file=db_fasta_name,
            title=args.blastdb_name,
        )

    # invoke the command
    stdout, stderr = cmd_makedb()

    # check the command was successfully exectured
    if len(stderr) != 0:
        logger.warning(f"Could not build non-CAZyme db.\nstdout={stdout}\nstderr={stderr}")

    return
