#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
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
"""A module for expanding the the local CAZY database beyond what is provided within CAZy."""


import logging
import json
import sys

from tqdm import tqdm


def get_cazy_dict(args):
    """Retrieve dict of CAZy family annotations of proteins."""
    logger = logging.getLogger(__name__)

    try:
        with open(args.dict, "r") as fh:
            cazy_dict = json.load(fh)

    except FileNotFoundError:
        logger.error(
            "Did not find the local CAZy dict (JSON) file.\n"
            "Check the path is correct.\n"
            "Terminating programme"
        )
        sys.exit(1)

    return cazy_dict


def get_qualifying_proteins(cazy_dict, config_data):
    """Identify proteins to retrieve sequences for, those that meet at least one config criteria.

    :param cazy_dict: dict of proteins catalogued in CAZy
        Keyed by protein GenBank accession, valued by list of CAZy family annotations
    :param config_data: dict of two sets: CAZy classes and families to retrieve seqs for

    Return list of proteins (1 protein = 1 GenBank acession) to retrieve seqs for.
    """
    proteins = set()  # proteins to retrieve sequences for, set prevents duplicates

    if (len(list(config_data["classes"])) == 0) and (len(list(config_data["families"])) != 0):
        # check only in family configuration data
        for protein in tqdm(cazy_dict, desc="Identifying CAZymes matching config data"):
            for fam in cazy_dict[protein]:
                if fam in config_data["families"]:
                    proteins.add(protein)
                    continue
    
    elif (len(list(config_data["classes"])) != 0) and (len(list(config_data["families"])) == 0):
        # check only class configuration data
        for protein in tqdm(cazy_dict, desc="Identifying CAZymes matching config data"):
            for cazy_class in config_data["classes"]:
                for fam in cazy_dict[protein]:
                    if fam.startswith(cazy_class):
                        proteins.add(protein)
                        continue
    
    else:
        # check both configuration data
        for protein in tqdm(cazy_dict, desc="Identifying CAZymes matching config data"):
            for fam in cazy_dict[protein]:
                if fam in config_data["families"]:
                    proteins.add(protein)
                    continue
            for cazy_class in config_data["classes"]:
                for fam in cazy_dict[protein]:
                    if fam.startswith(cazy_class):
                        proteins.add(protein)
                        continue
    
    return list(proteins)
