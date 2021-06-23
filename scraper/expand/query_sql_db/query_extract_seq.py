#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
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
"""Script contains queries to a local CAZyme SQL database"""


import logging

from tqdm import tqdm

from scraper.sql.sql_orm import (
    Cazyme,
    CazyFamily,
    Cazymes_Genbanks,
    Genbank,
    Kingdom,
    Taxonomy,
)


# Retrieve all records from specific CAZy classes and families


def get_prim_gnbk_acc_from_clss_fam_with_seq(session, config_dict):
    """Retrieve Primary GenBank accession and associated Taxonomy and EC number records of CAZymes.

    :param session: open SQL database session
    :param config_dict: dict defining CAZy classes and families to retrieved records from

    Return list of query results.
    """
    logger = logging.getLogger(__name__)

    genbank_query_class = None
    genbank_query_family = None

    if len(config_dict["classes"]) != 0:
        logger.info(f"Retrieve records for specific CAZy classes: {config_dict['classes']}")

        # retrieve list of CAZy classes to get sequences for
        cazy_classes = config_dict["classes"]

        for cazy_class in tqdm(cazy_classes, desc="Parsing CAZy classes"):
            # retrieve class name abbreviation
            cazy_class = cazy_class[((cazy_class.find("(")) + 1):((cazy_class.find(")")) - 1)]

            # get the CAZymes within the CAZy class
            class_subquery = session.query(Cazyme.cazyme_id).\
                join(CazyFamily, Cazyme.families).\
                filter(CazyFamily.family.regexp(rf"{cazy_class}\d+")).\
                subquery()

            genbank_query_class = session.query(
                Genbank,
                Cazymes_Genbanks,
                Cazyme,
                Taxonomy,
                Kingdom,
            ).\
                join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
                join(Cazyme, (Cazyme.taxonomy_id == Taxonomy.taxonomy_id)).\
                join(Cazymes_Genbanks, (Cazymes_Genbanks.cazyme_id == Cazyme.cazyme_id)).\
                join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
                filter(Cazyme.cazyme_id.in_(class_subquery)).\
                filter(Cazymes_Genbanks.primary == True).\
                filter(Genbank.sequence != None).\
                all()

    # Retrieve protein sequences for specified families
    for key in config_dict:
        if key == "classes":
            continue
        if config_dict[key] is None:
            continue  # no families to parse

        logger.info(f"Retrieving records from CAZy families in {key}")

        for family in tqdm(config_dict[key], desc=f"Parsing families in {key}"):

            if family.find("_") != -1:  # subfamily
                # Retrieve GenBank accessions catalogued under the subfamily
                family_subquery = session.query(Cazyme.cazyme_id).\
                    join(CazyFamily, Cazyme.families).\
                    filter(CazyFamily.subfamily == family).\
                    subquery()

            else:  # family
                # Retrieve GenBank accessions catalogued under the family
                family_subquery = session.query(Cazyme.cazyme_id).\
                    join(CazyFamily, Cazyme.families).\
                    filter(CazyFamily.family == family).\
                    subquery()

            genbank_query_family = session.query(
                Genbank,
                Cazymes_Genbanks,
                Cazyme,
                Taxonomy,
                Kingdom,
            ).\
                join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
                join(Cazyme, (Cazyme.taxonomy_id == Taxonomy.taxonomy_id)).\
                join(Cazymes_Genbanks, (Cazymes_Genbanks.cazyme_id == Cazyme.cazyme_id)).\
                join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
                filter(Cazyme.cazyme_id.in_(family_subquery)).\
                filter(Cazymes_Genbanks.primary == True).\
                filter(Genbank.sequence != None).\
                all()

    return genbank_query_class, genbank_query_family


def get_all_gnbk_acc_from_clss_fams_with_seq(session, config_dict):
    """Retrieve all GenBank records and associated taxonomic and biochemical (EC number records) of
    CAZymes from user specified CAZy classes and families.

    :param session: open SQL database session
    :param config_dict: dict defining CAZy classes and families to retrieved records from

    Return list of query results.
    """
    logger = logging.getLogger(__name__)

    genbank_query_class = None
    genbank_query_family = None

    if len(config_dict["classes"]) != 0:
        logger.info(f"Retrieve records from specific CAZy classes: {config_dict['classes']}")

        # retrieve list of CAZy classes to get sequences for
        cazy_classes = config_dict["classes"]

        for cazy_class in tqdm(cazy_classes, desc="Parsing CAZy classes"):
            # retrieve class name abbreviation
            cazy_class = cazy_class[((cazy_class.find("(")) + 1):((cazy_class.find(")")) - 1)]

            # get the CAZymes within the CAZy class
            class_subquery = session.query(Cazyme.cazyme_id).\
                join(CazyFamily, Cazyme.families).\
                filter(CazyFamily.family.regexp(rf"{cazy_class}\d+")).\
                subquery()

            genbank_query_class = session.query(
                Genbank,
                Cazymes_Genbanks,
                Cazyme,
                Taxonomy,
                Kingdom,
            ).\
                join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
                join(Cazyme, (Cazyme.taxonomy_id == Taxonomy.taxonomy_id)).\
                join(Cazymes_Genbanks, (Cazymes_Genbanks.cazyme_id == Cazyme.cazyme_id)).\
                join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
                filter(Cazyme.cazyme_id.in_(class_subquery)).\
                filter(Genbank.sequence != None).\
                all()

    # Retrieve protein sequences for specified families
    for key in config_dict:
        if key == "classes":
            continue
        if config_dict[key] is None:
            continue  # no families to parse
    
        logger.info(f"Retrieve records for specific CAZy families from: {key}")

        for family in tqdm(config_dict[key], desc=f"Parsing families in {key}"):

            if family.find("_") != -1:  # subfamily
                # Retrieve GenBank accessions catalogued under the subfamily
                family_subquery = session.query(Cazyme.cazyme_id).\
                    join(CazyFamily, Cazyme.families).\
                    filter(CazyFamily.subfamily == family).\
                    subquery()

            else:  # family
                # Retrieve GenBank accessions catalogued under the family
                family_subquery = session.query(Cazyme.cazyme_id).\
                    join(CazyFamily, Cazyme.families).\
                    filter(CazyFamily.family == family).\
                    subquery()

            genbank_query_family = session.query(
                Genbank,
                Cazymes_Genbanks,
                Cazyme,
                Taxonomy,
                Kingdom,
            ).\
                join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
                join(Cazyme, (Cazyme.taxonomy_id == Taxonomy.taxonomy_id)).\
                join(Cazymes_Genbanks, (Cazymes_Genbanks.cazyme_id == Cazyme.cazyme_id)).\
                join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
                filter(Cazyme.cazyme_id.in_(family_subquery)).\
                filter(Genbank.sequence != None).\
                all()
    
    return genbank_query_class, genbank_query_family


# Retrieve records from all CAZy classes and families


def get_prim_genbank_accessions_with_seq(session):
    """Retrieve all PRIMARY GenBank accessions that have a sequence in the database

    :param session: open SQL database session

    Return database query result.
    """
    genbank_query = session.query(Genbank, Cazymes_Genbanks, Cazyme, Taxonomy, Kingdom).\
        join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
        join(Cazyme, (Cazyme.taxonomy_id == Taxonomy.taxonomy_id)).\
        join(Cazymes_Genbanks, (Cazymes_Genbanks.cazyme_id == Cazyme.cazyme_id)).\
        join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
        filter(Cazymes_Genbanks.primary == True).\
        filter(Genbank.sequence != None).\
        all()
    return genbank_query


def get_genbank_accessions_with_seq(session):
    """Retrieve ALL GenBank accessions that have a sequence in the database

    :param session: open SQL database session

    Return database query result.
    """
    genbank_query = session.query(Genbank, Cazymes_Genbanks, Cazyme, Taxonomy, Kingdom).\
        join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
        join(Cazyme, (Cazyme.taxonomy_id == Taxonomy.taxonomy_id)).\
        join(Cazymes_Genbanks, (Cazymes_Genbanks.cazyme_id == Cazyme.cazyme_id)).\
        join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
        filter(Genbank.sequence != None).\
        all()
    return genbank_query


def get_user_accessions(accessions_list, session):
    """Retrieve GenBank records from the local CAZyme database that contain an accession from a user
    specified list of GenBank accessions.

    :param accessions_list: list of GenBank protein accessions
    :param session: open SQL db session

    Return list of query results.
    """
    logger = logging.getLogger(__name__)

    all_query_results = []

    for accession in tqdm(accessions_list, desc="Retrieving records for provided accessions"):
        genbank_query = session.query(Genbank, Cazymes_Genbanks, Cazyme, Taxonomy, Kingdom).\
            join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
            join(Cazyme, (Cazyme.taxonomy_id == Taxonomy.taxonomy_id)).\
            join(Cazymes_Genbanks, (Cazymes_Genbanks.cazyme_id == Cazyme.cazyme_id)).\
            join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
            filter(Genbank.genbank_accession == accession).\
            filter(Genbank.sequence != None).\
            all()
        if len(genbank_query) == 0:
            all_query_results.append(genbank_query)
    if len(all_query_results) == 0:
        logger.warning(
            "No GenBank records in the local CAZyme database were found for the provided accessions"
        )

    return genbank_query
