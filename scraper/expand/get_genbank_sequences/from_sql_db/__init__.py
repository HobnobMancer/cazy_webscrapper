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
"""Retrieve proteins sequences from GenBank and populate the local database and write to FASTA"""


import logging
import os
import re
import sys
import time

import pandas as pd

from datetime import datetime
from typing import List, Optional

from Bio import Entrez, SeqIO
from tqdm import tqdm

from scraper.sql.sql_orm import (
    Cazyme,
    CazyFamily,
    Cazymes_Genbanks,
    Genbank,
    Kingdom,
    Taxonomy,
    get_db_session,
)
from scraper.utilities import config_logger, file_io, parse_configuration
from scraper.utilities.parsers import build_genbank_sequences_parser


def sequences_for_proteins_from_db(date_today, args):
    """Coordinate retrievel of protein sequences for proteins in a SQL database.
    
    :param date_today: str, date script was invoked, used for naming files
    :param args: cmd-line args parser
    
    Return nothing.
    """
    logger = logging.getLogger(__name__)

    # get database session
    session = get_db_session(args)

    file_io.make_output_directory(args.fasta, args.force, args.nodelete)

    if args.blastdb is not None:  # build directory to store FASTA file for BLAST db
        file_io.make_output_directory(args.blastdb, args.force, args.nodelete)

    # retrieve configuration data, as a dict of CAZy classes and families to retrieve seqs for
    parse_configuration_path = parse_configuration.__file__
    (
        config_dict, taxonomy_filters, kingdoms, ec_filters,
    ) = parse_configuration.parse_configuration_for_cazy_database(
        parse_configuration_path,
        args,
    )

    # check if config dict contains any specific CAZy classes or families

    if config_dict is None:
        # retrieve all protein GenBank accessions, and apply taxonomic and EC number filters

        if args.update:
            # get sequence for everything without a sequence and those with newer remote sequence
            add_and_update_all_sequences(date_today, taxonomy_filters, kingdoms, session, args)

        else:
            # get sequences for everything without a sequence
            get_all_missing_sequences(
                date_today,
                taxonomy_filters,
                kingdoms,
                session,
                args,
            )

    else:
        # apply CAZy class and family filters, and apply taxonomic and EC number filters

        # get sequences for only specified classes/families
        if args.update:
            update_sequences_for_specific_records(
                date_today,
                config_dict,
                taxonomy_filters,
                kingdoms,
                session,
                args,
            )

        else:
            get_missing_sequences_for_specific_records(
                date_today,
                config_dict,
                taxonomy_filters,
                kingdoms,
                session,
                args,
            )

    return


def add_and_update_all_sequences(date_today, taxonomy_filters, kingdoms, session, args):
    """Retrieve sequences for all proteins in the database.

    For records with no sequences, add the retrieved sequence.
    For records with a sequence, check if the remove sequence is more recent than the existing
    sequence. It it is, update the local sequence.

    :param date_today: str, today's date, used for logging the date the seq is retrieved in the db
    :param taxonomy_filters: set of genera, species and strains to retrieve sequences for
    :param kingdoms: set of taxonomy Kingdoms to retrieve sequences for
    :param session: open SQLite db session
    :param args: cmd-line argument parser

    Return nothing.
    """
    logger = logging.getLogger(__name__)

    # retrieve only sequences for primary GenBank accessions, and those without sequences
    if args.primary is True:
        logger.warning(
            "Retrieving sequences for all primary GenBank accessions that do not have sequences\n"
            "and those whose sequences have been updated in NCBI "
            "since they were retrieved previously"
        )
        genbank_query = session.query(Genbank, Cazymes_Genbanks, Cazyme, Taxonomy, Kingdom).\
            join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
            join(Cazyme, (Cazyme.taxonomy_id == Taxonomy.taxonomy_id)).\
            join(Cazymes_Genbanks, (Cazymes_Genbanks.cazyme_id == Cazyme.cazyme_id)).\
            join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
            filter(Cazymes_Genbanks.primary == True).\
            all()

    # retrieve sequences for all GenBank accessions
    else:
        logger.warning(
            "Retrieving sequences for all GenBank accessions that do not have sequences\n"
            "and those whose sequences have been updated in NCBI "
            "since they were retrieved previously"
        )
        genbank_query = session.query(Genbank, Cazymes_Genbanks, Cazyme, Taxonomy, Kingdom).\
            join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
            join(Cazyme, (Cazyme.taxonomy_id == Taxonomy.taxonomy_id)).\
            join(Cazymes_Genbanks, (Cazymes_Genbanks.cazyme_id == Cazyme.cazyme_id)).\
            join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
            all()

    # create dictionary of {genbank_accession: 'sequence update date' (str)}
    accessions = extract_accessions_and_dates(genbank_query, taxonomy_filters)

    if len(accessions.keys()) == 0:
        logger.warning(
            "Did not retrieve any GenBank accessions from the local database.\n"
            "Not adding sequences to the local database."
        )
        return

    accessions = get_accessions_for_new_sequences(accessions)  # list of genkbank_accession

    if len(accessions) == 0:
        logger.warning(
            "Did not retrieve any GenBank accessions whose sequences need updating.\n"
            "Not adding sequences to the local database."
        )
        return

    # separate accesions in to separate lists of length args.epost, epost doesn't like more than 200
    accessions = get_accession_chunks(accessions, args.epost)  # args.epost = number per chunk
    for lst in accessions:
        get_sequences_add_to_db(lst, date_today, session, args)
    
    return


def get_all_missing_sequences(date_today, taxonomy_filters, kingdoms, session, args):
    """Retrieve protein sequences for all CAZymes in the local CAZy database that don't have seq.

    :param date_today: str, today's date, used for logging the date the seq is retrieved in the db
    :param taxonomy_filters: set of genera, species and strains to restrict sequence retrieval
    :param kingdoms: set of taxonomy Kingdoms to retrieve sequences for
    :param session: open SQLite db session
    :param args: cmd-line argument parser

    Return nothing.
    """
    logger = logging.getLogger(__name__)

    # retrieve only sequences for primary GenBank accessions, and those without sequences
    if args.primary is True:
        logger.warning(
            "Retrieving sequences for all primary GenBank accessions that do not have sequences"
        )
        genbank_query = session.query(Genbank, Cazymes_Genbanks, Cazyme, Taxonomy, Kingdom).\
            join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
            join(Cazyme, (Cazyme.taxonomy_id == Taxonomy.taxonomy_id)).\
            join(Cazymes_Genbanks, (Cazymes_Genbanks.cazyme_id == Cazyme.cazyme_id)).\
            join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
            filter(Cazymes_Genbanks.primary == True).\
            filter(Genbank.sequence == None).\
            all()

    # retrieve sequences for all GenBank accessions without sequences
    else:
        logger.warning(
            "Retrieving sequences for all GenBank accessions that do not have sequences"
        )
        genbank_query = session.query(Genbank, Cazymes_Genbanks, Cazyme, Taxonomy, Kingdom).\
            join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
            join(Cazyme, (Cazyme.taxonomy_id == Taxonomy.taxonomy_id)).\
            join(Cazymes_Genbanks, (Cazymes_Genbanks.cazyme_id == Cazyme.cazyme_id)).\
            join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
            filter(Genbank.sequence == None).\
            all()

    # retrieve the genbank_accessions
    accessions = extract_accessions(genbank_query, taxonomy_filters)

    if len(accessions) == 0:
        logger.warning(
            "Did not retrieve any GenBank accessions from the local database\n"
            "that have sequences missing. Not adding sequences to the local database."
        )
        return

    # separate accesions in to separate lists of length args.epost, epost doesn't like more than 200
    accessions = get_accession_chunks(accessions, args.epost)  # args.epost = number per chunk
    for lst in accessions:
        get_sequences_add_to_db(lst, date_today, session, args)
    return




def get_missing_sequences_for_specific_records(
    date_today,
    config_dict,
    taxonomy_filters,
    kingdoms,
    session,
    args,
):
    """Coordinate getting the sequences for specific CAZymes, not with seqs in the db.

    :param date_today: str, today's date, used for logging the date the seq is retrieved in the db
    :param config_dict: dict, defines CAZy classes and families to get sequences for
    :param taxonomy_filters: set of genera, species and strains to restrict sequence retrieval
    :param kingdoms: set of taxonomy Kingdoms to retrieve sequences for
    :param session: open SQL database session
    :param args: cmd-line args parser

    Return nothing.
    """
    logger = logging.getLogger(__name__)
    logger.warning(
        "Retrieving sequences for GenBank accessions that do not have a sequence in the database"
    )

    # start with the classes
    if len(config_dict["classes"]) != 0:
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

            # retrieve the GenBank accessions of the CAZymes in the CAZy class without seqs
            if args.primary:
                genbank_query = session.query(Genbank, Cazymes_Genbanks, Cazyme, Taxonomy, Kingdom).\
                    join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
                    join(Cazyme, (Cazyme.taxonomy_id == Taxonomy.taxonomy_id)).\
                    join(Cazymes_Genbanks, (Cazymes_Genbanks.cazyme_id == Cazyme.cazyme_id)).\
                    join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
                    filter(Cazyme.cazyme_id.in_(class_subquery)).\
                    filter(Cazymes_Genbanks.primary == True).\
                    filter(Genbank.sequence == None).\
                    all()
            else:
                genbank_query = session.query(Genbank, Cazymes_Genbanks, Cazyme, Taxonomy, Kingdom).\
                    join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
                    join(Cazyme, (Cazyme.taxonomy_id == Taxonomy.taxonomy_id)).\
                    join(Cazymes_Genbanks, (Cazymes_Genbanks.cazyme_id == Cazyme.cazyme_id)).\
                    join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
                    filter(Cazyme.cazyme_id.in_(class_subquery)).\
                    filter(Genbank.sequence == None).\
                    all()

            # retrieve the genbank_accessions from the sql collection object returned by the query
            accessions = extract_accessions(genbank_query, taxonomy_filters)

            if len(accessions) == 0:
                logger.warning(
                    f"Did not retrieve any GenBank accessions for the CAZy class {cazy_class}\n"
                    "that have missing sequences. Not adding sequences to the local database."
                )
                continue

            # separate accesions in to separate lists of length args.epost
            # epost doesn't like posting more than 200 at once
            accessions = get_accession_chunks(accessions, args.epost)  # args.epost = number/chunk
            for lst in accessions:
                get_sequences_add_to_db(lst, date_today, session, args)
            continue

    # Retrieve protein sequences for specified families
    for key in config_dict:
        if key == "classes":
            continue
        if config_dict[key] is None:
            continue  # no families to parse

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

            # get the GenBank accessions of thes CAZymes, without sequences
            if args.primary:
                genbank_query = session.query(Genbank, Cazymes_Genbanks, Cazyme, Taxonomy, Kingdom).\
                    join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
                    join(Cazyme, (Cazyme.taxonomy_id == Taxonomy.taxonomy_id)).\
                    join(Cazymes_Genbanks, (Cazymes_Genbanks.cazyme_id == Cazyme.cazyme_id)).\
                    join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
                    filter(Cazyme.cazyme_id.in_(family_subquery)).\
                    filter(Cazymes_Genbanks.primary == True).\
                    filter(Genbank.sequence == None).\
                    all()
            else:
                genbank_query = session.query(Genbank, Cazymes_Genbanks, Cazyme, Taxonomy, Kingdom).\
                    join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
                    join(Cazyme, (Cazyme.taxonomy_id == Taxonomy.taxonomy_id)).\
                    join(Cazymes_Genbanks, (Cazymes_Genbanks.cazyme_id == Cazyme.cazyme_id)).\
                    join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
                    filter(Cazyme.cazyme_id.in_(family_subquery)).\
                    filter(Genbank.sequence == None).\
                    all()

            # retrieve a list of GenBank accessions from the sql collection returned from the query
            accessions = extract_accessions(genbank_query, taxonomy_filters)

            if len(accessions) == 0:
                logger.warning(
                    f"Did not retrieve any GenBank accessions for the CAZy class {family}\n"
                    "that have missing sequences. Not adding sequences to the local database."
                )
                continue
            # separate accesions in to separate lists of length args.epost
            # epost doesn't like posting more than 200 at once
            accessions = get_accession_chunks(accessions, args.epost)  # args.epost = acc/chunk
            for lst in accessions:
                get_sequences_add_to_db(lst, date_today, session, args)

    return


def update_sequences_for_specific_records(
    date_today,
    config_dict,
    taxonomy_filters,
    kingdoms,
    session,
    args,
):
    """Coordinate getting the sequences for specific CAZymes, not with seqs in the db nad those
    whose seq in NCBI has been updated since the last retrieval.

    For records with no sequences, add the retrieved sequence.
    For records with a sequence, check if the remove sequence is more recent than the existing
    sequence. It it is, update the local sequence.

    :param date_today: str, today's date, used for logging the date the seq is retrieved in the db
    :param config_dict: dict, defines CAZy classes and families to get sequences for
    :param taxonomy_filters: set of genera, species and strains to restrict sequence retrieval
    :param kingdoms: set of taxonomy Kingdoms to retrieve sequences for
    :param session: open SQL database session
    :param args: cmd-line args parser

    Return nothing.
    """
    logger = logging.getLogger(__name__)
    logger.warning(
        "Retrieving sequences for GenBank accessions that do not have a sequence in the database,\n"
        "and those whose sequence in NCBI has been updated since they were previously retrieved."
    )

    # start with the classes
    if len(config_dict["classes"]) != 0:
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

            # retrieve the GenBank accessions of the CAZymes in the CAZy class without seqs
            if args.primary:
                genbank_query = session.query(Genbank, Cazymes_Genbanks, Cazyme, Taxonomy, Kingdom).\
                    join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
                    join(Cazyme, (Cazyme.taxonomy_id == Taxonomy.taxonomy_id)).\
                    join(Cazymes_Genbanks, (Cazymes_Genbanks.cazyme_id == Cazyme.cazyme_id)).\
                    join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
                    filter(Cazyme.cazyme_id.in_(class_subquery)).\
                    filter(Cazymes_Genbanks.primary == True).\
                    all()
            else:
                genbank_query = session.query(Genbank, Cazymes_Genbanks, Cazyme, Taxonomy, Kingdom).\
                    join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
                    join(Cazyme, (Cazyme.taxonomy_id == Taxonomy.taxonomy_id)).\
                    join(Cazymes_Genbanks, (Cazymes_Genbanks.cazyme_id == Cazyme.cazyme_id)).\
                    join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
                    filter(Cazyme.cazyme_id.in_(class_subquery)).\
                    all()

            # create dictionary of genbank_accession: 'sequence update date' (str)
            accessions = extract_accessions_and_dates(genbank_query, taxonomy_filters)

            if len(accessions.keys()) == 0:
                logger.warning(
                    f"Did not retrieve any GenBank accessions for the CAZy class {cazy_class}.\n"
                    "Not adding sequences to the local database."
                )
                continue

            accessions = get_accessions_for_new_sequences(accessions)  # list of genkbank_accession

            if len(accessions) == 0:
                logger.warning(
                    "Did not retrieve any GenBank accessions whose sequences need updating for "
                    f"the CAZy class {cazy_class}.\n"
                    "Not adding sequences to the local database."
                )
                continue
            # separate accesions in to separate lists of length args.epost
            # epost doesn't like posting more than 200 at once
            accessions = get_accession_chunks(accessions, args.epost)  # args.epost = acc/chunk
            for lst in accessions:
                get_sequences_add_to_db(lst, date_today, session, args)

    # Retrieve protein sequences for specified families
    for key in config_dict:
        if key == "classes":
            continue
        if config_dict[key] is None:
            continue  # no families to parse

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

            # get the GenBank accessions of thes CAZymes, without sequences
            if args.primary:
                genbank_query = session.query(Genbank, Cazymes_Genbanks, Cazyme, Taxonomy, Kingdom).\
                    join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
                    join(Cazyme, (Cazyme.taxonomy_id == Taxonomy.taxonomy_id)).\
                    join(Cazymes_Genbanks, (Cazymes_Genbanks.cazyme_id == Cazyme.cazyme_id)).\
                    join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
                    filter(Cazyme.cazyme_id.in_(family_subquery)).\
                    filter(Cazymes_Genbanks.primary == True).\
                    filter(Genbank.sequence == None).\
                    all()
            else:
                genbank_query = session.query(Genbank, Cazymes_Genbanks, Cazyme, Taxonomy, Kingdom).\
                    join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
                    join(Cazyme, (Cazyme.taxonomy_id == Taxonomy.taxonomy_id)).\
                    join(Cazymes_Genbanks, (Cazymes_Genbanks.cazyme_id == Cazyme.cazyme_id)).\
                    join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
                    filter(Cazyme.cazyme_id.in_(family_subquery)).\
                    filter(Genbank.sequence == None).\
                    all()

            # create dictionary of {genbank_accession: 'sequence update date' (str)}
            accessions = extract_accessions_and_dates(genbank_query, taxonomy_filters)

            if len(accessions.keys()) == 0:
                logger.warning(
                    f"Did not retrieve any GenBank accessions for the CAZy class {family}.\n"
                    "Not adding sequences to the local database."
                )
                continue

            accessions = get_accessions_for_new_sequences(accessions)  # list of genkbank_accession

            if len(accessions) == 0:
                logger.warning(
                    "Did not retrieve any GenBank accessions whose sequences need updating for "
                    f"the CAZy class {family}.\n"
                    "Not adding sequences to the local database."
                )
                continue
            # separate accesions in to separate lists of length args.epost
            # epost doesn't like posting more than 200 at once
            accessions = get_accession_chunks(accessions, args.epost)  # args.epost = acc/chunk
            for lst in accessions:
                get_sequences_add_to_db(lst, date_today, session, args)

    return


# The following functions are retrieving the list of Genbank accessions to retrieve sequences for #


def extract_accessions(genbank_query, taxonomy_filters):
    """The query contains GenBank accessions and Cazymes_Genbanks records, retrieve the accessions.

    :param genbank_query: sql collection
    :param taxonomy_filters: set of genera, species and strains to restrict retrieval of sequences

    Return a list of GenBank accessions. Each element is a string of a unique accession.
    """
    if taxonomy_filters is None:
        accessions = [item[0] for item in genbank_query]
        return [x for x in accessions if "NA" != x]

    else:
        accessions = []
        for item in genbank_query:
            if item[0] != "NA":  # if GenBank accession not stored as 'NA'
                source_organism = item[-1].genus + item[-1].species
                if any(filter in source_organism for filter in taxonomy_filters):
                    accessions.append(item[0])
    return accessions


def extract_accessions_and_dates(genbank_query, taxonomy_filters):
    """Retrieve the GenBank accessions and retrieval dates of existing sequences from the db query.

    :param genbank_query: sql collection
    :param taxonomy_filters: set of genera, species and strains to restrict retrieval of sequences

    Return a dict {GenBank_accession: retrieval_date}
    """
    accessions = {}

    if taxonomy_filters is None:
        for item in genbank_query:
            if item[0].genbank_accession == "NA":  # no GenBank accession stored in CAZy
                continue
            accessions[item[0].genbank_accession] = item[0].seq_update_date

    else:
        for item in genbank_query:
            if item[0].genbank_accession == "NA":  # no GenBank accession stored in CAZy
                continue
            source_organism = item[-1].genus + item[-1].species
            if any(filter in source_organism for filter in taxonomy_filters):
                accessions[item[0].genbank_accession] = item[0].seq_update_date

    return accessions


def get_accessions_for_new_sequences(accessions):
    """Get the GenBank accessions of sequences to be added to the local database.

    For records currently with no protein sequence, the retrieved protein sequence will be added
    to the record. For records with a sequence, the 'UpdateDate' for the sequence from NCBI will
    be compared against the  'seq_update_date' in the local database. The 'seq_update_date' is the
    'UpdateDate' previosuly retrieved from NCBI. If the NCBI sequence is newer,
    the local database will be updated with the new sequence.

    :param accessions: dict, {GenBank accessions (str):sequence retrieval data (str)}
    :param session: open SQL database session

    Return nothing.
    """
    logger = logging.getLogger(__name__)

    accessions_list = list(accessions.keys())
    accessions_string = ",".join(accessions_list)
    # perform batch query of Entrez
    epost_result = Entrez.read(
        entrez_retry(
            Entrez.epost, "Protein", id=accessions_string, retmode="text",
        )
    )
    # retrieve the web environment and query key from the Entrez post
    epost_webenv = epost_result["WebEnv"]
    epost_query_key = epost_result["QueryKey"]

    # retrieve summary docs to check the sequence 'UpdateDates' in NCBI
    with entrez_retry(
        Entrez.efetch,
        db="Protein",
        query_key=epost_query_key,
        WebEnv=epost_webenv,
        rettype="docsum",
        retmode="xml",
    ) as handle:
        summary_docs = Entrez.read(handle)

    for doc in summary_docs:
        try:
            temp_accession = doc["AccessionVersion"]  # accession of the current working protein
        except KeyError:
            logger.warning(
                f"Retrieved protein with accession {temp_accession} but this accession is not in "
                "the local database.\n"
                "Not retrieving a sequence for this accession."
            )
            continue
        previous_data = accessions[temp_accession]
        if previous_data is not None:
            # sequence retrieved previosuly, thus check if the NCBI seq has been updated since
            previous_data = previous_data.split("/")  # Y=[0], M=[1], D=[]
            update_date = doc["UpdateDate"]
            update_date = update_date.split("/")  # Y=[0], M=[1], D=[]
            if datetime.date(
                previous_data[0], previous_data[1], previous_data[2],
            ) < datetime.data(
                update_date[0], update_date[1], update_date[2],
            ) is False:
                # the sequence at NCBI has not been updated since the seq was retrieved
                # thus no need to retrieve it again
                accessions_list.remove(temp_accession)

    return accessions_list
