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

from scraper.expand.get_genbank_sequences.from_sql_db import query_sql_db
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
    try:
        session = get_db_session(args)
    except Exception as err:
        logger.error(
            "Could not connect to local CAZyme database.\nThe following error was raised:\n"
            f"{err}\nTerminating program\n"
        )
        sys.exit(1)

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

    genbank_accessions = get_genbank_accessions(
        args,
        session,
        config_dict,
        taxonomy_filters,
        kingdoms,
        ec_filters,
    )

    return


def get_genbank_accessions(args, session, config_dict, taxonomy_filters, kingdoms, ec_filters):
    """Retrieve the GenBank accessions for all proteins for which a sequence will be retrieved.

    :param args: cmd-line argument parser
    :param session: open SQLite db session
    :param config_dict: dict, defines CAZy classes and families to get sequences for
    :param taxonomy_filters: set of genera, species and strains to retrieve sequences for
    :param kingdoms: set of taxonomy Kingdoms to retrieve sequences for
    :param ec_filters: set of EC numbers annotations CAZymes must have at least one to retrieve
        a sequence

    Return a list of GenBank accessions, containing no duplicate GenBank accessions
    """
    logger = logging.getLogger(__name__)

    if config_dict:  # there are specific CAZy classes/families to retrieve sequences for
        genbank_query_class, genbank_query_family = get_cazy_class_fam_genbank_records(
            session,
            config_dict,
        )

        class_genbank_accessions = parse_genbank_query(
            genbank_query_class,
            taxonomy_filters,
            kingdoms,
            ec_filters,
        )

        family_genbank_accessions = parse_genbank_query(
            genbank_query_family,
            taxonomy_filters,
            kingdoms,
            ec_filters,
        )

        genbank_accessions = class_genbank_accessions + family_genbank_accessions

    else:
        if args.update:  # retrieve all GenBank accessions

            if args.primary:
                logger.warning(
                    "Retrieving sequences for all PRIMARY GenBank accessions that:\n"
                    "Do not have a sequence in the db OR the sequence has been updated in NCBI"
                )
                genbank_query = query_sql_db.get_all_prim_genbank_acc(session)

            else:
                logger.warning(
                    "Retrieving sequences for all PRIMARY GenBank accessions that\n"
                    "do not have a sequence in the db"
                )
                genbank_query = query_sql_db.get_all_genbank_acc(session)

        else:  # retrieve GenBank accesions of records that don't have a sequence
            if args.primary:
                logger.warning(
                    "Retrieving sequences for all PRIMARY GenBank accessions that\n"
                    "do not have a sequence in the db"
                )
                genbank_query = query_sql_db.get_prim_genbank_accessions_with_no_seq(session)

            else:
                logger.warning(
                    "Retrieving sequences for ALL GenBank accessions that\n"
                    "do not have a sequence in the db"
                )
                genbank_query = query_sql_db.get_genbank_accessions_with_no_seq(session)

        genbank_accessions = parse_genbank_query(
            genbank_query,
            taxonomy_filters,
            kingdoms,
            ec_filters,
        )

    return list(set(genbank_accessions))  # prevent quering the same accession multiple times


def get_cazy_class_fam_genbank_records():
    if args.update:  # retrieve all GenBank accessions
        if args.primary:
            logger.warning(
                "Retrieving sequences for PRIMARY GenBank accessions that:\n"
                "belong to specific CAZy classes/families AND\n"
                "do not have a sequence in the db OR the sequence has been updated in NCBI"
            )
        else:
            logger.warning(
                "Retrieving sequences for PRIMARY GenBank accessions that:\n"
                "belong to specific CAZy classes/families AND\n"
                "do not have a sequence in the db OR the sequence has been updated in NCBI"
            )
    else:  # retrieve GenBank accesions of records that don't have a sequence
        if args.primary:
            logger.warning(
                "Retrieving sequences for PRIMARY GenBank accessions that:\n"
                "belong to specific CAZy classes/families AND do not have a sequence in the db"
            )
        else:
            logger.warning(
                "Retrieving sequences for PRIMARY GenBank accessions that:\n"
                "belong to specific CAZy classes/families AND do not have a sequence in the db"
            )
    return


def parse_genbank_query(genbank_query, taxonomy_filters, kingdoms, ec_filters):
    """Parse SQL query result and retrieve GenBank accessions of CAZymes that meet the user cirteria

    :param:

    Return list of GenBank accessions.
    """
    return genbank_accessions





###
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
