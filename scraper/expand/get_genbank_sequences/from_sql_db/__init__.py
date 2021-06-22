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
import math
import sys

from tqdm import tqdm

from scraper.expand import get_accession_chunks
from scraper.expand.get_genbank_sequences import ncbi
from scraper.expand.get_genbank_sequences.from_sql_db import query_sql_db
from scraper.sql.sql_orm import get_db_session
from scraper.sql.sql_interface import log_scrape_in_db
from scraper.utilities import parse_configuration


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
            "Could not connect to local CAZyme database.\n"
            "The following error was raised:\n"
            f"{err}\nTerminating program\n"
        )
        sys.exit(1)

    # retrieve configuration data, as a dict of CAZy classes and families to retrieve seqs for
    (
        config_dict,
        taxonomy_filters,
        kingdoms,
        ec_filters,
    ) = parse_configuration.parse_configuration_for_cazy_database(args)

    if args.fasta_only is not None:
        logger.info("Adding log of sequence retrieval to the local CAZyme database")
        log_scrape_in_db(
            data_addition="GenBank sequences",
            time_stamp=f"<Seq addition>",
            config_dict=config_dict,
            taxonomy_filters=taxonomy_filters,
            kingdoms=kingdoms,
            ec_filters=ec_filters,
            session=session,
            args=args,
        )

    logger.info("Retrieving GenBank accessions that match provided criteria")
    genbank_accessions = get_genbank_accessions(
        args,
        session,
        date_today,
        config_dict,
        taxonomy_filters,
        kingdoms,
        ec_filters,
    )

    if genbank_accessions is None:
        logger.warning(
            "Retrieved no GenBank accessions matching provided criteria\n"
            "Therefore, retrieveing no protein sequences from NCBI"
        )
        return "EXIT"
    
    if len(genbank_accessions) == 0:
        logger.warning(
            "Retrieved 0 GenBank accessions matching provided criteria\n"
            "Therefore, retrieveing no protein sequences from NCBI"
        )
        return "EXIT"

    logger.warning(f"Retrieving sequences for {len(genbank_accessions)} proteins")

    # break up protein_list into multiple, smaller lists for batch querying Entrez
    # batches of greater than 200 can be rejected by Entrez during busy periods
    # args.epost=size of chunks

    accessions_lists_for_individual_queries = []

    for accession_list in tqdm(
        get_accession_chunks(genbank_accessions, args.epost),
        desc="Batch retrieving sequences from NCBI",
        total=(math.ceil(len(genbank_accessions) / args.epost)),
    ):
        try:
            accession_list.remove("NA")
        except ValueError:
            pass

        if args.fasta_only:   # Retrieve seqs and write to FASTA ONLY do not write to db
            try:
                ncbi.get_sequences(
                    accession_list,
                    args,
                )
            except RuntimeError as err:  # typically Some IDs have invalid value and were omitted.
                logger.warning(
                    "RuntimeError raised for accession list. "
                    "Will query accessions individualy after.\n"
                    f"The following error was raised:\n{err}"
                )
                accessions_lists_for_individual_queries.append(accession_list)

        else:  # write seqs to db (and FASTA if enabled)
            try:
                ncbi.get_sequences_add_to_db(
                    accession_list,
                    date_today,
                    session,
                    args,
                )
            except RuntimeError as err:  # typically Some IDs have invalid value and were omitted.
                logger.warning(
                    "RuntimeError raised for accession list. "
                    "Will query accessions individualy after.\n"
                    f"The following error was raised:\n{err}"
                )
                accessions_lists_for_individual_queries.append(accession_list)

    if len(accessions_lists_for_individual_queries) != 0:
        for accession_list in tqdm(
            accessions_lists_for_individual_queries,
            desc="Performing individual queries for records that previously raised errors",
        ):

            if args.fasta_only:  # Retrieve seqs and write to FASTA ONLY do not write to db
                for accession in tqdm(accession_list, desc="Retrieving individual sequences"):
                    try:
                        ncbi.get_sequences(
                            [accession],
                            args,
                        )

                    except RuntimeError as err:
                        logger.warning(
                            f"Queried NCBI for {accession} raised the following RuntimeError:\n"
                            f"{err}"
                        )

            else:  # write seqs to db (and FASTA if enabled)
                for accession in tqdm(accession_list, desc="Retrieving individual sequences"):
                    try:
                        ncbi.get_sequences_add_to_db(
                            [accession],
                            date_today,
                            session,
                            args,
                        )

                    except RuntimeError as err:
                        logger.warning(
                            f"Queried NCBI for {accession} raised the following RuntimeError:\n"
                            f"{err}"
                        )

    return


def get_genbank_accessions(
    args,
    session,
    date_today,
    config_dict,
    taxonomy_filters,
    kingdoms,
    ec_filters,
):
    """Retrieve the GenBank accessions for all proteins for which a sequence will be retrieved.

    :param args: cmd-line argument parser
    :param session: open SQLite db session
    :param date_today: str, date script was invoked
    :param config_dict: dict, defines CAZy classes and families to get sequences for
    :param taxonomy_filters: set of genera, species and strains to retrieve sequences for
    :param kingdoms: set of taxonomy Kingdoms to retrieve sequences for
    :param ec_filters: set of EC numbers annotations CAZymes must have at least one to retrieve
        a sequence

    Return a list of GenBank accessions, containing no duplicate GenBank accessions
    """
    logger = logging.getLogger(__name__)

    if config_dict:  # there are specific CAZy classes/families to retrieve sequences for

        if (args.update) or (args.fasta_only):
            logger.info("Enabled updating sequences in local CAZyme database")
            
            if args.primary:
                logger.warning(
                    "Retrieving sequences for PRIMARY GenBank accessions that:\n"
                    "belong to specific CAZy classes/families AND\n"
                    "do not have a sequence in the db OR the sequence has been updated in NCBI"
                )
                (
                    genbank_query_class,
                    genbank_query_family,
                ) = query_sql_db.get_prim_gnbk_acc_from_clss_fams(
                    session,
                    config_dict,
                )

            else:
                logger.warning(
                    "Retrieving sequences for ALL GenBank accessions that:\n"
                    "belong to specific CAZy classes/families AND\n"
                    "do not have a sequence in the db OR the sequence has been updated in NCBI"
                )
                (
                    genbank_query_class,
                    genbank_query_family,
                ) = query_sql_db.get_all_gnbk_acc_from_clss_fams(
                    session,
                    config_dict,
                )
        
        else:  # retrieve GenBank accesions of records that don't have a sequence
            logger.info(
                "Only retrieving sequences for proteins that do not have a seq in the CAZyme db"
            )

            if args.primary:
                logger.warning(
                    "Retrieving sequences for PRIMARY GenBank accessions that:\n"
                    "belong to specific CAZy classes/families AND do not have a sequence in the db"
                )
                (
                    genbank_query_class,
                    genbank_query_family,
                ) = query_sql_db.get_prim_gnbk_acc_from_clss_fams_no_seq(
                    session,
                    config_dict,
                )
            
            else:
                logger.warning(
                    "Retrieving sequences for ALL GenBank accessions that:\n"
                    "belong to specific CAZy classes/families AND do not have a sequence in the db"
                )
                (
                    genbank_query_class,
                    genbank_query_family,
                ) = query_sql_db.get_all_gnbk_acc_from_clss_fams_no_seq(
                    session,
                    config_dict,
                )

        query_results = genbank_query_class + genbank_query_family
        logger.info(
            f"Retrieved {len(query_results)} records from the specified CAZy classes and families\n"
            "prior to apply any taxonomy and EC number filters"
        )

    else:
        if (args.update) or (args.fasta_only):  # retrieve all GenBank accessions
            logger.info("Enabled updating sequences in local CAZyme database")

            if args.primary:
                logger.warning(
                    "Retrieving sequences for all PRIMARY GenBank accessions that:\n"
                    "Do not have a sequence in the db OR the sequence has been updated in NCBI"
                )
                query_results = query_sql_db.get_prim_genbank_acc_for_update(session)

            else:
                logger.warning(
                    "Retrieving sequences for ALL GenBank accessions that\n"
                    "do not have a sequence in the db OR the sequence has been updated in NCBI"
                )
                query_results = query_sql_db.get_all_genbank_acc_for_update(session)

        else:  # retrieve GenBank accesions of records that don't have a sequence
            logger.info(
                "Only retrieving sequences for proteins that do not have a seq in the CAZyme db"
            )

            if args.primary:
                logger.warning(
                    "Retrieving sequences for all PRIMARY GenBank accessions that\n"
                    "do not have a sequence in the db"
                )
                query_results = query_sql_db.get_prim_genbank_accessions_with_no_seq(session)

            else:
                logger.warning(
                    "Retrieving sequences for ALL GenBank accessions that\n"
                    "do not have a sequence in the db"
                )
                query_results = query_sql_db.get_genbank_accessions_with_no_seq(session)
            
            logger.info(
                f"Retrieved {len(query_results)} records from the local CAZyme database\n"
                "Prior to apply any taxonomy and EC number filters"
            )

    # check if any records were retrived from the querying of the local CAZyme database
    if len(query_results) == 0:
        logger.warning(
            "Retrieved no records from the local CAZyme database mathcing provided criteria"
        )
        return

    # apply taxonomic and EC number filters
    logger.info(
        "Applying any provided taxonomic and EC number filters to records"
        "retrieved from the local CAZyme database"
    )

    filtered_query_results = parse_genbank_query(
        query_results,
        taxonomy_filters,
        kingdoms,
        ec_filters,
        session,
    )

    if args.update:
        logger.info(
            f"Out of {len(filtered_query_results)} records, checking which have no seq in the "
            "local db and\n"
            "which have a seq to update"
        )
        
        genbank_accessions = check_if_to_update(
            filtered_query_results,
            args,
        )

    else:
        logger.info(
            "Extracting GenBank accessions from records retrieved from the local CAZyme database"
        )
        try:
            genbank_accessions = [
                query_result.genbank_accession for query_result in filtered_query_results
            ]
        except AttributeError:
            genbank_accessions = [
                query_result[0].genbank_accession for query_result in filtered_query_results
            ]

    return list(set(genbank_accessions))  # prevent quering the same accession multiple times


def parse_genbank_query(
    genbank_query_results,
    taxonomy_filters,
    kingdoms,
    ec_filters,
    session,
):
    """Parse SQL query result and retrieve GenBank accessions of CAZymes that meet the user cirteria

    :param genbank_query_results: list of query results from the local CAZyme database
    :param taxonomy_filters: list of genera, species and strains to retrict the retrival of seq to
    :param kingdoms: list of taxonomic kingdoms to restrict the retrieval of CAZymes to
    :param ec_filters: list of EC numbers, a CAZyme must be annotated with at least one EC number
        from the list to have its protein sequence retrived.

    Return list of GenBank records from the query results that meet the user's criteria.
    """
    logger = logging.getLogger(__name__)

    tax_filtered_genbank_accessions = []

    # apply no filters
    if (taxonomy_filters is None) and (kingdoms is None) and (ec_filters is None):
        logger.info("Apply no taxonomy or EC number filters")
        return genbank_query_results
    
    # apply only taxonomy filters
    elif (taxonomy_filters is not None) and (kingdoms is None):
        logger.info("Applying taxonomy filters")
        for result in genbank_query_results:
            if result[0].genbank_accession == "NA":
                continue
                
            # check if CAZyme records meets the taxonomy criteria
            source_organism = result[-2].genus + " " + result[-2].species
            if any(filter in source_organism for filter in taxonomy_filters):
                tax_filtered_genbank_accessions.append(result[0])
    
    # apply only kingdom filter
    elif (taxonomy_filters is None) and (kingdoms is not None):
        logger.info("Applying kingdom filters")
        for result in genbank_query_results:
            if result[0].genbank_accession == "NA":
                continue

            if result[-1].kingdom.upper() in (name.upper() for name in kingdoms):
                tax_filtered_genbank_accessions.append(result[0])
    
    # apply taxonomy and kingdom filter
    else:
        logger.info("Applying taxonomy and kingdom filters")
        for result in genbank_query_results:
            if result[0].genbank_accession == "NA":
                continue
                
            # check if CAZyme records meets the taxonomy criteria
            source_organism = result[-2].genus + " " + result[-2].species
            if any(filter in source_organism for filter in taxonomy_filters):
                
                # apply kingdom filter
                if result[-1].kingdom.upper() in (name.upper() for name in kingdoms):
                    tax_filtered_genbank_accessions.append(result[0])
                    continue

    if ec_filters is None:
        return tax_filtered_genbank_accessions
    
    logger.info("Applying EC number filters")
    final_filtered_genbank_accessions = [] 
    for query_result in tqdm(
        tax_filtered_genbank_accessions,
        desc="Applying EC number filter",
    ):
        ec_annotations = query_sql_db.query_ec_number(
            session,
            query_result.genbank_accession,
        )

        # check if any of the EC number annotations for the protein (identified by its genbank 
        # accession) are included in the EC numbers specificed by the user
        if (set(ec_annotations) and set(ec_filters)):
            final_filtered_genbank_accessions.append(query_result)

    return final_filtered_genbank_accessions


def check_if_to_update(genbank_records, args):
    """Coordinate checking if need to update sequences in the local CAZyme database.

    :param genbank_records: list of GenBank records retrieved from the local CAZyme db
    :param args: cmd-lines args parser

    Return list of GenBank accessions from GenBank records from the local CAZyme database that
    either do not have a sequence, or the sequence in the NCBI database has been updated since
    the sequence was last retrieved and added to the local CAZyme database.
    """
    logger = logging.getLogger(__name__)
    
    # separate records with and without sequences, and extract GenBank accession
    logger.warning("Separating GenBank records with and without sequences in the local CAZyme db")
    
    gbk_records_with_seq = {}  # {accession: db_genbank_record}
    gbk_records_without_seq = []

    for record in genbank_records:
        if record.sequence is None:
            gbk_records_without_seq.append(record.genbank_accession)
        else:
            gbk_records_with_seq[record.genbank_accession] = record

    logger.info("Checking which protein sequences are to be updated")
    genbank_seq_to_update = []

    accessions_lists_for_individual_queries = []

    for accession_list in tqdm(
        get_accession_chunks(
            list(gbk_records_with_seq.keys()),
            args.epost,
        ),
        desc="Batch retrieving NCBI to check if to update seq",
        total=(math.ceil(len(gbk_records_with_seq) / args.epost)),
    ):
        try:
            genbank_to_update = ncbi.check_ncbi_seq_data(
                accession_list,
                gbk_records_with_seq,
                args,
            )
            genbank_seq_to_update += genbank_to_update

        except RuntimeError as err:  # typically Some IDs have invalid value and were omitted.
            logger.warning(
                "RuntimeError raised for accession list. Will query accessions individualy after.\n"
                f"The following error was raised:\n{err}"
            )
            accessions_lists_for_individual_queries.append(accession_list)

    if len(accessions_lists_for_individual_queries) != 0:
        for accession_list in tqdm(
            accessions_lists_for_individual_queries,
            desc="Performing individual queries for records that previously raised errors",
        ):
            for accession in tqdm(accession_list, desc="Checking NCBI seq date"):

                try:
                    genbank_to_update = ncbi.check_ncbi_seq_data(
                        [accession],
                        gbk_records_with_seq,
                        args,
                    )
                    genbank_seq_to_update += genbank_to_update

                except RuntimeError as err:
                    logger.warning(
                        f"Queried NCBI for {accession} raised the following RuntimeError:\n"
                        f"{err}"
                    )

    genbank_accessions = gbk_records_without_seq + genbank_seq_to_update

    return genbank_accessions
