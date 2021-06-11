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
#
# Bio.PDB reference:
# Hamelryck, T., Manderick, B. (2003) PDB parser and structure class 
# implemented in Python. Bioinformatics 19: 2308–2310
"""Retrieve PDB structures from RSCB PDB and write to disk"""


import logging
import os
import sys

import pandas as pd

from datetime import datetime
from typing import List, Optional

import Bio.PDB

from tqdm import tqdm

from scraper.expand import get_accession_chunks
from scraper.sql.sql_orm import (
    Cazyme,
    CazyFamily,
    EC,
    Kingdom,
    Pdb,
    Taxonomy,
    get_db_session,
)
from scraper.utilities import config_logger, file_io, parse_configuration
from scraper.utilities.parsers import build_pdb_structures_parser


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Set up programme and initate run."""
    start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # used in terminating message
    start_time = pd.to_datetime(start_time)

    # parse cmd-line arguments
    if argv is None:
        parser = build_pdb_structures_parser()
    else:
        parser = build_pdb_structures_parser(argv)
    args = parser.parse_args()

    # build logger
    if logger is None:
        logger = logging.getLogger(__name__)
        config_logger(args)
    
    # validate PDB file choices
    validate_pdb_file_types(args)

    session = get_database_session(args)

    # create output directory
    if args.outdir is not None:
        file_io.make_output_directory(args.outdir, args.force, args.nodelete)
    # else write files to the CWD

    # get list of all PDB accessions to retrieve structure files for
    pdb_accessions = get_pdb_accessions(args=args, session=session)
    pdb_accessions = [acc.pdb_accession for acc in pdb_accessions]

    download_pdb_structures(pdb_accessions, args)

    end_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # used in terminating message
    end_time = pd.to_datetime(start_time)
    end_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    end_time = pd.to_datetime(end_time)
    total_time = end_time - start_time

    logger.info(
        "Finished dowloading protein structures from PDB."
        "Terminating program.\n"
        f"Scrape initated at {start_time}\n"
        f"Scrape finished at {end_time}\n"
        f"Total run time: {total_time}"
    )

    print(
        "=====================cazy_webscraper-expand-pdb_structures=====================\n"
        "Finished dowloading protein structures from PDB."
        "Terminating program.\n"
        f"Scrape initated at {start_time}\n"
        f"Scrape finished at {end_time}\n"
        f"Total run time: {total_time}\n"
    )


def validate_pdb_file_types(args):
    """Check valid file types for PDB were specified by the user.

    :param args: cmd-line args parser

    Return nothing.
    """
    logger = logging.getLogger(__name__)

    valid_choices = ['mmCif', 'pdb', 'xml', 'mmtf', 'bundle']

    invalid_choices = []

    user_choices = (args.pdb).split(",")

    for choice in user_choices:
        if choice not in valid_choices:
            invalid_choices.append(choice)

    if len(invalid_choices) != 0:
        logger.error(
            f"Invalid file option selected: {invalid_choices}.\n"
            f"The valid choices are: {valid_choices}.\n"
            "Terminating program."
        )
        sys.exit(1)

    return


def get_database_session(args):
    """Get session to local CAZyme database.

    :param args: cmd-line args parser
    
    Return open SQL database session
    """
    logger = logging.getLogger(__name__)

    # check path to the local db is valid
    if os.path.isfile(args.database) is False:
        logger.error(
            "Could not find local CAZy database.\n"
            "Check path is correct.\n"
            "Terminating programme."
        )
        sys.exit(1)

    # get database session
    try:
        session = get_db_session(args)
    except Exception as err:
        logger.error(
            "Could not connect to local CAZyme database.The following error was raised:\n"
            f"{err}\nTerminating program\n"
        )
        sys.exit(1)
    
    return session


def get_pdb_accessions(args, session):
    """Retrieve the PDB accessions of CAZymes matching user criteria.

    :param args: cmd-line args parser
    :param session: open SQL database session

    Return list of unique PDB accession.
    """
    # retrieve configuration data, as a dict of CAZy classes and families to retrieve seqs for
    parse_configuration_path = parse_configuration.__file__
    (
        config_dict,
        taxonomy_filters,
        kingdoms,
        ec_filters,
    ) = parse_configuration.parse_configuration_for_cazy_database(
        parse_configuration_path,
        args,
    )

    if config_dict is None:
        pdb_query_results = get_all_pdb_accessions(session)

    else:
        pdb_query_results = get_pdb_acc_from_clss_fams(session, config_dict)

    # parse the PDB query results to retain only those that match the user's criteria
    # object order Pdb, Cazyme, Taxonomy, Kingdom, EC

    pdb_accessions = []

    # apply no filters
    if (taxonomy_filters is None) and (kingdoms is None) and (ec_filters is None):
        accessions = [item[0] for item in pdb_query_results]
        pdb_accessions = [x for x in accessions if "NA" != x]

    # apply only taxonomy_filters
    elif (taxonomy_filters is not None) and (kingdoms is None) and (ec_filters is None):
        for result in pdb_query_results:
            if result[0] == "NA":
                continue
                
            # check if CAZyme records meets the taxonomy criteria
            source_organism = result[-3].genus + " " + result[-3].species
            if any(filter in source_organism for filter in taxonomy_filters):
                pdb_accessions.append(result[0])
                continue
    
    # apply only kingdoms filters
    elif (taxonomy_filters is None) and (kingdoms is not None) and (ec_filters is None):
        for result in pdb_query_results:
            if result[0] == "NA":
                continue

            if result[-2].kingdom in kingdoms:
                pdb_accessions.append(result[0])
                continue
        
    # apply only EC number filters
    elif (taxonomy_filters is None) and (kingdoms is None) and (ec_filters is not None):
        for result in pdb_query_results:
            if result[0] == "NA":
                continue

            # check if the CAZyme record meets the EC filter requirements
            if result[-1].ec_number in ec_filters:
                pdb_accessions.append(result[0])
                continue
    
    # apply taxonomy and kingdom filters
    elif (taxonomy_filters is not None) and (kingdoms is not None) and (ec_filters is None):
        for result in pdb_query_results:
            if result[0] == "NA":
                continue

            # check if CAZyme records meets the taxonomy criteria
            source_organism = result[-3].genus + " " + result[-3].species
            if any(filter in source_organism for filter in taxonomy_filters):

                # check if CAZyme meets kingdom criteria (it should do but worth to check)
                if result[-2].kingdom in kingdoms:
                    pdb_accessions.append(result[0])
                    continue

    # apply taxonomy and EC number filters
    elif (taxonomy_filters is not None) and (kingdoms is None) and (ec_filters is not None):
        for result in pdb_query_results:
            if result[0] == "NA":
                continue

            # check if CAZyme records meets the taxonomy criteria
            source_organism = result[-3].genus + " " + result[-3].species
            if any(filter in source_organism for filter in taxonomy_filters):

                # check if the CAZyme record meets the EC filter requirements
                if result[-1].ec_number in ec_filters:
                    pdb_accessions.append(result[0])
                    continue

    # apply kingdom and EC filters
    elif (taxonomy_filters is None) and (kingdoms is not None) and (ec_filters is not None):
        for result in pdb_query_results:
            if result[0] == "NA":
                continue

            # check kingdom requirements are met
            if result[-2].kingdom in kingdoms:

                # check if the CAZyme record meets the EC filter requirements
                if result[-1].ec_number in ec_filters:
                    pdb_accessions.append(result[0])
                    continue

    # apply taxonomy, kingdom and EC number filters
    else:
        for result in pdb_query_results:
            if result[0] == "NA":
                continue

            # check if CAZyme records meets the taxonomy criteria
            source_organism = result[-3].genus + " " + result[-3].species
            if any(filter in source_organism for filter in taxonomy_filters):

                # check kingdom requirements are met
                if result[-2].kingdom in kingdoms:

                    # check if the CAZyme record meets the EC filter requirements
                    if result[-1].ec_number in ec_filters:
                        pdb_accessions.append(result[0])
                        continue

    return list(set(pdb_accessions))


def get_pdb_acc_from_clss_fams(session, config_dict):
    """Get PDB accessions of proteins from specific CAZy classes and/or families.

    :param session: open SQL db session
    :param config_dict: dict of CAZy classes and families to retrieve structures for

    Return database query result.
    """
    pdb_query_class = []
    pdb_query_family = []

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

            # Retrieve PDB accessions for the selected CAZymes
            pdb_query_class = session.query(Pdb, Cazyme, Taxonomy, Kingdom, EC).\
                join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
                join(Cazyme, (Cazyme.taxonomy_id == Taxonomy.taxonomy_id)).\
                join(Cazyme.pdbs).\
                join(Cazyme.ecs).\
                filter(Cazyme.cazyme_id.in_(class_subquery)).\
                all()       

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

            pdb_query_family = session.query(Pdb, Cazyme, Taxonomy, Kingdom, EC).\
                join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
                join(Cazyme, (Cazyme.taxonomy_id == Taxonomy.taxonomy_id)).\
                join(Cazyme.pdbs).\
                join(Cazyme.ecs).\
                filter(Cazyme.cazyme_id.in_(family_subquery)).\
                all()

    pdb_query_results = pdb_query_class + pdb_query_family
    return pdb_query_results


def get_all_pdb_accessions(session):
    """Retrieve ALL PDB accessions in the database.

    :param session: open SQL database session

    Return database query result.
    """
    pdb_query = session.query(Pdb, Cazyme, Taxonomy, Kingdom, EC).\
        join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
        join(Cazyme, (Cazyme.taxonomy_id == Taxonomy.taxonomy_id)).\
        join(Cazyme.pdbs).\
        join(Cazyme.ecs).\
        all()

    return pdb_query


def download_pdb_structures(pdb_accessions, args):
    """Download protein structure from the RSCB PDB database

    :param pdb_accession: list of PDB accessions
    :param args: cmd-line args parser

    Return nothing.
    """
    pdbl = Bio.PDB.PDBList()

    logger = logging.getLogger(__name__)
    logger.warning("Starting downloading of structure files from PDB")

    # remove chain annotations from pdb accession
    formated_pdb_accessions = [acc.split("[")[0] for acc in pdb_accessions]

    if args.outdir is None:
        logger.warning("Downloading to current working directory")
        for accession_list in get_accession_chunks(formated_pdb_accessions, args.batch_limit):
            for file_type in tqdm((args.pdb).split(","), desc="Downloading"):
                pdbl.download_pdb_files(
                    pdb_codes=accession_list,
                    file_format=file_type,
                    overwrite=args.overwrite,
                )

    else:
        logger.warning(f"Downloading structures to {args.outdir}")
        for accession_list in get_accession_chunks(formated_pdb_accessions, args.batch_limit):
            for file_type in tqdm((args.pdb).split(","), desc="Downloading"):
                pdbl.download_pdb_files(
                    pdb_codes=accession_list,
                    file_format=file_type,
                    overwrite=args.overwrite,
                    pdir=args.outdir,
                )

    return


if __name__ == "__main__":
    main()
