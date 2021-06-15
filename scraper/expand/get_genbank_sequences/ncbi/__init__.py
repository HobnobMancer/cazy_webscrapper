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
"""Module for quering and retrieving data from NCBI.Entrez"""


import logging
import re
import time

from datetime import datetime

from Bio import Entrez, SeqIO

from scraper.sql.sql_orm import Genbank
from scraper.utilities import file_io


def get_sequences_for_dict(accessions, args):
    """Retrieve protein sequences from Entrez for proteins from CAZy dictionary (JSON file).

    :param accessions: list, GenBank accessions of proteins
    :param args: cmb-line args parser

    Return nothing.
    """
    logger = logging.getLogger(__name__)

    # perform batch query of Entrez
    accessions_string = ",".join(accessions)

    # Runtime error captured by try/except function call
    epost_result = Entrez.read(
        entrez_retry(
            Entrez.epost,
            args,
            db="Protein",
            id=accessions_string,
        )
    )

    # retrieve the web environment and query key from the Entrez post
    epost_webenv = epost_result["WebEnv"]
    epost_query_key = epost_result["QueryKey"]

    # retrieve the protein sequences
    with entrez_retry(
        Entrez.efetch,
        args,
        db="Protein",
        query_key=epost_query_key,
        WebEnv=epost_webenv,
        rettype="fasta",
        retmode="text",
    ) as seq_handle:
        for record in SeqIO.parse(seq_handle, "fasta"):
            # retrieve the accession of the record
            temp_accession = record.id  # accession of the current working protein record

            if temp_accession.find("|") != -1:  # sometimes multiple items are listed
                print("temp accession: ", temp_accession)
                success = False   # will be true if finds protein accession

                temp_accessions = temp_accession.split("|")
                for item in temp_accessions:
                    # check if a accession number
                    try:
                        re.match(
                            (
                                r"(\D{3}\d{5,7}\.\d+)|(\D\d(\D|\d){3}\d)|"
                                r"(\D\d(\D|\d){3}\d\D(\D|\d){2}\d)"
                            ),
                            item,
                        ).group()
                        temp_accession = item
                        success = True
                        break
                    except AttributeError:  # raised if not an accession
                        continue

            else:
                success = True  # have protein accession number

            if success is False:
                logger.error(
                    f"Could not retrieve single protein accession from {record.id}\n"
                    f"writing out record with complete record ID: {record.id}"
                )

            # check the retrieved protein accession is in the list of retrieved accession
            if temp_accession not in accessions:
                logger.warning(
                    f"Retrieved the accession {temp_accession} from the record id={record.id},\n"
                    "but this accession was not in the origina CAZy dict.\n"
                    f"Retrieving sequencing, writing the record ID ({record.id}), as the accession"
                )

            if args.fasta is not None:
                file_io.write_out_fasta(record, temp_accession, args)

            if args.blastdb is not None:
                # need all squences in a single FASTA file to create a BLASTdb
                file_io.write_fasta_for_db(record, args)

            # remove the accession from the list
            try:
                accessions.remove(temp_accession)
            except ValueError:
                logger.warning(
                    f"Tried to remove {temp_accession} from list of accessions, "
                    "but it was not in the list of accessions.\n"
                    "The returned accession and the one present in CAZy do not match."
                )

    if len(accessions) != 0:
        logger.warning("Protein sequences for the following CAZymes were not retrieved:")
        for acc in accessions:
            logger.warning(f"GenBank accession: {acc}")

    return


def check_ncbi_seq_data(genbank_accessions, gbk_records_without_seq, args):
    """Query NCBI to see if sequence has been updated since the last retrieval.

    :param genbank_accessions: list of GenBank accessions from the local database
    :param gbk_records_without_seq: dict keyed by GenBank accession and valued by GenBank record
        from the local CAZyme database with the GenBank accession
    :param args: cmd-line args parser

    Return list of GenBank records whose sequence needs updating.
    """
    logger = logging.getLogger(__name__)

    genbank_records_to_update = []

    accessions_list = ",".join(genbank_accessions)

    epost_result = Entrez.read(
        entrez_retry(
            Entrez.epost,
            args,
            db="Protein",
            id=accessions_list,
            retmode="text",
        )
    )

    # retrieve the web environment and query key from the Entrez post
    epost_webenv = epost_result["WebEnv"]
    epost_query_key = epost_result["QueryKey"]

    # retrieve summary docs to check the sequence 'UpdateDates' in NCBI
    with entrez_retry(
        Entrez.efetch,
        args,
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
        except KeyError as err:
            logger.warning(
                "Could not retrieve AccessionVersion from the current record.\n"
                f"{err}"
            )
            continue
        
        if temp_accession not in genbank_accessions:
            logger.warning(
                f"Retrieved protein with accession {temp_accession}\n"
                "but this accession is not in the local database.\n"
                "Not retrieving a sequence for this accession."
            )
            continue

        # convert strings to dates
        previous_retrieval_data = gbk_records_without_seq[temp_accession].seq_update_date
        previous_retrieval_data = previous_retrieval_data.split("/")  # Y=[0], M=[1], D=[]
        previous_retrieval_data = datetime(
            previous_retrieval_data[0],
            previous_retrieval_data[1],
            previous_retrieval_data[2],
        )

        ncbi_seq_date = doc["UpdateDate"]
        ncbi_seq_date = ncbi_seq_date.split("/")  # Y=[0], M=[1], D=[]
        ncbi_seq_date = datetime(
            ncbi_seq_date[0],
            ncbi_seq_date[1],
            ncbi_seq_date[2],
        )

        if ncbi_seq_date > previous_retrieval_data:
            # the sequence at NCBI has been updated since the seq was retrieved, need to update seq
            genbank_records_to_update.append(temp_accession)

    return genbank_records_to_update


def get_sequences_add_to_db(accessions, date_today, session, args):
    """Retrieve protein sequences from Entrez and add to the local database.

    :param accessions: list, GenBank accessions
    :param date_today: str, YYYY/MM/DD
    :param session: open SQL database session
    :param args: cmb-line args parser

    Return nothing.
    """
    logger = logging.getLogger(__name__)
    # perform batch query of Entrez
    accessions_string = ",".join(accessions)
    epost_result = Entrez.read(
        entrez_retry(
            Entrez.epost,
            args,
            db="Protein",
            id=accessions_string,
        )
    )
    # retrieve the web environment and query key from the Entrez post
    epost_webenv = epost_result["WebEnv"]
    epost_query_key = epost_result["QueryKey"]

    # retrieve the protein sequences
    with entrez_retry(
        Entrez.efetch,
        args,
        db="Protein",
        query_key=epost_query_key,
        WebEnv=epost_webenv,
        rettype="fasta",
        retmode="text",
    ) as seq_handle:
        for record in SeqIO.parse(seq_handle, "fasta"):
            # retrieve the accession of the record
            temp_accession = record.id

            success = False   # will be true if finds protein accession

            temp_accession = temp_accession.split("|")

            for item in temp_accession:
                try:  # check if item is an accession number
                    re.match(
                        (
                            r"(\D{3}\d{5,7}\.\d+)|(\D\d(\D|\d){3}\d)|"
                            r"(\D\d(\D|\d){3}\d\D(\D|\d){2}\d)"
                        ),
                        item,
                    ).group()
                    temp_accession = item
                    success = True
                    break
                except AttributeError:  # raised if not an accession
                    continue

            if success is False:
                logger.error(
                    f"Could not retrieve accession from {record.id}, therefore, "
                    "protein sequence not added to the database,\n"
                    "because cannot retrieve the necessary CAZyme record"
                )
                continue

            # check the retrieve protein accession is in the list of retrieved accession
            if temp_accession not in accessions:
                logger.warning(
                    f"Retrieved the accession {temp_accession} from the record id={record.id}, "
                    "but this accession is not in the database.\n"
                    "Therefore, not adding this protein seqence to the local database"
                )
                continue

            # retrieve the GenBank record from the local data base to add the seq to
            genbank_record = session.query(Genbank).\
                filter(Genbank.genbank_accession == temp_accession).first()

            retrieved_sequence = str(record.seq)  # convert to a string becuase SQL expects a string
            genbank_record.sequence = retrieved_sequence
            genbank_record.seq_update_date = date_today
            session.commit()

            if args.fasta is not None:
                file_io.write_out_fasta(record, temp_accession, args)

            if args.blastdb is not None:
                file_io.write_fasta_for_db(record, args)

            # remove the accession from the list
            accessions.remove(temp_accession)

    if len(accessions) != 0:
        logger.warning(
            "Protein sequences were not retrieved for the following CAZymes in the local database"
        )
        for acc in accessions:
            logger.warning(f"GenBank accession: {acc}")

    return


def entrez_retry(entrez_func, args, *func_args, **func_kwargs):
    """Call to NCBI using Entrez.

    Maximum number of retries is 10, retry initated when network error encountered.

    :param logger: logger object
    :param retries: parser argument, maximum number of retries excepted if network error encountered
    :param entrez_func: function, call method to NCBI

    :param *func_args: tuple, arguments passed to Entrez function
    :param ** func_kwargs: dictionary, keyword arguments passed to Entrez function

    Returns record.
    """
    logger = logging.getLogger(__name__)
    record, retries, tries = None, args.retries, 0

    while record is None and tries < retries:
        try:
            record = entrez_func(*func_args, **func_kwargs)

        except IOError:
            # log retry attempt
            if tries < retries:
                logger.warning(
                    f"Network error encountered during try no.{tries}.\nRetrying in 10s",
                    exc_info=1,
                )
                time.sleep(10)
            tries += 1

    if record is None:
        logger.error(
            "Network error encountered too many times. Exiting attempt to call to NCBI"
        )
        return

    return record
