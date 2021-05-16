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
"""Module for quering and retrieving data from NCBI.Entrez"""


import logging
import re
import time

from Bio import Entrez, SeqIO

from scraper.utilities import file_io


def get_sequences_for_dict(accessions, date_today, args):
    """Retrieve protein sequences from Entrez and write out to FASTA files.

    :param accessions: list, GenBank accessions of proteins
    :param date_today: str, YYYY/MM/DD
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
            temp_accession = record.id  # accession of the current working protein record
            print("temp accession: ", temp_accession)

            if temp_accession.find("|") != -1:  # sometimes multiple items are listed
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

            file_io.write_out_fasta(record, temp_accession, args)

            if args.blastdb is not None:
                # need all squences in a single FASTA file to create a BLASTdb
                file_io.write_fasta_for_db(record, args)

            # remove the accession from the list
            accessions.remove(temp_accession)

    if len(accessions) != 0:
        logger.warning("Protein sequences for the following CAZymes were not retrieved:")
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

