from argparse import Namespace
from pathlib import Path
import logging
import os
import time

import pandas as pd

from datetime import datetime
from typing import List, Optional

from Bio import Entrez, SeqIO
from tqdm import tqdm

from scraper import file_io
from scraper.sql.sql_orm import (
    Cazyme,
    CazyFamily,
    Cazymes_Genbanks,
    Genbank,
    get_db_session,
)
from scraper.utilities import build_logger, build_genbank_sequences_parser


db = Path("E:\cazy_webscraper\test_build_3\cazy_scrape_2021-01-27--10-19-04.db")
args = {
    "args": Namespace(
        database=db,
    )
}
logger=None
get_db_session(args["args"], logger)
cazy_class = "PL"
class_query = session.query(CazyFamily, Cazyme, Cazymes_Genbanks, Genbank).\
    join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
    join(Cazyme, (Cazyme.cazyme_id == Cazymes_Genbanks.cazyme_id)).\
    filter(CazyFamily.family.regexp(rf"{cazy_class}\d+")).\
    all()
print(len(class_query))
