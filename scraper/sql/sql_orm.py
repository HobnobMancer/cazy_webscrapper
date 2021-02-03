#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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
"""Submodule to build a local SQL database"""

import logging
import os
import re
import sys

import sqlite3

from sqlalchemy import (
    Boolean,
    Column,
    ForeignKey,
    Integer,
    PrimaryKeyConstraint,
    String,
    Table,
    create_engine,
    event,
    exc,
)
from sqlalchemy.engine import Engine
from sqlalchemy.ext.compiler import compiles
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, sessionmaker
from sqlalchemy.sql.expression import BinaryExpression, func, literal
from sqlalchemy.sql.operators import custom_op


# Use the declarative system
# Database structured in NF1
Base = declarative_base()
Session = sessionmaker()


# Enable regular expression searching of the database
class ReString(String):
    """Enchanced version of standard SQLAlchemy's :class:`String`.

    Supports additional operators that can be used while constructing filter expressions.
    """
    class comparator_factory(String.comparator_factory):
        """Contains implementation of :class:`String` operators related to regular expressions."""
        def regexp(self, other):
            return RegexMatchExpression(self.expr, literal(other), custom_op('~'))

        def iregexp(self, other):
            return RegexMatchExpression(self.expr, literal(other), custom_op('~*'))

        def not_regexp(self, other):
            return RegexMatchExpression(self.expr, literal(other), custom_op('!~'))

        def not_iregexp(self, other):
            return RegexMatchExpression(self.expr, literal(other), custom_op('!~*'))


class RegexMatchExpression(BinaryExpression):
    """Represents matching of a column againsts a regular expression."""

@compiles(RegexMatchExpression, 'sqlite')
def sqlite_regex_match(element, compiler, **kw):
    """Compile the SQL expression representing a regular expression match for the SQLite engine."""
    # determine the name of a custom SQLite function to use for the operator
    operator = element.operator.opstring
    try:
        func_name, _ = SQLITE_REGEX_FUNCTIONS[operator]
    except (KeyError, ValueError) as e:
        would_be_sql_string = ' '.join((compiler.process(element.left),
                                        operator,
                                        compiler.process(element.right)))
        raise exc.StatementError(
            "unknown regular expression match operator: %s" % operator,
            would_be_sql_string, None, e)

    # compile the expression as an invocation of the custom function
    regex_func = getattr(func, func_name)
    regex_func_call = regex_func(element.left, element.right)
    return compiler.process(regex_func_call)

@event.listens_for(Engine, 'connect')
def sqlite_engine_connect(dbapi_connection, connection_record):
    """Listener for the event of establishing connection to a SQLite database.

    Creates the functions handling regular expression operators
    within SQLite engine, pointing them to their Python implementations above.
    """
    if not isinstance(dbapi_connection, sqlite3.Connection):
        return

    for name, function in SQLITE_REGEX_FUNCTIONS.values():
        dbapi_connection.create_function(name, 2, function)

# Mapping from the regular expression matching operators
# to named Python functions that implement them for SQLite.
SQLITE_REGEX_FUNCTIONS = {
    '~': ('REGEXP',
          lambda value, regex: bool(re.match(regex, value))),
    '~*': ('IREGEXP',
           lambda value, regex: bool(re.match(regex, value, re.IGNORECASE))),
    '!~': ('NOT_REGEXP',
           lambda value, regex: not re.match(regex, value)),
    '!~*': ('NOT_IREGEXP',
            lambda value, regex: not re.match(regex, value, re.IGNORECASE)),
}


# define association/relationship tables


# linker table between cazymes and CAZy family and subfamilies
cazymes_families = Table(
    "cazymes_families",
    Base.metadata,
    Column("cazyme_id", Integer, ForeignKey("cazymes.cazyme_id")),
    Column("family_id", Integer, ForeignKey("families.family_id")),
    PrimaryKeyConstraint("cazyme_id", "family_id"),
)


# linker table between cazymes and ec numbers
cazymes_ecs = Table(
    "cazymes_ecs",
    Base.metadata,
    Column("cazyme_id", Integer, ForeignKey("cazymes.cazyme_id")),
    Column("ec_id", Integer, ForeignKey("ecs.ec_id")),
    PrimaryKeyConstraint("cazyme_id", "ec_id"),
)


# linker table between cazymes and UniProt accessions of CAZymes
cazymes_uniprots = Table(
    "cazymes_uniprots",
    Base.metadata,
    Column("cazyme_id", Integer, ForeignKey("cazymes.cazyme_id")),
    Column("uniprot_id", Integer, ForeignKey("uniprots.uniprot_id")),
    PrimaryKeyConstraint("cazyme_id", "uniprot_id"),
)


# linker table between CAZymes and PDB structures
cazymes_pdbs = Table(
    "cazymes_pdbs",
    Base.metadata,
    Column("cazyme_id", Integer, ForeignKey("cazymes.cazyme_id")),
    Column("pdb_id", Integer, ForeignKey("pdbs.pdb_id")),
    PrimaryKeyConstraint("cazyme_id", "pdb_id"),
)


# define models


class Cazyme(Base):
    """Describes a CAZyme, which is a protein single entry in CAZy.

    Every CAZyme will have a name, a source organism, at least one CAZy family, and at least
    a primary GenBank accession. A CAZyme may also have non-primary GenBank accessions, EC
    number annotations, UniProt accessions and PDB accessions.
    """
    __tablename__ = "cazymes"

    cazyme_id = Column(Integer, primary_key=True)
    cazyme_name = Column(String)
    taxonomy_id = Column(Integer, ForeignKey("taxs.taxonomy_id"))

    taxonomy = relationship("Taxonomy", back_populates="cazymes")

    families = relationship(
        "CazyFamily",
        secondary=cazymes_families,
        back_populates="cazymes",
        lazy="dynamic",
    )
    cazymes_genbanks = relationship(
        "Cazymes_Genbanks",
        back_populates="cazymes",
        lazy="dynamic",
    )

    # Not all CAZymes will have EC numbers, UniProt accessions or PDB accessions
    ecs = relationship(
        "EC",
        secondary=cazymes_ecs,
        back_populates="cazymes",
        lazy="dynamic",
    )
    uniprots = relationship(
        "Uniprot",
        secondary=cazymes_uniprots,
        back_populates="cazymes",
        lazy="dynamic",
    )
    pdbs = relationship(
        "Pdb",
        secondary=cazymes_pdbs,
        back_populates="cazymes",
        lazy="dynamic",
    )

    def __str__(self):
        return f"-CAZyme name={self.cazyme_name}, id={self.cazyme_id}-"

    def __repr__(self):
        return f"<Class Cazyme: name={self.cazyme_name}, id={self.cazyme_id}>"


class Taxonomy(Base):
    """Describes the source organism of CAZymes."""
    __tablename__ = "taxs"

    taxonomy_id = Column(Integer, primary_key=True)
    genus = Column(String)
    species = Column(String)

    cazymes = relationship("Cazyme", back_populates="taxonomy")

    def __str__(self):
        return f"-Source organism, Genus={self.genus}, Species={self.species}"

    def __repr__(self):
        return (
            f"<Class Taxonomy: genus={self.genus}, species={self.species}, id={self.taxonomy_id}>"
        )


class CazyFamily(Base):
    """Describes a CAZy family.

    Every unique CAZy family-subfamily pair will be given a unique family_id. For example, 
    if a CAZyme is catalogued under a subfamily, the parent CAZy family and the CAZy subfamily
    will be listed together, and given a single family_id. If another protein is catalogued
    under only the parent CAZy family, another entry with for the CAZy family will be made with
    a null value for the subfamily and a different family_id. """
    __tablename__ = "families"
    family_id = Column(Integer, primary_key=True)
    family = Column(ReString)
    subfamily = Column(String)

    cazymes = relationship(
        "Cazyme",
        secondary=cazymes_families,
        back_populates="families",
        lazy="dynamic",
    )

    def __str__(self):
        if self.subfamily is None:
            return f"-CAZy family {self.family}, id={self.family_id}-"
        else:
            return f"-CAZy subfamily {self.subfamily}, parent={self.family}, id={self.family_id}-"

    def __repr__(self):
        """Return string representation of source organism."""
        return(
            f"<Class Family, family={self.family}, subfamily={self.subfamily}, id={self.family_id}"
        )


class Genbank(Base):
    """Describe a GenBank accession number of protein sequences.

    The associated GenBank protein record is the source record from which CAZy retrieves the
    protein sequence for the CAZyme.
    """
    __tablename__ = "genbanks"

    genbank_id = Column(Integer, primary_key=True)
    genbank_accession = Column(String)
    sequence = Column(String)
    seq_update_date = Column(String)  # 'YYYY/MM/DD'

    cazymes_genbanks = relationship(
        "Cazymes_Genbanks",
        back_populates="genbanks",
        lazy="dynamic",
    )

    def __str__(self):
        return f"-Genbank accession={self.genbank_accession}-"

    def __repr__(self):
        return f"<Class GenBank acc={self.genbank_accession}>"


class Cazymes_Genbanks(Base):
    """Represent assoication between a CAZyme and its primary and non-primary GenBank accessions.

    The primary GenBank accession is the only hyperlinked
    GenBank accession for the protein, and believed to be used by CAZy to indicate the source
    GenBank protein record for the record in CAZy. It can not be guareenteed that a GenBank
    accession will only be recorded as a primary OR a non-primary accession. It may be possible
    that a GenBank accession is the primary accession for one CAZyme and a non-primary accession
    for another. This is believed to be possible becuase CAZy does not appear to ID unique proteins
    by the GenBank accession because duplicate entries for CAZyme can be found within CAZy.
    """
    __tablename__ = "cazymes_genbanks"

    link_id = Column(Integer, primary_key=True)  # unique ID of the CAZyme-GenBank relationship
    cazyme_id = Column(Integer, ForeignKey("cazymes.cazyme_id"))
    genbank_id = Column(Integer, ForeignKey("genbanks.genbank_id"))

    primary = Column(Boolean)  # indicate if primary or non-primary accession of the CAZyme

    cazymes = relationship("Cazyme", back_populates="cazymes_genbanks")
    genbanks = relationship("Genbank", back_populates="cazymes_genbanks")

    def __str__(self):
        return f"cazyme_id={self.cazyme_id}--genbank_id={self.genbank_id}-primary={self.primary}-"

    def __repr__(self):
        return(
            f"<Class Cazymes_GenBanks cazyme_id={self.cazyme_id}-"
            f"-genbank_id={self.genbank_id}-primary={self.primary}>"
        )


# Not all CAZymes will have EC numbers, UniProt accessions or PDB accessions


class EC(Base):
    """Describe EC numbers."""
    __tablename__ = "ecs"

    ec_id = Column(Integer, primary_key=True)
    ec_number = Column(String)

    cazymes = relationship("Cazyme", secondary=cazymes_ecs, back_populates="ecs", lazy="dynamic")

    def __str__(self):
        return f"-EC{self.ec_number}-ec_id={self.ec_number}-"

    def __repr__(self):
        return f"<Class EC, EC{self.ec_number}, ec_id={self.ec_number}>"


class Uniprot(Base):
    """Describe a UniProt accession number.

    The primary UniProt accession is the first UniProt accession that is lsited in UniProt for
    the CAZyme.
    """
    __tablename__ = "uniprots"

    uniprot_id = Column(Integer, primary_key=True)
    uniprot_accession = Column(String)
    primary = Column(Boolean)
    sequence = Column(String)

    cazymes = relationship(
        "Cazyme",
        secondary=cazymes_uniprots,
        back_populates="uniprots",
        lazy="dynamic",
    )

    def __str__(self):
        return(
            f"-UniProt accession={self.uniprot_accession}, "
            f"id={self.uniprot_id}, primary={self.primary}-"
        )

    def __repr__(self):
        return(
            f"<Class Uniprot accession={self.uniprot_accession}, "
            f"id={self.uniprot_id}, primary={self.primary}>"
        )


class Pdb(Base):
    """Describe a PDB accession number of protein structure.

    The primary PDB accession is the first PDB accession that is lsited in UniProt for
    the CAZyme.
    """
    __tablename__ = "pdbs"

    pdb_id = Column(Integer, primary_key=True)
    pdb_accession = Column(String)
    primary = Column(Boolean)

    cazymes = relationship("Cazyme", secondary=cazymes_pdbs, back_populates="pdbs", lazy="dynamic")

    def __str__(self):
        return(
            f"-PDB accession={self.pdb_accession}, "
            f"id={self.pdb_id}, primary={self.primary}-"
        )

    def __repr__(self):
        return(
            f"<Class Pdb accession={self.pdb_accession}, "
            f"id={self.pdb_id}, primary={self.primary}>"
        )


def build_db(time_stamp, args):
    """Build an empty SQL database and open a session.

    :param time_stamp: str, date and time stamp of when scrape was initated
    :param args: cmd args parser

    Return an open database session.
    """
    logger = logging.getLogger(__name__)
    logger.info("Building empty db to store data")

    # build database engine
    if args.database is None:
        if args.output is sys.stdout:
            # write to cwd, this is deleted in scrape is successful
            cwd = os.getcwd()
            db_path = cwd / f"cazy_scrape_temp_{time_stamp}.db"

        else:
            # write to specified output directory
            db_path = args.output / f"cazy_scrape_{time_stamp}.db"

    else:
        # user specificed an existing local CAZy SQL database
        db_path = args.database

    engine = create_engine(f"sqlite+pysqlite:///{db_path}", echo=False)
    Base.metadata.create_all(engine)
    Session.configure(bind=engine)

    return Session()


def get_db_session(args):
    """Create open session to local CAZy SQL database.

    :param args: cmd args parser

    Return an open database session.
    """
    logger = logging.getLogger(__name__)
    logger.info("Building empty db to store data")

    db_path = args.database

    engine = create_engine(f"sqlite+pysqlite:///{db_path}", echo=False)
    Base.metadata.create_all(engine)
    Session.configure(bind=engine)

    return Session()
