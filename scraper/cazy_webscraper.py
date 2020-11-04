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
"""
Web scraper to scrape CAZy website and retrieve all protein data.

:cmd_args --config: path to configruration file
:cmd args --data_split: [None, class, family] how data is to be split/separated into dataframes
:cmd_args --force: force overwriting content in exisiting output directory
:cmd_args --log: path to log file, enables writing out log messages to a log file
:cmd_args --nodelete: if true does not delete content in pre-existing output directory
:cmd_args --output: path to output directory
:cmd_args --retries: the maxmimum number of times to retry connection if network error raised
:cmd_args --subfamily: enable retrieval of subfamilies from CAZy
:cmd_args --verbose: change logger level from warning to info, verbose logging

:func main: coordinate scraping of CAZy database
:func get_cazy_class_urls: retrieve URLs to CAZy class summary pages from CAZy homepage
:func get_cazy_family_urls: retrieve URLs to families on CAZy class summary page
:func get_subfamily_links: retrieve URLs to subfamilies on CAZy class summary page
:func parse_family: build Family class object to represent CAZy family
:func parse_family_pages: retrieve all URLs to pages containing protein records for CAZy family
:func parse_proteins: retrieve protein records from protein table page
:func row_to_protein: parse the protein record to build a Protein class object
:func browser_decorator: decorator for get_page() to coordinate retrying failed connections
:func get_page: connect to webpage and retrieval page as BeautifulSoup4 object

:class Protein: A single protein from CAZy database
:class Family: A single family from CAZy containing proteins
"""

import logging
import re
import sys

from collections import defaultdict
from typing import List, Optional
from requests.exceptions import ConnectionError

import mechanicalsoup

from tqdm import tqdm

from scraper import utilities, file_io, parse


class Protein:
    """A single protein.

    Each protein has a name, source organism (source), and links to external databases. The links to
    external databases are stored in a dictionary, keyed by the external database name ('str') with
    'list' values becuase there may be multiple links per database.
    """

    def __init__(self, name, source, ec, links=None):
        self.name = name
        self.ec = ec
        self.source = source
        if links is None:
            self.links = defaultdict(list)
        else:
            self.links = links

    def __str__(self):
        """Create representative string of class object"""
        return f"{self.name} ({self.source}): links to {self.links.keys()}"

    def __repr__(self):
        """Create representative object"""
        return(
            (
                f"<Protein: {id(self)}: {self.name}, "
                f"({self.source}), {len(self.links)} to external databases>"
            )
        )

    def get_protein_dict(self):
        """Return a dictionary containing all the data of the protein."""
        protein_dict = {
            "Protein_name": [self.name],
            "EC#": [self.ec],
            "Source_organism": [self.source],
        }
        if type(self.links) is dict:
            for database in ["GenBank", "UniProt", "PDB/3D"]:
                try:
                    protein_dict[database] = self.links[database]
                except KeyError:
                    protein_dict[database] = []
        else:
            for database in ["GenBank", "UniProt", "PDB/3D"]:
                protein_dict[database] = []
        return protein_dict


class Family:
    """A single CAZy family."""

    members = set()  # holds Protein instances

    def __init__(self, name, cazy_class):
        self.name = name
        self.cazy_class = cazy_class

    def __str__(self):
        return f"CAZy family {self.name}: {len(self.members)} protein members"

    def __repr__(self):
        return f"<Family: {id(self)}: {self.name}, {len(self.members)} protein members"

    def get_proteins(self):
        """Return a list of all protein members of the CAZy family."""
        return self.members

    def get_family_name(self):
        """Return family name"""
        return self.name


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Set up parser, logger and coordinate overal scrapping of CAZy.

    The collected data can be stored as a singel dataframe containing (not split), split into
    separate dataframes by class or by family. Excluded classes are CAZy classes not specified in
    the configuration file and thus, will not be scraped. User_cazy_families is the list of CAZy
    families specified to be scraped in the configration file.
    """
    # Program preparation
    if argv is None:
        parser = utilities.build_parser()
        args = parser.parse_args()
    else:
        args = utilities.build_parser(argv).parse_args()

    if logger is None:
        logger = utilities.build_logger("cazy_webscraper", args)
    logger.info("Run initiated")

    if args.output is not sys.stdout:
        file_io.make_output_directory(args.output, logger, args.force, args.nodelete)

    # retrieve configuration data
    excluded_classes, config_dict, cazy_dict = file_io.parse_configuration(args, logger)

    # Crawl through and scrape CAZy website/database
    cazy_home = "http://www.cazy.org"  # the CAZy homepage URL

    # retrieve links to CAZy class pages
    class_pages = get_cazy_class_urls(cazy_home, excluded_classes, args, logger)

    all_data = []  # stores all Family class objects if not splitting the data

    # retrieve links to CAZy family pages
    for class_url in class_pages:

        # retrieve class name from url and convert to synonym used in configuration file
        class_name = class_url[20: -5]
        for key in cazy_dict:
            if class_name in cazy_dict[key]:
                class_name = key

        # retrieve URLs to families under current working CAZy class
        family_urls = get_cazy_family_urls(class_url, cazy_home, args, logger)

        families = []  # store Family class objects if splitting data be class

        if (config_dict is None) or (len(config_dict[class_name]) == 0):
            # no (sub)families specified. Scrape all families in CAZy class
            for family_url in family_urls:
                family_name = family_url[(len(cazy_home) + 1): -5]
                family = parse_family(family_url, family_name, cazy_home, args, logger)

                if args.data_split == "family":
                    parse.proteins_to_dataframe([family], args, logger)
                else:
                    families.append(family)

        else:
            # scrape only (sub)families specified in config file
            for family_url in family_urls:
                family_name = family_url[(len(cazy_home) + 1): -5]
                if family_name in config_dict[class_name]:
                    family = parse_family(family_url, family_name, cazy_home, args, logger)

                    if args.data_split == "family":
                        parse.proteins_to_dataframe([family], args, logger)
                    else:
                        families.append(family)

        if args.data_split == "class":
            # Write dataframe for CAZy class
            parse.proteins_to_dataframe(families, args, logger)
        else:
            all_data += families

    if args.data_split is None:
        # Write dataframe containing all data from CAZy
        parse.proteins_to_dataframe(all_data, args, logger)

    logger.info("Program finished")


def get_cazy_class_urls(cazy_home, excluded_classes, args, logger):
    """Returns a list of CAZy class main/home page URLs for each specified class as the CAZy site.

    :param cazy_url: str, URL to the CAZy home page.
    :param excluded_classes: list, list of CAZy classes not to be scraped
    :param args: cmd-line args parser
    :param logger: logger object

    Return list of URLs.
    """
    logger.info("Retrieving URLs to summary CAZy class pages")

    # define items to be excluded from returned class list, ALWAYS exlide links to genomes
    if excluded_classes is not None:
        exclusions = tuple(["<strong>Genomes</strong>"] + excluded_classes)
    else:
        exclusions = ("<strong>Genomes</strong>")

    # scrape the home page
    home_page = get_page(cazy_home, args.retries, logger)

    return [f"{cazy_home}/{_['href']}" for _ in home_page.soup.find_all("a", {"class": "spip_out"})
            if (not _["href"].startswith("http")) and (str(_.contents[0]) not in exclusions)]


def get_cazy_family_urls(class_url, cazy_home, args, logger):
    """Retrieve all protein members of each CAZy family within the given CAZy class.

    :param class_url: str, URL to the CAZy class
    :param cazy_home: str, URL to CAZy home page
    :param subfam_retrieval: bool, enables subfamily retrieval if true
    :param args: cmd-line args parser
    :param logger: logger object

    Returns list of URLs to family pages.
    """
    logger.info(f"Retrieving URLs to families at {class_url}")

    # scrape the class page
    class_page = get_page(class_url, args.retries, logger)

    # retrieve the <h3> element that titles the div section containing the tables of family links
    family_h3_element = [_ for _ in class_page.soup.find_all("h3", {"class": "spip"}) if str(_.contents[0]) == "Tables for Direct Access"][0]

    # retrieve all tables within the parent div section of the <h3> element
    tables = family_h3_element.parent.find_all("table")

    # tables[0] is the table containing links to CAZy families
    # tables[1] is the table containing the link to unclassified proteins

    family_urls = family_urls = [f"{cazy_home}/{_['href']}" for _ in tables[0].find_all("a")]
    family_urls.append(f"{cazy_home}/{tables[1].a['href']}")
    if args.subfamilies is True:
        family_urls.append(get_subfamily_links(family_h3_element, cazy_home, logger))

    return family_urls


def get_subfamily_links(family_h3_element, cazy_home, logger):
    """Retrieve URL links to CAZy subfamilies.

    :param family_h3_element: bs4.element.Tag, h3 element titling the page div
    :param cazy_home: str, URL to CAZy home_page
    :param logger: logger object

    Return list of URLs to subfamilies.
    """
    all_links = family_h3_element.find_all("a")

    pattern = re.compile(r"\D+?\d+?_\d+?\.html")

    urls = []  # empty list to store subfamily URLs

    for link in all_links:
        try:
            search_result = re.search(pattern, link["href"])
            urls.append(f"{cazy_home}/{search_result.group()}")
        except (KeyError, AttributeError) as error:
            # KeyError raised if link does not have ['href']
            # AttributeError error raised if search_result is None becuase not subfam link
            pass

    return urls


def parse_family(family_url, family_name, cazy_home, args, logger):
    """Returns a Family object with Protein members, scraped from CAZy.

    :param family_url: str, URL to CAZy family summary page
    :param family_name: str, name of CAZy family
    :param cazy_home: str, URL to CAZy home page
    :param args: cmd-line args parser
    :param logger: logger object

    Return Family object.
    """
    logger.info(f"Starting retrieval of proteins for {family_name}")
    # retrieve family name from URL
    family_name = family_url[(len(cazy_home) + 1): -5]

    # retrieve class from family name
    pattern = re.compile(r"\D+")
    search_result = re.match(pattern, family_name)
    cazy_class = search_result.group()

    family = Family(family_name, cazy_class)

    for protein in tqdm(
        parse_family_pages(family_url, cazy_home, args, logger), desc=f"Parsing {family_name}"
    ):
        family.members.add(protein)

    return family


def parse_family_pages(family_url, cazy_home, args, logger):
    """Retrieve all protein records for given CAZy family.
    Protein records are listed in a pagination method, with 1000 proteins per page.

    :param family_url: str, URL to CAZy family main page
    :param cazy_home: str, URL to CAZy home page
    :param args: cmd-line args parser
    :param logger: logger object

    Return list of protein records.
    """
    logger.info(f"Retrieving URLs to all pages containing proteins for {family_url}")
    # compile URL to first family page of protein records
    first_pagination_url = family_url.replace(".html", "_all.html")
    first_pagination_page = get_page(first_pagination_url, args.retries, logger)

    protein_page_urls = [first_pagination_url]

    # retrieve the URL to the final page of protein records in the pagination listing
    try:
        last_pagination_url = first_pagination_page.soup.find_all(
            "a", {"class": "lien_pagination", "rel": "nofollow"})[-1]
    except IndexError:  # there is no pagination; a single-query entry
        last_pagination_url = None

    if last_pagination_url is not None:
        url_prefix = last_pagination_url["href"].split("PRINC=")[0] + "PRINC="
        last_princ_no = int(last_pagination_url["href"].split("PRINC=")[-1].split("#pagination")[0])
        url_suffix = "#pagination" + last_pagination_url["href"].split("#pagination")[-1]

        # Build list of urls to all pages in the pagination listing, increasing the PRINC increment
        protein_page_urls.extend([f"{cazy_home}/{url_prefix}{_}{url_suffix}" for _ in
                                  range(1000, last_princ_no + 1000, 1000)])

    # Process all URLs into a single collection - a generator
    return (y for x in (parse_proteins(url, args, logger) for url in protein_page_urls) for y in x)


def parse_proteins(protein_page_url, args, logger):
    """Returns generator of Protein objects for all protein rows on a single CAZy family page.

    :param protein_page_url, str, URL to the CAZy family page containing protein records
    :param args: cmd-line args parser
    :param logger: logger object

    Return generator object.
    """
    logger.info(f"Retrieving proteins from {protein_page_url}")

    protein_page = get_page(protein_page_url, args.retries, logger)

    # retrieve protein record table
    protein_table = protein_page.soup.find_all("table", {"class": "listing"})[0]
    protein_rows = [_ for _ in protein_table.descendants if (_.name == "tr") and
                    ("id" not in _.attrs) and ("class" not in _.attrs)]

    for row in protein_rows:
        yield row_to_protein(row)


def row_to_protein(row):
    """Returns a Protein object representing a single protein row from a CAZy family protein page.

    Each row, in order, contains the protein name, EC number, source organism, GenBank ID(s),
    UniProt ID(s), and PDB accession(s).

    :param row: tr element from CAZy family protein page

    Return Protein instance.
    """
    # retrieve list of cells ('td' elements) in row
    tds = list(row.find_all("td"))

    protein_name = tds[0].contents[0].strip()
    source_organism = tds[2].a.get_text()
    links = {}
    try:
        ec_number = tds[1].contents[0].strip()
    except TypeError:  # raised if protein is not annotated with an EC number
        ec_number = "N/A"

    # test for len(tds[x].contents) in case there is no link,
    # the check of .name then ensures link is captured
    if len(tds[3].contents) and tds[3].contents[0].name == "a":
        links["GenBank"] = [f"{_.get_text()}: {_['href']}" for _ in tds[3].contents if
                            _.name == "a"]
    if len(tds[4].contents) and tds[4].contents[0].name == "a":
        links["UniProt"] = [f"{_.get_text()}: {_['href']}" for _ in tds[4].contents if
                            _.name == "a"]
    if len(tds[5].contents) and tds[5].contents[0].name == "a":
        links["PDB"] = [f"{_.get_text()}: {_['href']}" for _ in tds[5].contents if _.name == "a"]

    return Protein(protein_name, source_organism, ec_number, links)


def browser_decorator(func):
    """Decorator to retry the wrapped function up to 'retries' times."""

    def wrapper(*args, retries=10, **kwargs):
        tries, success = 0, False
        while not success and (tries < retries):
            try:
                response = func(*args, **kwargs)
            except ConnectionError:
                logger.warning(f"Failed to make connection to CAZy on try {tries} of {retries}")
                response = None
            if str(response) == "<Response [200]>":  # response was successful
                success = True
            # if response from webpage was not successful
            tries += 1
        if not success:
            logging.warning("Ran out of connection retries, and was unable to connect to CAZy")
            return None
        else:
            return response

    return wrapper


@browser_decorator
def get_page(url, retries, logger):
    """Create browser and use browser to retrieve page for given URL.

    :param url: str, url to webpage
    :param retries: int, number of times to retry connection
    :param logger: logger object

    Return browser response object (the page).
    """
    # create browser object
    browser = mechanicalsoup.Browser()

    # create response object
    page = browser.get(url)

    return page


if __name__ == "__main__":
    main()
