================================
Retrieving Structures from PDB
================================

The CAZy webscraper supports the automated retrieval of protein structure files from PDB for CAZymes 
contained within the local CAZyme database, created using ``cazy_webscraper``. ``cazy_webscraper`` 
retrieves the PDB accessions from CAZymes that match user specified criteria (include taxonomy, CAZy 
family, EC number etc). The retrived PDB accessions are then passed to the ` ``Biopython`` module 
``Bio.PDB`` <https://biopython.org/docs/latest/api/Bio.PDB.html>`_. ``Biopython.PDB`` then retrieves 
the structure files from PDB and writes them out the a specified directory.

For specific information of the ``Bio.PDB`` module please see the 
`Biopython documentation <https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ>`_.


.. warning::
    If many requests are going to be made in a series to PDB (for example a series of 100 
    requests), then it is expected practise to do this **weekend** or 
    **outside peak times**.


Intro and quick start
-----------------------------

``cazy_webscraper`` can retrieve PDB structure files from proteins catalogued within a local CAZyme 
database, created using ``cazy_webscraper``. Therefore, if you do not already have a local CAZyme database, 
you need to make one.

To retrieve the structure file for every PDB accession listed in your local CAZyme database (which is the 
default behaviour), the command has the following structure:

.. code-block:: bash
  cazy_webscraper_get_pdb_structures <path_to_your_local_cazyme_db.db> <files_types>

The path to the local CAZyme database needs to point specifically at the database file created using 
``cazy_webscraper``, which has the name format: ``cazy_scrape_YYYY-MM-DD--HH-MM-ss.db``. For example, 
``cazy_scrape_2021-06-10--22-44-49.db``.

*<file_types>* is a list of file types to download every structure in. To list multiple file types, separate 
them with a single comma (with no spaces). For example: ``pdb,xml``.

The accepted file formats are (*these are determined by ``Biopython.PDB``, and the following is taken from 
the ``Biopython.PDB`` `documentation <https://biopython.org/docs/1.75/api/Bio.PDB.PDBList.html>`_):
*	``mmCif`` (default, PDBx/mmCif file),
*	``pdb`` (format PDB),
*	``xml`` (PDBML/XML format),
*	``mmtf`` (highly compressed),
*	``bundle`` (PDB formatted archive for large structure}


Output directory
-------------------

The retrieved protein structures can be written to a specified directory by using the ``--output`` or `-o` flags, 
followed by the path to desired output directory. If the directory does not exist, ``cazy_webscraper`` 
will build the directory and all necessary parent directories.

The ``Biopython`` ``Bio.PDB`` module does not support writing out the protein structures to standard 
out. Therefore, if ``--output`` is not used, ``cazy_webscraper`` will write out the protein structures 
to the current working directory.


Customising the retrieval of structure files
--------------------------------------------------

You may not wish to retrieve the structure file for every PDB accession stored in your local CAZyme database. 
You may wish to only retrieve the structure files that are for CAZymes derived from different taxa, CAZy 
classes/(sub)families, and/or annotated with specified EC numbers.

The same taxonomy (kingdom, genera, species and strain), biochemical (EC number) and CAZy classification 
(CAZy (sub)family, and CAZy class) filters used for customising the scraping of CAZy, can be used for the 
retrieval of structure files for only CAZymes matching these critiera.

These filters are applied in combination, in the same way for customising the scraping of CAZy. For example, 
If a kingdom and EC number is specified, structure files will only be retrieved for CAZymes from species from 
the taxonomy kingdom **and** annotated with at least one of the EC numbers.

For looking at how these filters work in detail, check the documentation for customising the scraping of CAZy. At 
the end of this page are also some example commands for customised retrieval of structure files.


Command line arguments cheat sheet
----------------------------------------

Positional arguments
^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash
  database - Path to the local CAZyme database containing CAZymes to retrieve structure files for
  pdb - List of file types to retrieve each structure file in. Choice: "mmCif", "pdb", "xml", "mmtf", "bundle".
    Separate each file type with a single comma, with no spaces

Optional arguments
^^^^^^^^^^^^^^^^^^^^^^

Many of the flags are the same flags for customising the scraping of CAZy. For ``--classess``, ``ec``, 
``--families``, ``--kingdoms``, ``--genera``, ``--species`` and ``--strains``, the default is not to apply the filter. 
The filter is only applied if the flag is called.

``--batch_limit`` - Number of PDB accessions passed to Biopython.PDB at once, to retrieve from PDB. Default: 200

``--cazy_synonyms`` - Path to JSON file containing class synonoms to used instead of those built into cazy_webscraper

``-c``, ``--classes`` - List of CAZy classes to limit the retrieval of structure files to. Separate classes with a single comma

``--ec`` - List of EC numbers, separate EC numbers with a single comma. CAZymes with at least one of the listed EC numbers will 
be selected for retrieving structure files for

``-f``, ``--force`` - Force writing to output directory that already exists. Default: False - will not write out to an output directory 
that already exists.

``--families`` - List of CAZy families to limit the retrieval of structure files to. Separate families with a single comma

``--kingdoms`` - List of taxonomy kingdoms to limit the retrieval structure files for CAZymes that are derived from species from the specified 
kingdoms. Options: archaea, bacteria, eukaryota, viruses, unclassified (not case sensitive). 

``--genera`` - List of genera to limit the retrieval structure files for CAZymes that are derived from species from the specified 
genera.

``-l``, ``--log`` - Write out log to a file. Define a path to write out the log to. Default: None, does not write the log to a log file.

``-n``, ``--nodelete`` - Enable not deleting content in an output directory that already exists. Default: False – 
content in the output directory **is deleted**. **The output directory is nuked!**. If enabled/called 
the content already present in an existing output directory **is not deleted**.

``-o``, ``--output`` - Path to output directory. Default: write out structure files to the current working directory.

``--overwrite`` - Enable overwriting local structure files if already present. Default: False – do not overwrite a file if already present.

``--species`` - List of species to limit the retrieval of structure files for CAZymes that are derived from these species. List a species will identify 
CAZymes from *all strains* of each listed species.

``--strains`` - List of species *strains* to limit of structure files for CAZymes that are derived from these specific species strains.

``-v``, ``--verbose`` - Enable verbose logging. Default: False.



Configuration file
--------------------------


A configuration file (written in YAML) can also be used instead or in combination with command line flags, for customising 
the retrieval of structure files. The configruation file is structured the exact same way it is used for 
customising the scraping of CAZy.


Example commands
------------------------

In these examples the local CAZyme database is stored in the ``my_cazy_db`` directory.

**Example 1**
Retrieve structure files for all CAZymes from bacterial species from PL1, PL2 and PL3, and retrieve all 
structures in pdb and xml file formats.
.. code-block:: bash
    cazy_webscraper_get_pdb_structures my_cazy_db/cazy_scrape_2021-06-10--22-44-49.db pdb,xml --kingdoms bacteria --families PL1,PL2,PL3

**Example 3**
Retrieve structure files for all CAZymes from all Aspergillus species, from any GH or CBM family, and annotated with at least one of 
EC3.2.1.21 or EC3.2.1.23.
.. code-block:: bash
    cazy_webscraper_get_pdb_structures my_cazy_db/cazy_scrape_2021-06-10--22-44-49.db pdb --classes GH,CBM --genera Aspergillus --ec 3.2.1.21,3.2.1.23,3.2.1.*

The 'EC' prefix of EC numbers is not essential, therefore, the following comamnd is also accepted:
.. code-block:: bash
    cazy_webscraper_get_pdb_structures my_cazy_db/cazy_scrape_2021-06-10--22-44-49.db pdb --classes GH,CBM --genera Aspergillus --ec EC3.2.1.21,EC3.2.1.23,EC3.2.1.*

For missing digits, either a dash ('-') or an asterix ('*') can be used. Is using dashes, encapsulate the EC number list with quotation marks (double of single, as long as 
it is consistent).
