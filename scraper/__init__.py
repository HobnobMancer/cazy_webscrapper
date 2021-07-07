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
"""Web scraper to scrape the CAZy website."""


from scraper.utilities import termcolour


VERSION_INFO = f"cazy_webscraper version: {__version__}"


CITATION_INFO = [
    termcolour(
        "If you use cazy_webscraper in your work, please cite the following publication:",
        "green",
    ),
    termcolour(
        "\tHobbs, E. E. M., Pritchard, L., Chapman, S., Gloster, T. M.,",
        "blue",
    ),
    termcolour(
        "\t(2021) cazy_webscraper Microbiology Society Annual Conference 2021 poster. ",
        "blue",
    ),
    termcolour(
        "\tFigShare. Poster.",
        "blue",
    ),
    termcolour(
        "\thttps://doi.org/10.6084/m9.figshare.14370860.v7",
        "blue",
    ),
]
