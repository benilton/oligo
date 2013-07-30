THE `OLIGO` PACKAGE
=================

The `oligo` package provides tools to preprocess oligonucleotide
arrays. It is designed to support all Affymetrix and NimbleGen chips and
offers tools for reading in intensity files in their native format
(Affymetrix CEL files and NimbleGen XYS files).

Installation
------------

`oligo` is a BioConductor package. The recommended way to install it is
loading `R` and using `biocLite`, as shown below:

     source('http://www.bioconductor.org/biocLite.R')
     biocLite('oligo')

What does it offer?
-------------------

`oligo` offers a number of tools for preprocessing:

* Data Import: Affymetrix CEL and NimbleGen XYS files;
* Background correction methods;
* Normalization methods;
* Summarization methods;
* Visualization tools;
* Quality control tools via PLMs;

Contributing
------------

1. Fork it.
2. Create a branch (`git checkout -b my_contrib`)
3. Commit your changes (`git commit -am "My Contributions"`)
4. Push to the branch (`git push origin my_contrib`)
5. Create an [Issue][1] with a link to your branch


Contact
-------

Benilton Carvalho
<firstName.lastName> AT cancer <dot> org <dot> uk

[1]: http://github.com/benilton/oligo/issues
