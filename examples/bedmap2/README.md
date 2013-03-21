Antarctic model using Bedmap2 data
=================

This is just a draft script for converting the Bedmap2 product into NetCDF.

See (http://www.antarctica.ac.uk//bas_research/our_research/az/bedmap2/index.php)

Cite this publication if you use this example:

* Fretwell, P., and others (2013) *Bedmap2: improved ice bed, surface and thickness datasets for Antarctica*.  The Cryosphere 7, 375--393, doi:10.5194/tc-7-375-2013. (http://www.the-cryosphere.net/7/375/2013/tc-7-375-2013.pdf)

Running the example
=================

First download and extract (https://secure.antarctica.ac.uk/data/bedmap2/bedmap2_bin.zip)

The script `readgeom.py` has two editable parameters, the string `BM2PATH` and
the flag `showgeom`.  Choose the values you like.

Do:

    $ ./readgeom.py

