# xlinkprophet
Validator of XL-MS crosslinks obtained using cleavable crosslinkers

Download files XLinkProphet.pl and uniprot_conversion_table.txt, and then modify XLinkProphet.pl line 73 to include the 
full path location of uniprot_conversion_table.txt on your computer.  

Running XLinkProphet on a cleavable crosslink dataset requires:
1. Validated search results of the cleaved peptides, using PeptideProphet and optional subsequent iProphet.
2. A pairing file with the columns 'scan1', 'scan2', 'prec mass', and 'prec z' that contain for each crosslink in the dataset
the scan of the first released peptide, scan of the second released peptide, the crosslink precursor mass, and crosslink 
precursor charge.  This information is necessary for XLinkProphet to properly pair together the two search results corresponding to both released peptides of a crosslink.
3. ProteinProphet must be present on the computer.


The input file is the validated search results in pepXML format.  The output is validated crosslinks in Kojak pepXML format.
Note that the XLinkProphet probability is encoded as PeptideProphet in the output so that multiple output files can be
combined with iProphet.

Usage:   XLinkProphet.pl < PeptideProphet or iProphet pepXML file > (options)

Type:    XLinkProphet.pl to view all options
