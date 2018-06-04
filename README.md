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

Run XLinkProphet on the ReACT analysis of the Control data set, consisting of 9 search result pepXML files and 9 react2.xls pairing files that link up the search results of the crosslink's 2 released peptides:

Download all Control Data Set files in the ReACT directory to a directory on your computer, along with XLinkProphet.pl.  You will need the Trans-Proteomic Pipeline (https://sourceforge.net/projects/sashimi/) on your computer as well.

In a terminal, go into the directory with your data and type the following three commands in succession:

1. Combine the 9 search result files into a single interact.pep.xml file, and run PeptideProphet to assign probabilities that each search result is correct:  <p style='color:red'>xinteract -OEAdP -p0 -PPM -l6 -drev_ 082117_BSA_1mMBDP_2hr_1.pep.xml 101117_a-casein_BDP_2hr_1.pep.xml 101117_a-lactalbumin_BDP_2hr_1.pep.xml 101117_ADH_BDP_2hr_1.pep.xml 101117_b-casein_BDP_2hr_1.pep.xml 101117_b-lactoglobulin_BDP_2hr_1.pep.xml 101117_cytochromeC_BDP_2hr_1.pep.xml 101117_histone_BDP_2hr_1.pep.xml 101117_myoglobin_BDP_2hr_1.pep.xml</p>

2. Run iProphet to further validate the search results with additional models, assigning revised probabilities that the search results are correct in output file iprophet.pep.xml:
InterProphetParser interact.pep.xml iprophet.pep.xml

3. Run XLinkProphet using the LOCAL_REACT setting to not use the raw data file locations in the search results to find the react2.xls pairing files, but rather find them in the current directory:
XLinkProphet.pl iprophet.pep.xml LOCAL_REACT

Your final output will be a pepXML file iprophet-xl.pep.xml and a tab delimited file iprophet-xl.xls.

