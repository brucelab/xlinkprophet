#!/usr/bin/perl

use strict;
use POSIX;

my $SELF = 'XLinkProphet_devxxx.pl';

# these learned from data to discriminate correct from incorrect cross-links
my %MODELS = (	'intra' => [0, 1], 
				'total_charge' => [2, 3, 4, 5, 6], 
				'massdiff_bin' => [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8], 
				'joint_score' => [0 .. 3],
				'nrx' => [0, 1, 2],
				'nsx' => [0 .. 1],
				'homopeptide' => [0, 1], 
				);
my %MODEL_BINS = ( 	'massdiff_bin' => 5,
					'deltacn_bin' => 0.2, 
				); # for those with adjusted values to create bins
my %MODEL_CONSTRAINTS = (
					'joint_score' => 'right',
				#	'nrx' => 'left'
				);

# these are modified peptide masses indicative of the crosslink stump
my @CROSSLINK_MODS = (); #({'n' => 198, 'K' => 325}); # light values for our crosslinker
my @CROSSLINK_STUMP_MODS = ({'n' => 198.040247, 'K' => 325.127385}); # light values for our crosslinker
my $MIN_PEPTIDE_LENGTH = 6; # exclude peps with length less than this, set to 0 to include all
my $ENCODE_IN_MANGO_RUNNERUP_SCORE = 1; # whether to write the -1 and -2 special codes for Mango runner ups as scores rather than probabilities in pepXML
my $NUM_MODPEP_DECIMALS = 2; # how many are stored in modified crosslink peptide name
my $MAX_ITERS = 20;
my $REPORTERMASS = 751.40508; # for BHP
my $DECOY_PREFIX = 'rev_';
my $SPEC_LABEL = "";
my $LIGHT_HEAVY = 0; # combined results, for react, penalize by parent scan and use full path react2xls to match with correct crosslink in dataset
my $SET_GENES_FROM_UNIPROT = 0; # whether to bother getting gene information to write into output PEPXML for use by XLinkDB upon data upload
my %XLINKDB = (); #("M.musculus" => 1);
my $XLINKDB_XL = {}; #getXLinkDBCrosslinks();
my $DISPLAY_INTER_NO_INTRA = 1; # report inter-link fraction rather than intra-link
my %NSX_VALS = ();
my %prot_crosslinks = (); # 2 proteins
my %protpair_instances = (); # accumulate total for each protpair
my $custom_reporter = 0; # when user specified
my $REMOVE_SILAC_LYS_ARG_MODS = 1; # whether to convert heavy modifications at R and K to their silac light masses
my $DECOY_FRACTION = 1;
my @SILAC_MODS = ({'K' => 136.11, 'R' => 162.12 } );

#$MODELS{'aa_score'} = [0, 1, 2, 3]; # set this

if(@ARGV == 0) {
	printf STDERR "\n";
	printf STDERR " usage:   $SELF < PeptideProphet or iProphet pepXML file > (options)\n\n";
	printf STDERR " options: MAX_FDR=xx (Filter output for results with computed FDR less than or equal to xx)\n";
	printf STDERR "          NO_NSX (Do not use NSX [number of sibling cross-links] distribution to compute probabilities)\n";
	printf STDERR "          OUTPUT_DECOYS (Output along with targets and include when computing FDRs)\n";
	printf STDERR "          TARGET_NEG (Do not use decoys to set negative distributions)\n";
	printf STDERR "          SPEC_CONTAINS=xxx (Filter input data and use only spectra matching xxx)\n";
	printf STDERR "          OUTPUT_SUFF=xxx (Add xxx onto typical output file names)\n";
	printf STDERR "          LOCAL_REACT (Use React2.xls files in current directory rather than those referenced in the inputfile)\n";
	printf STDERR "          PAIRFILE_SUFF=xxx (Suffix appended to search result pepXML file names for pairing files containing spectrum scan numbers of 2 peptides plus mass and charge of crosslink parent)\n";
	printf STDERR "          CROSSLINK_MODS=x:202.13343,y:329.23223,xxx (supplement default values of n:198.040247 and K:325.127385 with specified values for n/c termini or AA x,y,z...[e.g. `silac heavy K:333.141584,n:199.040247])\n";
	printf STDERR "          EMAIL=user\@org (send email at end of analysis with link to upload output to XLinkDB)\n";
	printf STDERR "          MIN_PEPLEN=xx (exclude peptides of length below xx [default: $MIN_PEPTIDE_LENGTH])\n";
	printf STDERR "          MAX_ITERS=xx (set maximum number of iterations for EM ccc [default: $MAX_ITERS])\n";
	printf STDERR "          REPORTERMASS=xxx (neutral mass of crosslink reporter [default: $REPORTERMASS for BHP])\n";
	printf STDERR "          DECOY_PREFIX=xxx (prefix of protein name indicating a decoy [default: $DECOY_PREFIX])\n";
	printf STDERR "          SPEC_LABEL=xxx (label spectrum tags so multiple instances of same spectrum are maintained)\n";
	printf STDERR "          LIGHT_HEAVY (label spectrum tags so multiple instances of same spectrum are maintained)\n";
	printf STDERR "          IQPIR_REPORTERMASSES=xxx:yyy,www:zzz,... (specify reporter masses yyy and zzz corresponding to lysine stump modification masses xxx and yyy, respectively, [with at least $NUM_MODPEP_DECIMALS decimal places, e.g. 325.13], e.g. IQPIR_REPORTERMASSES=325.13:811.455947,327.13:807.442528)\n";
	printf STDERR "          NON_UNIPROT (not all proteins in fasta database are Uniprot)\n";
	printf STDERR "          DECOY_FRACTION=xx (Accept only random xx fraction of all decoys [number between 0 and 1])\n";
	printf STDERR "          SILAC_MODS=K:xxx,R:yyy (Specify silac heavy mass difference xxx for lysine and yyy for arginine [default values %s:%s and %s:%s])\n", "K", $SILAC_MODS[0]->{'K'}, "R", $SILAC_MODS[0]->{'R'};
	printf STDERR "\n";
	printf STDERR " models:  %s\n", join(", ", sort {$a cmp $b} keys %MODELS);
	printf STDERR "\n";
	exit(1);
 
}

my $OUTPUT_SUFF = "";
my $CROSSLINK_PRIOR = 0.5;
my $OUTPUT_ADJUSTED = 1;
my $OUTPUT_DECOYS = 0;
my $TARGET_NEG = 0;
my $NO_NSX = 0;
my $NO_NSX = 0;
my $MAX_FDR = 1;
my $NO_DECOY_PRIOR_ADJ = 1;
my $NO_DECOY_NEGPRIOR_ADJ = 1;
my $MIN_PEP_LEN = 0; # one peptide must be at least this long
my $SPEC_TEXT = "";
my @NSX_BIN_BOUNDS = ();
my $LOCAL_REACT = 0;
my %LOCAL_REACT_BASES = ();
my $CURRENT_NSP = 0;
my %NSP_VALS = ();
my $INIT_UPDATE = 1;
my $USE_STRIPPED_PEP_NRX = 0;
my %PEP_PAIRS = (); # spectrum to peppair for nrx
my %CROSSLINKS = (); # peptide pair to hash of spectra
my %CROSSLINK_FDRs = (); # peptide pair to product of 1 - p of each instance
my $COMPUTE_COMPOSITE_FDRs = exists $MODELS{'nrx'};
my $PROTPROPH_INFO = {}; # read in proteins associated with each peptide
my $USE_PROTPROPH_PROTS = 1; # whether to read proteinprophet protein information and prioritize 
my $MIN_XL_PROB_4_PROTPROPH = 0; # only substitute with protein prophet proteins if XLinkProphet probability at least this great
my %GENE_NAMES = ();
my $UNIPROT_CONVERSION_FILE = "uniprot_conversion_table.txt"; # full path to uniprot_conversion_table.txt file
my %PEP_MAXPROBS = (); # holds the max xlink prob with each peptide, so know whether to bother pursuing GENENAMES
my $PAIRFILE_SUFF = ""; # used only for universal acceptance
my $PAIRING = "";
my %MANGO_SPECS = (); # spectrum_labels to make sure iProphet can accommodate multiple results for the same spectrum
my $SET_DECOY_PROBS_2_ZERO = 0; # whether to set decoy probs to 0 when computing distributions (doesn't have much effect on control data sample)
my $EMAIL_ADDRESS = "";
my @CROSSLINK_MODMASSES = (); # to give these to XLinkDB upload page....
my %ORGANISMS = (); # to give these to XLinkDB upload page....
my $QUERY_UNIPROT_GENENAMES = 1; # whether to to to uniprot web page to get proper gene names for insertion into output pepXML
my %IQPIR_REPORTERMASSES = (); # hashed from stump mass to corresponding reporter mass, likely used with LIGHT_HEAVY option IQPIR_REPORTERMASSES=325.13:811.455947,327.13:807.442528
my %FULL_STUMP_MASSES = ();
my $DIR_DIVISOR_REGX = $^O =~ /Win/ ? "\\\\" : "\/";
my $DIR_DIVISOR = $^O =~ /Win/ ? "\\" : "/";

my $file = $ARGV[0];
if($file !~ /^$DIR_DIVISOR_REGX/) {
	$file =  getcwd() . $DIR_DIVISOR . $file; # must be full length
}
for(my $k = 1; $k < @ARGV; $k++) {
	if($ARGV[$k] =~ /^MAX\_FDR\=(\S+)/) {
		$MAX_FDR = $1;
		printf STDERR "Setting max fdr to %s\n", $MAX_FDR;
	}
	elsif($ARGV[$k] eq 'NO_NSX') {
		$NO_NSX = 1;
		printf STDERR "Setting to NO_NSX to 1 (Not using NSX to compute probabilities)\n";
	}
	elsif($ARGV[$k] eq 'OUTPUT_DECOYS') {
		$OUTPUT_DECOYS = 1;
		printf STDERR "Setting to OUTPUT_DECOYS to 1\n";
	}
	elsif($ARGV[$k] eq 'TARGET_NEG') {
		$TARGET_NEG = 1;
		printf STDERR "Setting to TARGET_NEG to 1\n";
	}
	elsif($ARGV[$k] =~ /^SPEC\_CONTAINS\=(\S+)/) {
		$SPEC_TEXT = $1;
		printf STDERR "Setting spectrum filter to %s\n", $SPEC_TEXT;
	}
	elsif($ARGV[$k] =~ /^OUTPUT\_SUFF\=(\S+)/) {
		$OUTPUT_SUFF = $1;
		$OUTPUT_SUFF = "-" . $OUTPUT_SUFF if($OUTPUT_SUFF !~ /^\-/);
		printf STDERR "Setting output suffix to %s\n", $OUTPUT_SUFF;
	}
	# tab delimited text file with header including columns 'scan1','scan2','prec mass','prec z','prec scan' with individual scan numbers for peptides 1 and 2 in run, crosslink precursor mass charge
	elsif($ARGV[$k] =~ /^PAIRFILE\_SUFF\=(\S+)/) {
		$PAIRFILE_SUFF = $1;
		printf STDERR "Setting pairfile suffix to %s\n", $PAIRFILE_SUFF;
		$PAIRING = "Universal";
	}
	elsif($ARGV[$k] =~ /^CROSSLINK\_MODS\=(\S+)/) {
		my @next = split(",", $1);
		for(my $j = 0; $j < @next; $j++) {
			my $next_crosslink_mods = {};
			my @next_mods = split(":", $next[$j]);
			if(@next_mods != 2) {
				printf STDERR "Error: format of CROSSLINK_MODS $1 is not correct\n";
				sendEmail("Error: format of CROSSLINK_MODS $1 is not correct.", "The required format is the modified amino acid single letter abbreviation followed by colon and the modified amino acid mass in the pepXML files\n.  Multiple modifications are allowed, separated by commas.  For example: 'CROSSLINK_MODS=n:198.040247,K:325.127385'.");
				exit(1);
			}
			$next_crosslink_mods->{$next_mods[0]} = $next_mods[1];
			printf "Setting mod for $next_mods[0] to $next_mods[1]...\n";
			if(scalar keys %$next_crosslink_mods > 0) {
				push(@CROSSLINK_STUMP_MODS, $next_crosslink_mods);
				printf STDERR "Added additional crosslink mods %s\n", $1;
				for(my $j = 0; $j < @CROSSLINK_STUMP_MODS; $j++) {
					my @next = keys %{$CROSSLINK_STUMP_MODS[$j]};
					for(my $i = 0; $i < @next; $i++) {
						printf "number %d: %s => %s\n", $j, $next[$i], $CROSSLINK_STUMP_MODS[$j]->{$next[$i]};
					}
				}
			}
			else {
				die "Error, have no CROSSLINK_MODS $1 specified\n";
			}
		}
	}
	elsif($ARGV[$k] eq 'LOCAL_REACT') {
		$LOCAL_REACT = 1;
		printf STDERR "Setting to LOCAL_REACT to 1\n";
	}
	elsif($ARGV[$k] =~ /^EMAIL\=(\S+)/) {
		if($DIR_DIVISOR eq "/") {
			$EMAIL_ADDRESS = $1;
			printf STDERR "Setting email address to %s\n", $EMAIL_ADDRESS;
		}
	}
	elsif($ARGV[$k] =~ /^MIN\_PEPLEN\=(\d+)/) {
		$MIN_PEPTIDE_LENGTH = $1;
		printf STDERR "Setting MIN_PEPTIDE_LENGTH to %d\n", $MIN_PEPTIDE_LENGTH;
	}
	elsif($ARGV[$k] =~ /^MAX\_ITERS\=(\d+)/) {
		$MAX_ITERS = $1;
		printf STDERR "Setting MAX_ITERS to %d\n", $MAX_ITERS;
	}
	elsif($ARGV[$k] =~ /^REPORTERMASS\=(\S+)/) {
		$REPORTERMASS = $1;
		$custom_reporter = 1;
		printf STDERR "Setting REPORTERMASS to %s\n", $REPORTERMASS; #exit(1);
	}
	elsif($ARGV[$k] =~ /^DECOY_PREFIX\=(\S+)/) {
		$DECOY_PREFIX = $1;
		printf STDERR "Setting DECOY_PREFIX to %s\n", $DECOY_PREFIX; #exit(1);
	}
	elsif($ARGV[$k] =~ /^SPEC\_LABEL\=(\S+)/) {
		$SPEC_LABEL = $1;
		printf STDERR "Setting SPEC_LABEL to %s\n", $SPEC_LABEL; #exit(1);
	}
	elsif($ARGV[$k] eq 'LIGHT_HEAVY') {
		$LIGHT_HEAVY = 1;
		printf STDERR "Setting to LIGHT_HEAVY to 1\n";
	}
	elsif($ARGV[$k] eq 'NON_UNIPROT') {
		$SET_GENES_FROM_UNIPROT = 0;
		printf STDERR "Setting to SET_GENES_FROM_UNIPROT to 0\n";
	}
	elsif($ARGV[$k] =~ /^IQPIR\_REPORTERMASSES\=(\S+)/) {
		my @next = split(",", $1);
		for(my $j = 0; $j < @next; $j++) {
			my @div = split(":", $next[$j]);
			$div[0] = sprintf "%0.".$NUM_MODPEP_DECIMALS."f", $div[0];
			$IQPIR_REPORTERMASSES{$div[0]} = $div[1];
			printf "Setting reporter mass for $div[0] to $div[1]...\n";
		}
	}
	elsif($ARGV[$k] =~ /^DECOY\_FRACTION\=(\S+)/) {
		$DECOY_FRACTION = $1;
		if($DECOY_FRACTION < 0 || $DECOY_FRACTION > 1) {
			printf "Error: decoy fraction must be between 0 and 1, not $DECOY_FRACTION\n";
			exit(1);
		}
		printf "Setting decoy fraction to $DECOY_FRACTION ...\n"; #exit(1);
	}
	elsif($ARGV[$k] =~ /^SILAC\_MODS\=(\S+)/) {
 		my @mods = split(",", $1);
 		if(@mods != 2) {
			printf STDERR "Error: you must specify heavy SILAC masses for both K and R\n";
			exit(1);
  		}
 		for(my $m = 0; $m < @mods; $m++) {
 			my @next = split(":", $mods[$m]);
 			if($next[0] != 'K' && $next[0] != 'R') {
 				printf STDERR "Error: you must specify heavy SILAC masses for both K and R\n";
 				exit(1);
 			}
 			$SILAC_MODS[0]->{$next[0]} = $next[1];
			printf "Setting SILAC mass for $next[0] to $next[1] ...\n"; #exit(1);
 		}
	}
	else {
		die "Error: $ARGV[$k] is not valid option\n";
	}
}
$DECOY_PREFIX =~ s/\_/\\\_/g;

my $REMOVE_STATIC_MODS = 0; #getRemoveStaticMods($file); # look for mod amino acid with variable mod


# now get the crosslinkstump mods without decimal places that appear in the search results modified peptide sequence
my %seen = ();
for(my $z = 0; $z < @CROSSLINK_STUMP_MODS; $z++) {
	my %next = ();
	my @nextstumps = keys %{$CROSSLINK_STUMP_MODS[$z]};
	my $next = "";
	my $total = "";
	for(my $k = 0; $k < @nextstumps; $k++) {
		$next = $nextstumps[$k] . ":" . sprintf("%0.2f", $CROSSLINK_STUMP_MODS[$z]{$nextstumps[$k]}); 
		if(exists $seen{$next}) {
			@nextstumps = @nextstumps[0 .. $k-1, $k+1 .. $#nextstumps];
			$k--;
			next;
		}
		$seen{$next}++;
		$total .= $nextstumps[$k] . ":" . $CROSSLINK_STUMP_MODS[$z]{$nextstumps[$k]} . " ";
	}
	for(my $j = 0; $j < @nextstumps; $j++) {
		$next{$nextstumps[$j]} = sprintf "%0.".$NUM_MODPEP_DECIMALS."f", $CROSSLINK_STUMP_MODS[$z]->{$nextstumps[$j]};
	}
	push(@CROSSLINK_MODS, \%next);
}
my @stump_modmasses = ();
for(my $z = 0; $z < @CROSSLINK_MODS; $z++) {
	my @nextstumps = keys %{$CROSSLINK_MODS[$z]};
	for(my $j = 0; $j < @nextstumps; $j++) {
		push(@stump_modmasses, $nextstumps[$j] . ":" . $CROSSLINK_MODS[$z]->{$nextstumps[$j]});
	}
}

my $outfile = $file =~ /^(\S+)\.pep.xml$/ ? $1 . $OUTPUT_SUFF. $SPEC_TEXT. '-xl.xls' : die "Error: $file is not the required .pep.xml type\n"; #$ARGV[1];

my $COMPUTE_COMPOSITE_FDRs = exists $MODELS{'nrx'}; # compute 

my $REACT_UNFILT_PRIOR = '';
my $REACT_FILT_PRIOR = '';
my $REACT_USE_ADJ_PRIOR = 0;

# determine what type of pairing used to associate spectra with common crosslink precursor ion
if($PAIRING eq '') {
	if(isReact($file)) {
		$PAIRING = "ReACT";
	}
	else {
		$PAIRING = "Mango";
	}
}
printf STDERR "Pairing Algorithm: %s\n", $PAIRING;

my $OUTPUT_PEPXML = $outfile =~ /^(\S+)\.[t,x]/ ? $1 . '.pep' : die "Error parsing $outfile\n"; # .xml or .xls


parseResultsForDistributions($file, $outfile); 



sub getPriorAdjustedProbability {
(my $prob, my $prior, my $adjusted_prior) = @_;
my $adjusted_prob = ($prob * (1 - $prior) * $adjusted_prior) / (($adjusted_prior * ($prob - $prior)) + ($prior * (1 - $prob)));

}

sub getPeptideLength {
(my $modpep) = @_;
return length stripPeptide($modpep);
}

# remove modifications
sub stripPeptide {
(my $pep) = @_;
my $copy = $pep;
while($copy =~ /^([^[]+)\[.*?\](\S*)$/) {
	$copy = $1 . $2;
}
return $copy;

}

# reads precursor mass and charge as well as partner scans in each run
# general usage, specify suffix following search results pep.xml name...
sub readCombinedPairfileTexts {
(my $fileptr, my $pairfilesuff) = @_;
my %scan_info = (); # hashed by spectrum to peptide1 and peptide2 results
my $tot = 0;
for(my $k = 0; $k < @{$fileptr}; $k++) {
	my $filename = '';
	my $spectrum_base = '';
	if($fileptr->[$k] =~ /$DIR_DIVISOR_REGX([^$DIR_DIVISOR_REGX]+)$/) {
		$filename = $fileptr->[$k];
		if($filename =~ /^(\S+)\.pep\.xml\.$pairfilesuff$/) {
			$spectrum_base = $1;
		}
		else {
			die "1. here with $filename\n";
		}
	}
	else {
		die "Cannot parse filename from $fileptr->[$k]\n";
	}

	my %headers = ('scan1' => -1, 'scan2' => -1, 'prec mass' => -1, 'prec z' => -1); 
	my @headers = keys %headers;
	my $first = 1;
	open(FILE, $fileptr->[$k]) or die "cannot read $fileptr->[$k] $!\n";
	while(<FILE>) {
		chomp();
		my @parsed = split("\t");
		if($first) {
			$first = 0;
			for(my $k = 0; $k < @parsed; $k++) {
				$headers{$parsed[$k]} = $k if(exists $headers{$parsed[$k]});
			}
			for(my $k = 0; $k < @headers; $k++) {
				die "Could not find $headers[$k] in $_\n" if($headers{$headers[$k]} == -1);
			}
		}
		else {
			$scan_info{$spectrum_base} = {} if(! exists $scan_info{$spectrum_base});
			my $scan1 = $parsed[$headers{'scan1'}];
			while(length $scan1 < 5) {
				$scan1 = '0' . $scan1;
			}
			my $scan2 = $parsed[$headers{'scan2'}];
			while(length $scan2 < 5) {
				$scan2 = '0' . $scan2;
			}
			$scan_info{$spectrum_base}->{$scan1} = {};
			$scan_info{$spectrum_base}->{$scan1}->{'charge'} = $parsed[$headers{'prec z'}];
			$scan_info{$spectrum_base}->{$scan1}->{'mass'} = $parsed[$headers{'prec mass'}];
			$scan_info{$spectrum_base}->{$scan2}->{'partner'} = $spectrum_base . '_000_' . $scan1 . '.' . $scan1;
			$tot++;
		}
	}
	close(FILE);
}
printf STDERR "Stored %d base files with a total of %d values from %s\n", scalar keys %scan_info, $tot, join(",", @{$fileptr}); 
return \%scan_info;
}


# reads precursor mass and charge as well as partner scans in each run from react2.xls files
sub readCombinedReact2Texts {
(my $fileptr) = @_;
my %scan_info = (); # hashed by spectrum to peptide1 and peptide2 results
my $tot = 0;
for(my $k = 0; $k < @{$fileptr}; $k++) {
	my $filename = '';
	my $spectrum_base = '';
	if($fileptr->[$k] =~ /$DIR_DIVISOR_REGX([^$DIR_DIVISOR_REGX]+)$/) {
		$filename = $fileptr->[$k];
		if($filename =~ /^(\S+)\.pep\.xml\.react2\.xls$/) {
			$spectrum_base = $1;
		}
		else {
			die "2. here with $filename\n";
		}
	}
	else {
		die "Cannot parse filename from $fileptr->[$k]\n";
	}

	my %headers = ('ms3 scan1' => -1, 'ms3 scan2' => -1, 'ms2 prec mass' => -1, 'ms2 prec z' => -1, 'ms2 scan' => -1); 
	my @headers = keys %headers;
	my $first = 1;
	open(FILE, $fileptr->[$k]) or die "cannot read $fileptr->[$k] $!\n";
	while(<FILE>) {
		chomp();
		my @parsed = split("\t");
		if($first) {
			$first = 0;
			for(my $k = 0; $k < @parsed; $k++) {
				$headers{$parsed[$k]} = $k if(exists $headers{$parsed[$k]});
			}
			for(my $k = 0; $k < @headers; $k++) {
				die "Could not find $headers[$k] in $_\n" if($headers{$headers[$k]} == -1);
			}
		}
		else {
			$scan_info{$spectrum_base} = {} if(! exists $scan_info{$spectrum_base});
			my $scan1 = $parsed[$headers{'ms3 scan1'}];
			while(length $scan1 < 5) {
				$scan1 = '0' . $scan1;
			}
			my $scan2 = $parsed[$headers{'ms3 scan2'}];
			while(length $scan2 < 5) {
				$scan2 = '0' . $scan2;
			}
			$scan_info{$spectrum_base}->{$scan1} = {};
			$scan_info{$spectrum_base}->{$scan1}->{'charge'} = $parsed[$headers{'ms2 prec z'}];
			$scan_info{$spectrum_base}->{$scan1}->{'mass'} = $parsed[$headers{'ms2 prec mass'}];
			$scan_info{$spectrum_base}->{$scan2}->{'partner'} = $spectrum_base . '_000_' . $scan1 . '.' . $scan1;
			$scan_info{$spectrum_base}->{$scan1}->{'parent_scan'} = $parsed[$headers{'ms2 scan'}]; 
			$tot++;
		}
	}
	close(FILE);
}
printf STDERR "Stored %d base files with a total of %d values from %s\n", scalar keys %scan_info, $tot, join(",", @{$fileptr}); 
return \%scan_info;
}


# for precursor mass and charge based on spectrum scan for mango results
sub readPeaksText {
(my $file, my $compute_pepmassdiff) = @_;
my %headers = ('scan' => -1, 'intact_mass' => -1, 'intact_charge' => -1, 'pep1_mass' => -1, 'pep2_mass' => -1);
my @headers = keys %headers;
my $first = 1;
my %scan_info = (); # hashed by spectrum to peptide1 and peptide2 results
open(FILE, $file) or die "cannot read $file $!\n";
while(<FILE>) {
	chomp();
	my @parsed = split("\t");
	if($first) {
		$first = 0;
		for(my $k = 0; $k < @parsed; $k++) {
			$headers{$parsed[$k]} = $k if(exists $headers{$parsed[$k]});
		}
		for(my $k = 0; $k < @headers; $k++) {
			die "Could not find $headers[$k] in $_\n" if($headers{$headers[$k]} == -1);
		}
	}
	elsif(! exists $scan_info{$file . $parsed[$headers{'scan'}]}) {
	
		$scan_info{$parsed[$headers{'scan'}]} = {};
		$scan_info{$parsed[$headers{'scan'}]}->{'charge'} = $parsed[$headers{'intact_charge'}];
		$scan_info{$parsed[$headers{'scan'}]}->{'mass'} = $parsed[$headers{'intact_mass'}];
		# this is deprecated
		if($compute_pepmassdiff) {
			$scan_info{$parsed[$headers{'scan'}]}->{'pepmass_diff'} = [$parsed[$headers{'intact_mass'}] - $parsed[$headers{'pep1_mass'}] - $parsed[$headers{'pep2_mass'}] - $REPORTERMASS];
		}
	}
	# this is deprecated
	elsif($compute_pepmassdiff) {
		push(@{$scan_info{$parsed[$headers{'scan'}]}->{'pepmass_diff'}}, $parsed[$headers{'intact_mass'}] - $parsed[$headers{'pep1_mass'}] - $parsed[$headers{'pep2_mass'}] - $REPORTERMASS);
	}
}
printf STDERR "Stored %d values from %s\n", scalar keys %scan_info, $file;
return \%scan_info;
}

# util
sub stats {
(my $dataptr) = @_;
if(@{$dataptr} == 0) {
        return ("NA", "NA");
}
my $mean = 0;
my $meansq = 0;
foreach(@{$dataptr}) {
        $mean += $_;
        $meansq += $_ * $_;
}
$mean /= @{$dataptr};
$meansq /= @{$dataptr};
$meansq -= $mean * $mean;

$meansq = 0 if($meansq < 0);
$meansq = sqrt($meansq);
return($mean, $meansq, scalar @{$dataptr});
}

# weighted
sub weighted_stats {
(my $dataptr, my $wtptr) = @_;
if(@{$dataptr} == 0) {
        return ("NA", "NA");
}
die "problem with ", scalar @{$dataptr}, " data and ", scalar @{$wtptr}, " weights\n" if(@{$dataptr} != @{$wtptr});
my $mean = 0;
my $meansq = 0;
my $wt = 0;
for(my $k = 0; $k < @{$dataptr}; $k++) {
        $mean += $dataptr->[$k] * $wtptr->[$k];
        $meansq += $dataptr->[$k] * $dataptr->[$k] * $wtptr->[$k];
        $wt += $wtptr->[$k];
}
$mean /= $wt;
$meansq /= $wt;
$meansq -= $mean * $mean;

$meansq = 0 if($meansq < 0);
$meansq = sqrt($meansq);
return($mean, $meansq, (sprintf "%0.0f", $wt));
}

# here $data is spectrum -> relevant score...need one for each score
# do all cross-link distributions, exempting the original peptideprophet score
sub updateDistribution {
(my $posdist, my $negdist, my $data, my $probptr, my $maxdiff, my $prior, my $decoyptr, my $use_decoy_set_neg, my $ordered) = @_;
my %nextpos = ();
my %nextneg = ();
my $totprob = 0;
my $totneg = 0;
my $output = 0;
printf STDERR "setting decoys to neg (%d data total)\n", scalar keys %{$data} if($use_decoy_set_neg->[0] && ! $use_decoy_set_neg->[1]);
foreach(keys %{$data}) {
	next if($data->{$_} eq '*');

	my $nextprob = $probptr->{$_} < 0 ? 0 : $probptr->{$_};

	$nextpos{$data->{$_}} += $nextprob;
	$totprob += $nextprob;
	if($use_decoy_set_neg->[0]) {
		if(! $use_decoy_set_neg->[1] && $decoyptr->{$_}) {
			$nextneg{$data->{$_}}++;
			$totneg++;
		}
	}
	else {
		$nextneg{$data->{$_}} += 1 - $nextprob;
	}
}
foreach(keys %{$posdist}) {
	if($nextpos{$_} < $totprob * $prior) {
		$totprob += ($totprob * $prior - $nextpos{$_});
		$nextpos{$_} += $totprob * $prior - $nextpos{$_};
	}
}
if($totprob == 0) {
	printf "Warning: have totprob of 0\n";
	sendEmail("Warning: No correct crosslinks found.", "Please check your run parameters and number of correct search results in your input files.");
	exit(1);
}
if($use_decoy_set_neg->[0] && ! $use_decoy_set_neg->[1]) {
	foreach(keys %{$negdist}) {
	 	if($nextneg{$_} < $totneg * $prior) {
			$totneg += ($totneg * $prior - $nextneg{$_});
			$nextneg{$_} += $totneg * $prior - $nextneg{$_};
		}
	}
}

# now normalize and compare with previous
foreach(keys %nextpos) {
	$nextpos{$_} /= $totprob;
	$nextpos{$_} = $prior if($nextpos{$_} < $prior);
	$output = 1 if(abs($nextpos{$_} - $posdist->{$_}) > $maxdiff);
}
my $num_data = scalar keys %{$data};
foreach(keys %nextneg) {
	if($use_decoy_set_neg->[0]) {
		if(! $use_decoy_set_neg->[1]) {
			$nextneg{$_} /= $totneg;
			$nextneg{$_} = $prior if($nextneg{$_} < $prior);
		}
	}
	else {
		$nextneg{$_} /= ($num_data - $totprob);
		$nextneg{$_} = $prior if($nextneg{$_} < $prior);
		$output = 1 if(abs($nextneg{$_} - $negdist->{$_}) > $maxdiff);
	}
}
if($output) {
	foreach(keys %{$posdist}) {
		$posdist->{$_} = $nextpos{$_};
	}
	foreach(keys %{$negdist}) {
		if(! $use_decoy_set_neg->[0] || ! $use_decoy_set_neg->[1]) {
			$negdist->{$_} = $nextneg{$_};
		}
	}
}
# distribution should increase as the pos values increase
if($ordered eq 'right') {
	my @pos_values = sort {$a <=> $b} keys %{$posdist};
	my $tot = $posdist->{$pos_values[$#pos_values]};
	for(my $k = $#pos_values - 1; $k >= 0; $k--) {
		if($posdist->{$pos_values[$k]} >= $posdist->{$pos_values[$k+1]}) {
			$posdist->{$pos_values[$k]} = $posdist->{$pos_values[$k+1]};
		}
		$tot += $posdist->{$pos_values[$k]};
	}
	if($tot > 0 && $tot < 1) {
		foreach(keys %{$posdist}) {
			$posdist->{$_} /= $tot;
		}
	}
}
# distribution should increase as the pos values decrease
elsif($ordered eq 'left') {
	my @pos_values = sort {$a <=> $b} keys %{$posdist};
	my $tot = $posdist->{$pos_values[0]};
	for(my $k = 1; $k < @pos_values; $k++) {
		if($posdist->{$pos_values[$k]} >= $posdist->{$pos_values[$k-1]}) {
			$posdist->{$pos_values[$k]} = $posdist->{$pos_values[$k-1]};
		}
		$tot += $posdist->{$pos_values[$k]};
	}
	if($tot > 0 && $tot < 1) {
		foreach(keys %{$posdist}) {
			$posdist->{$_} /= $tot;
		}
	}
}

return $output;
}


sub updateDistributions {
(my $posdist, my $negdist, my $probptr, my $data, my $maxdiff, my $prior, my $use_decoy_set_neg) = @_;
my $output = 0;
foreach(keys %{$posdist}) {
	if(updateDistribution($posdist->{$_}, $negdist->{$_}, $data->{$_}, $probptr, $maxdiff, $prior, $data->{'decoy'}, $use_decoy_set_neg, 
		exists $MODEL_CONSTRAINTS{$_} ? $MODEL_CONSTRAINTS{$_} : "")) {
		$output = 1 ;
		printf STDERR "Exceeded max diff %s with %s distribution\n", $maxdiff, $_;
	}
}
$use_decoy_set_neg->[1] = 1 if($use_decoy_set_neg->[0] && ! $use_decoy_set_neg->[1]);
return $output;
}


# list of posdistributions and negdistributions, list of data values for those distributions, original pepprophet score probs, and current probs
sub updateProbabilities {
(my $posdisthash, my $negdisthash, my $datahash, my $probptr, my $maxdiff, my $decoyptr) = @_;
my %probs = ();
my $output = 0;
my @dists = keys %{$posdisthash};
my $newprior = 0;
my $tot = 0;
foreach(keys %{$probptr}) {
	my $prob = $INIT_UPDATE ? $datahash->{'product_probability'}->{$_} : 0.5;
	if($SET_DECOY_PROBS_2_ZERO && $decoyptr->{$_}) {
		$prob = 0;
	}
	# can add here to set $prob to zero if it is a decoy.....
	# Only update distributions for results with product probability greater than 0
	if(exists $datahash->{'product_probability'}->{$_} && $datahash->{'product_probability'}->{$_} > 0) {
		
		my $next_pos = $CROSSLINK_PRIOR; 
		my $next_neg = 1 - $CROSSLINK_PRIOR; 
		for(my $k = 0; $k < @dists; $k++) {
			next if($datahash->{$dists[$k]}->{$_} eq '*'); # ignore these spectra since not valid to participate in distribution

			$next_pos *= $posdisthash->{$dists[$k]}->{$datahash->{$dists[$k]}->{$_}};
			$next_neg *= $negdisthash->{$dists[$k]}->{$datahash->{$dists[$k]}->{$_}};
			printf "Update with $next_pos and $next_neg for $dists[$k]\n" if(0 && ! $decoyptr->{$_});
		}
		if(! $OUTPUT_DECOYS) {
			$next_pos *= $NO_DECOY_PRIOR_ADJ;
			$next_neg *= $NO_DECOY_NEGPRIOR_ADJ;
		}
		$probs{$_} = $prob * $next_pos / (($prob * $next_pos) + ((1 - $prob) * $next_neg));
		$output = 1 if(! exists $probptr->{$_} || ($probptr->{$_} >= 0 && abs($probptr->{$_} - $probs{$_}) > $maxdiff));
	}
	if(exists $probptr->{$_} && ! ($probptr->{$_} eq '') && $probptr->{$_} > 0) {
		$tot++;
		$newprior += $probs{$_};
		printf "Adding $probs{$_} for this guy....\n" if(0 && ! $decoyptr->{$_});
	}
}
if($output) {
	foreach(keys %{$probptr}) {
		$probptr->{$_} = $probs{$_} if(! exists $probptr->{$_} || $probptr->{$_} >= 0);
	}
}
printf "FINAL NEWPORIOR %0.0f for %d .....\n", $newprior, $tot;
$CROSSLINK_PRIOR = $newprior / $tot if($tot > 0);
printf "New PRIOR: %0.3f for %d probs (new prior %0.1f over tot $tot)\n", $CROSSLINK_PRIOR, scalar keys %{$probptr}, $newprior; 
$INIT_UPDATE = 0;

return $output;
}

# first step load the distribution lists with their defined keys
# second is to parse data and enter the values for each, making sure not to exceed the min or max vals of distribution
# update distrivbutions
sub parseResultsForDistributions {
(my $pprophfile, my $outfile) = @_;

my %posdists = %MODELS;
		
my $MIN_NUM_LOWS = 50; # need this many to compute decoy based fdr
my @RUNS = ();
my @RUN_SPECTRA = ();
				
my %negdists = ();
# these are all pieces of information associated with each cross-link result and included in the final output
my %distvals = ('decoy' => {}, 'decoy-decoy' => {}, 'massdiff_ppm' => {}, 'product_probability' => {}, 'homopeptide' => {}, 'peptide_pair' => {}, 
	'max_expect' => {}, 'peptide1_len' => {}, 'peptide2_len' => {}); # dist name => spectrum => score val

$distvals{'spectrum2'} = {} if($PAIRING eq "ReACT" || $PAIRING eq "Universal"); #if($REACT);
delete $posdists{'nsx'} if($NO_NSX);

my $compute_decoy_fdr = 1;
if($compute_decoy_fdr) {
	$distvals{'decoy_expect_fdr'} = {}; # use decoys to compute fdrs for results sorted by expect
	$distvals{'decoy_prob_fdr'} = {}; # use decoys to compute fdrs for results sorted by computed probability
}
my $num_low_decoys = 0;
my $num_lows = 0;
my $min_low_expect = 10;

my %sorted_distribution_values = ();
foreach(keys %posdists) {
	$distvals{$_} = {};
	# now copy values for negdistr
	$negdists{$_} = {};
	my @next = @{$posdists{$_}};
	$posdists{$_} = {}; # change
	for(my $k = 0; $k < @next; $k++) {
		$posdists{$_}->{$next[$k]} = 1 / @next; # do it even
		$negdists{$_}->{$next[$k]} = $posdists{$_}->{$next[$k]};
	}
	$sorted_distribution_values{$_} = [sort {$a <=> $b} keys %{$posdists{$_}}] if(exists $posdists{$_});
}

my $search_conditions_alert = 0;

# with mango, can have multiple queries for a spectrum, each based on different parent masses that sum to precursor and reporter
# all queries will be assigned probabilities independently unless one of the following two options is selected
my $penalize_alternative_xlinks = 0; # adjust probabilities so don't sum to more than unity
my $top1_alternative_xlinks = 0; # 12.27.17 turned this off since all hits can be correct! # only report query with highest probability with preference for non-homo peptides
die "Error: cannot choose both to penalize alternative xlinks and choose only top1\n" if($penalize_alternative_xlinks && $top1_alternative_xlinks);

my %scoreprobs = ();
# need separate scan_info for each input file, know which to use by spectrum name.....
my %scan_info = ();
my $scan_dir = getcwd() . $DIR_DIVISOR; #'/';
if($pprophfile =~ /^(\S+$DIR_DIVISOR_REGX)[^$DIR_DIVISOR_REGX]+$/) {
	$scan_dir = $1;
}
else {
	printf STDERR "Warning: using scan directory $scan_dir since not found in $pprophfile name\n";
}
my %headers = ('spectrum' => -1, 'probability' => -1, 'protein' => -1, 'peptide' => -1, 'assumed_charge' => -1, 'calc_neutral_pep_mass' => -1, 'expect' => -1);
my %result_groups = (); # by spectrum prefix to collect together alternative parent pair results for same spectrum
if($PAIRING eq 'ReACT') {

	my @react_files = ();
	
	if($LOCAL_REACT) {
		@react_files = glob($scan_dir . '*pep.xml.react2.xls');
		if(@react_files == 0) {
			open GREP, "grep '<inputfile' $file | grep -v 'name=\"int' | grep -v 'name=\"iproph' |" or die "cannot grep $!\n";
			while(<GREP>) {
				chomp();
				if(/\<inputfile name\=\"(\S+?)\" directory\=\"(\S+?)\"/) {
					push(@react_files, $2 . $DIR_DIVISOR . $1 . ".react2.xls");
					if(! -e $react_files[-1]) {
						printf STDERR "Warning: react2 file $react_files[-1] not found\n";
						sendEmail("Warning: react2.xls file $react_files[-1] not found.", "Place your react2.xls files corresponding to each search result pepXML file in the $scan_dir.");
						exit(1);
					}
					if($react_files[-1] =~ /^(\S+$DIR_DIVISOR_REGX)([^$DIR_DIVISOR_REGX]+)\.pep.xml.react2.xls$/) {
						if($LIGHT_HEAVY) {
							$LOCAL_REACT_BASES{$1 . $2} = $1;
						}
						else {
							$LOCAL_REACT_BASES{$2} = $1;
						}
					}
					else {
						die "3. here with $react_files[-1] for entry $_\n";
					}
				}
			}
			close(GREP);
		}
	}
	else {
		open GREP, "grep '<search_summary base_name' $file |" or die "cannot grep $!\n";
		while(<GREP>) {
			chomp();
			if(/\<search\_summary base\_name\=\"(\S+?)\"/) {
				my $next = $1 . ".pep.xml.react2.xls";
				if(! -e $next) {
					printf STDERR "Warning: react2 file $next not found\n";
					sendEmail("Warning: react2.xls file $next not found.", "Place your react2.xls files corresponding to each search result pepXML file in the $scan_dir.");
					exit(1);
				}
				push(@react_files, $next);
			}
			else {
				die "Problem parsing $_ for grep\n";
			}
		}
		close(GREP);
	}

	if(@react_files == 0) {
		printf STDERR "Error: No react2.xls files found in %s\n", $scan_dir;
		sendEmail("Error: No react2.xls files found in $scan_dir.", "Place your react2.xls files in the same directory with your search result pepXML files.");
		exit(1);
	}
	%scan_info = %{readCombinedReact2Texts([@react_files])};
}
elsif($PAIRING eq 'Universal') {

	my @pairing_files = ();
	
	if($LOCAL_REACT) {
		@pairing_files = glob($scan_dir . "*$PAIRFILE_SUFF");
	}
	else {
		open GREP, "grep '<search_summary base_name' $file |" or die "cannot grep $!\n";
		while(<GREP>) {
			chomp();
			if(/\<search\_summary base\_name\=\"(\S+?)\"/) {
				my $next = $1 . ".pep.xml$PAIRFILE_SUFF";
				if(! -e $next) {
					printf STDERR "Warning: pairing file $next not found\n";
					sendEmail("Warning: pairing file $next not found.", "Place your pairing files corresponding to each search result pepXML file in the $scan_dir.");
					exit(1);
				}
				push(@pairing_files, $next);
			}
			else {
				die "Problem parsing $_ for grep\n";
			}
		}
		close(GREP);
	}

	if(@pairing_files == 0) {
		printf STDERR "Error: No pairing files '*$PAIRFILE_SUFF' found in %s\n", $scan_dir;
		sendEmail("Error: No pairing files '*$PAIRFILE_SUFF' found in $scan_dir.", "Place your pairing files in the same directory with your search result pepXML files.");
		exit(1);
	}

	%scan_info = %{readCombinedPairfileTexts([@pairing_files], $PAIRFILE_SUFF)};
}
my $use_decoy_set_negdist = [! $TARGET_NEG, 0]; # 0/1 whether to use decoys, 0/1 whether they were used and the distributions are set
			
my $SPEC_TOT = 0;
			
my @headers = keys %headers;
my $first = 1;
my %results = (); # hashed by spectrum to peptide1 and peptide2 results
my @pepxml_header = ();
my $header_on = 0;
(my $sec,my $min,my $hour,my $mday,my $mon,my $year,my $wday,my $yday,my $isdst) = localtime();
my $timestamp = sprintf "%s-%s-%sT%s:%s:%s", $year+1900, $mon, $mday, $hour, $min, $sec; 
my @INPUT_FILES = ();
my @next_hit = ();
my $hit_on = 0;
my $interprophet = 0;
my $prob_type_known = 0;
my $base = "";
my $nomod = 0;
my $current_run = "";
my $naked_peptide = ""; # without mods
open(FILE, $pprophfile) or die "cannot read $pprophfile $!\n";
if($pprophfile =~ /xml$/) {

	printf STDERR "Parsing pepXML file %s\n", $pprophfile;
	my $spectrum = '';
	my $assumed_charge = -1;
	my $probability = -1;
	my $expect = '';
	my $pepmass = '';
	my $deltacn = '';
	my $protein = '';
	my $peptide = '';
	my $prev = '';
	my $foll = '';
	my $prec_mass = '';
	my $orig_spec = "";
	while(<FILE>) {
		chomp();
		if(! ($OUTPUT_PEPXML eq '')) {
			if(! $prob_type_known && /\<analysis\_summary analysis=\"interprophet\"/) {
				$interprophet = 1;
				$prob_type_known = 1;
				printf STDERR "Using iProphet probabilities in $pprophfile\n";
			}
			if(/\<inputfile name.*directory\=/) {
				;
			}
			elsif(/\<inputfile name/) {
				push(@INPUT_FILES, $_);
			}
			elsif(/\<msms\_run\_summary base\_name\=/) {
				push(@RUNS, [$_]);
				push(@RUN_SPECTRA, []);
				$header_on = 1;
				if($LIGHT_HEAVY && /base\_name\=\"(\S+?)\"/) {
					$current_run = $1;
				}
			}
			elsif($header_on) {
				if(/\<spectrum\_query/) {
					$header_on = 0; # done
				}
				elsif(/^(.*?\<search\_summary base\_name\=.*?search\_engine\=\")Comet(\".*)$/) {
					push(@{$RUNS[-1]}, $1 . 'Kojak' . $2);
				}
				elsif(/(\<analysis\_timestamp analysis\=\"peptideprophet\" time=\").*?(\" id=\"1\"\/\>)/) {
					push(@{$RUNS[-1]}, $1 . $timestamp . $2);
					$header_on = 0;
				}
				else {
					push(@{$RUNS[-1]}, $_);
				}
			}
		} 

		if(/\<spectrum_query spectrum\=\"(\S+?)\".*?precursor\_neutral\_mass=\"(\S+)\".*?assumed\_charge\=\"(\d+?)\"/) {
			$spectrum = $1;
			$prec_mass = $2;
			$assumed_charge = $3;
			$prob_type_known = 1; # past the header info
			$nomod = 0; # reset
			$orig_spec = $spectrum;
		}
		elsif(/\<search\_summary base\_name\=\"(\S+?)\"/) {
			$base = $1;
		}
		elsif(/\<search\_hit hit\_rank\=\"1\" peptide\=\"(\S+?)\".*?protein\=\"(\S+?)\".*?calc\_neutral\_pep\_mass\=\"(\S+?)\"/) {
			$peptide = $1;

			$protein = $2;
			$pepmass = $3;
			if(! ($OUTPUT_PEPXML eq '') && /\<search_hit hit\_rank\=\"1\"(.*)$/) {
				@next_hit = ('<linked_peptide' . $1);
				$hit_on = 1;
			}
			# 11.22.19
			if(1 || $REMOVE_STATIC_MODS) {
				$naked_peptide = $peptide; # used to assess mods
			}
			
			
		}
		elsif(/\<alternative\_protein protein\=\"(\S+?)\"/) {
			$protein .= ',' . $1;
		}
		elsif(/\<modification\_info modified\_peptide\=\"(\S+?)\"/) {
			$peptide = $1;
			
			$peptide = removeStaticMods($peptide);
		}
		elsif(! $REMOVE_STATIC_MODS && /\<mod\_aminoacid\_mass position\=\"(\d+)\" mass\=\"(\S+?)\".*?variable/) {
			$peptide = addModpepDecimalPlaces($peptide, $1, $2);
			# 11.22.19
			my $next_abbrev = substr($naked_peptide, $1 - 1, 1) . ":" . sprintf("%0.".$NUM_MODPEP_DECIMALS."f", $2);
			my $next_full = substr($naked_peptide, $1 - 1, 1) . ":" . $2;
			if(exists $FULL_STUMP_MASSES{$next_abbrev}) {
				if(! ($FULL_STUMP_MASSES{$next_abbrev} eq $next_full)) {
					printf "Error: have two full length modification masses in results: %s and %s\n", $FULL_STUMP_MASSES{$next_abbrev}, $next_full;
					exit(1);
				}
			}
			else {
				$FULL_STUMP_MASSES{$next_abbrev} = $next_full;
			}
		}
		elsif($REMOVE_STATIC_MODS && /\<mod\_aminoacid\_mass position\=\"(\d+)\" mass\=\"(\S+?)\"/) {
			my $pos = $1 - 1;
			my $mass = $2;
			my $aa = substr($naked_peptide, $pos, 1);
			if((! ($aa eq 'C') || $mass != 160.030649) && 
				(! ($aa eq 'K') || $mass != 136.109162) && (! ($aa eq 'R') || $mass != 162.121240)) {
			
				$peptide = addModpepDecimalPlaces($peptide, $1, $2);
			}
		}
		elsif(/\<search\_score name\=\"expect\" value\=\"(\S+?)\"/) {
			$expect = $1;
		}
		elsif(exists $distvals{'deltacn_bin'} && /\<search\_score name\=\"deltacn\" value\=\"(\S+?)\"/) {
			$deltacn = $1;
		}
		elsif(! $interprophet && /\<peptideprophet\_result probability\=\"(\S+?)\"/) {
			$probability = $nomod ? 0 : $1;
		}
		elsif($interprophet && /\<interprophet\_result probability\=\"(\S+?)\"/) {
			$probability = $nomod ? 0 : $1;
		}
		elsif(/\<\/spectrum_query\>/) {
			next if(! ($SPEC_TEXT eq '') && $spectrum !~ /$SPEC_TEXT/);
			
			if($DECOY_FRACTION < 1) {
				my @prots = split(",", $protein);
				# check for decoy
				my $decoy = 1;
				for(my $p = 0; $p < @prots; $p++) {
					if($prots[$p] !~ /rev\_/) {
						$decoy = 0;
						$p = @prots;
					}
				}
				next if($decoy && rand() < $DECOY_FRACTION);
			}
			
			next if($MIN_PEPTIDE_LENGTH > 0 && getPeptideLength($peptide) < $MIN_PEPTIDE_LENGTH);

			$SPEC_TOT++;
			
  			if($PAIRING eq 'Mango' && $spectrum =~ /^(\S+)(\_\d\d\d\_)([A,B])\.(\d+)\.(\d+)(\.\d)$/) {
				my $pref = $base . $2;
				my $index = $3 eq 'A' ? 1 : $3 eq 'B' ? 2 : die "Error: $3 is not A or B, cannot determine index\n";
				my $scan1 = $4;
				my $scan2 = $5;
				my $suff = $6;
				$scan_info{$base}++;
					$spectrum = $pref . $scan1 . '.' . $scan1;
				$orig_spec = "";
				$results{$spectrum} = {} if(! exists $results{$spectrum});
				$results{$spectrum}->{$index} = [$probability, $protein, 
					$peptide, $assumed_charge, $pepmass, $expect];	

				push(@{$results{$spectrum}->{$index}}, $deltacn) if(exists $distvals{'deltacn_bin'});
				if(! ($OUTPUT_PEPXML eq '')) {
					if(@next_hit > 0 && $next_hit[0] =~ /(<linked.*?)\>/) {
						$next_hit[0] = $1 . (sprintf " complement_mass=\"%s\" precursor_neutral_mass=\"%s\" designation=\"%s\" assumed_charge=\"%d\">",
						$pepmass, $prec_mass, $index == 1 ? "alpha" : "beta", $assumed_charge);
					}
					else {
						die "problem parsing $next_hit[0] of next_hit\n";
					}
					push(@{$results{$spectrum}->{$index}}, [@next_hit]);
				}

				$pref = $1 if($pref =~ /^(\S+)\_\d\d\d\_$/);
				$pref .= '_' . $scan1;
				$result_groups{$pref} = {} if(! exists $result_groups{$pref});
				$result_groups{$pref}->{$spectrum}++;
				
				push(@{$RUN_SPECTRA[-1]}, $spectrum) if(! ($OUTPUT_PEPXML eq '') && $index == 1);
			}
			elsif(! ($PAIRING eq 'Mango') && $spectrum =~ /^(\S+)(\.)(\d+)\.(\d+)(\.\d)$/) {
				if($LOCAL_REACT) {
	
					if(scalar keys %LOCAL_REACT_BASES > 0) {
						my $found = 0;
						if($LIGHT_HEAVY && exists $LOCAL_REACT_BASES{$current_run}) { # have to substitute correct base position
							$base = $LOCAL_REACT_BASES{$current_run} . $1;
							$found = 1;
						}
						elsif(! $LIGHT_HEAVY && exists $LOCAL_REACT_BASES{$1}) {
							$base = $LOCAL_REACT_BASES{$1} . $1;
							$found = 1;
						}
						if(! $found) {
							printf STDERR "Warning: no local react react2.xls file found for $1\n";
							printf STDERR "Make sure your search result file runs match with current pepxml file locations\n";
							exit(1);
						}
					}
					else {
						$base = $scan_dir . $1;
					}
				}
				my $pref = $base . '_000_'; 
				my $scan1 = $3;
				my $scan2 = $4;
				my $suff = $5;
				my $index = 1;
				# check whether this is peptide1 or 2 for a pair
				if(exists $scan_info{$base} && exists $scan_info{$base}->{$scan1}) {
				
					if($REACT_USE_ADJ_PRIOR) {
						$probability = sprintf "%0.3f", getPriorAdjustedProbability($probability, $REACT_UNFILT_PRIOR, $REACT_FILT_PRIOR);
					}
					my $orig_spectrum = $spectrum;
					my $index = exists $scan_info{$base}->{$scan1}->{'mass'} ? 1 :
						exists $scan_info{$base}->{$scan1}->{'partner'} ? 2 : die "invalid scan type for $base and $scan1\n";
					$spectrum = $index == 1 ? $pref . $scan1 . '.' . $scan1 : $scan_info{$base}->{$scan1}->{'partner'};
					
					if($peptide =~ /K\[13/) {
						printf "Here with peptide $peptide\n"; exit(1);
					}
					$results{$spectrum} = {} if(! exists $results{$spectrum});
					$results{$spectrum}->{$index} = [$probability, $protein, 
					$peptide, $assumed_charge, $pepmass, $expect, $scan1];
					push(@{$results{$spectrum}->{$index}}, $deltacn) if(exists $distvals{'deltacn_bin'});
					if(! ($OUTPUT_PEPXML eq '')) {
						if(@next_hit > 0 && $next_hit[0] =~ /(<linked.*?)\>/) {
							$next_hit[0] = $1 . (sprintf " complement_mass=\"%s\" precursor_neutral_mass=\"%s\" designation=\"%s\" ms3_scan=\"%d\" assumed_charge=\"%d\">",
							$pepmass, $prec_mass, $index == 1 ? "alpha" : "beta", $scan1, $assumed_charge);
						}
						else {
							die "problem parsing $next_hit[0] of next_hit\n";
						}
						push(@{$results{$spectrum}->{$index}}, [@next_hit]);
						push(@{$results{$spectrum}->{$index}->[-1]}, $scan_info{$base}->{$scan1}->{'parent_scan'}) if($PAIRING eq 'ReACT' && $index==1);
					}
			
					$pref = $1 if($pref =~ /^(\S+)\_\d\d\d\_$/);
					if($index==1 && $pref =~ /$DIR_DIVISOR_REGX([^$DIR_DIVISOR_REGX]+)$/) {
						
						$pref = $1 . '_' . $scan_info{$base}->{$scan1}->{'parent_scan'};
						my @next = keys %{$scan_info{$base}};
						my @next2 = keys %{$scan_info{$base}->{$scan1}};
						$result_groups{$pref} = {} if(! exists $result_groups{$pref});
						$result_groups{$pref}->{$spectrum}++;
					}
					$distvals{'spectrum2'}->{$spectrum} = $base . '.' . $scan1 . '.' . $scan1 if($index == 2);
					push(@{$RUN_SPECTRA[-1]}, $spectrum) if(! ($OUTPUT_PEPXML eq '') && $index == 1);
				}
				elsif(! exists $scan_info{$base}) {
					next;
				}
				else {
					next;
				}
			}
			else {
				die "--> problem parsing $spectrum\n";
			}
		}
		if($hit_on) {
			if(/\<search_score/) {
				printf STDERR "Warning: have unmodified peptide: check search conditions!!!\n$spectrum with peptide $peptide\n" if(! $search_conditions_alert);
				$search_conditions_alert = 1;
				$hit_on = 0;
				$nomod = 1; # record no mod
			}
			if(/\<\/modification_info\>/) {
				$hit_on = 0;
			}
			push(@next_hit, $_) if(! /\<search_hit(.*)$/ && ! /\<search_score/);
		}
	}
}
else {
	printf STDERR "Parsing XLS file %s\n", $pprophfile;
	while(<FILE>) {
		chomp();
		my @parsed = split("\t");
		if($first) {
			$first = 0;
			for(my $k = 0; $k < @parsed; $k++) {
				$headers{$parsed[$k]} = $k if(exists $headers{$parsed[$k]});
			}
			for(my $k = 0; $k < @headers; $k++) {
				die "Could not find $headers[$k] in $_\n" if($headers{$headers[$k]} == -1);
			}
		}
		else {
			if($PAIRING eq 'Mango' && $parsed[$headers{'spectrum'}] =~ /^(\S+)(\_\d\d\d\.)(\d+)\.(\d+)(\.\d)$/) {
				my $base = $1;
				printf "next base: $base\n";
				my $pref = $base . $2;
				my $scan1 = $3;
				my $scan2 = $4;
				my $suff = $5;
				my $index = 1;
				if($scan1 =~ /^1(\d\d\d\d\d)$/) {
					$scan1 = $1;
					$index = 2;
				}
				elsif($scan1 !~ /\d\d\d\d\d/) {
					die "problem parsing $parsed[$headers{'spectrum'}]\n";
				}
				$scan_info{$base}++;
			
				$parsed[$headers{'spectrum'}] = $pref . $scan1 . '.' . $scan1;
				$results{$parsed[$headers{'spectrum'}]} = {} if(! exists $results{$parsed[$headers{'spectrum'}]});
				$results{$parsed[$headers{'spectrum'}]}->{$index} = [$parsed[$headers{'probability'}], $parsed[$headers{'protein'}], 
					$parsed[$headers{'peptide'}], $parsed[$headers{'assumed_charge'}], $parsed[$headers{'calc_neutral_pep_mass'}],
					$parsed[$headers{'expect'}]];	

				push(@{$results{$parsed[$headers{'spectrum'}]}->{$index}}, $parsed[$headers{'deltacn'}]) if(exists $distvals{'deltacn_bin'});

				$pref = $1 if($pref =~ /^(\S+)\_\d\d\d\.$/);
				$pref .= '_' . $scan1;
				$result_groups{$pref} = {} if(! exists $result_groups{$pref});
				$result_groups{$pref}->{$parsed[$headers{'spectrum'}]}++;
			}
			elsif(! ($PAIRING eq 'Mango') && $parsed[$headers{'spectrum'}] =~ /^(\S+)(\.)(\d+)\.(\d+)(\.\d)$/) {
				my $base = $1;
				my $pref = $base . '_000_' . $2;
				my $scan1 = $3;
				my $scan2 = $4;
				my $suff = $5;
				my $index = 1;
				# check whether this is peptide1 or 2 for a pair
				if(exists $scan_info{$base} && exists $scan_info{$base}->{$scan1}) {
					my $index = exists $scan_info{$base}->{$scan1}->{'mass'} ? 1 :
						exists $scan_info{$base}->{$scan1}->{'partner'} ? 2 : die "invalid scan type for $base and $scan1\n";
					my $spectrum = $index == 1 ? $pref . $scan1 . '.' . $scan1 : $scan_info{$base}->{$scan1}->{'partner'};
					$results{$spectrum} = {} if(! exists $results{$spectrum});
					$results{$spectrum}->{$index} = [$parsed[$headers{'probability'}], $parsed[$headers{'protein'}], 
						$parsed[$headers{'peptide'}], $parsed[$headers{'assumed_charge'}], $parsed[$headers{'calc_neutral_pep_mass'}],
						$parsed[$headers{'expect'}]];
						
					push(@{$results{$spectrum}->{$index}}, $parsed[$headers{'deltacn'}]) if(exists $distvals{'deltacn_bin'});
					$pref = $1 if($pref =~ /^(\S+)\_\d\d\d\.$/);
					$pref .= '_' . $scan1;
					$result_groups{$pref} = {} if(! exists $result_groups{$pref});
					$result_groups{$pref}->{$spectrum}++;
					$distvals{'spectrum2'}->{$spectrum} = $base . '.' . $scan1 . '.' . $scan1 if($index == 2);
				}
				elsif(! exists $scan_info{$base}) {
					next;
				}
				else {
					next;
				}
			}
			else {
				die "--> problem parsing $parsed[$headers{'spectrum'}]\n";
			}
		}
	}
} # pepxls
close(FILE); 
printf STDERR "Read in results for %d spectra\n", scalar keys %results;
printf STDERR "Read in %d unique input file bases: %s\n", scalar keys %scan_info, join(",", keys %scan_info);
printf "SPEC TOT: $SPEC_TOT\n";
if($PAIRING eq 'Mango') {
	foreach my $next (keys %scan_info) {
	#printf "Reading next peak file $next.peaks\n";
		$scan_info{$next} = readPeaksText($next.".peaks", exists $posdists{'massdiff_offset'});
	}
}
my $tot_bal = 0;
my $tot_unbal = 0;
my %seen_mango_scan_pps = ();
foreach(keys %results) {
	if(! exists $results{$_}->{1} || ! exists $results{$_}->{2} || @{$results{$_}->{1}} == 0 || @{$results{$_}->{2}} == 0) {
		delete $results{$_};
		$tot_unbal++;
		next;
	}

	$tot_bal++;
	$results{$_}->{1}->[2] = $1 if($results{$_}->{1}->[2] =~ /^\S\.(\S+)\.\S$/);
	$results{$_}->{2}->[2] = $1 if($results{$_}->{2}->[2] =~ /^\S\.(\S+)\.\S$/);
	$distvals{'peptide_pair'}->{$_} = $results{$_}->{1}->[2]  . '_' . $results{$_}->{2}->[2];
	if(exists $distvals{'nrx'}) {
		my $next_nrx_value = $LIGHT_HEAVY || $USE_STRIPPED_PEP_NRX ? getPeptidePairForNrx($results{$_}->{1}->[2], $results{$_}->{2}->[2], 
			stripPeptide($results{$_}->{1}->[2]), stripPeptide($results{$_}->{2}->[2])) : $distvals{'peptide_pair'}->{$_};
		printf "Have nrx %s for %s, l-h: %d versus %s\n", $next_nrx_value, $_, $LIGHT_HEAVY ? 1 : 0, getPeptidePairForNrx($results{$_}->{1}->[2], $results{$_}->{2}->[2], 
			stripPeptide($results{$_}->{1}->[2]), stripPeptide($results{$_}->{2}->[2])) if(0);
		$PEP_PAIRS{$_} = $next_nrx_value;

		if($PAIRING eq 'Mango' && /^(\S+)\_(\d\d\d)\_(\S+)$/) {
			my $next_entry = $1  . '.' . $3 . ":" . $next_nrx_value;
			next if(exists $seen_mango_scan_pps{$next_entry});  # only keep the first (highest prob) instance cross-link for each mango spectrum assigned to a peptide parir
			$seen_mango_scan_pps{$next_entry}++;
		}

		$CROSSLINKS{$next_nrx_value} = {} if(! exists $CROSSLINKS{$next_nrx_value});
		$CROSSLINKS{$next_nrx_value}->{$_}++;
	}
	recordProteinPairs($results{$_}->{1}->[1], $results{$_}->{2}->[1], \%prot_crosslinks, $results{$_}->{1}->[0] * $results{$_}->{2}->[0], $distvals{'peptide_pair'}->{$_});
		$distvals{'spectrum2'}->{$_} .= '.' . $results{$_}->{2}->[3] if(! ($PAIRING eq 'Mango')); 

}
printf "Have a total of %d peptide pairs, having excluded %d without partners\n", $tot_bal, $tot_unbal;

my %probabilities = ();
my %parent_massdiffs = ();
my %parent_massdiffs_ppm = ();
my %parent_chargediffs = ();

my $prior_with_decoys = 0;
my $num_decoys = 0;
my $prior_without_decoys = 0;
my $tot_num = 0;
my @NSX_VALS = ();
my $MIN_PENALTY = 0.99; 
my @IQPIR_reportermasses = keys %IQPIR_REPORTERMASSES;
foreach(sort {$a cmp $b} keys %results) {


	next if(scalar keys %{$results{$_}} != 2);
	$distvals{'product_probability'}->{$_} = sprintf "%0.4f", $results{$_}->{1}->[0] * $results{$_}->{2}->[0] * $MIN_PENALTY;
    $probabilities{$_} = $distvals{'product_probability'}->{$_};
	# now have option of using more info...
	my $base = '';
	my $query_no = -1;
	if(/^(\S+)\_(\d\d\d)\_/) {
		$base = $1;
		$query_no = sprintf "%d", $2;
	}
	else {
		die "Error, no base in $_\n";
	}
	die "no scan file available for $base of $_\n" if(! exists $scan_info{$base});
	my $scan = -1;
	if(/\.0*([^\.]+)$/) {
		$scan = $1;
		while(! ($PAIRING eq 'Mango') && length $scan < 5) {
			$scan = '0' . $scan;
		}
		# compute the new scores
		if(exists $scan_info{$base}->{$scan}) {
			$distvals{'decoy'}->{$_} = getDecoy($results{$_}->{1}->[1], $results{$_}->{2}->[1]) || getNonCrosslinkedPeptidePair($results{$_}->{1}->[2], $results{$_}->{2}->[2]);
			$distvals{'decoy-decoy'}->{$_} = getDecoyDecoy($results{$_}->{1}->[1], $results{$_}->{2}->[1]);
		 
			$distvals{'max_expect'}->{$_} = $results{$_}->{1}->[5] > $results{$_}->{2}->[5] ? $results{$_}->{1}->[5] : $results{$_}->{2}->[5];
			$distvals{'peptide1_len'}->{$_} = getPeptideLength($results{$_}->{1}->[2]);
			$distvals{'peptide2_len'}->{$_} = getPeptideLength($results{$_}->{2}->[2]);
			
			if($compute_decoy_fdr && $distvals{'max_expect'}->{$_} > $min_low_expect) {
				$num_low_decoys++ if($distvals{'decoy'}->{$_});
				$num_lows++;
			}

			if(exists $distvals{'reporter_charge'}) {
				$distvals{'reporter_charge'}->{$_} = $scan_info{$base}->{$scan}->{'charge'} - $results{$_}->{1}->[3] - $results{$_}->{2}->[3];
				$distvals{'reporter_charge'}->{$_} = geValueWithinDistributionBounds($distvals{'reporter_charge'}->{$_}, $sorted_distribution_values{'reporter_charge'});
			}
			# here have to change REPORTER MASS depending on the peptides $results{$_}->{1}->[2], $results{$_}->{2}->[2]
			if(@IQPIR_reportermasses > 0) {
				$REPORTERMASS = -1; 
				my @reportermasses = (0, 0);
				for(my $z = 0; $z < @IQPIR_reportermasses; $z++) {
					if($results{$_}->{1}->[2] =~ /K\[$IQPIR_reportermasses[$z]\]/ && $results{$_}->{2}->[2] =~ /K\[$IQPIR_reportermasses[$z]\]/) {
						$REPORTERMASS = $IQPIR_REPORTERMASSES{$IQPIR_reportermasses[$z]};
						$z = @IQPIR_reportermasses;
					}
					else {
						if($results{$_}->{1}->[2] =~ /\[$IQPIR_reportermasses[$z]\]/) {
							$reportermasses[0] = $IQPIR_REPORTERMASSES{$IQPIR_reportermasses[$z]};
						}
						if($results{$_}->{2}->[2] =~ /\[$IQPIR_reportermasses[$z]\]/) {
							$reportermasses[1] = $IQPIR_REPORTERMASSES{$IQPIR_reportermasses[$z]};
						}
					}
				}
				if($REPORTERMASS == -1) {
					if($reportermasses[0] != 0 && $reportermasses[1] != 0 && $reportermasses[0] == $reportermasses[1]) {
						$REPORTERMASS = $reportermasses[0];
					}
					else {
						die "Error: no reporter mass matching " . join(",", @IQPIR_reportermasses). " found consistent with mods on $results{$_}->{1}->[2] and/or $results{$_}->{2}->[2]\n";
					}
				}
			}
			
			$distvals{'massdiff_ppm'}->{$_} = $scan_info{$base}->{$scan}->{'mass'} - $results{$_}->{1}->[4] - $results{$_}->{2}->[4] - $REPORTERMASS; 
			# can remove the offsets and convert to ppm
			if($distvals{'massdiff_ppm'}->{$_} < -4) {
				$distvals{'massdiff_ppm'}->{$_} += 4;
			}
			elsif($distvals{'massdiff_ppm'}->{$_} < -3) {
				$distvals{'massdiff_ppm'}->{$_} += 3;
			}
			elsif($distvals{'massdiff_ppm'}->{$_} < -2) {
				$distvals{'massdiff_ppm'}->{$_} += 2;
			}
			elsif($distvals{'massdiff_ppm'}->{$_} < -1) {
				$distvals{'massdiff_ppm'}->{$_} += 1;
			}
			if($distvals{'massdiff_ppm'}->{$_} > 4) {
				$distvals{'massdiff_ppm'}->{$_} -= 4;
			}
			elsif($distvals{'massdiff_ppm'}->{$_} > 3) {
				$distvals{'massdiff_ppm'}->{$_} -= 3;
			}
			elsif($distvals{'massdiff_ppm'}->{$_} > 2) {
				$distvals{'massdiff_ppm'}->{$_} -= 2;
			}
			elsif($distvals{'massdiff_ppm'}->{$_} > 1) {
				$distvals{'massdiff_ppm'}->{$_} -= 1;
			}
			if($distvals{'massdiff_ppm'}->{$_} > 0.5) {
				$distvals{'massdiff_ppm'}->{$_} = 1 - $distvals{'massdiff_ppm'}->{$_};
			}
			$distvals{'massdiff_ppm'}->{$_} = sprintf "%0.4f", $distvals{'massdiff_ppm'}->{$_} * 1000000 / $scan_info{$base}->{$scan}->{'mass'};
			

			if(exists $distvals{'massdiff_offset'}) {
				die "problem with query no $query_no and ", scalar @{$scan_info{$base}->{$scan}->{'pepmass_diff'}}, "\n" if($query_no > $#{$scan_info{$base}->{$scan}->{'pepmass_diff'}});
				$distvals{'massdiff_offset'}->{$_} = sprintf "%0.0f", $scan_info{$base}->{$scan}->{'pepmass_diff'}->[$query_no];
				$distvals{'massdiff_offset'}->{$_} = '0' if($distvals{'massdiff_offset'}->{$_} eq '-0');
				$distvals{'massdiff_offset'}->{$_} = geValueWithinDistributionBounds($distvals{'massdiff_offset'}->{$_}, $sorted_distribution_values{'massdiff_offset'});
			}
			if(exists $distvals{'xlinkdb'}) {
				$distvals{'xlinkdb'}->{$_} = exists $PEP_PAIRS{$_} && exists $XLINKDB_XL->{$PEP_PAIRS{$_}} ? $XLINKDB_XL->{$PEP_PAIRS{$_}} : 0;
				if($distvals{'xlinkdb'}->{$_} > 1) {
					if($distvals{'xlinkdb'}->{$_} < 3) {
						$distvals{'xlinkdb'}->{$_} = 1;
					}
					elsif($distvals{'xlinkdb'}->{$_} < 10) {
						$distvals{'xlinkdb'}->{$_} = 2;
					}
					else {
						$distvals{'xlinkdb'}->{$_} = 3;
					}
				}
			}
			if(exists $distvals{'aa_score'}) {
				$distvals{'aa_score'}->{$_} = exists $PEP_PAIRS{$_} ? getAAScore($PEP_PAIRS{$_}) : 0;
				if($distvals{'aa_score'}->{$_} < 0.2) {
					$distvals{'aa_score'}->{$_} = 0;
				}
				elsif($distvals{'aa_score'}->{$_} < 0.3) {
					$distvals{'aa_score'}->{$_} = 1;
				}
				elsif($distvals{'aa_score'}->{$_} < 0.4) {
					$distvals{'aa_score'}->{$_} = 2;
				}
				else {
					$distvals{'aa_score'}->{$_} = 3;
				}
			}
			if(! $OUTPUT_DECOYS) {
				if(! $distvals{'decoy'}->{$_}) {
					$prior_without_decoys += $distvals{'product_probability'}->{$_};
				}
				else {
					$num_decoys++;
				}
				$prior_with_decoys += $distvals{'product_probability'}->{$_};
				$tot_num++;
			}			
		}
		else {
			die "Error: no scaninfo available for $scan of $_ with base $base\n";
		}
	}
	else {
		die "Error: now scan for spectrum $_\n";
	}
	
	if(exists $distvals{'massdiff_bin'}) {
		$distvals{'massdiff_bin'}->{$_} = getBinValue($distvals{'massdiff_ppm'}->{$_}, 'massdiff_bin');
		$distvals{'massdiff_bin'}->{$_} = geValueWithinDistributionBounds($distvals{'massdiff_bin'}->{$_}, $sorted_distribution_values{'massdiff_bin'});
	}
		
	if(exists $distvals{'total_charge'}) {
		$distvals{'total_charge'}->{$_} = $results{$_}->{1}->[3] + $results{$_}->{2}->[3];
		$distvals{'total_charge'}->{$_} = geValueWithinDistributionBounds($distvals{'total_charge'}->{$_}, $sorted_distribution_values{'total_charge'});
	}

	$distvals{'homopeptide'}->{$_} = $results{$_}->{1}->[2] eq $results{$_}->{2}->[2] ? 1 : 0;

	if(exists $distvals{'intra'}) {
		$distvals{'intra'}->{$_} = $distvals{'homopeptide'}->{$_} ? '*' : getIntra($results{$_}->{1}->[1], $results{$_}->{2}->[1]);
	}
	if(exists $distvals{'nsx'}) {
		$distvals{'nsx'}->{$_} = $distvals{'homopeptide'}->{$_} ? '*' : getNsx($results{$_}->{1}->[1], $results{$_}->{2}->[1], \%prot_crosslinks, $distvals{'peptide_pair'}->{$_}); 
		if(! $distvals{'homopeptide'}->{$_}) {
			push(@NSX_VALS, $distvals{'nsx'}->{$_}) if($distvals{'decoy'}->{$_} && $distvals{'nsx'}->{$_} > 0);
			$distvals{'nsx'}->{$_} = geValueWithinDistributionBounds($distvals{'nsx'}->{$_}, $sorted_distribution_values{'nsx'});
			$NSP_VALS{$_} = $CURRENT_NSP; # get the current global variable value
		}
	}
	if(exists $distvals{'nrx'}) {
		$distvals{'nrx'}->{$_} = getNrx($_, \%probabilities);
	}
		
	if(exists $distvals{'joint_score'}) {
		if($distvals{'product_probability'}->{$_} < 0.05) {
			$distvals{'joint_score'}->{$_} = 0;
		}
		elsif($distvals{'product_probability'}->{$_} < 0.5) {
			$distvals{'joint_score'}->{$_} = 1;
		}
		elsif($distvals{'product_probability'}->{$_} < 0.9) {
			$distvals{'joint_score'}->{$_} = 2;
		}
		else {
			$distvals{'joint_score'}->{$_} = 3;
		}
	}
}

printf "HAVE %d PROBS\n", scalar keys %probabilities; 
printf "TOT: $tot_num with decoys: $num_decoys\n";# exit(1);
if(! $OUTPUT_DECOYS) {
	$prior_without_decoys /= ($tot_num - $num_decoys) if($tot_num - $num_decoys > 0);
	$prior_with_decoys /= $tot_num if($tot_num > 0);
	if($prior_with_decoys > 0 && $prior_with_decoys < 1) {
		$NO_DECOY_PRIOR_ADJ = $prior_without_decoys / $prior_with_decoys;
		$NO_DECOY_NEGPRIOR_ADJ = (1 - $prior_without_decoys) / (1 - $prior_with_decoys);
		printf STDERR "Using NO_OUTPUT_DECOYS prior adjustment of %0.2f and negprior adjustment of %0.3f\n", 
			$NO_DECOY_PRIOR_ADJ, $NO_DECOY_NEGPRIOR_ADJ;
	}
	else {
		printf STDERR "Warning: Have prior $prior_with_decoys\n";
	}
}
if($compute_decoy_fdr) {
	if($num_lows < $MIN_NUM_LOWS) {
		printf STDERR "Warning: %d low scoring results with e-score above %s, fewer than $MIN_NUM_LOWS required to compute decoy fdr values\n", $num_lows, $min_low_expect;
		$compute_decoy_fdr = 0;
	}
	elsif($num_low_decoys == 0) {
		printf STDERR "Warning: no low decoys with e-score above %s, cannot compute decoy fdr values\n", $min_low_expect;
		$compute_decoy_fdr = 0;
	}
}

if($compute_decoy_fdr) {
	$num_lows /= $num_low_decoys;
	$num_lows -= 1 if(! $OUTPUT_DECOYS);
	printf STDERR "Using decoy factor %0.2f to compute decoy fdrs\n", $num_lows;
	my @sorted_expects = sort {$distvals{'max_expect'}->{$a} <=> $distvals{'max_expect'}->{$b}} keys %results;
	my $tot = 0;
	my $num_decoys = 0;
	my $num_decoy_decoys = 0;
	my $prev_fdr = 0;
	for(my $k = 0; $k < @sorted_expects; $k++) {
		$tot++ if($OUTPUT_DECOYS || ! $distvals{'decoy'}->{$sorted_expects[$k]});
		$num_decoys++ if($distvals{'decoy'}->{$sorted_expects[$k]});
		$num_decoy_decoys++ if($distvals{'decoy-decoy'}->{$sorted_expects[$k]});
		my $init = $k;
		while($k < @sorted_expects - 1 && $distvals{'max_expect'}->{$sorted_expects[$k+1]} == $distvals{'max_expect'}->{$sorted_expects[$init]}) {
			$k++;
			$tot++ if($OUTPUT_DECOYS || ! $distvals{'decoy'}->{$sorted_expects[$k]});
			$num_decoys++ if($distvals{'decoy'}->{$sorted_expects[$k]});
			$num_decoy_decoys++ if($distvals{'decoy-decoy'}->{$sorted_expects[$k]});
		}
		my $fdr = $tot > 0 && $num_decoys > 2 * $num_decoy_decoys ? ($num_decoys - 2 * $num_decoy_decoys) / $tot : 0;
		$fdr = $prev_fdr if($fdr < $prev_fdr);
		for(my $j = $init; $j <= $k; $j++) {
			$distvals{'decoy_expect_fdr'}->{$sorted_expects[$j]} = sprintf "%0.4f", $fdr;
		}
		$prev_fdr = $fdr;
	}
}

my %adjusted = ();
if($MIN_PEP_LEN > 0) {
	foreach(sort {$a cmp $b} keys %probabilities) {
		if($probabilities{$_} > 0 && $distvals{'peptide1_len'}->{$_} < $MIN_PEP_LEN && $distvals{'peptide2_len'}->{$_} < $MIN_PEP_LEN) {
			printf STDERR "Setting probability of %s (%s) to -1 since peptide 1 and 2 lengths %d and %d smaller than %d\n", 
				$_, $probabilities{$_}, $distvals{'peptide1_len'}->{$_}, $distvals{'peptide2_len'}->{$_}, $MIN_PEP_LEN 
					if($probabilities{$_} > 0.5);
			$probabilities{$_} = -2;
			$adjusted{$_} = '' if($OUTPUT_ADJUSTED);
		}
	}
}

# have the option of pre-filtering mango's result grouops---- exclude all runner-ups that are less than 0.5 the original probability of the top
my $prefilter_runnersup = 1;
if($PAIRING eq 'Mango' && $prefilter_runnersup) {
	my $min_factor = 0.5;
	my $min_top_prob = 0.75;
	my $tot_rejected = 0;
	foreach(sort {$a cmp $b} keys %result_groups) {
		my @ranked_spectra = ();
		my @next = keys %{$result_groups{$_}};
		for(my $k = 0; $k < @next; $k++) {
			push(@ranked_spectra, [$next[$k], $probabilities{$next[$k]}, $distvals{'peptide_pair'}->{$next[$k]}]);
		}
		@ranked_spectra = sort {$b->[1] <=> $a->[1]} @ranked_spectra; # sort by prob
		if($ranked_spectra[0]->[1] > $min_top_prob) {
			for(my $k = 1; $k < @ranked_spectra; $k++) {
				if($ranked_spectra[$k]->[1] < $min_factor * $ranked_spectra[0]->[1]) { # done and discard
					for(my $j = $k; $j < @ranked_spectra; $j++) {
						$probabilities{$ranked_spectra[$j]->[0]} = -1;
						printf STDERR "Setting probability of %s (%s) to -1 since less than %s that of top ranking hit %s, %0.4f (%s)\n", 
							$ranked_spectra[$j]->[0], $ranked_spectra[$j]->[2], $min_factor, 
							$ranked_spectra[0]->[0], $ranked_spectra[0]->[1], $ranked_spectra[0]->[2] if(0);
						$tot_rejected++;
					}
					$k = @ranked_spectra;
				}
			}
		}
	}
	printf STDERR "Set probabilities of %d results to -1 since they were less than %s times the top ranking hit with probability %s or more\n", 
		$tot_rejected, $min_factor, $min_top_prob;
}
# update distributions and probabilities
my $maxdiff = 0.001;
my $prior = 0.01; 
my $max_iters = $MAX_ITERS; #20;
my $iter = 1;
while($iter < $max_iters && updateDistributions(\%posdists, \%negdists, \%probabilities, \%distvals, $maxdiff, $prior, $use_decoy_set_negdist)) {
	printf STDERR "Iteration %d....\n", $iter;
	printDistributions(\%posdists, \%negdists, *STDERR);
	updateProbabilities(\%posdists, \%negdists, \%distvals, \%probabilities, $maxdiff, $distvals{'decoy'});
	$iter++;
}

# now check the probabilities within each spectrum group
if($PAIRING eq 'Mango' && $penalize_alternative_xlinks) {
	foreach(sort {$a cmp $b} keys %result_groups) {
		my $tot = 0;
		my @peps = ();
		my @ranked_spectra = ();
		my @next = keys %{$result_groups{$_}};
		for(my $k = 0; $k < @next; $k++) {
			$tot += $probabilities{$next[$k]} if($probabilities{$next[$k]} > 0);
			push(@peps, $distvals{'peptide_pair'}->{$next[$k]}) if($probabilities{$next[$k]} >= 0.8);
			push(@ranked_spectra, [$next[$k], $probabilities{$next[$k]}, $distvals{'peptide_pair'}->{$next[$k]}]);
		}
		if($tot > 1.2) {
			# these can be reduced in the future
			@ranked_spectra = sort {$b->[1] <=> $a->[1]} @ranked_spectra; # sort by prob
			# For each peptide pair, keep only the top one	
			my $total_prob = $probabilities{$ranked_spectra[0]->[0]};
			for(my $k = 1; $k < @ranked_spectra; $k++) {
				if($ranked_spectra[$k]->[2] eq $ranked_spectra[0]->[2]) {
					$probabilities{$ranked_spectra[$k]->[0]} = -1;
					@ranked_spectra = @ranked_spectra[0 .. $k-1, $k+1 .. $#ranked_spectra];
					$k--;
				}
				else {
					$total_prob += $probabilities{$ranked_spectra[$k]->[0]} if($probabilities{$ranked_spectra[$k]->[0]} > 0);
				}
			}
			my $found = 0;
			if($total_prob > 1.2) {
				for(my $k = 0; $k < @ranked_spectra; $k++) {
					next if($probabilities{$ranked_spectra[$k]->[0]} <= 0);
					printf STDERR "Adjusting prob of $ranked_spectra[$k]->[0] from %0.4f to %0.4f with %s by factor of %0.2f....\n",
						$probabilities{$ranked_spectra[$k]->[0]}, $probabilities{$ranked_spectra[$k]->[0]} / $total_prob, $ranked_spectra[$k]->[2], $total_prob;
					$probabilities{$ranked_spectra[$k]->[0]} = sprintf "%0.4f", $probabilities{$ranked_spectra[$k]->[0]} / $total_prob;
					$found = 1;
					$adjusted{$ranked_spectra[$k]->[0]} = $probabilities{$ranked_spectra[$k]->[0]} if($OUTPUT_ADJUSTED);
				}
				printf STDERR "\n" if($found);
			}
		}
	}  
} # if penalize
elsif($PAIRING eq 'Mango' && $top1_alternative_xlinks) {
	foreach(sort {$a cmp $b} keys %result_groups) {
		my $tot = 0;
		my @peps = ();
		my @ranked_spectra = ();
		my @next = keys %{$result_groups{$_}};
		for(my $k = 0; $k < @next; $k++) {
			$tot += $probabilities{$next[$k]} if($probabilities{$next[$k]} > 0);
			push(@peps, $distvals{'peptide_pair'}->{$next[$k]}) if($probabilities{$next[$k]} >= 0.8);
			push(@ranked_spectra, [$next[$k], $probabilities{$next[$k]}, $distvals{'peptide_pair'}->{$next[$k]}, $distvals{'homopeptide'}->{$next[$k]}]);
		}
		if($tot > 1.2) {
			# these can be reduced in the future
			@ranked_spectra = sort {$b->[1] <=> $a->[1] or $a->[3] <=> $b->[3]} @ranked_spectra; # sort by prob and prefer non-homo
			# For each peptide pair, keep only the top one	
			my $total_prob = $probabilities{$ranked_spectra[0]->[0]};
			for(my $k = 1; $k < @ranked_spectra; $k++) {
				if($ranked_spectra[$k]->[2] eq $ranked_spectra[0]->[2]) {
					$probabilities{$ranked_spectra[$k]->[0]} = -1;
					@ranked_spectra = @ranked_spectra[0 .. $k-1, $k+1 .. $#ranked_spectra];
					$k--;
				}
			}
			my @found = (sprintf "Top rank prob of $ranked_spectra[0]->[0] %0.4f with %s", $probabilities{$ranked_spectra[0]->[0]},
				$ranked_spectra[0]->[2]);
			for(my $k = 1; $k < @ranked_spectra; $k++) {
				next if($probabilities{$ranked_spectra[$k]->[0]} <= 0.01);
				push(@found, sprintf "Adjusting prob of $ranked_spectra[$k]->[0] from %0.4f to 0 with %s....",
					$probabilities{$ranked_spectra[$k]->[0]}, $ranked_spectra[$k]->[2]);
				$probabilities{$ranked_spectra[$k]->[0]} = -2;
				$adjusted{$ranked_spectra[$k]->[0]} = '' if($OUTPUT_ADJUSTED);
			}
			printf STDERR "%s\n\n", join("\n", @found) if(@found > 1);
		}
	}
} # if penalize
elsif($PAIRING eq 'ReACT' && $LIGHT_HEAVY) {
	foreach(sort {$a cmp $b} keys %result_groups) {
		my $tot = 0;
		my @peps = ();
		my @ranked_spectra = ();
		my @next = keys %{$result_groups{$_}};
		for(my $k = 0; $k < @next; $k++) {
			$tot += $probabilities{$next[$k]} if($probabilities{$next[$k]} > 0);
			push(@peps, $distvals{'peptide_pair'}->{$next[$k]}) if($probabilities{$next[$k]} >= 0.8);
			push(@ranked_spectra, [$next[$k], $probabilities{$next[$k]}, $distvals{'peptide_pair'}->{$next[$k]}]);
		}
		if($tot > 1.2) {
			# these can be reduced in the future
			@ranked_spectra = sort {$b->[1] <=> $a->[1]} @ranked_spectra; # sort by prob
			# For each peptide pair, keep only the top one	
			my $total_prob = $probabilities{$ranked_spectra[0]->[0]};
			for(my $k = 1; $k < @ranked_spectra; $k++) {
				if($ranked_spectra[$k]->[2] eq $ranked_spectra[0]->[2]) {
					$probabilities{$ranked_spectra[$k]->[0]} = -1;
					@ranked_spectra = @ranked_spectra[0 .. $k-1, $k+1 .. $#ranked_spectra];
					$k--;
				}
				else {
					$total_prob += $probabilities{$ranked_spectra[$k]->[0]} if($probabilities{$ranked_spectra[$k]->[0]} > 0);
				}
			}
			my $found = 0;
			if($total_prob > 1.2) {
				for(my $k = 0; $k < @ranked_spectra; $k++) {
					next if($probabilities{$ranked_spectra[$k]->[0]} <= 0);
					printf STDERR "Adjusting prob of $ranked_spectra[$k]->[0] from %0.4f to %0.4f with %s by factor of %0.2f....\n",
						$probabilities{$ranked_spectra[$k]->[0]}, $probabilities{$ranked_spectra[$k]->[0]} / $total_prob, $ranked_spectra[$k]->[2], $total_prob;
					$probabilities{$ranked_spectra[$k]->[0]} = sprintf "%0.4f", $probabilities{$ranked_spectra[$k]->[0]} / $total_prob;
					$found = 1;
					$adjusted{$ranked_spectra[$k]->[0]} = $probabilities{$ranked_spectra[$k]->[0]} if($OUTPUT_ADJUSTED);
				}
				printf STDERR "\n" if($found);
			}
		}
	}
} # if penalize
if(! ($PAIRING eq "ReACT") && $CROSSLINK_PRIOR==1) {
	printf "Error: Found no incorrect results in dataset.\n";
	if($PAIRING eq "ReACT") {
		printf "Make sure your react2.xls files were generated using the react2csv -c9999 -F options.\n";
	}
	sendEmail("Error: Found no incorrect results in dataset.", "Make sure your react2.xls files were generated using the react2csv -c9999 -F options.");
	exit(1);
}

my @sorted = sort {$probabilities{$b} <=> $probabilities{$a}} keys %probabilities;

# a future group strategy might be to apply toponly when the runner up is a homo version of one of the top ranking peptides

my @ERRORPT_FDRS = (0, 0.0001, 0.0002, 0.0003, 0.0004, 0.0005, 0.0006, 0.0007, 0.0008, 0.0009, 0.001, 0.0015, 0.002, 0.0025, 0.003, 0.004, 0.005, 
	0.006, 0.007, 0.008, 0.009, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05, 0.075, 0.1);
my $ERRORPT_IND = 0;
my %ERRORPTS = ();
my %COMPOSITE_ERRORPTS = ();
my $COMPOSITE_ERRORPT_IND = 0;
my $COMPOSITE_FDR = 0;
my %COMPOSITE_ROCS = ();
my $COMPOSITE_ROC_IND = 0;
my @ROC_MINPROBS = (0.9999, 0.999, 0.990, 0.98, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05, 0);
my $ROC_IND = 0;
my %ROCS = ();
my %DECOY_ROCS = ();
my %DECOY_ERRORPTS = ();
my $prev_fdr = 0;
my %COMPOSITE_REFS = (); # peptide pair to the reference spectrum with highest prob

%seen_mango_scan_pps = ();
if($COMPUTE_COMPOSITE_FDRs) {
	foreach(@sorted) {
		if($PAIRING eq 'Mango' && /^(\S+)\_(\d\d\d)\_(\S+)$/) {
			my $next_entry = $1  . '.' . $3 . ":" . $PEP_PAIRS{$_};
			next if(exists $seen_mango_scan_pps{$next_entry});  # only keep the first (highest prob) instance cross-link for each mango spectrum assigned to a peptide parir
			$seen_mango_scan_pps{$next_entry}++;
		}
		my $nextprob = $probabilities{$_} > 0 ? $probabilities{$_} : 0;
	
 		if($OUTPUT_DECOYS || ! $distvals{'decoy'}->{$_}) {
			next if(! exists $PEP_PAIRS{$_});
			if(! exists $CROSSLINK_FDRs{$PEP_PAIRS{$_}}) {
				$CROSSLINK_FDRs{$PEP_PAIRS{$_}} = 1 - $nextprob;
			}
			else {
				$CROSSLINK_FDRs{$PEP_PAIRS{$_}} *= (1 - $nextprob); 
			}
		}
	}
	my @composite_sorted = sort {$CROSSLINK_FDRs{$a} <=> $CROSSLINK_FDRs{$b}} keys %CROSSLINK_FDRs;
	printf "HAVE %d composite sorted crosslinks with outer values %0.3f and %0.3f\n", scalar @composite_sorted,
		$CROSSLINK_FDRs{$composite_sorted[0]}, $CROSSLINK_FDRs{$composite_sorted[-1]};
	$COMPOSITE_FDR = 0;
	my $lastprob = 1;
	for(my $k = 0; $k < @composite_sorted; $k++) {
		die "Here with blank\n" if($composite_sorted[$k] eq '');
	
		my $nextprob = 1 - $CROSSLINK_FDRs{$composite_sorted[$k]};
		$COMPOSITE_FDR += $CROSSLINK_FDRs{$composite_sorted[$k]};
		while($k < @composite_sorted - 1 && 
				$CROSSLINK_FDRs{$composite_sorted[$k+1]} == $CROSSLINK_FDRs{$composite_sorted[$k]}) {
			$k++;
			$COMPOSITE_FDR += $CROSSLINK_FDRs{$composite_sorted[$k]};
		}
		my $next_comp_incorr = sprintf "%0.0f", $COMPOSITE_FDR;
		my $next_fdr = $COMPOSITE_FDR / ($k+1);
		$next_fdr = 0 if(! $OUTPUT_DECOYS && $next_fdr < 0); # could happen if have high prob decoy near top....
		$next_fdr = $prev_fdr if($next_fdr < $prev_fdr);
		while($COMPOSITE_ROC_IND < @ROC_MINPROBS && $nextprob < $ROC_MINPROBS[$COMPOSITE_ROC_IND]) {
			printf "Adding for ROC minprob thresh $ROC_MINPROBS[$COMPOSITE_ROC_IND] at $lastprob with FDR $prev_fdr....\n" if(0);

			$COMPOSITE_ROCS{$ROC_MINPROBS[$COMPOSITE_ROC_IND]} = [$k + 1 - $next_comp_incorr, $next_comp_incorr, $prev_fdr];
			$COMPOSITE_ROC_IND++;
		}
		while($COMPOSITE_ERRORPT_IND < @ERRORPT_FDRS && $COMPOSITE_FDR / ($k+1) > $ERRORPT_FDRS[$COMPOSITE_ERRORPT_IND]) {
			printf "HERE with composite FDR %0.3f for threshold %0.2f and index $k for error point %s\n", 
			$COMPOSITE_FDR / ($k+1), $lastprob, $ERRORPT_FDRS[$COMPOSITE_ERRORPT_IND] if(0);
			
			$COMPOSITE_ERRORPTS{$ERRORPT_FDRS[$COMPOSITE_ERRORPT_IND]} = [$lastprob, $k + 1 - $next_comp_incorr, $next_comp_incorr];
			$COMPOSITE_ERRORPT_IND++;
		}	
		$lastprob = $nextprob;
		$prev_fdr = $next_fdr;

	}
}
if($compute_decoy_fdr) {
	my $lastprob = 1;
	my $tot = 0;
	my $num_decoys = 0;
	my $num_decoy_decoys = 0;
	my $prev_fdr = 0;
	for(my $k = 0; $k < @sorted; $k++) {
		$tot++ if($OUTPUT_DECOYS || ! $distvals{'decoy'}->{$sorted[$k]});
		$num_decoys++ if($distvals{'decoy'}->{$sorted[$k]});
		$num_decoy_decoys++ if($distvals{'decoy-decoy'}->{$sorted[$k]});
		my $init = $k;
		while($k < @sorted - 1 && $probabilities{$sorted[$k+1]} == $probabilities{$sorted[$init]}) {
			$k++;
			$tot++ if($OUTPUT_DECOYS || ! $distvals{'decoy'}->{$sorted[$k]});
			$num_decoys++ if($distvals{'decoy'}->{$sorted[$k]});
			$num_decoy_decoys++ if($distvals{'decoy-decoy'}->{$sorted[$k]});
		}
		my $fdr = $num_decoys > 2 * $num_decoy_decoys ? ($num_decoys - 2 * $num_decoy_decoys) / $tot : 0;
		$fdr = $prev_fdr if($fdr < $prev_fdr);
		for(my $j = $init; $j <= $k; $j++) {
			$distvals{'decoy_prob_fdr'}->{$sorted[$j]} = sprintf "%0.4f", $fdr;
		}
		my $nextprob = $probabilities{$sorted[$init]} > 0 ? $probabilities{$sorted[$init]} : 0;
		if(! ($OUTPUT_PEPXML eq '')) {
			my $tot_incorr = $fdr * $tot;
			my $next_corr = $tot - $tot_incorr - $nextprob;
			$next_corr = 0 if($next_corr < 0);
			while($ERRORPT_IND < @ERRORPT_FDRS && $fdr > $ERRORPT_FDRS[$ERRORPT_IND]) {
				$DECOY_ERRORPTS{$ERRORPT_FDRS[$ERRORPT_IND]} = [$lastprob, $next_corr, $tot_incorr - 1 + $nextprob];
				$ERRORPT_IND++;
			}
			while($ROC_IND < @ROC_MINPROBS && $nextprob < $ROC_MINPROBS[$ROC_IND]) {
				$DECOY_ROCS{$ROC_MINPROBS[$ROC_IND]} = [$next_corr, $tot_incorr - 1 + $nextprob, $prev_fdr];
				$ROC_IND++;
			}
		}
		$prev_fdr = $fdr;
		$lastprob = $nextprob;
	} # next sorted
	if(! ($OUTPUT_PEPXML eq '')) {
		my $tot_incorr = $prev_fdr * $tot;
		my $next_corr = $tot - $tot_incorr - $lastprob;
		$next_corr = 0 if($next_corr < 0);
		if($ROC_IND == $#ROC_MINPROBS && $lastprob == 0) {
			my $tot_incorr = $prev_fdr * $tot;
			my $next_corr = $tot - $tot_incorr ;
			$next_corr = 0 if($next_corr < 0);
			$DECOY_ROCS{$ROC_MINPROBS[$ROC_IND]} = [$next_corr, $tot_incorr, $prev_fdr];
		}
	}
	
	
}

my %FDRs = ();
my $tot = 0;
my $tot_incorr = 0;
$ROC_IND = 0;
$ERRORPT_IND = 0;
my %SEEN_COMPOSITES = ();

open(OUT, ">$outfile") or die "cannot write to $outfile $!\n";

printf OUT "%s\t%s%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%s%s", 'probability', 'spectrum', (! ($PAIRING eq 'Mango') ? '' : "\tmango_query"), 'peptide1', 'peptide2', 
	'protein1', 'protein2', 'charge1', 'charge2', 'probability1', 'probability2', 'peptide1_mass', 'peptide2_mass', 'parent_charge', 'parent_mass', 
	$PAIRING eq 'ReACT' ? "\tparent_scan" : "", $COMPUTE_COMPOSITE_FDRs ? "\tcomposite_probability\tcomposite_id" : "";

my @distributions_used = sort {exists $posdists{$a} <=> exists $posdists{$b} or $a cmp $b} keys %distvals;	
foreach(@distributions_used) {
	printf OUT "\t%s", $_;
}
printf OUT "\t%s\t%s\n", 'fdr', 'orig_nsx';
my $prev_fdr = 0;
my $lastprob = 1;
my $tot_non_redundant = 0; # how many composite correct 
my $tot_non_redundant_intra = 0; # how many of those are intra-links
my %observed_stumps = (); # for LIGHT_HEAVY, make sure have at least 2
foreach(@sorted) {
	next if(! $OUTPUT_DECOYS && $distvals{'decoy'}->{$_});
	
	my $next_parent = undef;
	my $base = '';
	
	if(/^(\S+)\_\d\d\d\_/) {
		$base = $1;
	}
	else {
		die "Error, no base in $_\n";
	}
	
	
	my $scan = -1;
	if(/\.0*([^\.]+)$/) {
		$scan = $1;
		while(! ($PAIRING eq 'Mango') && length $scan < 5) {
			$scan = '0' . $scan;
		}
		# compute the new scores
		if(exists $scan_info{$base}->{$scan}) {
			$next_parent = $scan_info{$base}->{$scan};
		}
		else {
			die "problem finding parent for $scan and $_\n";
		}
	}
	next if($probabilities{$_} eq '');
	my $nextprob = $probabilities{$_} > 0 ? $probabilities{$_} : 0;
	
 	if($OUTPUT_DECOYS || ! $distvals{'decoy'}->{$_}) {
		$tot++;
		$tot_incorr += (1 - $nextprob);
		
		my $next_fdr = $tot_incorr / $tot if($tot > 0);
		$next_fdr = 0 if(! $OUTPUT_DECOYS && $next_fdr < 0); # could happen if have high prob decoy near top....

		my $next_corr = $tot - $tot_incorr - $nextprob;
		$next_corr = 0 if($next_corr < 0);
		if(! ($OUTPUT_PEPXML eq '')) {
			while($ERRORPT_IND < @ERRORPT_FDRS && $next_fdr > $ERRORPT_FDRS[$ERRORPT_IND]) {
				$ERRORPTS{$ERRORPT_FDRS[$ERRORPT_IND]} = [$lastprob, $next_corr, $tot_incorr - 1 + $nextprob];
				$ERRORPT_IND++;
			}
			while($ROC_IND < @ROC_MINPROBS && $nextprob < $ROC_MINPROBS[$ROC_IND]) {
				$ROCS{$ROC_MINPROBS[$ROC_IND]} = [$next_corr, $tot_incorr - 1 + $nextprob, $prev_fdr];
				$ROC_IND++;
			}
		}
		
		$lastprob = $nextprob;
		$next_fdr = $prev_fdr if($next_fdr < $prev_fdr);

		printf OUT "%0.4f\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%0.4f\t%0.4f\t%4f\t%4f\t%d\t%0.4f%s%s%s", #\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", 
		$probabilities{$_}, (! ($PAIRING eq 'Mango') && /^(\S+)\_000\_(\S+)$/ ? $1 . '.' . $2 . '.' . $results{$_}->{1}->[3] : 
			$PAIRING eq 'Mango' &&  /^(\S+)\_(\d\d\d)\_(\S+)$/ ? $1  . '.' . $3 . '.' . $results{$_}->{1}->[3] . "\t" . ($2+1) : $_), 
		$results{$_}->{1}->[2], $results{$_}->{2}->[2],
		$results{$_}->{1}->[1], $results{$_}->{2}->[1], $results{$_}->{1}->[3], $results{$_}->{2}->[3], 
		$results{$_}->{1}->[0], $results{$_}->{2}->[0], 
		$results{$_}->{1}->[4], $results{$_}->{2}->[4], 
		$next_parent->{'charge'}, $next_parent->{'mass'}, $PAIRING eq 'ReACT' ? "\t" . pop @{$results{$_}->{1}->[-1]} : "",
		($COMPUTE_COMPOSITE_FDRs ? (sprintf "\t%s", exists $PEP_PAIRS{$_} && ! exists $SEEN_COMPOSITES{$PEP_PAIRS{$_}} ? (sprintf "%0.2f", 1 - $CROSSLINK_FDRs{$PEP_PAIRS{$_}}) : "") : ""),
		($COMPUTE_COMPOSITE_FDRs ? (sprintf "\t%s", exists $PEP_PAIRS{$_} ? $PEP_PAIRS{$_} : "") : "");
			$tot_non_redundant += 1 - $CROSSLINK_FDRs{$PEP_PAIRS{$_}} if(exists $PEP_PAIRS{$_} && ! exists $SEEN_COMPOSITES{$PEP_PAIRS{$_}}); 
			$tot_non_redundant_intra += 1 - $CROSSLINK_FDRs{$PEP_PAIRS{$_}} if(exists $PEP_PAIRS{$_} && ! exists $SEEN_COMPOSITES{$PEP_PAIRS{$_}} && $distvals{'intra'}->{$_} eq '1');
		$SEEN_COMPOSITES{$PEP_PAIRS{$_}}++ if($COMPUTE_COMPOSITE_FDRs && exists $PEP_PAIRS{$_});
		for(my $k = 0; $k < @distributions_used; $k++) {
			printf OUT "\t%s", $distvals{$distributions_used[$k]}->{$_};
		}
		printf OUT "\t%0.4f\t%s\n", $next_fdr, $NSP_VALS{$_};
		$prev_fdr = $next_fdr;
		$FDRs{$_} = sprintf "%0.3f", $next_fdr;
		if($nextprob >= $MIN_XL_PROB_4_PROTPROPH) {
			my $stripped1 = stripPeptide($results{$_}->{1}->[2]);
			$PROTPROPH_INFO->{$stripped1} = {} if(! exists $PROTPROPH_INFO->{$stripped1});
			my $stripped2 = stripPeptide($results{$_}->{2}->[2]);
			$PROTPROPH_INFO->{$stripped2} = {} if(! exists $PROTPROPH_INFO->{$stripped2});
			$PEP_MAXPROBS{$stripped1} = $nextprob if(! exists $PEP_MAXPROBS{$stripped1} || $nextprob > $PEP_MAXPROBS{$stripped1});
			$PEP_MAXPROBS{$stripped2} = $nextprob if(! exists $PEP_MAXPROBS{$stripped2} || $nextprob > $PEP_MAXPROBS{$stripped2});
		}
		
		if($LIGHT_HEAVY && $probabilities{$_} >= 0.9) {
			for(my $k = 0; $k < @CROSSLINK_MODS; $k++) {
				my @next = keys %{$CROSSLINK_MODS[$k]};
				for(my $j = 0; $j < @next; $j++) {
					if($next[$j] eq 'K' && ! exists $observed_stumps{$CROSSLINK_MODS[$k]->{$next[$j]}} && 
						$results{$_}->{1}->[2]=~ /$next[$j]\[$CROSSLINK_MODS[$k]->{$next[$j]}\]/) {
						$observed_stumps{$CROSSLINK_MODS[$k]->{$next[$j]}}++;
						#printf "Added stump $CROSSLINK_MODS[$k]->{$next[$j]} for cross-link\n";
						
					}
					#printf "$next[$j]: $CROSSLINK_MODS[$k]->{$next[$j]}\n";
				}
			}
		}
	}
	
	delete $adjusted{$_} if(exists $adjusted{$_});

	last if($MAX_FDR < 1 && $tot_incorr / $tot > $MAX_FDR);
	
}
if($LIGHT_HEAVY && scalar keys %observed_stumps < 2) {
	printf "Error: cross-linked peptides only with stump modification ".join(",", keys %observed_stumps)." observed\n";
	printf "Please check heavy and light mod settings or seek help\n";
	exit(1);
}
if($ROC_IND == $#ROC_MINPROBS && $lastprob == 0) {
	my $tot_incorr = $prev_fdr * $tot;
	my $next_corr = $tot - $tot_incorr ;
	$next_corr = 0 if($next_corr < 0);
	$ROCS{$ROC_MINPROBS[$ROC_IND]} = [$next_corr, $tot_incorr, $prev_fdr];
}

my $total_correct = $tot - $tot_incorr;
foreach(keys %prot_crosslinks) {
	my @next = keys %{$prot_crosslinks{$_}};
	my $tot = 0;
	for(my $k = 0; $k < @next; $k++) {
		$tot += $prot_crosslinks{$_}->{$next[$k]};
	}
	$protpair_instances{$_} = $tot;
}
printf "Computed protpair instances for %d prot pairs\n", scalar keys %protpair_instances;

# now can print out remaining adjusted guys....
if(scalar keys %adjusted > 0) {
	my @sorted = keys %adjusted;
	foreach(@sorted) {
		my $next_parent = undef;
		my $base = '';
		next if($adjusted{$_} eq '');
		if(/^(\S+)\_\d\d\d\_/) {
			$base = $1;
		}
		else {
			die "Error, no base in $_\n";
		}
		my $scan = -1;
		if(/\.0*([^\.]+)$/) {
			$scan = $1;
			while(! ($PAIRING eq 'Mango') && length $scan < 5) {
				$scan = '0' . $scan;
			}
			# compute the new scores
			if(exists $scan_info{$base}->{$scan}) {
				$next_parent = $scan_info{$base}->{$scan};
			}
			else {
				die "problem finding parent for $scan and $_\n";
			}
		}
		my $nextprob = $probabilities{$_} > 0 ? $probabilities{$_} : 0;

		$distvals{'decoy_prob_fdr'}->{$_} = ""; # don't want to report this

		printf OUT "%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%0.4f\t%0.4f\t%4f\t%4f\t%d\t%0.4f\t%s%s", #\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", 
			$adjusted{$_}, (! ($PAIRING eq 'Mango') && /^(\S+)\_000\_(\S+)$/ ? $1 . '.' . $2 . '.' . $results{$_}->{1}->[3] : 
			$PAIRING eq 'Mango' &&  /^(\S+)\_(\d\d\d)\_(\S+)$/ ? $1  . '.' . $3 . '.' . $results{$_}->{1}->[3] . "\t" . ($2+1) : $_), 
		$results{$_}->{1}->[2], $results{$_}->{2}->[2],
		$results{$_}->{1}->[1], $results{$_}->{2}->[1], $results{$_}->{1}->[3], $results{$_}->{2}->[3], 
		$results{$_}->{1}->[0], $results{$_}->{2}->[0], 
		$results{$_}->{1}->[4], $results{$_}->{2}->[4], 
		$next_parent->{'charge'}, $next_parent->{'mass'}, $PAIRING eq 'ReACT' ? "\t" . pop @{$results{$_}->{1}->[-1]} : "",

		$COMPUTE_COMPOSITE_FDRs ? (sprintf "\t%s", exists $PEP_PAIRS{$_} && ! exists $SEEN_COMPOSITES{$PEP_PAIRS{$_}} ? (sprintf "%0.2f\t%s", 1 - $CROSSLINK_FDRs{$PEP_PAIRS{$_}}, $PEP_PAIRS{$_}) : "\t") : "";
		$tot_non_redundant += 1 - $CROSSLINK_FDRs{$PEP_PAIRS{$_}} if(exists $PEP_PAIRS{$_} && ! exists $SEEN_COMPOSITES{$PEP_PAIRS{$_}}); 

		$SEEN_COMPOSITES{$PEP_PAIRS{$_}}++ if($COMPUTE_COMPOSITE_FDRs && exists $PEP_PAIRS{$_});

		for(my $k = 0; $k < @distributions_used; $k++) {
			printf OUT "\t%s", $distvals{$distributions_used[$k]}->{$_};
		}
		printf OUT "\t%s\n", '' if($OUTPUT_DECOYS || ! $distvals{'decoy'}->{$_});
		if($USE_PROTPROPH_PROTS && $nextprob >= $MIN_XL_PROB_4_PROTPROPH) {
			my $stripped1 = stripPeptide($results{$_}->{1}->[2]);
			$PROTPROPH_INFO->{$stripped1} = {} if(! exists $PROTPROPH_INFO->{$stripped1});
			my $stripped2 = stripPeptide($results{$_}->{2}->[2]);
			$PROTPROPH_INFO->{$stripped2} = {} if(! exists $PROTPROPH_INFO->{$stripped2});
		}	
	}
}
close(OUT);



if(! ($OUTPUT_PEPXML eq '')) {
	if(scalar keys %{$PROTPROPH_INFO} > 0) {
		if($file =~ /^(\S+)\.pep.xml$/) {
			my $protproph = $1 . ".minprob0.01.prot.xml";
			if(! -e $protproph) {
				my $command = "ProteinProphet $file $protproph MINPROB0.01";
				printf STDERR "%s\n", $command;
				system($command);
				if(! -e $protproph) {
					die "Error running ProteinProphet with $command: $protproph no present\n";
				}
			}
			attributePeptides2Proteins($protproph);
			setGeneNamesFromUniprot() if($SET_GENES_FROM_UNIPROT);
		}
		else {
			die "Error with input file $file\n";
		}
	} # if pursuing proteinprophet results


	$tot = 0;
	
	# print out crosslink stump mod masses
	my @stump_modmasses = ();
	for(my $z = 0; $z < @CROSSLINK_STUMP_MODS; $z++) {
		my @nextstumps = keys %{$CROSSLINK_STUMP_MODS[$z]};
		for(my $j = 0; $j < @nextstumps; $j++) {
			if(exists $FULL_STUMP_MASSES{$nextstumps[$j] . ":" . $CROSSLINK_STUMP_MODS[$z]->{$nextstumps[$j]}}) {
				printf "Subsituting stump mod %s for %s in pepXML\n",  $FULL_STUMP_MASSES{$nextstumps[$j] . ":" . $CROSSLINK_STUMP_MODS[$z]->{$nextstumps[$j]}},
					$nextstumps[$j] . ":" . $CROSSLINK_STUMP_MODS[$z]->{$nextstumps[$j]};
				push(@stump_modmasses, $FULL_STUMP_MASSES{$nextstumps[$j] . ":" . $CROSSLINK_STUMP_MODS[$z]->{$nextstumps[$j]}});
			}
			else {
				push(@stump_modmasses, $nextstumps[$j] . ":" . $CROSSLINK_STUMP_MODS[$z]->{$nextstumps[$j]});
			}
		}
	}
	
	open(PEPXML, ">$OUTPUT_PEPXML.xml") or die "cannot write to file $OUTPUT_PEPXML.xml $!\n";
	printf PEPXML "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
	printf PEPXML "<?xml-stylesheet type=\"text/xsl\" href=\"%s.xsl\"?>\n", $OUTPUT_PEPXML; # must take off extension first
	printf PEPXML "<msms_pipeline_analysis xsi:schemaLocation=\"http://regis-web.systemsbiology.net/pepXML /usr/local/tppschema/pepXML_v120.xsd\" summary_xml=\"%s.xml\" date=\"%s\">\n",
		$OUTPUT_PEPXML, $timestamp;
	printf PEPXML "<analysis_summary analysis=\"peptideprophet\" time=\"%s\">\n", $timestamp;
	printf PEPXML "<peptideprophet_summary version=\"XLinkProphet\" author=\"AKeller\@UW\" type=\"unlinked\" min_prob=\"0.00\" options=\"%s\" xlinker_stump_modmasses=\"%s\" est_tot_num_correct=\"%0.1f\">\n",
		join(" ", @ARGV[1 .. $#ARGV]), join(",", @stump_modmasses), $total_correct;
	foreach(@INPUT_FILES) {
		printf PEPXML "%s\n", $_;
	}
	# now the roc stuff
	printf PEPXML "<roc_error_data charge=\"all\">\n";
	foreach(sort {$b <=> $a} keys %ROCS) {
		printf PEPXML "<roc_data_point min_prob=\"%0.4f\" sensitivity=\"%0.4f\" error=\"%0.4f\" num_corr=\"%0.0f\" num_incorr=\"%0.0f\"/>\n",
		$_, $ROCS{$_}->[0]/$total_correct, $ROCS{$_}->[2], $ROCS{$_}->[0], $ROCS{$_}->[1] || 0;
	}
	foreach(sort {$a <=> $b} keys %ERRORPTS) {
		printf PEPXML "<error_point error=\"%0.4f\" min_prob=\"%0.4f\" num_corr=\"%0.0f\" num_incorr=\"%0.0f\"/>\n",
			$_, $ERRORPTS{$_}->[0], $ERRORPTS{$_}->[1], $ERRORPTS{$_}->[2];
	}
	printf PEPXML "</roc_error_data>\n";
	if($compute_decoy_fdr) {
		printf PEPXML "<roc_error_data charge=\"decoy\">\n";
		foreach(sort {$b <=> $a} keys %DECOY_ROCS) {
			printf PEPXML "<roc_data_point min_prob=\"%0.4f\" sensitivity=\"%0.4f\" error=\"%0.4f\" num_corr=\"%0.0f\" num_incorr=\"%0.0f\"/>\n",
			$_, $DECOY_ROCS{$_}->[0]/$total_correct, $DECOY_ROCS{$_}->[2], $DECOY_ROCS{$_}->[0], $DECOY_ROCS{$_}->[1] || 0;
		}
		foreach(sort {$a <=> $b} keys %DECOY_ERRORPTS) {
			printf PEPXML "<error_point error=\"%0.4f\" min_prob=\"%0.4f\" num_corr=\"%0.0f\" num_incorr=\"%0.0f\"/>\n",
				$_, $DECOY_ERRORPTS{$_}->[0], $DECOY_ERRORPTS{$_}->[1], $DECOY_ERRORPTS{$_}->[2];
		}
		printf PEPXML "</roc_error_data>\n";
	}
	if($COMPUTE_COMPOSITE_FDRs) {
		my $COMP_TOT_CORR = sprintf "%0.0f", (scalar keys %CROSSLINK_FDRs) - $COMPOSITE_FDR;
		printf PEPXML "<roc_error_data charge=\"composite\">\n";
		foreach(sort {$b <=> $a} keys %COMPOSITE_ROCS) {
			printf PEPXML "<roc_data_point min_prob=\"%0.4f\" sensitivity=\"%0.4f\" error=\"%0.4f\" num_corr=\"%0.0f\" num_incorr=\"%0.0f\"/>\n",
			$_, $COMPOSITE_ROCS{$_}->[0]/$COMP_TOT_CORR, $COMPOSITE_ROCS{$_}->[2], $COMPOSITE_ROCS{$_}->[0], $COMPOSITE_ROCS{$_}->[1] || 0;
		}
		foreach(sort {$a <=> $b} keys %COMPOSITE_ERRORPTS) {
			printf PEPXML "<error_point error=\"%0.4f\" min_prob=\"%0.4f\" num_corr=\"%0.0f\" num_incorr=\"%0.0f\"/>\n",
				$_, $COMPOSITE_ERRORPTS{$_}->[0], $COMPOSITE_ERRORPTS{$_}->[1], $COMPOSITE_ERRORPTS{$_}->[2];
		}
		printf PEPXML "</roc_error_data>\n";
	}
	printf PEPXML "<mixture_model precursor_ion_charge=\"1\" prior_probability=\"%0.3f\" est_tot_correct=\"%0.1f\" tot_num_spectra=\"%d\" num_iterations=\"%d\">\n",
		$total_correct / (scalar keys %results), $total_correct, scalar keys %results, $iter - 1;
	foreach(keys %posdists) {
		printf PEPXML "<mixturemodel_distribution name=\"%s\">\n", $_;
		printf PEPXML "<posmodel_distribution>\n";
		my @next = sort {$a <=> $b} keys %{$posdists{$_}};
		for(my $k = 0; $k < @next; $k++) {
			printf PEPXML "<parameter name=\"%s\" value=\"%0.3f\"/>\n", exists $MODEL_BINS{$_} ? $next[$k] * $MODEL_BINS{$_} : $next[$k], 
				$posdists{$_}->{$next[$k]};
		}
		printf PEPXML "</posmodel_distribution>\n";
		printf PEPXML "<negmodel_distribution>\n";
		for(my $k = 0; $k < @next; $k++) {
			printf PEPXML "<parameter name=\"%s\" value=\"%0.3f\"/>\n", exists $MODEL_BINS{$_} ? $next[$k] * $MODEL_BINS{$_} : $next[$k], 
				$negdists{$_}->{$next[$k]};
		}
		printf PEPXML "</negmodel_distribution>\n";
		printf PEPXML "</mixturemodel_distribution>\n";
	}
	printf PEPXML "</mixture_model>\n";
	
	
	printf PEPXML "</peptideprophet_summary>\n";
	printf PEPXML "</analysis_summary>\n";
	printf PEPXML "<dataset_derivation generation_no=\"0\"/>\n";
	foreach(@pepxml_header) {
		printf PEPXML "%s\n", $_;
	}
	for(my $k = 0; $k < @RUNS; $k++) {
		foreach(@{$RUNS[$k]}) {
			printf PEPXML "%s\n", $_;
		}
		foreach(@{$RUN_SPECTRA[$k]}) {
			next if(! $OUTPUT_DECOYS && $distvals{'decoy'}->{$_});
			next if(! exists $results{$_} || ! (ref($results{$_}->{1}->[-1]) eq 'ARRAY') || ! (ref($results{$_}->{2}->[-1]) eq 'ARRAY'));
			next if($probabilities{$_} eq '');
		
			my $next_parent = undef;
			my $base = '';
	
			my $query_no = -1;
			if(/^(\S+)\_(\d\d\d)\_/) {
				$base = $1;
				$query_no = sprintf "%d", $2 if($PAIRING eq 'Mango');
			}
			else {
				die "Error, no base in $_\n";
			}
			my $scan = -1;
			my $ms2scan = "";
			if(/\.0*([^\.]+)$/) {
				$scan = $1;
				while(! ($PAIRING eq 'Mango') && length $scan < 5) {
					$scan = '0' . $scan;
				}
				# compute the new scores
				if(exists $scan_info{$base}->{$scan}) {
					$next_parent = $scan_info{$base}->{$scan};
				}
				else {
					die "problem finding parent for $scan and $_\n";
				}
			}
			if($PAIRING eq 'ReACT') {
				if(! exists $next_parent->{'parent_scan'}) {
					exit("ERROR, no parent scan for $scan\n");
				}
				$ms2scan = " ms2scan=\"".$next_parent->{'parent_scan'}."\"";
			}

			my $nextprob = $probabilities{$_} > 0 ? $probabilities{$_} : 0;
 			if($OUTPUT_DECOYS || ! $distvals{'decoy'}->{$_}) {

				$tot++;
				my $spectrum = $_;
				$spectrum = $1 if($spectrum =~ /$DIR_DIVISOR_REGX([^$DIR_DIVISOR_REGX]+)$/); # clear off full path (added to distinguish among same raw files searched multiple times)
				$spectrum = $1 . '.' . $3 . '.' . $results{$_}->{1}->[3] if($spectrum =~ /(\S+)(\_\d\d\d\_)(\d+\.\d+)/);
				# when have multiple results with same spectrum, need to add label so iprophet can distinguish
				my $spectrum_label = "";
				if($PAIRING eq 'Mango') {
					$MANGO_SPECS{$spectrum}++;
					$spectrum_label = sprintf " spectrum_label=\"%s%d\"", $SPEC_LABEL, $MANGO_SPECS{$spectrum};
				}
				elsif(! ($SPEC_LABEL eq '')) {
					$spectrum_label = sprintf " spectrum_label=\"%s\"", $SPEC_LABEL;
				}

				printf PEPXML "<spectrum_query spectrum=\"%s\" start_scan=\"%d\" end_scan=\"%d\"%s precursor_neutral_mass=\"%s\" assumed_charge=\"%d\" index=\"%d\">\n",
					(! ($PAIRING eq 'Mango') && /^(\S+)\_000(\.\S+)$/ ? $1 . $2 . '.' . $results{$_}->{1}->[3] : $spectrum), $scan, $scan, $spectrum_label,
						$next_parent->{'mass'}, $next_parent->{'charge'}, $tot;
				printf PEPXML " <search_result>\n";
				
				
				printf PEPXML " <search_hit hit_rank=\"1\" peptide=\"-\" peptide_prev_aa=\"-\" peptide_next_aa=\"-\" protein=\"-\" num_tot_proteins=\"1\" calc_neutral_pep_mass=\"1000\" massdiff=\"0\"%s xlink_type=\"%s\">\n",
						$ms2scan, "xl";  # need mass and massdiff here
				
				my $strip1 = $results{$_}->{1}->[2];
				my $strip2 = $results{$_}->{2}->[2];
				if($nextprob >= $MIN_XL_PROB_4_PROTPROPH) {
					$strip1 = stripPeptide($results{$_}->{1}->[2]);
					$strip2 = stripPeptide($results{$_}->{2}->[2]);
				}
				
				# here use overall protein's nsp to estimate instead of saying 0, for proteins with prob above a minimum?
				if(exists $PROTPROPH_INFO->{$strip1} && scalar keys %{$PROTPROPH_INFO->{$strip1}} == 0) {
					my $used = 0;
					for(my $k = 0; $k < @{$results{$_}->{1}->[-1]}-1; $k++) {
						if($results{$_}->{1}->[-1]->[$k] =~ /\<linked\_peptide.*?protein=\"(\S+)\".*$/) {
							my $nextprot = $1;
							if($OUTPUT_DECOYS || $nextprot !~ /^$DECOY_PREFIX/) {
								$PROTPROPH_INFO->{$strip1}->{$1} = [1, @{$results{$_}->{1}->[-1]} == 2 ? "Y" : "N", 0];
								$used++;
							}
						}
						elsif($results{$_}->{1}->[-1]->[$k] =~ /\<alternative\_protein protein=\"(\S+)\" protein_descr=\".*?peptide_next_aa=\"\S\"/) {
							my $nextprot = $1;
							if($OUTPUT_DECOYS || $nextprot !~ /^$DECOY_PREFIX/) {
								$PROTPROPH_INFO->{$strip1}->{$1} = [1, "N", 0];
								$used++;
							}
						}
					}
					my @used = keys %{$PROTPROPH_INFO->{$strip1}};
					for(my $k = 0; $k < @used; $k++) {
						$PROTPROPH_INFO->{$strip1}->{$used[$k]}->[0]/=@used;
					}
				}
				if(exists $PROTPROPH_INFO->{$strip2} && scalar keys %{$PROTPROPH_INFO->{$strip2}} == 0) {
					my $used = 0;
					for(my $k = 0; $k < @{$results{$_}->{2}->[-1]}-1; $k++) {
						if($results{$_}->{2}->[-1]->[$k] =~ /\<linked\_peptide.*?protein=\"(\S+)\".*$/) {
							my $nextprot = $1;
							if($OUTPUT_DECOYS || $nextprot !~ /^$DECOY_PREFIX/) {
								$PROTPROPH_INFO->{$strip2}->{$1} = [1, @{$results{$_}->{2}->[-1]} == 2 ? "Y" : "N", 0];
								$used++;
							}
						}
						elsif($results{$_}->{2}->[-1]->[$k] =~ /\<alternative\_protein protein=\"(\S+)\" protein_descr=\".*?peptide_next_aa=\"\S\"/) {
							my $nextprot = $1;
							if($OUTPUT_DECOYS || $nextprot !~ /^$DECOY_PREFIX/) {
								$PROTPROPH_INFO->{$strip2}->{$1} = [1, "N", 0];
								$used++;
							}
						}
					}
					my @used = keys %{$PROTPROPH_INFO->{$strip2}};
					for(my $k = 0; $k < @used; $k++) {
						$PROTPROPH_INFO->{$strip2}->{$used[$k]}->[0]/=@used;
					}
				}
				
				my $protInfo = $nextprob >= $MIN_XL_PROB_4_PROTPROPH && exists $PROTPROPH_INFO->{$strip1} && scalar keys %{$PROTPROPH_INFO->{$strip1}} > 0 && 
					exists $PROTPROPH_INFO->{$strip2} && scalar keys %{$PROTPROPH_INFO->{$strip2}} > 0
					? getProtInfo($strip1, 
					$strip2, $PROTPROPH_INFO) : [];
				# if has data, have to substitute into the original results........
				if(@{$protInfo} > 0) {
					my %protData = (); # store description. prevaa, nextaa, num tol term
					for(my $k = 0; $k < @{$results{$_}->{1}->[-1]}; $k++) {
						if($results{$_}->{1}->[-1]->[$k] =~ /\<linked\_peptide.*?peptide\_prev\_aa\=\"(\S)\" peptide\_next\_aa\=\"(\S)\" protein=\"(\S+)\".*?num\_tol\_term\=\"(\d)\".*?protein\_descr\=\"(.*?)\"/) {
							$protData{$3} = [$5, $1, $2, $4];
						}
						elsif($results{$_}->{1}->[-1]->[$k] =~ /\<alternative\_protein protein\=\"(\S+)\" protein\_descr=\"(.*?)\" num\_tol\_term\=\"(\d)\" peptide\_prev\_aa\=\"(\S)\" peptide\_next\_aa\=\"(\S)\"/) {
							$protData{$1} = [$2, $4, $5, $3];
						}
					}
					my $index = 0;
					for(my $k = 0; $k < @{$results{$_}->{1}->[-1]}-1; $k++) {
						if($results{$_}->{1}->[-1]->[$k] =~ /(\<linked\_peptide.*?peptide\_prev\_aa\=\")\S(\" peptide\_next\_aa\=\")\S(\" protein=\")\S+(\".*?num\_tol\_term\=\")\d(\".*?protein\_descr\=\").*?(\".*)$/) {
							$results{$_}->{1}->[-1]->[$k] = $1 . $protData{$protInfo->[0]}->[1] . $2 . $protData{$protInfo->[0]}->[2] . $3 . $protInfo->[0] . $4 . 
								$protData{$protInfo->[0]}->[3] . $5 . $protData{$protInfo->[0]}->[0] . (sprintf "\" gene=\"%s\" weight=\"%0.2f\" nsp=\"%0.1f", $GENE_NAMES{$protInfo->[0]}, $PROTPROPH_INFO->{$strip1}->{$protInfo->[0]}->[0], $PROTPROPH_INFO->{$strip1}->{$protInfo->[0]}->[2]) . $6;
						}
						elsif($results{$_}->{1}->[-1]->[$k] =~ /\<alternative\_protein/) {
							if($index > $#{$protInfo->[1]}) {
								@{$results{$_}->{1}->[-1]} = @{$results{$_}->{1}->[-1]}[0 .. $k-1, $k+1 .. $#{$results{$_}->{1}->[-1]}];
								$k--;
							}
							else {
								my @next = split(":", $protInfo->[1]->[$index]);
								$results{$_}->{1}->[-1]->[$k] = sprintf "<alternative_protein protein=\"%s\" protein_descr=\"%s\" num_tol_term=\"%d\" peptide_prev_aa=\"%s\" peptide_next_aa=\"%s\" gene=\"%s\" weight=\"%0.2f\" nsp=\"%0.1f\"/>",
									$next[0], $protData{$next[0]}->[0], $protData{$next[0]}->[3], $protData{$next[0]}->[1], $protData{$next[0]}->[2], $GENE_NAMES{$next[0]}, $next[1],
									$PROTPROPH_INFO->{$strip1}->{$next[0]}->[2];
								$index++;
							}
						}
						elsif(/<modification_info/) {
							while($k > 0 && $index < @{$protInfo->[1]}) {
								my @next = split(":", $protInfo->[1]->[$index]);
								$results{$_}->{1}->[-1]->[$k-1] .= sprintf "\n<alternative_protein protein=\"%s\" protein_descr=\"%s\" num_tol_term=\"%d\" peptide_prev_aa=\"%s\" peptide_next_aa=\"%s\" gene=\"%s\" weight=\"%0.2f\" nsp=\"%0.1f\"/>",
									$next[0], $protData{$next[0]}->[0], $protData{$next[0]}->[3], $protData{$next[0]}->[1], $protData{$next[0]}->[2], $GENE_NAMES{$next[0]}, $next[1],
									$PROTPROPH_INFO->{$strip1}->{$next[0]}->[2];
								$index++;							
							}
						}	
					}
					%protData = (); # store description. prevaa, nextaa, num tol term
					for(my $k = 0; $k < @{$results{$_}->{2}->[-1]}; $k++) {
						if($results{$_}->{2}->[-1]->[$k] =~ /\<linked\_peptide.*?peptide\_prev\_aa\=\"(\S)\" peptide\_next\_aa\=\"(\S)\" protein=\"(\S+)\".*?num\_tol\_term\=\"(\d)\".*?protein\_descr\=\"(.*?)\"/) {
							$protData{$3} = [$5, $1, $2, $4];
						}
						elsif($results{$_}->{2}->[-1]->[$k] =~ /\<alternative\_protein protein\=\"(\S+)\" protein\_descr=\"(.*?)\" num\_tol\_term\=\"(\d)\" peptide\_prev\_aa\=\"(\S)\" peptide\_next\_aa\=\"(\S)\"/) {
							$protData{$1} = [$2, $4, $5, $3];
						}
					}
					my $index = 0;
					for(my $k = 0; $k < @{$results{$_}->{2}->[-1]}; $k++) {
						if($results{$_}->{2}->[-1]->[$k] =~ /(\<linked\_peptide.*?peptide\_prev\_aa\=\")\S(\" peptide\_next\_aa\=\")\S(\" protein=\")\S+(\".*?num\_tol\_term\=\")\d(\".*?protein\_descr\=\").*?(\".*)$/) {
							$results{$_}->{2}->[-1]->[$k] = $1 . $protData{$protInfo->[2]}->[1] . $2 . $protData{$protInfo->[2]}->[2] . $3 . $protInfo->[2] . $4 . 
								$protData{$protInfo->[2]}->[3] . $5 . $protData{$protInfo->[2]}->[0] . (sprintf "\" gene=\"%s\" weight=\"%0.2f\" nsp=\"%0.1f", $GENE_NAMES{$protInfo->[2]}, $PROTPROPH_INFO->{$strip2}->{$protInfo->[2]}->[0], $PROTPROPH_INFO->{$strip2}->{$protInfo->[2]}->[2]) . $6;
						}
						elsif($results{$_}->{2}->[-1]->[$k] =~ /\<alternative\_protein/) {
							if($index > $#{$protInfo->[3]}) {
								@{$results{$_}->{2}->[-1]} = @{$results{$_}->{2}->[-1]}[0 .. $k-1, $k+1 .. $#{$results{$_}->{2}->[-1]}];
								$k--;
							}
							else {
								my @next = split(":", $protInfo->[3]->[$index]);
								$results{$_}->{2}->[-1]->[$k] = sprintf "<alternative_protein protein=\"%s\" protein_descr=\"%s\" num_tol_term=\"%d\" peptide_prev_aa=\"%s\" peptide_next_aa=\"%s\" gene=\"%s\" weight=\"%0.2f\" nsp=\"%0.1f\"/>",
									$next[0], $protData{$next[0]}->[0], $protData{$next[0]}->[3], $protData{$next[0]}->[1], $protData{$next[0]}->[2], $GENE_NAMES{$next[0]}, $next[1],
									$PROTPROPH_INFO->{$strip2}->{$next[0]}->[2];
								$index++;
							}
						}
						elsif(/<modification_info/) {
							while($k > 0 && $index < @{$protInfo->[3]}) {
								my @next = split(":", $protInfo->[3]->[$index]);
								$results{$_}->{2}->[-1]->[$k-1] .= sprintf "\n<alternative_protein protein=\"%s\" protein_descr=\"%s\" num_tol_term=\"%d\" peptide_prev_aa=\"%s\" peptide_next_aa=\"%s\" gene=\"%s\" weight=\"%0.2f\" nsp=\"%0.1f\"/>",
									$next[0], $protData{$next[0]}->[0], $protData{$next[0]}->[3], $protData{$next[0]}->[1], $protData{$next[0]}->[2], $GENE_NAMES{$next[0]}, $next[1],
									$PROTPROPH_INFO->{$strip2}->{$next[0]}->[2];
								$index++;							
							}
						}	
					}
				}

				# get the stripped peptide so can look up the protein information for the pair
				printf PEPXML "   <xlink identifier=\"BDP-NHP\" mass=\"200.00\">\n";
				for(my $k = 0; $k < @{$results{$_}->{1}->[-1]}; $k++) {
					printf PEPXML "%s\n", $results{$_}->{1}->[-1]->[$k];
				}
				#printf PEPXML " <xlink_score name=\"link\" value=\"%d\"/>\n", (length $results{$_}->{1}->[2])+1;
				printf PEPXML " <xlink_score name=\"score\" value=\"%s\"/>\n", (sprintf "%0.4f", $results{$_}->{1}->[0]);
				printf PEPXML "</linked_peptide>\n";
				
				for(my $k = 0; $k < @{$results{$_}->{2}->[-1]}; $k++) {
					printf PEPXML "%s\n", $results{$_}->{2}->[-1]->[$k];
				}
				#printf PEPXML " <xlink_score name=\"link\" value=\"%d\"/>\n", (length $results{$_}->{2}->[2])+1;
				printf PEPXML "<xlink_score name=\"score\" value=\"%s\"/>\n", (sprintf "%0.4f", $results{$_}->{2}->[0]);
				printf PEPXML "</linked_peptide>\n";
     			printf PEPXML "</xlink>\n";
				my $delta_score = '';
				if($distvals{'homopeptide'}->{$_}) {
					$delta_score = 'homopeptide:1';
				}
				else {
					$delta_score = 'intra:';
					$delta_score .= exists $distvals{'intra'} ? $distvals{'intra'}->{$_} : '';
					$delta_score .= '_nsx:';
					$delta_score .= exists $distvals{'nsx'} ? $distvals{'nsx'}->{$_} : '';
				}
				printf PEPXML "   <search_score name=\"expect\" value=\"%s\"/>\n", $distvals{'max_expect'}->{$_};
				printf PEPXML "   <search_score name=\"kojak_score\" value=\"%s\"/>\n", $distvals{'max_expect'}->{$_}; 
				printf PEPXML "   <search_score name=\"delta_score\" value=\"%s\"/>\n", $delta_score;
				printf PEPXML "   <search_score name=\"ppm_error\" value=\"%s\"/>\n", $distvals{'massdiff_ppm'}->{$_};
				if($PAIRING eq 'Mango' && $ENCODE_IN_MANGO_RUNNERUP_SCORE && $probabilities{$_} < 0) {	
					printf PEPXML "   <search_score name=\"mango_runnerup\" value=\"%s\"/>\n", $probabilities{$_};
				}			
				# add composite_probability in case
				if($COMPUTE_COMPOSITE_FDRs) {
					if(exists $PEP_PAIRS{$_} && exists $CROSSLINK_FDRs{$PEP_PAIRS{$_}}) {
						printf PEPXML "   <search_score name=\"composite_probability\" value=\"%0.4f\"/>\n", 1 - $CROSSLINK_FDRs{$PEP_PAIRS{$_}};
						printf PEPXML "   <search_score name=\"composite_id\" value=\"%s\"/>\n", $PEP_PAIRS{$_};
					}
					else {
						printf PEPXML "   <search_score name=\"composite_probability\" value=\"0\"/>\n";
					}
				}
				for(my $k = 0; $k < @distributions_used; $k++) {
					next if(exists $posdists{$distributions_used[$k]});
     				printf PEPXML "   <search_score name=\"%s\" value=\"%s\"/>\n", $distributions_used[$k], $distvals{$distributions_used[$k]}->{$_};
				}
		
				printf PEPXML " <analysis_result analysis=\"peptideprophet\">\n";
				if($PAIRING eq 'Mango' && $ENCODE_IN_MANGO_RUNNERUP_SCORE && $probabilities{$_} < 0) {
					printf PEPXML " <peptideprophet_result probability=\"%0.4f\" all_ntt_prob=\"(%0.4f,%0.4f,%0.4f)\">\n", 0, 0, 0, 0;
				}
				else {
					printf PEPXML " <peptideprophet_result probability=\"%0.4f\" all_ntt_prob=\"(%0.4f,%0.4f,%0.4f)\">\n", 
						$probabilities{$_}, $probabilities{$_}, $probabilities{$_}, $probabilities{$_};
				}
				printf PEPXML " <search_score_summary>\n";
				my $fval = ! ($PAIRING eq 'Mango') ? $distvals{'spectrum2'}->{$_} : $query_no + 1;
				printf PEPXML "    <parameter name=\"%s\" value=\"%s\"/>\n", "fval", $fval;  
				for(my $k = 0; $k < @distributions_used; $k++) {
					next if(! exists $posdists{$distributions_used[$k]});
     				printf PEPXML "   <parameter name=\"%s\" value=\"%s\"/>\n", $distributions_used[$k], $distvals{$distributions_used[$k]}->{$_};
				}
				printf PEPXML " </search_score_summary>\n";
				printf PEPXML " </peptideprophet_result>\n";
				printf PEPXML "</analysis_result>\n";
    			printf PEPXML "  </search_hit>\n";
   				printf PEPXML " </search_result>\n";
  				printf PEPXML "</spectrum_query>\n";
  			} # only if output
  		} # next spectrum
		printf PEPXML "</msms_run_summary>\n";
	} # next run
	printf PEPXML "</msms_pipeline_analysis>\n";
	close(PEPXML);
} # output pepXML


# now print out the distributions
printf STDERR "Distributions after %d iterations\n", $iter - 1;
printDistributions(\%posdists, \%negdists, *STDERR);
printf STDERR "decoy fdr multiplier: %0.2f\n\n", $num_lows if($compute_decoy_fdr);
my $fraction_intra = "";
my $fraction_intra_nonred = "";
if($DISPLAY_INTER_NO_INTRA) {
	$fraction_intra = exists $posdists{'intra'} ? (sprintf " (%0.0f%% inter-protein),", (1 - $posdists{'intra'}->{1}) * 100) : "";
	$fraction_intra_nonred = $tot_non_redundant > 0 ? sprintf " (%0.0f%% inter-protein)", (1 - $tot_non_redundant_intra /$tot_non_redundant)*100 : "";
}
else {
	$fraction_intra = exists $posdists{'intra'} ? (sprintf " (%0.0f%% intra-protein),", $posdists{'intra'}->{1} * 100) : "";
	$fraction_intra_nonred = $tot_non_redundant > 0 ? sprintf " (%0.0f%% intra-protein)", $tot_non_redundant_intra * 100/$tot_non_redundant : "";
}
printf STDERR "Computed probabilities with estimated %0.0f correct cross-links%s including %0.0f non-redundant%s, written to output file%s %s%s\n\n", 
	$total_correct, $fraction_intra, $tot_non_redundant, $fraction_intra_nonred, 
	$OUTPUT_PEPXML eq '' ? '' : 's', $outfile, $OUTPUT_PEPXML eq '' ? '' : ' and ' . $OUTPUT_PEPXML . '.xml';

unlink "$OUTPUT_PEPXML-MODELS.html" if(! ($OUTPUT_PEPXML eq '') && -e "$OUTPUT_PEPXML-MODELS.html");

if(! ($OUTPUT_PEPXML eq '') && ! ($EMAIL_ADDRESS eq '')) {

	my %orgs = ("HUMAN" => "H.sapiens", "YEAST" => "S.cerevisiae", "ARATH" => "A.thaliana", "MOUSE" => "M.musculus", "BOVIN" => "B.taurus", "ECOLI" => "E.coli");
	my $url_extn = ""; #&org=H.sapiens";
	my @organism = keys %ORGANISMS;
	if(0 && @organism == 1 && exists $orgs{$organism[0]}) {
		$url_extn = "&org=$orgs{$organism[0]}";
	}
	
	if($custom_reporter) {
		$url_extn = "&reporter_mass={$REPORTERMASS}";
	}
	
	my $subject = "XLinkProphet $file analysis complete and ready for upload to XLinkDB";
	my $message = sprintf "Computed probabilities with estimated %0.0f correct cross-links%s including %0.0f non-redundant, written to output file%s %s%s.\n\n", $total_correct, $fraction_intra, $tot_non_redundant, 
		$OUTPUT_PEPXML eq '' ? '' : 's', $outfile, $OUTPUT_PEPXML eq '' ? '' : ' and ' . $OUTPUT_PEPXML . '.xml';

	$message .= "View crosslink results on tephra at: https://proteomicsresource.washington.edu".$OUTPUT_PEPXML.".xml.\n\n";

	$message .= "To upload your data to XLinkDB, please click on the link: http://xlinkdb.gs.washington.edu/xlinkdb/new_index.php?email=$EMAIL_ADDRESS&input=$OUTPUT_PEPXML.xml$url_extn and on the right Upload New Local Data section, choose organism, write experiment name and description, verify all other settings, and upload.";

	sendEmail($subject, $message);

}


}

sub printDistributions {
(my $posdist, my $negdist, my $out) = @_;
foreach(sort {$a cmp $b} keys %{$posdist}) {
	printf $out "%s distribution\n", $_;
	printf $out "%s\t%s\t%s\n", "value", "pos", "neg";
	my @next = sort {$a <=> $b} keys %{$posdist->{$_}};
	for(my $k = 0; $k < @next; $k++) {
		printf $out "%s\t%0.3f\t%0.3f\n", exists $MODEL_BINS{$_} ? $next[$k] * $MODEL_BINS{$_} : $next[$k], 
			$posdist->{$_}->{$next[$k]}, $negdist->{$_}->{$next[$k]};
	}
	printf $out "\n";
}
}

# records the peptide cross link attaching the two proteins (for later calculation of nsp)
sub recordProteinPairs {
(my $prot1, my $prot2, my $prot_pair_ptr, my $prob, my $pep_pair) = @_;
my @prots1 = split(",", $prot1);
my @prots2 = split(",", $prot2);

for(my $k = 0; $k < @prots1; $k++) {
	for(my $j = 0; $j < @prots2; $j++) {
		if($prots1[$k] cmp $prots2[$j]) {
			$prot_pair_ptr->{$prots1[$k] . '_' . $prots2[$j]} = {} if(! exists $prot_pair_ptr->{$prots1[$k] . '_' . $prots2[$j]});
			$prot_pair_ptr->{$prots1[$k] . '_' . $prots2[$j]}->{$pep_pair} = $prob if($prot_pair_ptr->{$prots1[$k] . '_' . $prots2[$j]}->{$pep_pair} < $prob); # take the max
		}
		else {
			$prot_pair_ptr->{$prots2[$j] . '_' . $prots1[$k]} = {} if(! exists $prot_pair_ptr->{$prots2[$j] . '_' . $prots1[$k]});
			$prot_pair_ptr->{$prots2[$k] . '_' . $prots1[$j]}->{$pep_pair} = $prob if($prot_pair_ptr->{$prots2[$k] . '_' . $prots1[$j]}->{$pep_pair} < $prob); # take the max
		}
	}
}
}

# number of other instances of crosslink
sub getNrx {
(my $spectrum, my $probptr) = @_;
my $nrx = 0;
my $nrx_value = exists $PEP_PAIRS{$spectrum} ? $PEP_PAIRS{$spectrum} : die "No Nrx pep pair for spectrum $spectrum\n";
return $nrx if(scalar keys %{$CROSSLINKS{$nrx_value}} < 2);
foreach(keys %{$CROSSLINKS{$nrx_value}}) {
	next if($_ eq $spectrum);
	$nrx += $probptr->{$_}; 
}
if($nrx < 0.5) {
	return 0;
}
elsif($nrx < 3) {
	return 1;
}
return 2;


}

# number of other crosslinks spanning same 2 proteins
sub getNsx {
(my $prot1, my $prot2, my $prot_pair_ptr, my $pep_pair) = @_;
my @prots1 = split(",", $prot1);
my @prots2 = split(",", $prot2);
my $nsp = 0;
my $MIN_PROB = 0.5;
my $NSX_BIN_DIVIDER = 2.5; 
for(my $k = 0; $k < @prots1; $k++) {
	for(my $j = 0; $j < @prots2; $j++) {
		if($prots1[$k] cmp $prots2[$j]) {
			my $next = 0;
			my @peps = keys %{$prot_pair_ptr->{$prots1[$k] . '_' . $prots2[$j]}};
			for(my $p = 0; $p < @peps; $p++) {
				next if($prot_pair_ptr->{$prots1[$k] . '_' . $prots2[$j]}->{$peps[$p]} < $MIN_PROB);
				$next += $prot_pair_ptr->{$prots1[$k] . '_' . $prots2[$j]}->{$peps[$p]} if(! ($peps[$p] eq $pep_pair));
			}
			
			$nsp = $next if($nsp < $next);
		}
		else {
			my $next = 0;
			my @peps = keys %{$prot_pair_ptr->{$prots2[$j] . '_' . $prots1[$k]}};
			for(my $p = 0; $p < @peps; $p++) {
				next if($prot_pair_ptr->{$prots2[$k] . '_' . $prots1[$j]}->{$peps[$p]} < $MIN_PROB);
				$next += $prot_pair_ptr->{$prots2[$j] . '_' . $prots1[$k]}->{$peps[$p]} if(! ($peps[$p] eq $pep_pair));
			}
			$nsp = $next if($nsp < $next);
		}
	}
}
$CURRENT_NSP = sprintf "%0.1f", $nsp;

my $output = (sprintf "%0.0f", $nsp);
$output = 0 if($output eq '-0');
if(@NSX_BIN_BOUNDS > 0) {
	for(my $k = 0; $k < @NSX_BIN_BOUNDS; $k++) {
		if($output < $NSX_BIN_BOUNDS[$k]) {
			return $k;
		}
	}
	return scalar @NSX_BIN_BOUNDS;
}
return $nsp > $NSX_BIN_DIVIDER ? 1 : 0;
}


# whether the cross-link is possibly between 2 peptides of same protein
sub getIntra {
(my $prot1, my $prot2) = @_;
my @prots1 = split(",", $prot1);
my @prots2 = split(",", $prot2);
for(my $k = 0; $k < @prots1; $k++) {
	for(my $j = 0; $j < @prots2; $j++) {
		return 1 if($prots1[$k] eq $prots2[$j]);
	}
}
return 0;
}


sub getPeptideNrxValueForLightHeavyResults {
(my $modpep) = @_;
die "Problem with null value in peptidepair for Nrx\n" if($modpep eq '');
my $copy = $modpep;
my $verbose = 0;
if(0 && $copy =~ /K\[136/) {
	printf "HERE with $copy\n";
	$verbose = 1;
}
for(my $k = 0; $k < @CROSSLINK_MODS; $k++) {
	my @next = keys %{$CROSSLINK_MODS[$k]};
	for(my $j = 0; $j < @next; $j++) {
		my $merge = $next[$j] . "x"; 
		my $next_aa = $next[$j];
		$copy =~ s/$next_aa\[$CROSSLINK_MODS[$k]->{$next_aa}\]/$merge/g;
	}
}
for(my $k = 0; $k < @SILAC_MODS; $k++) {
	my @next = keys %{$SILAC_MODS[$k]};
	for(my $j = 0; $j < @next; $j++) {
		my $merge = $next[$j]; 
		my $next_aa = $next[$j];
		$copy =~ s/$next_aa\[$SILAC_MODS[$k]->{$next_aa}\]/$merge/g;
	}
}

if($REMOVE_SILAC_LYS_ARG_MODS) {		
	$copy =~ s/K\[178.12\]/K\[170.11\]/g; # acetyl Lys
	$copy =~ s/K\[178.16\]/K\[170.14\]/g; # trimethyl Lys
	$copy =~ s/K\[164.14\]/K\[156.13\]/g; # dimethyl Lys
	$copy =~ s/K\[150.12\]/K\[142.11\]/g; # methyl Lys
	$copy =~ s/R\[176.14\]/R\[170.12\]/g; # methyl Arg
	$copy =~ s/R\[190.15\]/R\[184.13\]/g; # dimethyl Arg
}

if($verbose) {
	printf "Have final $copy\n"; exit(1);
}
return $copy;
}

# unique peptide pair
sub getPeptidePairForNrx {
(my $modpep1, my $modpep2, my $stripped1, my $stripped2) = @_;
my $first =  $LIGHT_HEAVY ? getPeptideNrxValueForLightHeavyResults($modpep1) : getPeptideNrxValue($modpep1);
my $second =  $LIGHT_HEAVY ? getPeptideNrxValueForLightHeavyResults($modpep2) : getPeptideNrxValue($modpep2);
if(($stripped1 cmp $stripped2) > 0) {
	return $second . "_" . $first;
}
return $first . "_" . $second;
}

# have to maintain only the crosslink position in peptide
sub getPeptideNrxValue {
(my $modpep) = @_;
die "Problem with null value in peptidepair for Nrx\n" if($modpep eq '');
my $copy = $modpep;
for(my $k = 0; $k < @CROSSLINK_MODS; $k++) {
	my @next = keys %{$CROSSLINK_MODS[$k]};
	for(my $j = 0; $j < @next; $j++) {
		my $merge = $next[$j] . "x"; 
		my $next_aa = $next[$j];
		$copy =~ s/$next_aa\[$CROSSLINK_MODS[$k]->{$next_aa}\]/$merge/g;
	}
}
$copy =~ s/\[.*?\]//g;
return $copy;
}

# each peptide must have one and only one legal modification either at the n-terminus or at non-C-terminal lysine
sub getNonCrosslinkedPeptidePair {
(my $modpep1, my $modpep2) = @_;
my $found = 0;
my $crosslink_index = -1;
my $verbose = 0; #$modpep1 eq 'AK[214.08]SIVFHR' && $modpep2 eq 'SNYNFEK[214.08]PFLWLAR';
for(my $k = 0; $k < @CROSSLINK_MODS; $k++) {
	my @next = keys %{$CROSSLINK_MODS[$k]};
	for(my $j = 0; $j < @next; $j++) {
		my $next_mod = $next[$j].'\['.$CROSSLINK_MODS[$k]->{$next[$j]}.'\]';
		if($modpep1 =~ /$next_mod$/) { # can't be C-terminal
			return 1;
		}
		if($modpep1 =~ /$next_mod/) {
			if($found) {
				return 1; # cannot have two mods in one peptide
			}
			if($crosslink_index > -1 && $crosslink_index != $k) {
				return 1;
			}
			$found = 1;
			$crosslink_index = $k;
		}
	}
}
return 1 if(! $found);
$found = 0;
my @next = keys %{$CROSSLINK_MODS[$crosslink_index]};
for(my $j = 0; $j < @next; $j++) {
	my $next_mod = $next[$j].'\['.$CROSSLINK_MODS[$crosslink_index]->{$next[$j]}.'\]';
	if($modpep2 =~ /$next_mod$/) { # can't be C-terminal
		return 1;
	}
	if($modpep2 =~ /$next_mod/) {
		if($found) {
		return 1; # cannot have two mods in one peptide
		}
		$found = 1;
	}
}
if($verbose) {
	printf "Returning with %d\n", $found ? " found" : " not found";
}
return 1 if(! $found);
return 0; # ok

}

# a decoy has no target in either prot1 or prot2
sub getDecoy {
(my $prot1, my $prot2) = @_;
my @prots1 = split(",", $prot1);
my @prots2 = split(",", $prot2);
my @result = (1, 1); # assume have decoy in both unless proven false
for(my $k = 0; $k < @prots1; $k++) {
	$result[0] = 0 if($prots1[$k] !~ /^$DECOY_PREFIX/);
}
for(my $k = 0; $k < @prots2; $k++) {
	$result[1] = 0 if($prots2[$k] !~ /^$DECOY_PREFIX/);
}
return $result[0] || $result[1];
}

sub getPeptideLength {
(my $pep) = @_;
while($pep =~ /^(\S+)\[.*?\](\S*)$/) {
	$pep = $1 . $2;
}
return length $pep;
}

# a decoydecoy has no target in both prot1 and prot2
sub getDecoyDecoy {
(my $prot1, my $prot2) = @_;
my @prots1 = split(",", $prot1);
my @prots2 = split(",", $prot2);
my @result = (1, 1); # assume have decoy in both unless proven false
for(my $k = 0; $k < @prots1; $k++) {
	$result[0] = 0 if($prots1[$k] !~ /^$DECOY_PREFIX/);
}
for(my $k = 0; $k < @prots2; $k++) {
	$result[1] = 0 if($prots2[$k] !~ /^$DECOY_PREFIX/);
}
return $result[0] && $result[1];
}

sub getPeptideLength {
(my $pep) = @_;
while($pep =~ /^(\S+)\[.*?\](\S*)$/) {
	$pep = $1 . $2;
}
return length $pep;
}

# set any overflow values within a bin
sub geValueWithinDistributionBounds {
(my $value, my $sortedptr) = @_;
return $sortedptr->[0] if($value < $sortedptr->[0]);
return $sortedptr->[-1] if($value > $sortedptr->[-1]);
return $value;
}

# ReACT (Comet) or Mango
sub isReact {
(my $pprophfile) = @_;
if($pprophfile =~ /\.xml$/) {
	open GREP, "grep search_engine $pprophfile | head -1 |" or die "cannot grep $pprophfile $!\n";
	my @results = <GREP>;
	close(GREP);
	if(@results == 1) {
		chomp $results[0];
		if($results[0] =~ /search\_engine\_version\=\"Mango/) {
			return 0;
		}
	}
	else {
		die "error with ", scalar @results, " hits rather than 1\n";
	}
}
elsif($pprophfile =~ /\.xls$/) {
	open HEAD, "head -2 $pprophfile | tail -1 |" or die "cannot head $pprophfile $!\n";
	my @results = <HEAD>;
	close(HEAD);
	if(@results == 1) {
		chomp $results[0];
		if($results[0] =~ /\_\d\d\d(.\d\d\d\d\d)$1\.\d/) {
			return 0;
		}
	}
	else {
		die "error with ", scalar @results, " hits rather than 1\n";
	}
}
else {
	die "Illegal file format for $pprophfile\n";
}
return 1;
}

sub getBinValue {
(my $value, my $distr) = @_;
if(exists $MODEL_BINS{$distr}) {
	$value = sprintf "%0.0f", $value / $MODEL_BINS{$distr};
	$value = '0' if($value eq '-0');
}
return $value;
}

# reads in ProteinProphet output
sub attributePeptides2Proteins {
(my $protxml) = @_;
my %pepProts = ();
my @prots = ();
open(FILE, $protxml) or die "cannot read $protxml $!\n";
while(<FILE>) {
	chomp();
	if(/\<protein protein\_name=\"(\S+?)\"/) {
		@prots = ($1);
	}
	elsif(/\<indistinguishable\_protein protein\_name\=\"(\S+?)\"\>/) {
		push(@prots, $1);
		
	}
	elsif(/\<peptide peptide\_sequence\=\"(\S+?)\".*?weight\=\"(\S+?)\" is\_nondegenerate_evidence\=\"(\S)\".*?n\_sibling\_peptides\=\"(\S+?)\"/) {
		my $pep = $1;
		my $wt = $2;
		my $degen = $3;
		my $nsp = $4;
		my $verbose = $pep eq 'SPLAKK';
		for(my $k = 0; $k < @prots; $k++) {
			my $prot = $prots[$k];
			next if($prot =~ /^rev\_/);
			if(exists $PROTPROPH_INFO->{$pep} && 
				! exists $PROTPROPH_INFO->{$pep}->{$prot} &&
				($OUTPUT_DECOYS || $prot !~ /^$DECOY_PREFIX/) &&
				$wt > 0) {
				$PROTPROPH_INFO->{$pep}->{$prot} = [$wt/@prots, $degen, $nsp];
				$GENE_NAMES{$prot} = '' if(exists $PEP_MAXPROBS{$pep} && $PEP_MAXPROBS{$pep} >= 0.05);
			}
		}
	}
}
close(FILE);
my $tot = 0;
foreach(keys %{$PROTPROPH_INFO}) {
	$tot++ if(scalar keys %{$PROTPROPH_INFO->{$_}} > 0);
}
printf STDERR "Read in protein information for %d peptides out of %d total results from file %s\n", $tot, scalar keys %{$PROTPROPH_INFO}, $protxml; 
}


# returns [prot1, [runnerups:wt], prot2, [runnerups:wt]]
sub getProtInfo {
(my $pep1, my $pep2, my $pepProtptr) = @_;
my %pairs = ();
my $verbose = 0; #$pep1 eq 'VAPEEHPVLLTEAPLNPKANR' && $pep2 eq 'GILTLKYPIEHGIITNWDDMEK';
return [] if(! exists $pepProtptr->{$pep1} || ! exists $pepProtptr->{$pep2});
my @nextProt1 = keys %{$pepProtptr->{$pep1}};
my @nextProt2 = keys %{$pepProtptr->{$pep2}};

for(my $k = 0; $k < @nextProt1; $k++) {
	for(my $j = 0; $j < @nextProt2; $j++) {
		if($verbose) {
			printf "$k and $j: ($pepProtptr->{$pep1}->{$nextProt1[$k]}->[0], $pepProtptr->{$pep2}->{$nextProt2[$j]}->[0]), ($pepProtptr->{$pep1}->{$nextProt1[$k]}->[2], $pepProtptr->{$pep2}->{$nextProt2[$j]}->[2])\n";
		}
		$pairs{$k."_".$j} = [$pepProtptr->{$pep1}->{$nextProt1[$k]}->[0]*$pepProtptr->{$pep2}->{$nextProt2[$j]}->[0], ($nextProt1[$k] eq $nextProt2[$j] ? 1 : 0), 
			$pepProtptr->{$pep1}->{$nextProt1[$k]}->[2]+$pepProtptr->{$pep2}->{$nextProt2[$j]}->[2], substr($nextProt1[$k],10) . '_' . substr($nextProt2[$j],10), 
			(exists $protpair_instances{$nextProt1[$k] . "_" . $nextProt2[$j]} ? $protpair_instances{$nextProt1[$k] . "_" . $nextProt2[$j]} : 0)];
	}
} # next k
		# now sort them by wt
# orig
my @sorted = sort {$pairs{$b}->[0] <=> $pairs{$a}->[0] or $pairs{$b}->[1] <=> $pairs{$a}->[1] or $pairs{$b}->[2] <=> $pairs{$a}->[2] or $pairs{$a}->[3] cmp $pairs{$b}->[3] } keys %pairs;

if($verbose) {
	printf "Have %d sorted: %s for $pep1 and $pep2\n", scalar @sorted, join(",", @sorted);
	printf "PROts1: %s, PROts2: %s\n", join(",", @nextProt1), join(",", @nextProt2);
}
# will have to order prots1 and prots2 also so can assign the alternatives in order
my %sorted_prots1 = ();
my %sorted_prots2 = ();
my $winner = 0;
my $min_wtsq = 0.05;
my $max_wtsq_diff = 0.3;
for(my $k = 0; $k < @sorted; $k++) {
	if($sorted[$k] =~ /^(\d+)\_(\d+)$/) {
		printf "%s\t%s\t%s\t%s\t%0.2f\t%d\t%0.1f\n", $pep1, $nextProt1[$1], $pep2, $nextProt2[$2], $pairs{$sorted[$k]}->[0], $pairs{$sorted[$k]}->[1], 
			$pairs{$sorted[$k]}->[2] if(0);
		$sorted_prots1{$nextProt1[$1]} = (scalar keys %sorted_prots1) if(! exists $sorted_prots1{$nextProt1[$1]});
		$sorted_prots2{$nextProt2[$2]} = (scalar keys %sorted_prots2) if(! exists $sorted_prots2{$nextProt2[$2]});
	}
}
		# rules for finding the winner
		
# wtsq, intra, nsp_total sort by these
if(@sorted > 1 && ! $pairs{$sorted[0]}->[1]) {
	for(my $k = 1; $k < @sorted; $k++) {
		if($pairs{$sorted[$k]}->[1] && $pairs{$sorted[$k]}->[0] >= $min_wtsq &&
			$pairs{$sorted[0]}->[0] - $pairs{$sorted[$k]}->[0] <= $max_wtsq_diff) {
			$winner = $k;
			$k = @sorted;
		}
	}
}
if($verbose) {
	printf "Winner %d with value %s\n", $winner, $sorted[$winner];
}
if($sorted[$winner] =~ /^(\d+)\_(\d+)$/) {
		
	delete $sorted_prots1{$nextProt1[$1]};
	delete $sorted_prots2{$nextProt2[$2]};
	my @nextAlt1 = sort {$sorted_prots1{$a} <=> $sorted_prots1{$b}} keys %sorted_prots1;
	my @nextAlt2 = sort {$sorted_prots2{$a} <=> $sorted_prots2{$b}} keys %sorted_prots2;
			
	for(my $k = 0; $k < @nextAlt1; $k++) {
		if(! exists $pepProtptr->{$pep1}->{$nextAlt1[$k]}) {
			
			printf "Here without value $pep1 with $nextAlt1[$k]\n";
			return [];
		}
		if($pepProtptr->{$pep1}->{$nextAlt1[$k]}->[0] > 0) {
			$nextAlt1[$k] .= (sprintf ":%0.2f",  $pepProtptr->{$pep1}->{$nextAlt1[$k]}->[0]);
		}
		else {
			$nextAlt1[$k] .= ":0";
		}
	}
	for(my $k = 0; $k < @nextAlt2; $k++) {
		if(! exists $pepProtptr->{$pep2}->{$nextAlt2[$k]}) {
			
			printf "Here without value $pep2 with $nextAlt2[$k]\n";
			return [];
		}
		if($pepProtptr->{$pep2}->{$nextAlt2[$k]}->[0] > 0) {
			$nextAlt2[$k] .= (sprintf ":%0.2f", $pepProtptr->{$pep2}->{$nextAlt2[$k]}->[0]);
		}
		else {
			$nextAlt2[$k] .= ":0";
		}
	}
	printf "Returniung $nextProt1[$1] as winner\n" if($verbose);
	return [$nextProt1[$1], \@nextAlt1, $nextProt2[$2], \@nextAlt2];

}
else {
	die "Error with winner $winner\n";
}
}

# this can be bypassed if desired
sub setGeneNamesFromUniprot {
my @genes = keys %GENE_NAMES;
my $tot = 0;
my %APPENDED_NAMES = (); # record the new ones and append to file
my %orgs = ("HUMAN" => "Homo sapiens (Human)", "YEAST" => "Saccharomyces cerevisiae (strain ATCC 204508 / S288c) (Baker's yeast)", 
			"ARATH" => "Arabidopsis thaliana (Mouse-ear cress)", "MOUSE" => "Mus musculus (Mouse)",
			"BOVIN" => "Bos taurus (Bovine)", "ECOLI" => "Escherichia coli (strain K12)");

if(-e $UNIPROT_CONVERSION_FILE) {
	open(UNI, $UNIPROT_CONVERSION_FILE) or die "cannot read $UNIPROT_CONVERSION_FILE $!\n";
	my $first = 1;
	while(<UNI>) {
		if($first) {
			$first = 0;
		}
		else {
			chomp();
			my @parsed = split("\t");
			if(exists $GENE_NAMES{"sp|".$parsed[0]."|".$parsed[1]}) {
				$GENE_NAMES{"sp|".$parsed[0]."|".$parsed[1]} = $parsed[4] eq '' ? '-' : (split(" ", $parsed[4]))[0];
				$tot++;
			}
			elsif(exists $GENE_NAMES{"tr|".$parsed[0]."|".$parsed[1]}) {
				$GENE_NAMES{"tr|".$parsed[0]."|".$parsed[1]} = $parsed[4] eq '' ? '-' : (split(" ", $parsed[4]))[0];
				$tot++;
			}
		}
	}
	close(UNI);
	printf STDERR "Read in %d genes from %s\n", $tot, $UNIPROT_CONVERSION_FILE;
}
for(my $k = 0; $k < @genes; $k++) {
	if($GENE_NAMES{$genes[$k]} eq '') {
		$GENE_NAMES{$genes[$k]} = getGeneNameFromUniprot($genes[$k]);
		if($genes[$k] =~ /^sp\|(\S+?)\|(\S+\_)(\S+)$/) {
			$APPENDED_NAMES{$1} = sprintf "%s\t%s\tunknown\t\t%s\t%s", $1, $2 . $3, 
				($GENE_NAMES{$genes[$k]} eq '-' ? "" : $GENE_NAMES{$genes[$k]}), (exists $orgs{$3} ? $orgs{$3} : "");
			$ORGANISMS{$3}++;
		}
		elsif($genes[$k] =~ /^tr\|(\S+?)\|(\S+\_)(\S+)$/) {
			$APPENDED_NAMES{$1} = sprintf "%s\t%s\tunknown\t\t%s\t%s", $1, $2 . $3, 
				($GENE_NAMES{$genes[$k]} eq '-' ? "" : $GENE_NAMES{$genes[$k]}), (exists $orgs{$3} ? $orgs{$3} : "");
			$ORGANISMS{$3}++;
		}
		$tot++ if(! ($GENE_NAMES{$genes[$k]} eq ''));
	}
	elsif($GENE_NAMES{$genes[$k]} eq '-') {
		$GENE_NAMES{$genes[$k]} = ""; # revert
	}
}
printf STDERR "Obtained gene names for $tot entries\n";
if(scalar keys %APPENDED_NAMES > 0) {
	open(CONV, ">>$UNIPROT_CONVERSION_FILE") or die "cannot append to $UNIPROT_CONVERSION_FILE $!\n";
	foreach(keys %APPENDED_NAMES) {
		printf CONV "%s\n", $APPENDED_NAMES{$_};
	}
	close(CONV);
	printf STDERR "Appended %d new entries to %s\n", scalar keys %APPENDED_NAMES, $UNIPROT_CONVERSION_FILE;
}
close(GLOG);
}

# assumes the original search database contained protein names in uniprot format
sub getGeneNameFromUniprot {
(my $prot) = @_;
return '' if(! $QUERY_UNIPROT_GENENAMES);
if($prot =~ /^sp|tr\|(\S+?)\|/) {
	my $uni = $1;
	my $output = "uniprot.txt";
	my $command = "wget -O $output \"http://www.uniprot.org/uniprot/$uni\" > tmp.out";
	system($command);
	open(UNI, $output) or die "cannot open $output for $command\n";
	while(<UNI>) {
		if(/Gene\<\/div\>\<div id\=\"content\-gene\" class\=\"entry\-overview\-content\"\>\<h2\>(\S+?)\<\/h2\>/) {
			my $gene = $1;
			printf STDERR "Found gene $gene for $uni!\n";
			close(UNI);
			
			unlink $output if(-e $output);
			return $gene;
		}
	}
	close(UNI);
	unlink $output if(-e $output);
	return '';
}
elsif(0 && $prot =~ /^tr\|(\S+?)\|/) {
	my $uni = $1;
	# still here, go to uniprot manually
	my $output = "uniprot.txt";
	my $command = "wget -O $output \"http://www.uniprot.org/uniprot/$uni\" > tmp.out";
	system($command);
	open(UNI, $output) or die "cannot open $output for $command\n";
	while(<UNI>) {
		if(/Gene\<\/div\>\<div id\=\"content\-gene\" class\=\"entry\-overview\-content\"\>\<h2\>(\S+?)\<\/h2\>/) {
			my $gene = $1;
			printf STDERR "Found gene $gene for $uni!\n";
			close(UNI);
			
			unlink $output if(-e $output);
			return $gene;
		}
	}
	close(UNI);
	unlink $output if(-e $output);
	return '';
}
else {
	die "Problem with protein name $prot\n";
}
die "Illegal protein format without uniprot identifier (sp|xxxx|...) or (tr|xxxx|...): $prot\n";
return '';


}

sub sendEmail {
(my $subject, my $message) = @_;
	

}


sub removeStaticMods {
(my $modpep) = @_;
$modpep =~ s/C\[160\]/C/g;
$modpep =~ s/K\[136\]/K/g;
$modpep =~ s/R\[162\]/R/g;
return $modpep;
}

sub addModpepDecimalPlaces {
(my $modpep, my $modpos, my $modmass) = @_;
return $modpep if($NUM_MODPEP_DECIMALS==0);
my $pos = 0;
my $output = "";
for(my $k = 0; $k < (length $modpep); $k++) {
	my $next = substr($modpep, $k, 1);
	if($next =~ /[A-Z]/) {
		$pos++;
		$output .= $next;
		if($pos==$modpos) {
			$output .= '[' .sprintf("%0.".$NUM_MODPEP_DECIMALS."f", $modmass) . ']';
		}
	}
	else {
		if($pos!=$modpos) {
			$output .= $next;
		}
	}	
}
return $output;
}


sub getAAScore {
(my $id) = @_;
my %distr = (
"A" =>0.094,
"P" =>0.046,
"D" =>0.050,
"Q" =>0.049,
"L" =>0.119,
"C" =>0.011,
"N" =>0.044,
"Y" =>0.031,
"T" =>0.054,
"E" =>0.082,
"V" =>0.081,
"H" =>0.020,
"F" =>0.044,
"M" =>0.026,
"G" =>0.091,
"I" =>0.062,
"S" =>0.067,
"W" =>0.018,
"R" =>0.005,
"K" =>0.008
);
$id = stripPeptide($id);
my @next = split("_", $id);
my $tot = 0;
for(my $k = 0; $k < @next; $k++) {
	my @stump_env = split("Kx", $next[$k]);
	if($stump_env[0] eq '') {
		$tot += $distr{substr($stump_env[1], 0, 1)} + $distr{substr($stump_env[1], 1, 1)};
	}
	elsif($stump_env[1] eq '') {
		$tot += $distr{substr($stump_env[0], (length $stump_env[0]) - 1, 1)} + $distr{substr($stump_env[0], (length $stump_env[0]) - 2, 1)};
	}
	else {
		$tot += $distr{substr($stump_env[0], (length $stump_env[0]) - 1, 1)} + $distr{substr($stump_env[1], 0, 1)};
	}
}
return $tot;
}

# reads pepXML to detect whether modification masses have variable mod
sub getRemoveStaticMods {
(my $file, my $max_num_mods = 5) = @_;
my $tot = 0;
open(FILE, $file) or die "cannot read $file $!\n";
while(<FILE>) {
	if(/^<mod_aminoacid_mass/) {
		if(/variable/) {
			close(FILE);
			return 1;
		}
		$tot++;
		if($tot >= $max_num_mods) {
			close(FILE);
			return 0;
		}
	}
}
close(FILE);

return 0;
}
