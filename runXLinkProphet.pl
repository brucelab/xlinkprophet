#!/usr/bin/perl

use strict;
use POSIX;

my $SELF = 'runXLinkProphet_devxxx.pl';
my $XLINKPROPHET = "c:\\Users\\fullpathto\\XLinkProphet_devxxx.pl";


if(@ARGV == 0) {
	my $ARG_DELIMITER = $^O =~ /Win/ ? "\"" : "'";
	printf STDERR "\n";
	printf STDERR " usage:   $SELF < PepXML Files > ( options )\n";
	printf STDERR "\n";
	printf STDERR " e.g.:    \033[1mrunXLinkProphet.pl %s*.pep.xml%s\033[0m\n", $ARG_DELIMITER, $ARG_DELIMITER;
	printf STDERR " options: TEST (to print all the commands of the analysis, including the input pepXML files, without running the analysis)\n";
	printf STDERR "          PARAMS (to view xlinkprophet.params file sample contents)\n";
	printf STDERR "          WRITEPARAMS (to write generic xlinkprophet.params file in current directory)\n";
	printf STDERR "          OUTFILESUFF=xxx (add '-xxx' to the names of all files created during analysis, e.g. 'interact-xxx.pep.xml', 'iproph-xxx.pep.xml', 'iproph-xxx-xl.pep.xml')\n";
	printf STDERR "          SAMPLEMAP=xxx (Filter input pepxml files for those present in the MasschroQ samplemap file)\n";
	printf STDERR "          SPEC_LABEL=xxx (label spectrum tags so multiple instances of same spectrum are maintained [e.g. light and heavy searches])\n";
	printf STDERR "          REPORTERMASS=xxx (neutral mass of crosslink reporter [default: 751.40508 for BHP])\n";
	printf STDERR " \033[1;32mor just released for combined light and heavy search results:\033[0;39m\n";
	printf STDERR "          $SELF LIGHT=< Light/Shortarm PepXML Files > HEAVY=< Heavy/Longarm PepXML Files >( options )\n";
	printf STDERR " e.g.:    \033[1mrunXLinkProphet.pl %sLIGHT=light_dir/*.pep.xml%s %sHEAVY=heavy_dir/*.pep.xml%s\033[0m\n", $ARG_DELIMITER, $ARG_DELIMITER, $ARG_DELIMITER, $ARG_DELIMITER;
	printf STDERR "          IQPIR_REPORTERMASSES=xxx:yyy,www:zzz,... (specify reporter masses yyy and zzz corresponding to lysine stump modification masses xxx and yyy, respectively, [with at least 2 decimal places, e.g. 325.13], e.g. IQPIR_REPORTERMASSES=325.13:811.455947,327.13:807.442528)\n";
	printf STDERR "          DSSO [equivalent to IQPIR_REPORTERMASSES=214.077603:-13.96152,182.105523:49.98264]\n"; # long arm, short arm  214.037603
	printf STDERR "          SILAC_MODS=K:xxx,R:yyy (Specify silac heavy mass difference xxx for lysine and yyy for arginine [default values K:136.11,R:162.12])\n";
	printf STDERR "\n";
	printf STDERR "          ReACT: react2.xls files generated using react2csv options -c9999 -F must be present alongside input pepXML files\n";
	printf STDERR "          Mango: peaks files must be present alongside input pepXML files\n\n";
	printf STDERR "          Uses \033[1mxlinkprophet.params\033[0m file in current directory to read in run parameters for analysis\n\n";
	printf STDERR " \n";
	exit(1);

}

# ---> runXLinkProphet.pl

my $REQUIRE_PARAMSFILE = 1;
my $PARAMSFILE = "xlinkprophet.params";
my $PEP_PROPH_OPTIONS = "-OEAdP -p0 -PPM -l6 -drev_";
my $RUN_IPROPH = 1;
my $NO_NSX = 0;
my $OUTPUT_DECOYS = 0;
my $TARGET_NEG = 0;
my $SPEC_CONTAINS = "";
my $LOCAL_REACT = 0;
my $PAIRIFLE_SUFF = "";
my $CROSSLINK_MODS = "";
my $EMAIL = "";
my $XL_PROPH_OPTIONS = "";
my $TEST = 0;
my $pairfile = 0;
my $SAMPLEMAPFILES = {}; #getSamplemapFiles("/net/gs/vol4/shared/brucelab/search/jdchavez/Jimmy_linux_MCQ_test/SUBSET2/sample_map.txt"); #{};
my $FILTER_PEPXML_WITHOUT_PAIRFILE = 1;
my $MIN_PEP_LEN = "";
my $FORCE_NTT2 = 0;
my $SEND_EMAIL = 1;
my $STUMP_MOD = "";
my $DIR_DIVISOR_REGX = $^O =~ /Win/ ? "\\\\" : "\/";
my $DIR_DIVISOR = $^O =~ /Win/ ? "\\" : "/";
my $IS_WINDOWS = $^O =~ /Win/;
my @COMMANDS = ();

if($ARGV[0] eq 'PARAMS' || (@ARGV > 1 && $ARGV[1] eq 'PARAMS')) {
	displayParamsFile(\*STDOUT, "");
	exit(1);
}
if($ARGV[0] eq 'WRITEPARAMS' || (@ARGV > 1 && $ARGV[1] eq 'WRITEPARAMS')) {
	writeParamsFile();
	exit(1);
}
my $option_startind = 1;
my $light_heavy = 0;
my @light_files = ();
my @heavy_files = ();
my @files = ();
if($ARGV[0] =~ /^LIGHT\=(\S+)/) {
	@light_files = glob($1);
	if(@ARGV > 1 && $ARGV[1] =~ /^HEAVY\=(\S+)/) {
		@heavy_files = glob($1);
	}
	else {
		die "Error: no heavy files submitted as second parameter\n";
	}
	$option_startind = 2;
	$light_heavy = 1;
}
else {
	@files = glob($ARGV[0]);
}

my $suffix = "";

my @DSSO_STUMP_MASSES = (0, 0); # short and long arm

my $paramsfile = getcwd() . $DIR_DIVISOR . $PARAMSFILE;
if(-e $paramsfile) {
	$pairfile = readParams($paramsfile);
}
else {
	if($REQUIRE_PARAMSFILE) {
		printf STDERR "Error: required paramsfile $paramsfile not found. Please copy one into the current directory and edit to select desired settings.\n";
		exit(1);
	}
	printf STDERR "Warning: using default run parameters in the absence of params file $paramsfile\n";
}

$XLINKPROPHET = "\"" . $XLINKPROPHET . "\"" if($IS_WINDOWS && $XLINKPROPHET =~ /\s/);
my $exec = 1;
my $all = 1;
my $REPORTERMASS = "";
my $IQPIR_REPORTERMASSES = "";
for(my $k = $option_startind; $k < @ARGV; $k++) {
	if($ARGV[$k] eq 'TEST') {
		$exec = 0;
	}
	elsif($ARGV[$k] =~ /^OUTFILESUFF\=(\S+)/) {
		printf STDERR "Setting outfile suffix to $1\n";
		$suffix = "-" . $1;
	}
	elsif($ARGV[$k] eq 'PARAMS') {
		displayParamsFile(\*STDERR, "");
		exit(1);
	}
	elsif($ARGV[$k] =~ /^SAMPLEMAP\=(\S+)/) {
		printf STDERR "Filtering for files present in samplemap file $1\n";
		$SAMPLEMAPFILES = getSamplemapFiles($1);
	}
	elsif($ARGV[$k] =~ /^SPEC\_LABEL\=(\S+)/) {
		printf STDERR "Setting SPEC_LABEL to %s\n", $1; #exit(1);
		$XL_PROPH_OPTIONS .= "SPEC_LABEL=$1 "; 
	}
	elsif($ARGV[$k] =~ /^REPORTERMASS\=(\S+)/) {
		$REPORTERMASS = $1;
		printf STDERR "Setting REPORTERMASS to %s\n", $REPORTERMASS; #exit(1);
	}
	elsif($ARGV[$k] =~ /^IQPIR\_REPORTERMASSES\=(\S+)/) {
		if(0 && $STUMP_MOD) {
			printf "Error: comment out xlinker_stump_mod_masses in xlinkprophet.params file before running with the IQPIR_REPORTERMASSES option\n";
			exit(1);
		}
		$IQPIR_REPORTERMASSES = $1;
		printf STDERR "Setting IQPIR_REPORTERMASSES to %s\n", $IQPIR_REPORTERMASSES; #exit(1);
		my @reporters = split(",", $IQPIR_REPORTERMASSES);
		my @new_mods = ();
		for(my $k = 0; $k < @reporters; $k++) {
			my @next_mass = split(":", $reporters[$k]);
			push(@new_mods, "K:" . $next_mass[0]);
		}
		$XL_PROPH_OPTIONS .= "CROSSLINK_MODS=" . join(",", @new_mods) . " ";
	}
	elsif($ARGV[$k] eq 'DSSO') {
		$IQPIR_REPORTERMASSES = "182.105523:49.98264,214.077603:-13.96152";
		printf STDERR "Setting DSSO REPORTERMASSES to %s\n", $IQPIR_REPORTERMASSES; #exit(1);
		if(! $DSSO_STUMP_MASSES[0] && ! $DSSO_STUMP_MASSES[1]) {
			$XL_PROPH_OPTIONS .= "CROSSLINK_MODS=K:182.105523,K:214.077603 ";
		}
		elsif(! $DSSO_STUMP_MASSES[0]) {
			$XL_PROPH_OPTIONS .= "CROSSLINK_MODS=K:182.105523 ";
		}
		elsif(! $DSSO_STUMP_MASSES[1]) {
			$XL_PROPH_OPTIONS .= "CROSSLINK_MODS=K:214.077603 ";
		}
	}
	elsif($ARGV[$k] =~ /^SILAC\_MODS\=(\S+)/) {
		printf STDERR "Setting SILAC_MODS to %s\n", $1; #exit(1);
		$XL_PROPH_OPTIONS .= $ARGV[$k] . " ";
	}
}

$XL_PROPH_OPTIONS .= $STUMP_MOD if($IQPIR_REPORTERMASSES eq '' && ! ($STUMP_MOD eq ''));

if($light_heavy) {
	for(my $k = 0; $k < @light_files; $k++) {
		if($light_files[$k] =~ /interact/ || $light_files[$k] =~ /iproph/) {
			printf STDERR "Removing $light_files[$k]\n";
			@light_files = @light_files[0 .. $k-1, $k+1 .. $#light_files];
			$k--;
		}
		elsif(scalar keys %{$SAMPLEMAPFILES} > 0 && $light_files[$k] =~ /$DIR_DIVISOR_REGX([^$DIR_DIVISOR_REGX]+)\.pep.xml$/ && ! exists $SAMPLEMAPFILES->{$1}) {
			printf STDERR "Removing $light_files[$k] since not in samplemap file\n";
			@light_files = @light_files[0 .. $k-1, $k+1 .. $#light_files];
			$k--;
		}
	}
	checkForReact2OrPeaks(\@light_files) if(! $pairfile);
	for(my $k = 0; $k < @heavy_files; $k++) {
		if($heavy_files[$k] =~ /interact/ || $heavy_files[$k] =~ /iproph/) {
			printf STDERR "Removing $heavy_files[$k]\n";
			@heavy_files = @heavy_files[0 .. $k-1, $k+1 .. $#heavy_files];
			$k--;
		}
		elsif(scalar keys %{$SAMPLEMAPFILES} > 0 && $heavy_files[$k] =~ /$DIR_DIVISOR_REGX([^$DIR_DIVISOR_REGX]+)\.pep.xml$/ && ! exists $SAMPLEMAPFILES->{$1}) {
			printf STDERR "Removing $heavy_files[$k] since not in samplemap file\n";
			@heavy_files = @heavy_files[0 .. $k-1, $k+1 .. $#heavy_files];
			$k--;
		}
	}
	checkForReact2OrPeaks(\@heavy_files) if(! $pairfile);
	my $light = "interact";
	my $heavy = "interact";
	if(! ($suffix eq '')) {
		$light .= $suffix;
		$heavy .= $suffix;
	}
	$light .= "-light.pep.xml";
	$heavy .= "-heavy.pep.xml";
	
	
	my $light_command = "xinteract $PEP_PROPH_OPTIONS -nP -N$light " . join(" ", @light_files);
	my $heavy_command = "xinteract $PEP_PROPH_OPTIONS -nP -Eheavy -N$heavy " . join(" ", @heavy_files);

	die "No light files\n" if(@light_files == 0);
	die "No heavy files\n" if(@heavy_files == 0);
	
	my $command = $light_command . "; " . ($FORCE_NTT2 ? addRefreshForNtt2($light, $light_files[0])."; " : "") . 
		$heavy_command . "; " . ($FORCE_NTT2 ? addRefreshForNtt2($heavy, $heavy_files[0])."; " : "") . "InteractParser interact";

	if($IS_WINDOWS) {
		push(@COMMANDS, $light_command);
		push(@COMMANDS, addRefreshForNtt2($light, $light_files[0])) if($FORCE_NTT2);
		push(@COMMANDS, $heavy_command);
		push(@COMMANDS, addRefreshForNtt2($heavy, $heavy_files[0])) if($FORCE_NTT2);
		push(@COMMANDS, "InteractParser interact");
	}
	my $windows_suffix = "";
	if(! ($suffix eq '')) {
		$command .= $suffix;
		$windows_suffix .= $suffix;
	}
	$command .=".pep.xml $light $heavy";
	$windows_suffix .=".pep.xml $light $heavy";
	if(! ($MIN_PEP_LEN eq '')) {
		$command .= " -L$MIN_PEP_LEN";
		$windows_suffix .= " -L$MIN_PEP_LEN";
	}
	$COMMANDS[-1] .= $windows_suffix if($IS_WINDOWS && ! ($windows_suffix eq ''));
	my $inputfile = "interact.pep.xml"; #"iprophet.pep.xml";

	$command .= "; PeptideProphetParser interact".$suffix.".pep.xml EXPECTSCORE ACCMASS DECOYPROBS NONPARAM MINPROB=0 PPM DECOY=rev_";
	push(@COMMANDS, "PeptideProphetParser interact".$suffix.".pep.xml EXPECTSCORE ACCMASS DECOYPROBS NONPARAM MINPROB=0 PPM DECOY=rev_") if($IS_WINDOWS);

	if($RUN_IPROPH) {
		my $command2 = "InterProphetParser interact".$suffix.".pep.xml iprophet".$suffix.".pep.xml";
		$inputfile = "iprophet".$suffix.".pep.xml";
		$command .= "; " . $command2 if($all);
		push(@COMMANDS, $command2) if($IS_WINDOWS && $all);
	}
	$XL_PROPH_OPTIONS .= "LIGHT_HEAVY";
	if(! ($REPORTERMASS eq '')) {
		$XL_PROPH_OPTIONS .= " REPORTERMASS=$REPORTERMASS";
	}
	if(! ($IQPIR_REPORTERMASSES eq '')) {
		$XL_PROPH_OPTIONS .= " IQPIR_REPORTERMASSES=$IQPIR_REPORTERMASSES";
	}
	my $command3 = "$XLINKPROPHET $inputfile $XL_PROPH_OPTIONS";
	$command .= "; " . $command3 if($all);
	push(@COMMANDS, $command3) if($IS_WINDOWS && $all);
	if(! $exec) {
		my $clean_command = $command;
		if($IS_WINDOWS) {
			$clean_command = join("\n\n", @COMMANDS);
		}
		else {
			$clean_command =~ s/\; /\n\n/g;
		}
		printf "\n%s\n", $clean_command;
	}
	else {
		if($IS_WINDOWS) {
			foreach(@COMMANDS) {
				printf "%s\n", $_;
				system($_);
			}
		}
		else {	
			printf "%s\n", $command;
			system($command)
		}
	}
	exit(1);

}

for(my $k = 0; $k < @files; $k++) {
	if($files[$k] =~ /interact/ || $files[$k] =~ /iproph/) {
		printf STDERR "Removing $files[$k]\n";
		@files = @files[0 .. $k-1, $k+1 .. $#files];
		$k--;
	}
	elsif(scalar keys %{$SAMPLEMAPFILES} > 0 && $files[$k] =~ /$DIR_DIVISOR_REGX([^$DIR_DIVISOR_REGX]+)\.pep.xml$/ && ! exists $SAMPLEMAPFILES->{$1}) {
		printf STDERR "Removing $files[$k] since not in samplemap file\n";
		@files = @files[0 .. $k-1, $k+1 .. $#files];
		$k--;
	}
}

checkForReact2OrPeaks(\@files) if(! $pairfile);

die "No input files\n" if(@files == 0);

my $command = "xinteract $PEP_PROPH_OPTIONS ";
if(! ($suffix eq '')) {
	$command .= "-Ninteract".$suffix.".pep.xml ";
}

$command .= join(' ', @files);
push(@COMMANDS, $command);
if($FORCE_NTT2) {
	$command .= '; ' . addRefreshForNtt2("interact".$suffix.".pep.xml", $files[0]);
}


my $inputfile = "interact.pep.xml"; #"iprophet.pep.xml";
if($RUN_IPROPH) {
	my $command2 = "InterProphetParser interact".$suffix.".pep.xml iprophet".$suffix.".pep.xml";
	$inputfile = "iprophet".$suffix.".pep.xml";
	$command .= '; ' . $command2 if($all);
	push(@COMMANDS, $command2) if($IS_WINDOWS && $all);
}
if(! ($REPORTERMASS eq '')) {
	$XL_PROPH_OPTIONS .= " REPORTERMASS=$REPORTERMASS";
}
if(! ($IQPIR_REPORTERMASSES eq '')) {
	$XL_PROPH_OPTIONS .= " IQPIR_REPORTERMASSES=$IQPIR_REPORTERMASSES";
}

my $command3 = "$XLINKPROPHET $inputfile $XL_PROPH_OPTIONS";
$command .=  '; ' . $command3 if($all);
push(@COMMANDS, $command3) if($IS_WINDOWS && $all);

if(! $exec) {
	my $clean_command = $command;
	if($IS_WINDOWS) {
		$clean_command = join("\n\n", @COMMANDS);
	}
	else {
		$clean_command =~ s/\; /\n\n/g;
	}
	printf "\n%s\n", $clean_command;
}
else {
	if($IS_WINDOWS) {
		foreach(@COMMANDS) {
			printf "%s\n", $_;
			system($_);
		}
	}
	else {	
		printf "%s\n", $command;
		system($command)
	}
}

# returns whether pairfile used versus having ReACT or Mango results with react2.xls and peaks files, respectively
sub readParams {
(my $paramsfile) = @_;
open(PARAMS, $paramsfile) or die "cannot read $paramsfile $!\n";
my $parifile = 0;
while(<PARAMS>) {
	s/\#.*//; # strip off comments at right of each line
	if($SEND_EMAIL && /^email\_address \= (\S+)/) {
		$XL_PROPH_OPTIONS .= "EMAIL=$1 ";
		if($1 eq 'youremail@u.washington.edu') {
			printf STDERR "Error: you must put your own email address in your xlinkprophet params file, not $1!\n";
			exit(1);
		}
	}
	elsif(/^peptideprophet\_options \= (\S+.*\S)/) {
		$PEP_PROPH_OPTIONS = $1;
		if(/\-l(\d+)/ && $1 != 6) {
			$MIN_PEP_LEN = $1;
			$XL_PROPH_OPTIONS .= "MIN_PEPLEN=$MIN_PEP_LEN "; # pass to XLinkProphet
		}
	}
	elsif(/^run\_iprophet \= (\d)/) {
		$RUN_IPROPH = $1;
	}
	elsif(/^no\_nsx \= (\d)/ && $NO_NSX != $1) {
		$XL_PROPH_OPTIONS .= "NO_NSX ";
	}
	elsif(/^output\_decoys \= (\d)/ && $OUTPUT_DECOYS != $1) {
		$XL_PROPH_OPTIONS .= "OUTPUT_DECOYS ";
	}
	elsif(/^target\_neg \= (\d)/ && $TARGET_NEG != $1) {
		$XL_PROPH_OPTIONS .= "TARGET_NEG ";
	}
	elsif(/^spec\_contains \=\s*(\S+)/) {
		$XL_PROPH_OPTIONS .= "SPEC_CONTAINS=$1 ";
	}
	elsif(/^local\_react \= (\d)/ && $LOCAL_REACT != $1) {
		$XL_PROPH_OPTIONS .= "LOCAL_REACT ";
	}
	elsif(/^pairfile\_suff \=\s*(\S+)/) {
		$XL_PROPH_OPTIONS .= "PAIRFILE_SUFF=$1 ";
		$pairfile = 1;
	}
	elsif(/^xlinker\_stump\_mod\_masses \=\s*(\S+)/) {
		my $masses = $1;
		$STUMP_MOD = "CROSSLINK_MODS=$masses ";
		printf STDERR  "----> Setting $XL_PROPH_OPTIONS\n";
		if($masses !~ /\:/) {
			die "Error: xlinker_stump_mod_masses is in wrong format.  Need A:235.343\n";
		}
		my @stump_masses = split(",", $masses);
		for(my $k = 0;  $k < @stump_masses; $k++) {
			if(index($stump_masses[$k], "K:182.1")==0) {
				$DSSO_STUMP_MASSES[0] = 1; # found light
			}
			elsif(index($stump_masses[$k], "K:214.0")==0) {
				$DSSO_STUMP_MASSES[1] = 1; # found heavy
			}
		}
	}
	elsif(/^crosslink\_mods \=\s*(\S+)/) {
		$XL_PROPH_OPTIONS .= "CROSSLINK_MODS=$1 ";
		if($1 !~ /\:/) {
			die "Error: crosslnk_mods is in wrong format.  Need A:235.343\n";
		}
	}

}
close(PARAMS);
return $pairfile;
}


# makes sure there are either react2.xls or peaks files for each analyzed pepxml
sub checkForReact2OrPeaks {
(my $fileptr) = @_;
my $type = "unknown";
my $filter = 0; # whether to remove the pepxml files without react rather than exit
for(my $k = 0; $k < @{$fileptr}; $k++) {
	if($fileptr->[$k]=~ /^(\S+\.)pep.xml$/) {
		if($type eq 'unknown') {
			if(-e $fileptr->[$k] . '.react2.xls') {
				$type = "react";
			}
			elsif(-e $1 . 'peaks') {
				$type = "mango";
			}
			else {
				die "Error, no react2.xls or peaks file for $fileptr->[$k]\n";
			}
		}
		elsif($type eq 'react' && ! -e $fileptr->[$k] . '.react2.xls') {
			if($FILTER_PEPXML_WITHOUT_PAIRFILE) {
				printf "Error, no react2.xls $fileptr->[$k].react2.xls found file for $fileptr->[$k], so skipping....\n";
				@{$fileptr} = @{$fileptr}[0 .. $k-1, $k+1 .. $#{$fileptr}];
				$k--;
			}
			else {
				die "Error, no react2.xls $fileptr->[$k].react2.xls found file for $fileptr->[$k]\n";
			}
		}
		elsif($type eq 'mango' && ! -e $1 . 'peaks') {
			if($FILTER_PEPXML_WITHOUT_PAIRFILE) {
				printf "Error, no peaks file found for $fileptr->[$k], so skipping....\n";
				@{$fileptr} = @{$fileptr}[0 .. $k-1, $k+1 .. $#{$fileptr}];
				$k--;
			}
			else {
				die "Error, no peaks file $1peaks found for $fileptr->[$k]\n";
			}
		}
	}
	else {
		die "Unknown file format $fileptr->[$k]\n";
	}
}



}

sub getSamplemapFiles {
(my $file) = @_;
open(FILE, $file) or die "cannot read $file $!\n";
my %output = ();
while(<FILE>) {
	chomp();
	my @parsed = split("\t");
	$output{$parsed[1]}++;
}

close(FILE);
printf STDERR "Read in %d files from sample map file %s\n", scalar keys %output, $file;
return \%output;
}


sub writeParamsFile {
	my $outfile = getcwd() . $DIR_DIVISOR . "xlinkprophet.params";
	if(-e $outfile) {
		printf STDERR "xlinkprophet params file $outfile already exists.  Delete it and re-run runXLinkProphet.pl WRITEPARAMS if you want to write new one.\n";
		exit(1);
	}
	open(OUT, ">$outfile") or die "cannot write to $outfile $!\n";
	displayParamsFile(*OUT, "");
	close(OUT);
	printf STDERR "xlinkprophet.params file $outfile written.  Make sure to edit the file with your correct email address\n";

}

sub addRefreshForNtt2 {
(my $interact, my $inputpepxml) = @_;
# first get database
open GREP, "grep 'search_database local_path' $inputpepxml |";
my @results = <GREP>;
close(GREP);
die "Error: no database found in $inputpepxml\n" if(@results == 0);
chomp $results[0];
my $database = "";
if($results[0] =~ /local\_path\=\"(\S+?)\"/) {
	$database = $1;
	die "Error: database $database not found\n"if(! -e $database);
}
else {
	die "Error: no database found in $results[0]\n"
}
return "RefreshParser $interact $database 2";


}


sub displayParamsFile {
(my $stream, my $tab) = @_;
printf $stream $tab . "email_address = youremail\@u.washington.edu\n";
printf $stream "\n";
printf $stream "peptideprophet_options = -OEAdP -p0 -PPM -l5 -drev_             # recommended use of Comet expect score with nonparametric modeling using decoys with protein names beginning with rev_ and excluding results with peptides of length less than 6
run_iprophet = 1        # whether to run iProphet after PeptideProphet to improve discrimination\n";
printf $stream "\n";
printf $stream $tab . "no_nsx = 0              #  NO_NSX (Do not use NSX [number of sibling cross-links] distribution to compute probabilities)\n";
printf $stream $tab . "output_decoys = 0       # OUTPUT_DECOYS (Output along with targets and include when computing FDRs)\n";
printf $stream $tab . "target_neg = 0          # (Do not use decoys to set negative distributions)\n";
printf $stream $tab . "spec_contains =         #  (Filter input data and use only spectra matching xxx)\n";
printf $stream $tab . "local_react = 1         # (Use React2.xls files in current directory rather than those referenced in the inputfile)\n";
printf $stream $tab . "pairfile_suff =         # PAIRFILE_SUFF=xxx (Suffix appended to search result pepXML file names for pairing files containing spectrum scan numbers of 2 peptides plus mass and charge of crosslink parent)\n";
printf $stream $tab . "xlinker_stump_mod_masses =        #  CROSSLINK_MODS=x:202,y:329,xxx (supplement default values of n:198 and K325 with specified values for n/c termini or AA x,y,z...)\n";
}

# example contents of xlinkprophet.params file
# email_address = xxx@u.washington.edu

# peptideprophet_options = -OEAdP -p0 -PPM -l6 -drev_
# run_iprophet = 1
# 
# no_nsx = 0  #  NO_NSX (Do not use NSX [number of sibling cross-links] distribution to compute probabilities)\n";
# output_decoys = 0   # OUTPUT_DECOYS (Output along with targets and include when computing FDRs)\n";
# target_neg = 0  # (Do not use decoys to set negative distributions)\n";
# spec_contains =   #  (Filter input data and use only spectra matching xxx)\n";
# output_suff =    #  (Add xxx onto typical output file names)\n";
# local_react = 0  # (Use React2.xls files in current directory rather than those referenced in the inputfile)\n";
# pairfile_suff =    # PAIRFILE_SUFF=xxx (Suffix appended to search result pepXML file names for pairing files containing spectrum scan numbers of 2 peptides plus mass and charge of crosslink parent)\n";
# crosslink_mods =   #  CROSSLINK_MODS=x:202,y:329,xxx (supplement default values of n:198 and K325 with specified values for n/c termini or AA x,y,z...)\n";
