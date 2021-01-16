#!/usr/bin/perl


use strict;


my $SELF = 'convertXLinkProphetXLStoUpload7col.pl';

if(@ARGV == 0) {
	printf STDERR "\n";
	printf STDERR " usage:   $SELF < XLinkProphet output XLS file > < minimum probability > < xlinkprophet.params file > (options)\n";
	printf STDERR " options: USE_PROBABILITY (Use probability rather than composite_probability [non-redundant] as filter criteria)\n";
	printf STDERR "          FASTA=/full/path/to/searchdatabase.fasta (Get crosslink modification position in protein sequence and output to 7th and 8th columns for first and second peptide, respectively)\n\n";
	printf STDERR "          Use to filter XLinkProphet output XLS file based on minimum probability to produce 7-column file [ pepA, proA, kposA, pepB, proB, kposB, non-redundant probability ] ready for upload to XLinkDB\n\n";
	printf STDERR " example: $SELF iprophet-xl.xls 0.9 xlinkprophet.params\n";
	printf STDERR "          $SELF iprophet-xl.xls 0.9 xlinkprophet.params FASTA=C:\\Databasses\\database.fasta\n";
	printf STDERR "\n";
	exit(1);
 
}
my @MODMASSES = ();
(my $STUMP_MODAAS, my $STUMP_MODMASSES) = readParams($ARGV[2]);
printf "Using stump modmasses %s and amino acids %s\n", join(",", @$STUMP_MODMASSES), join(",", @$STUMP_MODAAS); #exit(1);

my $USE_PROBABILITY = 0;
my $FASTA = "";
my $proteinSeqs = {};
for(my $k = 3; $k < @ARGV; $k++) {
	if($ARGV[$k] eq 'USE_PROBABILITY') {
		printf "Using probability to filter results\n";
		$USE_PROBABILITY = 1;
	}
	elsif($ARGV[$k] =~ /^FASTA\=(\S+)/) {
		$FASTA = $1;
		$proteinSeqs = readFasta($FASTA);
		printf "Using fasta $FASTA to compute peptide cross-link protein residue positions and output as 2 additional columns\n";
	}
}

parseXLinkProphetXls($ARGV[0], $ARGV[1]); exit(1);

sub readFasta {
(my $fasta) = @_;
my %protSeqs = ();
my $seq = "";
my $prot = "";
#printf "Here with $peptide and $protein for $FASTA\n";
open(FASTA, $FASTA) or die "cannot read fasta $FASTA $!\n";
while(<FASTA>) {
	chomp();
	if(/^\>(\S+)/) { 
		$protSeqs{$prot} = $seq if($prot !~ /^rev\_/ && ! ($seq eq ''));
		$prot = $1;
		$seq = ""; # signal start
		#printf "FOUnd protein at $_\n";
	}
	else {
		$seq .= $_;
	}

}

close(FASTA);
$protSeqs{$prot} = $seq if(! ($seq eq ''));
printf STDERR "Read in sequences for %d proteins from fasta %s\n", scalar keys %protSeqs, $fasta;
return \%protSeqs;

}

sub getPeptideCrosslinkModPosition {
(my $peptide, my $kpos, my $protein) = @_;
$peptide = stripMods($peptide);

return index($proteinSeqs->{$protein}, $peptide)+ $kpos if(exists $proteinSeqs->{$protein});
printf "Error: protein $protein for peptide $peptide with kpos$kpos not found in $FASTA\n";
return "";
}


sub readParams {
(my $paramsfile) = @_;
my @stump_masses = ();
my @stump_aas = ();
open(PARAMS, $paramsfile) or die "cannot read $paramsfile $!\n";
while(<PARAMS>) {

	if(/^xlinker\_stump\_mod\_masses \=\s*(\S+)/ && $1 !~ /\#/) {
		@stump_masses = split(",", $1);
	}
}
close(PARAMS);
for(my $k = 0; $k < @stump_masses; $k++) {
	my @next = split(":", $stump_masses[$k]);
	push(@stump_aas, $next[0]);
	$stump_masses[$k] = sprintf("%0.2f", $next[1]); # just in case
}
push(@stump_aas, "K");
push(@stump_masses, "325.13");
return (\@stump_aas, \@stump_masses);
}

sub stripMods {
(my $modpep) = @_;
	$modpep =~ s/\[.*?\]//g;
	return $modpep;
}

sub addPepKposX {
(my $pep) = @_;
	for(my $k = 0;  $k < @$STUMP_MODAAS; $k++) {
		$pep =~ s/$STUMP_MODAAS->[$k]\[$STUMP_MODMASSES->[$k]\]/$STUMP_MODAAS->[$k]x/g;
	}
	return $pep;
}

sub getKpos {
(my $pepWithX) = @_;
return index($pepWithX, "x")-1;
}

sub getUniprot {
(my $pro) = @_;
my @next = split('\|', $pro);
return $next[1];
}

sub parseXLinkProphetXls {
(my $file, my $filter_thresh) = @_;
my $outfile = $file;
if($file =~ /(\S+)\.xls/) {
	$outfile = $1 . "_7col.txt";
}
else {
	$outfile = $file . "_7col.txt";
}
my %headers = ("probability" => 0, "protein1" => 5, "protein2" => 6, "composite_probability" => 15, "composite_id" => 16, "peptide1" => 3, "peptide2" => 4);
my $first = 1;
my $tot = 0;
my %seen = (); # in case by prob
open(FILE, $file) or die "Error: cannot read $file $!\n";
open(OUT, ">$outfile") or die "Error: cannot write to $outfile $!\n";
my %seen = (); # don't re-write the same cross-link, for example from light and heavy
while(<FILE>) {
	chomp();
	my @parsed = split("\t");
	if($first) {
		$first = 0;
		foreach my $head (keys %headers) {
			if(! $parsed[$headers{$head}] eq $head) {
				printf STDERR "Error: expecting column %d header to be %s, not %s\n", $headers{$head}+1, $head, $parsed[$headers{$head}];
				exit(1);
			}
		}
	}
	else {
		next if($USE_PROBABILITY && (exists $seen{$parsed[$headers{"composite_id"}]} || $parsed[$headers{"probability"}] < $filter_thresh));
		next if(! $USE_PROBABILITY && ($parsed[$headers{"composite_probability"}] eq '' || $parsed[$headers{"composite_probability"}] < $filter_thresh));
		$seen{$parsed[$headers{"composite_id"}]}++;
		my $pepA = addPepKposX($parsed[$headers{"peptide1"}]);
		my $pepB = addPepKposX($parsed[$headers{"peptide2"}]);
		

		for(my $k = 0; $k < @$STUMP_MODAAS; $k++) {
			$parsed[$headers{"peptide1"}] =~ s/$STUMP_MODAAS->[$k]\[$STUMP_MODMASSES->[$k]\]/$STUMP_MODAAS->[$k]/g;
			$parsed[$headers{"peptide2"}] =~ s/$STUMP_MODAAS->[$k]\[$STUMP_MODMASSES->[$k]\]/$STUMP_MODAAS->[$k]/g;
		}
		my $kposA = getKpos(stripMods($pepA));
		die "Error: no stump modification found in $pepA, check xlinkprophet.params to make sure contains all stump modifications used\n" if($kposA < 0);
		my $kposB = getKpos(stripMods($pepB));
		die "Error: no stump modification found in $pepB, check xlinkprophet.params to make sure contains all stump modifications used\n" if($kposB < 0);
		
		$parsed[$headers{"protein1"}] = $1 if($parsed[$headers{"protein1"}] =~ /^([^,]+)\,/);
		$parsed[$headers{"protein2"}] = $1 if($parsed[$headers{"protein2"}] =~ /^([^,]+)\,/);
		
		
		my $extra_cols = "";
		if(! ($FASTA eq '')) {
			$extra_cols .= "\t";
			$extra_cols .= getPeptideCrosslinkModPosition($parsed[$headers{"peptide1"}], $kposA, $parsed[$headers{"protein1"}]);
			$extra_cols .= "\t";
			$extra_cols .= getPeptideCrosslinkModPosition($parsed[$headers{"peptide2"}], $kposB, $parsed[$headers{"protein2"}]);
		}
		printf OUT "%s\t%s\t%d\t%s\t%s\t%d\t%s%s\n", $parsed[$headers{"peptide1"}], getUniprot($parsed[$headers{"protein1"}]), $kposA, $parsed[$headers{"peptide2"}], getUniprot($parsed[$headers{"protein2"}]), 
			$kposB, $parsed[$headers{"composite_probability"}], $extra_cols; #, $parsed[$headers{"composite_id"}];

		$tot++;
	}
}
close(FILE);
close(OUT);
printf "%d cross-links in %s with %s %s or greater written to 7-column file %s ready for upload to XLinkDB\n", $tot, $file, $USE_PROBABILITY ? "probability" : "composite_probability", $filter_thresh, $outfile;

}


