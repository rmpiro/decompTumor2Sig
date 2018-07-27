#!/usr/bin/perl -w

use strict;
use Getopt::Std;

# turn off buffering of STDOUT
$| = 1;

my $progname = "extractSpecColumns.pl";
my $delimiter = "\t";  # default delimiter between columns is <tab>
my $columnNoString = "";
my $verboseOutput = 0;

# -------------------
# getting options ...
# -------------------

my %opts;
$opts{h} = "";
$opts{d} = "";
$opts{c} = "";
$opts{v} = "";

getopts('hc:d:v', \%opts);

if ( $opts{h} eq "1" ) {
    &printUsageInfo();
}

if ( $opts{d} ne "" ) {
    $delimiter = $opts{d};
}

if ( $opts{v} eq "1" ) {
    $verboseOutput = 1;
}

#if ( $opts{c} ne "" && $opts{c}=~/^[\d,\-]*$/) {
if ( $opts{c} ne "") {
    $columnNoString = $opts{c};
}

if ($columnNoString eq "") {
    print STDERR "Error: column numbers not specified, or not specified correctly (-h for help).\n";
    exit(1);
}


my @colIndices = ();

my $needHeader = 0; # if we have also non numeric column names!
my $haveHeader = 0; # if we have also non numeric column names!
my %colName2Index = (); # if we have also non numeric column names!
my @requiredColNames = (); # if we have also non numeric column names!

my @indexTokens = split (/,/, $columnNoString);


my $maxColIndex = 0;

# parse column string!
foreach my $iTok (@indexTokens) {

    if ($iTok =~ /^(\d+)$/) {  # one integer -> single column

	my $cI = $1;
	push(@colIndices, $cI);

	if ($verboseOutput) {
	    print "# extracting column with index: $cI\n";
	}
	if ($cI > $maxColIndex) {
	    $maxColIndex = $cI;
	}
    }
    elsif ($iTok =~ /^(\d+)-(\d+)$/) {  # range -> multiple columns
	my $first = $1;
	my $last = $2;

	if ($verboseOutput) {
	    print "# extracting column range with indices: ";
	}

	my $incr = 1;
	if ($last < $first) {
	    $incr = -1;
	}

	my $cI = $first;
	for ($cI = $first; $incr*$cI <= $incr*$last; $cI = $cI+$incr) {
	    push(@colIndices, $cI);

	    if ($verboseOutput) {
		if ($cI != $first) {
		    print ", ";
		}
		print "$cI";
	    }
	    if ($cI > $maxColIndex) {
		$maxColIndex = $cI;
	    }
	}

	if ($verboseOutput) {
	    print "\n";
	}
	
    } else {
	# this isn't numeric; interpret it as column name to be confornted with header!
	$needHeader = 1;
	push(@colIndices, $iTok);
	push(@requiredColNames, $iTok);
    }

}



if ($verboseOutput) {
    print "# maximum column index: $maxColIndex\n";
}

my $lineNum = 0;

while (<STDIN>) {
    my $line = $_;
    chomp($line);

    $lineNum++;

    next if ($line =~ /^\#/);

    if ($verboseOutput) {
	print "# line '$line'\n";
    }

    my @cols = split (/$delimiter/, $line);

    if ($needHeader && !$haveHeader) {
	# we need a header and this is the first non-commented line!
	if ($verboseOutput) {
	    print "# we have at least one column name specified; parsing the header for mapping to indices:\n";
	}

	my $hii;
	for ($hii = 0; $hii < scalar(@cols); $hii++) {
	    $colName2Index{$cols[$hii]} = $hii+1;
	    if ($verboseOutput) {
		print "# $colName2Index{$cols[$hii]}=<$cols[$hii]>\n";
	    }
	}

	$haveHeader = 1;
    }


    #if (scalar(@cols) < $maxColIndex) {
	#print STDERR "Error: input line $lineNum has less than $maxColIndex columns! Quitting ...\n";
	#exit(2);
    #}

    # if we have column names specified; check if we have indices for them!
    foreach my $rCN (@requiredColNames) {
	if (!exists($colName2Index{$rCN})) {
	    print STDERR "Error: column name '$rCN' not found in the header!\n";
	    exit(77);
	}
    }
    

    my $outLine = "";

    my $first = 1;

    foreach my $cI (@colIndices) {

	next if ($cI eq "" || ($cI =~ /^\d+$/ && $cI == 0));

	if (!$first) {
	    $outLine .= $delimiter; # add another delimiter
	}

	if ($cI !~ /^\d+$/) {
	    # not an integer; must be a column name; convert it to the column index!
	    $cI = $colName2Index{$cI};
	}

	if (scalar(@cols) >= $cI) {
	    $outLine .= $cols[$cI-1]; # add column content
	}

	$first = 0;
    }

    print "$outLine\n";
}

exit 0;


# ---------------------------------------------------------------------------
#
# Print usage information (help).
# arguments: none
#
# ---------------------------------------------------------------------------
sub printUsageInfo() {
    print "Author: Rosario M. Piro (piro\@to.infn.it)\n\n";
    print "Usage: cat <inFile> | $progname [OPTIONS]\n\n";
    print "  Where OPTIONS are:\n";
    print "  -h           Display this help message and exit.\n";
    print "  -d <delim>   Delimiter between columns of both output and input file.\n";
    print "               The default delimiter is <tab>.\n";
    print "  -c <cols>    Comma-serarated string with column numbers to extract (the\n";
    print "               first column of the input file corresponds to column number 1).\n";
    print "               Also column ranges can be specified. Example: \"-c 1,3-5,9-7,2\"\n";
    print "               will extract the columns in the order 1,3,4,5,9,8,7,2.\n";
    print "               Instead a column numbers also a column name as encoutnered in the\n";
    print "               header can be specified. In this case no ranges will be considered,\n";
    print "               because a minus symbol will be interpreted as part of the column name!\n";
    print "               Example: \"-c '1-3,5,Gene ID,Chromosome Name'\"\n";
    print "  -v           Verbose standard output.\n\n";
    print "Mandatory parameters: -c\n";

    exit 0;
}
