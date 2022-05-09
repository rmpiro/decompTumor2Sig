#!/usr/bin/perl -w

use strict;
use Getopt::Std;

# turn off buffering of STDOUT
$| = 1;

my $FT_GT = 1; # filter type: greater than
my $FT_GE = 2; # filter type: greater or equal (lower limit)
my $FT_LT = 3; # filter type: lower than
my $FT_LE = 4; # filter type: lower or equal (upper limit)
my $FT_EQ = 5; # filter type: equal
my $FT_NE = 6; # filter type: not euqal


my $progname = "filterLinesNumeric.pl";
my $delimiter = "\t";  # default delimiter between columns is <tab>
my $columnNo = -1;
my $filterType = $FT_GE; # default: greater or equal (lower limit)
my $filterValue = 0;
my $printNonNumericLines = 0;
my $printLineNum = 0;
my $verboseOutput = 0;

# -------------------
# getting options ...
# -------------------

my %opts;
$opts{h} = "";
$opts{d} = "";
$opts{c} = "";
$opts{T} = "";
$opts{V} = "";
$opts{p} = "";
$opts{n} = "";
$opts{v} = "";

getopts('hd:c:T:V:npv', \%opts);

if ( $opts{h} eq "1" ) {
    &printUsageInfo();
}

if ( $opts{d} ne "" ) {
    $delimiter = $opts{d};
}

if ( $opts{T} ne "" ) {
    if ($opts{T} =~ /^(gt)$/i) {
	$filterType = $FT_GT;
    }
    elsif ($opts{T} =~ /^(ge)$/i) {
	$filterType = $FT_GE;
    }
    elsif ($opts{T} =~ /^(lt)$/i) {
	$filterType = $FT_LT;
    }
    elsif ($opts{T} =~ /^(le)$/i) {
	$filterType = $FT_LE;
    }
    elsif ($opts{T} =~ /^(eq)$/i) {
	$filterType = $FT_EQ;
    }
    elsif ($opts{T} =~ /^(ne)$/i) {
	$filterType = $FT_NE;
    }
    else {
	print STDERR "Error: unknown filter type '$opts{T}' specified (-h for help)!\n";
	exit(1);
    }
}

if ( $opts{V} ne "" ) {
    if ($opts{V} =~ /^-?[\d]+(\.[\d]+)?([eE]-?[\d]+)?$/) {
	$filterValue = $opts{V};
    }
    else {
	print STDERR "Error: the threshold/filter value must be numeric (-h for help)!\n";
	exit(1);
    }
}

if ( $opts{p} eq "1" ) {
    $printNonNumericLines = 1;
}

if ( $opts{n} eq "1" ) {
    $printLineNum = 1;
}

if ( $opts{v} eq "1" ) {
    $verboseOutput = 1;
}

if ( $opts{c} ne "" && $opts{c}=~/^([\d]*)$/) {
    $columnNo = $1;
}

if ($columnNo <= 0) {
    print STDERR "Error: column number not specified, or not specified correctly (-h for help)!\n";
    exit(1);
}

if ($opts{V} eq "") {
    print STDERR "Error: the threshold/filter value must be specified (-h for help)!\n";
    exit(2);
}


# now filter input

my $line = "";
my $lineNum = 0;

while ($line = <STDIN>) {

    chomp($line);
    $lineNum++;

    if ($line =~ /^\#/ || $line =~ /^!/) { # commented lines
	if ($printNonNumericLines) {
	    if ($printLineNum) {
		print "$lineNum:";
	    }
	    print "$line\n";
	}
	next;
    }

    if ($verboseOutput) {
	print "# line '$line'\n ";
    }

    my @cols = split (/$delimiter/, $line);

    if (scalar(@cols) < $columnNo) {
	if ($printNonNumericLines) {
	    if ($printLineNum) {
		print "$lineNum:";
	    }
	    print "$line\n"; # just print it, no filtering!
	}
	else {
	    print "# Ignoring short line $lineNum (less than $columnNo columns).\n";
	}
	next;
    }


    my $matches = 0;

    if ($cols[$columnNo-1] =~ /^-?[\d]+(\.[\d]+)?([eE]-?[\d]+)?$/) {
	# this is a numeric field, we can check if for the threshold!

	if ( ($filterType == $FT_GT && $cols[$columnNo-1] > $filterValue)
	     || ($filterType == $FT_GE && $cols[$columnNo-1] >= $filterValue)
	     || ($filterType == $FT_LT && $cols[$columnNo-1] < $filterValue)
	     || ($filterType == $FT_LE && $cols[$columnNo-1] <= $filterValue)
	     || ($filterType == $FT_EQ && $cols[$columnNo-1] == $filterValue)
	     || ($filterType == $FT_NE && $cols[$columnNo-1] != $filterValue)
	     ) {

	    if ($printLineNum) {
		print "$lineNum:";
	    }
	    print "$line\n"; # satisfies filter/threshold!
	}
    }
    else {
	if ($printNonNumericLines) {
	    if ($printLineNum) {
		print "$lineNum:";
	    }
	    print "$line\n"; # user wants to print non-numeric content
	}
	else {
	    if ($verboseOutput) {
		print "# Ignoring non-numeric content of column $columnNo in line $lineNum.\n";
	    }
	}
    }

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
    print "  -c <col>     Number of the column that should be checked against the filter\n";
    print "               value (threshold) <thres> using the filter type (comparison\n";
    print "               oprator) <type>.\n";
    print "  -V <thres>   Filter value (threshold).\n";
    print "  -T <type>    Filter type (comparison operator); can be one of:\n";
    print "                 \"GT\" (greater than; ok if column content > threshold)\n";
    print "                 \"GE\" (greater or equal; ok if column content >= threshold)\n";
    print "                 \"LT\" (less than; ok if column content < threshold)\n";
    print "                 \"LE\" (less or equal; ok if column content <= threshold)\n";
    print "                 \"EQ\" (equal; ok if column content == threshold)\n";
    print "                 \"NE\" (not equal; ok if column content <> threshold)\n";
    print "               Default: \"GE\".\n";
    print "  -n           Print line numbers in grep -n style (i.e. \"<num>:<line>\").\n";
    print "  -p           Print lines also if they have a non-numeric content in\n";
    print "               column <col>, if they contain less columns, or if they are\n";
    print "               commented. Default: ignore such lines.\n";
    print "  -v           Verbose standard output.\n\n";
    print "Mandatory parameters: -c, -V.\n\n";

    exit 0;
}
