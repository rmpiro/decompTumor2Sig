#!/usr/bin/perl -w

use strict;
use Getopt::Std;

# turn off buffering of STDOUT
$| = 1;

my $progname = "filterLines.pl";
my $delimiter = "\t";  # default delimiter between columns is <tab>
my $columnNo = -1;
my $positiveFilters = 1; # 1 = one of the filter terms must be present(default)
                         # 0 = NONE of the filter terms must be present
my $partialFilters = 0; # default: complete filters
my $filterTermsFile = "";
my $verboseOutput = 0;

# -------------------
# getting options ...
# -------------------

my %opts;
$opts{h} = "";
$opts{d} = "";
$opts{c} = "";
$opts{f} = "";
$opts{N} = "";
$opts{P} = "";
$opts{v} = "";

getopts('hc:f:NPd:v', \%opts);

if ( $opts{h} eq "1" ) {
    &printUsageInfo();
}

if ( $opts{d} ne "" ) {
    $delimiter = $opts{d};
}

if ( $opts{N} eq "1" ) {
    $positiveFilters = 0; # use as negative filter
}

if ( $opts{f} ne "" ) {
    $filterTermsFile = $opts{f};
}

if ( $opts{P} eq "1" ) {
    $partialFilters = 1; # allow partial filters
}

if ( $opts{v} eq "1" ) {
    $verboseOutput = 1;
}

if ( $opts{c} ne "" && $opts{c}=~/^([\d]*)$/) {
    $columnNo = $1;
}

if ($columnNo <= 0) {
    print STDERR "Error: column number not specified, or not specified correctly (-h for help).\n";
    exit(1);
}

if ($filterTermsFile eq "") {
    print STDERR "Error: File containing the filter terms not specified (-h for help).\n";
    exit(2);
}


# first read the filter terms:
my %filterTerms = ();

if ($verboseOutput) {
    print "# Reading filter terms from $filterTermsFile\n";
}

unless ( open (FTFILE, "< $filterTermsFile") ) {
    print STDERR "Error: cannot open input file $filterTermsFile!\n";
    exit(3);
}

my $line = "";
my $lineNum = 0;
while ($line = <FTFILE>) {

    chomp($line);
    $lineNum++;

    next if ($line =~ /^\#/ || $line =~ /^!/);

    if ($verboseOutput) {
	print "# Filter $lineNum: '$line'\n";
    }

    $filterTerms{$line} = 1;
}

close (FTFILE);


# now filter input

$lineNum = 0;

while ($line = <STDIN>) {

    chomp($line);
    $lineNum++;

    if ($line =~ /^\#/ || $line =~ /^!/) {
	print "$line\n";
	next; # don't filter comments, just print them
    }

    if ($verboseOutput) {
	print "# line '$line'\n ";
    }

    my @cols = split (/$delimiter/, $line);

    if (scalar(@cols) < $columnNo) {
	print "$line\n"; # just print it, no filtering!
	next;
    }

    my $matches = 0;

    if (!$partialFilters) {
	if (exists($filterTerms{$cols[$columnNo-1]})) {
	    if ($verboseOutput) {
		print "# line $lineNum  with column $columnNo='$cols[$columnNo-1]' matches filter $filterTerms{$cols[$columnNo-1]}\n";
	    }

	    $matches = 1;
	}
    }
    else {
	# partial filers, this is slower!
	foreach my $filter (keys(%filterTerms)) {

	    if($cols[$columnNo-1] =~ /$filter/) {
		# this column content satisfies this filter.

		if ($verboseOutput) {
		    print "# line $lineNum  with column $columnNo='$cols[$columnNo-1]' matches filter $filter\n";
		}

		$matches = 1;
		last;
	    }

	}
    }

    if (!$matches) {
	if ($verboseOutput) {
	    print "# no matching filter for line $lineNum  with column $columnNo='$cols[$columnNo-1]'\n";
	}
    }

    if ( ($positiveFilters && $matches)
	 || (!$positiveFilters && !$matches)) {
	print "$line\n";
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
    print "  -c <col>     Number of the column that should be checked for the presence of\n";
    print "               the filter terms specified in <filters>. Each line having less\n";
    print "               then <col> fields, will simply be printed without filtering.\n";
    print "  -f <filters> File containing the filters for which to check column <col>.\n";
    print "               Must contain one filter per line. Can contain Perl-compatible\n";
    print "               regular expressions (e.g. \"ENSG[0-9]*\"). Filters must span\n";
    print "               the entire column content (e.g. filter \"ENSG\" will not\n";
    print "               be recognized as part of column content \"ENSG00000101391\").\n";
    print "  -P           Allow for partial filters. If not specified filters must span\n";
    print "               the entire column content (e.g. partial filter \"ENSG\" will\n";
    print "               not be recognized as part of column content \"ENSG00000101391\").\n";
    print "               Important: this can be much slower!\n";
    print "  -N           If specified only lines NOT containing any of the filter terms\n";
    print "               at the field of column <col> will be printed (negative filters).\n";
    print "               Otherwise only lines containing one of the filter terms will be\n";
    print "               printed (positive filters; default). Commented lines (staring\n";
    print "               with \# or !) will always be printed instead of being filtered.\n";
    print "  -v           Verbose standard output.\n\n";
    print "Mandatory parameters: -c\n";

    exit 0;
}
