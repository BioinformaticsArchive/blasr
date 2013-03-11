#!/usr/bin/env perl

use RunConfiguration;
use MappingStudy;
require ParameterNames;

if ($#ARGV < 0) {
		print "usage: RunMappingStudy.pl config_file.txt [jumploc]\n";
	  print "jumploc: \n";
	  print "  1,sa       Start with suffix array construction.\n";
	  print "  2,setup    Skip to read file generation and error simulation.\n";
    print "  3,align    Skip to alignment of reads.\n";
    print "  4,tabulate Skip to parsing alignments to gather results.\n";
		exit(1);
}
$configFileName = shift @ARGV;

		
&RunConfiguration::ReadConfigFile($configFileName, \%config);

#
# Check for required values.
#
&RunConfiguration::CheckRequiredParameter(\%config, $referenceParam);
&RunConfiguration::CheckRequiredParameter(\%config, $errorRatesParam);
&RunConfiguration::CheckRequiredParameter(\%config, $lengthsParam);
&RunConfiguration::CheckRequiredParameter(\%config, $mapperParam);
&RunConfiguration::CheckRequiredParameter(\%config, $runParam);
&RunConfiguration::CheckRequiredParameter(\%config, $setupParam);
&RunConfiguration::CheckRequiredParameter(\%config, $summarizeParam);

$setupScript     = $config{$setupParam};
$runScript       = $config{$runParam};
$summarizeScript = $config{$summarizeParam};

$reference   = $config{$referenceParam};
$suffixArray = $config{$suffixArrayParam};
$ctab        = $config{$ctabParam};
$tophitOut   = "$reference.tophit";
$sumOut      = "$reference.sum";

if (exists $config{$suffixArrayParam}) {
		$suffixArray = $config{$suffixArrayParam};
}


$jumploc = 0;
if ($#ARGV >= 0) {
		$jumploc = shift @ARGV;
    if ($jumploc == 1 or $jumploc eq "sa") {
        goto STEP1;
    }
    elsif ($jumploc == 2 or $jumploc eq "setup") {
        goto STEP2;
    }
    elsif ($jumploc == 3 or $jumploc eq "align") {
        goto STEP3;
    }
    elsif ($jumploc == 4 or $jumploc eq "tabulate") {
        goto STEP4;
    }
    else {
        print "ERROR! The step to start on must be 1,2,3,or 4, or name 'sa', \n'setup', 'align', or 'tabulate'\n";
        exit(1);
    }
}
#
# Build the suffix array if it does not exist
#
STEP1:
if (! -e $ctab) {
		&RunConfiguration::CheckRequiredParameter(\%config, $buildSuffixArrayParam);
    $ctabCmd= $config{$buildCountTableParam};
    $buildCtabCommand = "$ctabCmd $ctab 8 $reference";
    system($buildCtabCommand);
}
if (! -e $suffixArray) {
		$sawriter = $config{$buildSuffixArrayParam};
		$buildSuffixArrayCommand = "$sawriter $suffixArray $reference -larsson -blt 8";
		system($buildSuffixArrayCommand);
}

#
# Build the count table if it does not exist
#

if (! -e $suffixArray ) {
		print "ERROR, could not locate or build the suffix array.\n";
		print "The command to build the suffix array was: '$buildSuffixArrayCommand'\n";
		exit(1);
}


#
# Set up the files.
#
if ($jumploc > 0 || $jumploc eq "sa") {
    # only do the step specified
    exit(1)
}


STEP2:
system("$setupScript $configFileName");
if ($? != 0) {
		print "ERROR, running $setupScript\n";
		exit(1);
}

#
# Run the mapping
#

if ($jumploc > 0 || $jumploc eq "setup")  {
    # only do the step specified
    exit(1)
}
STEP3:
system("$runScript $configFileName");
if ($? != 0) {
		print "ERROR, running $runScript\n";
		exit(1);
}

if ($jumploc > 0 || $jumploc eq "align") {
    # only do the step specified
    exit(1)
}
STEP4:
print "computing summaries \n";
if (exists $config{$tophitResultParam} ) {
		$tophitOut = $config{$tophitResultParam};
}
if (exists $config{$sumResultParam} ) {
		$sumOut = $config{$sumResultParam};
}

$summarizeTopHitsCmd = "$summarizeScript $configFileName tophit > $tophitOut";
#$summarizeSumCmd     = "$summarizeScript $configFileName sum > $sumOut";
print "$summarizeTopHitsCmd\n";
system($summarizeTopHitsCmd);
#system($summarizeSumCmd);
