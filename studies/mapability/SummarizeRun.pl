#!/usr/bin/env perl

use MappingStudy;
use RunConfiguration;
require ParameterNames;
if ($#ARGV < 0) {
		print "ERROR, you must provide a configuration file.\n";
		exit(1);
}

$configFileName = shift @ARGV; 
$config = {};

# 
# Set up some defaults.
#

$config{$rawReadsDirParam} = "raw";
$config{$errorPrefixParam} = "Error";
$config{$lengthPrefixParam} = "Length";
$config{$readFileSuffixParam} = ".fasta";
$config{$errorSimParam} = "errorSim";
$config{$nProcParam}     = 1;
$config{$summaryParam} = "tophit";
$config{$alignSuffixParam} = "";

&RunConfiguration::ReadConfigFile($configFileName, \%config);

if ($#ARGV >= 0) {
		$config{$summaryParam} = shift @ARGV;
}
#
# Check for required values.
#
&RunConfiguration::CheckRequiredParameter(\%config, $errorRatesParam);
&RunConfiguration::CheckRequiredParameter(\%config, $lengthsParam);
&RunConfiguration::CheckRequiredParameter(\%config, $referenceParam);
&RunConfiguration::CheckRequiredParameter(\%config, $summaryParam);

#
# parse sample line: 0|33271|33571 0  -1453 33250 33577 0 303
#

@errors = split(/\s+/, $config{$errorRatesParam});
@lengths = split(/\s+/, $config{$lengthsParam});

@numMappedReads = ();
@numTopOverlaps = ();
@numNotTopOverlaps = ();
@numReads    = ();
$i = 0;
print "e l ratio\n";
for ($e = 0; $e < scalar @errors; $e++ ){
		@numMappedReads[$e] = ();
		$numTopOverlaps[$e] = ();
		$numNotTopOverlaps[$e] = ();
		@numReads[$e] = ();
		for ($l = 0; $l < scalar @lengths; $l++) {
				$dirName   = $config{$errorPrefixParam} . "_" . "$errors[$e]" . "_" . $config{$lengthPrefixParam} . "_" . "$lengths[$l]";
				$readsFile = "$dirName/reads$config{$readFileSuffixParam}";
				$alignFile = "$readsFile.$config{$alignSuffixParam}";
				$res = `/home/UNIXHOME/mchaisson/scripts/PrintMissing.py $readsFile $alignFile -summary`;
				@resarray = split(/\s+/, $res);
        $ratio = @resarray[1] / @resarray[0];
        print "$i $errors[$e] $lengths[$l] $ratio\n";
        $i++;
#				$numReads[$e][$l] = @resarray[0];
#				$numIncorrect[$e][$l] = @resarray[1];
		}
}

#@summary = ();
#
#
#		for ($e = 0; $e < scalar @errors; $e++ ){
#				$summary[$e] = ();
#				for ($l = 0; $l < scalar @lengths; $l++) {
#						if ($numReads[$e][$l] > 0) {
#								$summary[$e][$l] = $numIncorrect[$e][$l] / $numReads[$e][$l];
#						}
#						else {
#								$summary[$e][$l] = "NaN";
#						}
#				}
#		}
#
#
#
#
#print "e l ratio\n";
#$i = 0;
#for ($e = 0; $e < scalar @errors; $e++ ){
#		for ($l = 0; $l < scalar @lengths; $l++) {
#				print "$i $errors[$e] $lengths[$l] $summary[$e][$l]\n";
#				$i++;
#		}
#}

#print "e l ratio\n";
#$i = 0;
#for ($e = 0; $e < scalar @errors; $e++ ){
#		for ($l = 0; $l < scalar @lengths; $l++) {
#				print "$i $errors[$e] $lengths[$l] $summaryMis[$e][$l]\n";
#				$i++;
#		}
#}
#
#print "e l ratio\n";
#$i = 0;
#for ($e = 0; $e < scalar @errors; $e++ ){
#		for ($l = 0; $l < scalar @lengths; $l++) {
#				print "$i $errors[$e] $lengths[$l] $summaryInc[$e][$l]\n";
#				$i++;
#		}
#}


#print " l: ";
#for ($l = 0; $l < scalar @lengths; $l++) {
#		printf("%8d ", @lengths[$l]);
#}
#print "\n";
#
#for ($e = 0; $e < scalar @errors; $e++ ){
#		print "@errors[$e]: ";
#		for ($l = 0; $l < scalar @lengths; $l++) {
#				if ($numReads[$e][$l] > 0) {
#						$ratio = $numNotTopOverlaps[$e][$l] / $numReads[$e][$l];
#						printf("%3.6f ", $ratio);
#				}
#				else {
#						print "NaN ";
#				}
#
#		}
#		print "\n";
#}
