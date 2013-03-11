#!/usr/bin/env perl

use MappingStudy;
use RunConfiguration;
#
# Configure the options list for the 
#

require ParameterNames;

if ($#ARGV < 0) {
		print "usag: SetupMapabilityStudy.pl ConfigurationFile.txt\n";
		print "nprocparam: $nProcParam\n";
		exit(1);
}
$configFileName = shift @ARGV;
%config = {};

# 
# Set up some defaults.
#

$config{$rawReadsDirParam} = "raw";
$config{$errorPrefixParam} = "Error";
$config{$lengthPrefixParam} = "Length";
$config{$readFileSuffixParam} = ".fasta";
$config{$errorSimParam} = "evolve";
$config{$nProcParam}     = 1;

&RunConfiguration::ReadConfigFile($configFileName, \%config);

$config{$nProcParam}     = 1;

#
# Check for required values.
#
&RunConfiguration::CheckRequiredParameter(\%config, $errorRatesParam);
&RunConfiguration::CheckRequiredParameter(\%config, $lengthsParam);
&RunConfiguration::CheckRequiredParameter(\%config, $referenceParam);
#&RunConfiguration::CheckRequiredParameter(\%config, $coverageParam);

# Build the directory structure.

if (exists $config{$workingDirParam}) {
		chdir($config{$workingDirParam});
}

@errors = split(/\s+/, $config{$errorRatesParam});
@lengths = split(/\s+/, $config{$lengthsParam});


#
# Configure custom error rates if specified.
#
$useDefaultErrorProfile = 1;
@errorBreakdown = ();
if (exists $config{$errorBreakdownParam}) {
		$useDefaultErrorProfile = 0;
		@errorBreakdown = split(/\s+/, $config{$errorBreakdownParam});
}

# Make the directory that will hold the reads without errors.
`mkdir -p $config{$rawReadsDirParam}`;
$nLengths = scalar @lengths;
print "creating $nLengths lengths\n";
for ($l = 0; $l < scalar @lengths; $l++) {
		# compute the coverage
		if (exists $config{$coverageParam}) {
				$samplingStr = " -coverage $config{$coverageParam} ";
		}
		elsif (exists $config{$nReadsParam}) {
				$samplingStr = " -nReads $config{$nReadsParam} ";
		}

		$shredCmd = "$config{$simpleShredderParam} -inFile $config{$referenceParam} -readLength $lengths[$l] $samplingStr -readsFile $config{$rawReadsDirParam}/$config{$lengthPrefixParam}". "_" . "$lengths[$l]$config{$readFileSuffixParam}";
		print "$shredCmd\n";
		`$shredCmd`;
}

@errorSimCommands = ();
$numESC = 0;
for ($e = 0; $e < scalar @errors; $e++ ){ 
		for ($l = 0; $l < scalar @lengths; $l++) {
				
				# Create the output directory.
				$dirName = $config{$errorPrefixParam} . "_" . "$errors[$e]" . "_" . $config{$lengthPrefixParam} . "_" . "$lengths[$l]";
				`mkdir $dirName`;
				# simulate errors.
				$rawReadsFile = "$config{$rawReadsDirParam}/$config{$lengthPrefixParam}" . "_" . "$lengths[$l]$config{$readFileSuffixParam}";
				if ($errors[$e] != 0) {
						$errorPct = $errors[$e]/100;
						$insRate = $errorPct*0.621;
						$delRate = $errorPct*0.276;
						$mutRate = $errorPct*0.103;
						$errorSimCmd = $config{$errorSimParam} . " -refGenome $rawReadsFile -mutGenome $dirName/reads$config{$readFileSuffixParam} -i $insRate -d $delRate -m $mutRate";
#						if ($useDefaultErrorProfile == 1) {
#						}
#						else {
#								$branchRate =  ($errors[$e] / 100) * @errorBreakdown[0];
#								$darkRate  =   ($errors[$e] / 100) * @errorBreakdown[1];
#								$miscallRate = ($errors[$e] / 100) * @errorBreakdown[2];
#								$stickRate =   ($errors[$e] / 100) * @errorBreakdown[3];
#								$errorSimCmd = $config{$errorSimParam} . " $rawReadsFile -s $stickRate -m $miscallRate -b $branchRate -d $darkRate > $dirName/reads$config{$readFileSuffixParam}";
#						}
				}
				else {
						$errorSimCmd = "cp $rawReadsFile $dirName/reads$config{$readFileSuffixParam}";
				}
				@errorSimCommands[$numESC] = $errorSimCmd;
				++$numESC;
				print "esc: $errorSimCmd\n";
#				`$errorSimCmd`;
		}
}

$esc  = 0;
$esci = 0;

while ($esci < $numESC) {
		@commands = ();
		for ($e = 0; $e < $config{$nProcParam}; $e++ ){
				@commands[$e] = @errorSimCommands[$esci];
				++$esci;
		}
		&MappingStudy::LaunchJobs(\@commands, $config{$nProcParam});
}
