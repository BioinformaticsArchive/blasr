#!/usr/bin/env perl
use POSIX ":sys_wait_h";
use MappingStudy;
use RunConfiguration;




#
# move all this into the MappingStudy when time permits.
#
require ParameterNames;

#this is different.

$configFileName = shift @ARGV;
%config = {};

# 
# Set up some defaults.
#

$config{$rawReadsDirParam}    = "raw";
$config{$errorPrefixParam}    = "Error";
$config{$lengthPrefixParam}   = "Length";
$config{$readFileSuffixParam} = ".fasta";
$config{$errorSimParam}       = "errorSim";
$config{$maxExpandParam}      = "";
$config{$minMatchParam}       = "8";
$config{$nProcParam}          = 1;
$config{$bestNParam}          = 5;
$config{$alignSuffixParam}    = "align";


&RunConfiguration::ReadConfigFile($configFileName, \%config);

#
# Check for required values.
#
&RunConfiguration::CheckRequiredParameter(\%config, $referenceParam);
&RunConfiguration::CheckRequiredParameter(\%config, $errorRatesParam);
&RunConfiguration::CheckRequiredParameter(\%config, $lengthsParam);
&RunConfiguration::CheckRequiredParameter(\%config, $suffixArrayParam);
&RunConfiguration::CheckRequiredParameter(\%config, $mapperParam);

if (exists $config{$workingDirParam}) {
		chdir($config{$workingDirParam});
}
else {
		$config{$workingDirParam} = $ENV{"PWD"};
}

@errors  = split(/\s+/, $config{$errorRatesParam});
@lengths = split(/\s+/, $config{$lengthsParam});

@pids = ();

$gridCmdDir = "/mnt/secondary/Share/mchaisson/grid_commands/";
if (exists $config{$gridCommandDirParam}) {
    $gridCmdDir = $config{$gridCommandDirParam};
}

$cmdBaseName = "$gridCmdDir/cmd.$$.";
for ($e = 0; $e < scalar @errors; $e++ ){ 
		for ($l = 0; $l < scalar @lengths; $l++) {
				# Create the output directory.
				$dirName = $config{$errorPrefixParam} . "_" . "$errors[$e]" . "_" . $config{$lengthPrefixParam} . "_" . "$lengths[$l]";
 				# map reads
				$readsFile = "$dirName/reads$config{$readFileSuffixParam}";
        $readsBase = "reads$config{$readFileSuffixParam}";
				# launch the multiple jobs
				@pids = ();
				@commands = ();
				
				$ctabFileParamString = "";
				if (exists $config{$ctabParam}) {
						$ctabFileParamString = " -ctab $config{$ctabParam} ";
				}
        $cmdId = "$e.$l";
#				$command = "cd $config{$workingDirParam}; $config{$mapperParam} $readsFile $config{$referenceParam} -sa $config{$suffixArrayParam} $config{$maxExpandParam} -minMatch $config{$minMatchParam}  -bestn $config{$bestNParam} -nproc $config{$nProcParam} -out /scratch/$cmdId.$readsBase.$config{$alignSuffixParam} -nCandidates 10 -advanceExactMatches 5 -m 4 $ctabFileParamString -maxScore -200 $config{$extraOptsParam}; cat /scratch/$cmdId.$readsBase.$config{$alignSuffixParam}.* > $readsFile.$config{$alignSuffixParam} ";
				$command = "cd $config{$workingDirParam}; $config{$mapperParam} $readsFile $config{$referenceParam} -sa $config{$suffixArrayParam} $config{$maxExpandParam} -minMatch $config{$minMatchParam}  -bestn $config{$bestNParam} -nproc $config{$nProcParam} -out $readsFile.$config{$alignSuffixParam} -nCandidates 10 -advanceExactMatches 5 -m 4 $ctabFileParamString -maxScore -200 $config{$extraOptsParam}; ";
				if ($config{$runDistributedParam} eq "true") {
						$cmdName = $cmdBaseName . "$e.$l.sh";
						open(CMD, ">$cmdName") or die "cannot open cmd $cmdName\n";
						print CMD "$command\n";
						close(CMD);
						$startTime= time();
						print "submitting $command\n";
						system("qsub -V -o /tmp -pe smp $config{$nProcParam} $cmdName");
				}
				else {
						print "running $command\n";
						system($command);
				}


				if ($? != 0) {
						print "Error submitting job $command\n";
						exit(1);
				}
				$endTime = time();

				$outputFiles = "";

				$elapsedTime = $endTime - $startTime;
				system("echo $elapsedTime > $dirName/elapsed_time.txt");
		}
}
				
