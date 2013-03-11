#!/usr/bin/env perl

use AlignmentParser;

$in = shift @ARGV;
open(IN, "$in") or die "cannot open $ARGV[0]\n";

$nTotStick =0;
$nTotBranch = 0;
$nTotMiscall = 0;
$nTotDark = 0;
$nTotBase = 0;

while(<IN>) {
		%alignment = {};

		AlignmentParser::ParseCompareSequencesLine($_, \%alignment);

		$qString = $alignment{"qString"};
		$tString = $alignment{"tString"};

		# First compute the easy parameter - miscall
		$nMiscall = 0;
		if (length($qString) != length($tString)) {
				print "ERROR, the alignment lines must be the same length\n";
				print "$qString\n";
				print "$qLine\n";
				exit(1);
		}
		$lstr = length($qString);

		for ($p = 0; $p < length($qString); $p++ ){
				$q = substr($qString, $p, 1);
				$t = substr($tString, $p, 1);
				if ($q ne "-" && $t ne "-" and $q ne $t) {
						$nMiscall++;
				}
		}

		$nBranch = 0;
		$nDark = 0;
		# Next compute the dark fraction, another easy one.
		for ($p = 0; $p < length($qString); $p++) {
				$q = substr($qString, $p, 1);
				if ($q eq "-") {
						$nDark++;
				}
		}

		# Next compute the branching and sticking.  Branching is when 
		# the bases adjacent to a gap are the same, sticking is when they are different.
		$nBase = 0;
		$nStick = 0;
		$nBranch = 0;
		for ($p = 0; $p < length($qString); $p++ ){
				$t = substr($tString, $p, 1);
				if ($t eq "-") {
						$p2 = $p + 1;

						# look for the end of this gap.
						while ($p2 < length($qString)) {
								$t2 = substr($tString, $p2, 1);
								if ($t2 ne "-") {
										last;
								}
								else {
										$p2++;
								}
						}
						$sameOnLeft = 0;
						$qInsertedChar = substr($qString, $p, 1);
						if ($p > 0) {
								$leftChar = substr($tString, $p-1,1);
								if ($leftChar eq $qInsertedChar) {
										$sameOnLeft = 1;
								}
						}
						$sameOnRight = 0;
						if ($p2 < length($tString)-1) {
								$rightChar = substr($tString, $p2+1, 1);
								if ($rightChar eq $qInsertedChar) {
										$sameOnRight = 1;
								}
						}
						if ($sameOnLeft == 1 || $sameOnRight == 1) {
								$nBranch+= $p2 - $p;
						}
						else {
								$nStick++;
						}
						$p = $p2 - 1;
				}
				else {
						$nBase++;
				}
		}
		$nTotStick += $nStick;
		$nTotBranch += $nBranch;
		$nTotDark += $nDark;
		$nTotMiscall += $nMiscall;
		$nTotBase += $nBase;
}
						
$stickRate = $nTotStick / $nTotBase;
$branchRate = $nTotBranch / $nTotBase;
$darkRate = $nTotDark / $nTotBase;
$miscallRate = $nTotMiscall / $nTotBase;
$tot = $stickRate + $branchRate + $darkRate + $miscallRate;
print "$stickRate $branchRate $darkRate $miscallRate $tot\n";
$normStickRate = $stickRate / $tot;
$normBranchRate = $branchRate / $tot;
$normDarkRate = $darkRate / $tot;
$normMiscallRate = $miscallRate / $tot;
print "$normBranchRate $normDarkRate $normMiscallRate $normStickRate\n";
