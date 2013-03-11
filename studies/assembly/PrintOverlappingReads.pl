#!/usr/bin/env perl
$usage = "usage: PrintOverlappingReads.pl reads [-minIdent i] [-minAlignLegth a] [-minReadLength l] [-adjtable] [-equiv equivalency table]\n";
if ($#ARGV < 0) {
	print $usage;
	exit(1);
}

$in = shift @ARGV;
$minIdent = 0;
$minAlignLength = 0;
$minReadLength = 0;
$wiggle = 10;
$adjTableName = "";
$equivTableName = "";
while ($#ARGV >= 0) {
		$opt = shift @ARGV;
		if ($opt eq "-minIdent") {
				$minIdent = shift @ARGV;
		}
		elsif ($opt eq "-minAlignLength") {
				$minAlignLength = shift @ARGV;
		}
		elsif ($opt eq "-minReadLength") {
				$minReadLength = shift @ARGV;
		}
		elsif ($opt eq "-adj") {
				$adjTableName = shift @ARGV;
		}
		elsif ($opt eq "-equiv") {
				$equivTableName = shift @ARGV;
		}
		else {
				print "ERROR, bad option $opt\n";
				print $usage;
				exit(1);
		}
}

open(IN, "$in") or die "cannot open $in\n";

%leftToRight = {};
%rightToLeft = {};
if ($adjTableName) {
		open(ADJ, $adjTableName) or die "cannot open $adjTableName\n";
		while(<ADJ>) {
				chomp;
				@vals = split(/\s+/, $_);
				$leftToRight{@vals[0]} = [@vals];
				$rightToLeft{@vals[1]} = [@vals];

		}
}
%equivalencies = {};
if ($equivTableName) {
	open(EQ, "$equivTableName") or die "cannot open $equivTableName\n";
  while(<EQ>) {
	  @eqvals = split(/\s+/,$_);
	  $edgeIndex = shift @eqvals;
	  $neqvals =scalar @eqvals;
		if ($neqvals > 0) {
			$equivalencies[$edgeIndex] = [@eqvals];
    }
  }
} 

$prevReadTitle = "";
$curReadTitle = "";
%overlaps = ();
#print "source dest sourceStrand destStrand sourceAlignScore destAlignScore\n";
while(<IN>) {
		@vals = split(/\s+/, $_);
		$curReadTitle = @vals[0];
		$readLength   = @vals[9];
#		print "read length: $readLength\n";
		$alignLength = abs(@vals[8] - @vals[7]);
		if (@vals[0] eq $prevReadTitle and $alignLength > $minAlignLength && $readLength > $minReadLength) {
				#
				# Check for an overlap
				#
				# @vals[11] = refStart @vals[12] = refEnd @vals[13] refLength
				$overlap = 0;
				if (@vals[10] == 0 && @vals[11] < $wiggle) {
						$overlap = 1;
				}
				elsif (@vals[10] == 0  && @vals[13] - @vals[12] < $wiggle) {
						$overlap = 1;
				}
				elsif (@vals[10] == 1 && @vals[13] - @vals[11] < $wiggle) {
						$overlap = 1;
				}
				elsif (@vals[10] == 1 && @vals[12] < $wiggle) {
						$overlap = 1;
				}
				if (@vals[5] > $minIdent && @vals[4] < $maxScore && $overlap == 1) {
						if (not exists $equivalencies[@vals[1]]) {
						  push @readOverlaps, $_;
						}
						else{
#							print "skipped equiv: @vals[1]\n";
						}
				}
		}
		else {
				# 
				# process the previous set of overlaps.
				# 
				$nOverlaps = scalar @readOverlaps;
				if ($nOverlaps > 1) {
						@ovpMatrix = ();
						$ovpPrinted = 0;
						for ($i = 0; $i < $nOverlaps ;$i++ ){
#0                                              1     2    3        4      5       6 7   8   9   10  11  12  13
#x25_y74_2100223-0011_m100109_122815_Jan_p2_b15 22778 -863 -108.785 -28594 78.9474 0 151 418 418 0   0   254 1230
								@ovp = split(/\s+/, @readOverlaps[$i]);
								# 
								# Fix strandedness of overlap.
								#
								if (@ovp[6] == 1) {
										# 
										# make read in forward orientation.
										#
										$readAlignStart = @ovp[9] - @ovp[8];
										$readAlignEnd   = @ovp[9] - @ovp[7];
										@ovp[7] = $readAlignStart; 
										@ovp[8] = $readAlignEnd;
										$targetAlignStart = @ovp[13] - @ovp[11];
										$targetAlignEnd   = @ovp[13] - @ovp[12];
										@ovp[12] = $targetAlignStart;
										@ovp[11] = $targetAlignEnd;
										if (@ovp[12] == -1) { @ovp[12] = 0;}
										if (@ovp[11] == -1) { @ovp[11] = 0;}

										@ovp[6] = 0;
										if (@ovp[10] == 0) {
												@ovp[10] = 1;
										}
										else {
												@ovp[10] = 0;
										}
								}
								push @ovpMatrix, [@ovp];
						}
						$queryMapsAsRepeat = 0;
						for ($i = 0; $i < $#ovpMatrix; $i++ ){
								for ($j = $i + 1; $j <= $#ovpMatrix; $j++ ){
										if ($ovpMatrix[$i][1] eq $ovpMatrix[$j][1]) {
												$queryMapsAsRepeat = 1;
										}
								}
						}

						$queryStartPositionsAreUnique = 1;
						@queryStartPositions = ();
						$nOvp = $#ovpMatrix;
						for ($i = 0; $i <= $#ovpMatrix; $i++ ){
								for ($j = $i + 1; $j <= $#ovpMatrix; $j++ ){
										if ($ovpMatrix[$i][7] == $ovpMatrix[$j][7]) {
												$queryStartPositionsAreUnique = 0;
										}
								}
								push @queryStartPositions, $ovpMatrix[$i][7];
						}
						@sortedQueryStartPositions = sort {$a <=> $b} @queryStartPositions;
						if ($queryMapsAsRepeat == 0 && $queryStartPositionsAreUnique == 1) {
								$nsp = scalar @sortedQueryStartPositions;
								$nqsp = scalar @queryStartPositions;
								for ($i =0 ; $i < $#sortedQueryStartPositions; $i++) {
										$src = -1;
										$dest = -1;
										$srcIndex = 0;
										$destIndex = 0;
										for ($j = 0; $j <= $#ovpMatrix; $j++ ){
												if ($ovpMatrix[$j][7] == $sortedQueryStartPositions[$i]) {
														$src = $ovpMatrix[$j][1];
														$srcIndex = $j;
														last;
												}
										}
										for ($j = 0; $j <= $#ovpMatrix; $j++ ){
												if ($ovpMatrix[$j][7] == $sortedQueryStartPositions[$i+1]){
														$dest = $ovpMatrix[$j][1];
														$destIndex = $j;
														last;
												}
										}
										if ($src != -1 && $dest != -1) {
												$leftToRight = 1;
												$validated = 0;
												if (exists $leftToRight{$src}) {
														if ($leftToRight{$src}[1] eq $dest) {
																$validated = 1;
														}
												}
												if (exists $rightToLeft{$src}) {
														if ($rightToLeft{$src}[0] eq $dest) {
																$validated = 1;
														}
												}

												if ($adjTableName eq "" || $validated == 1) {
														print "@{$ovpMatrix[$srcIndex]}\n";
														print "@{$ovpMatrix[$destIndex]}\n";
												}
												$ovpPrinted = 1;
										}
										else {
												print "ERROR! Something went bad with looking up coordinates.\n";
												exit(1);
										}
								}
								$curStartPos = -1;
						}
				}
				$prevReadTitle = $curReadTitle;
				@readOverlaps = ($_);
		}
}
