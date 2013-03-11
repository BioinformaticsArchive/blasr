#!/usr/bin/env perl 

$inName = shift @ARGV;
open(RTAB, $inName) or die "cannot open $inName\n";
# get rid of header
<RTAB>;
%grid = ();
%xkeys = ();
%ykeys = ();
while(<RTAB>) {
  $_ =~ /(\d+)\s+(\d+)\s+(\d+)\s+(\S+)/;
  $index = $1;
  $x = $2;
	$y = $3;
  $score=$4;
	$grid{$x}{$y} = $score;
  $xkeys{$x} = 1;
  $ykeys{$y} = 1;
}

@xkeysa = sort{$a <=> $b} keys %xkeys;
@ykeysa = sort{$a <=> $b} keys %ykeys;
print "\"x/y\", ";
for ($yk = 0; $yk <= $#ykeysa; $yk++) {
  print "$ykeysa[$yk], ";
}
print "\n"; 
for ($xk = 0; $xk <= $#xkeysa; $xk++) {
  print "$xkeysa[$xk], ";
  for ($yk = 0; $yk <= $#ykeysa; $yk++) {
	 print "$grid{$xkeysa[$xk]}{$ykeysa[$yk]}, ";
  }
  print "\n";
}	
