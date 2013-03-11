#!/usr/bin/env perl
$inFile = shift @ARGV;

open(IN, "$inFile") or die "cannot open $inFile\n";
$cur = "";
while(<IN>) {
  $next = $_; 
	if ($cur eq "") {
	 	print $next;
  }
  else {
 		$cur =~ /(\S+) .*/;
	  $curTitle = $1;
 	 	$next =~ /(\S+) .*/;
		$nextTitle = $1;
  	if ($curTitle ne $nextTitle) {
			print "$next";
  	}
  } 
  $cur = $next;
}
