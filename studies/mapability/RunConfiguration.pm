#!/usr/bin/env perl

package RunConfiguration;

sub ReadConfigFile {
		my($configFileName, $config) = @_;
		open(CONF, "$configFileName") or die "cannot open $configFileName\n";
		while(<CONF>) {
				chomp;
				$opt = $_;
				if ($opt =~ /#/) {
						# there is a comment in this line, skip everything past it.
						$opt =~ s/([^#]*).*/$1/;
				}
				if ($opt =~ /(\S+)\s+(.*)/) {
						# otherwise, the option is a keyword followed by a space, then everything else left on the line.
						$key = $1;
						$value = $2;
						# trim whitespace from the value
						$value =~ s/^\s+|\s+$//g;
						# Check to make sure the key parameter is not ordered
						if ($key =~ /^@([^@]+)/) {
								# The key is in the format @key, which means that
								# it is a list-value key which the calling program knows 
								# how to deal with.
								$key = $1;
								push @{${$config}{$key}}, $value;
						}
						else {
								# standard keyword/value pair.
								${$config}{$key} = $value;
						}
				}
		}
}

sub CheckRequiredParameter {
		my ($config, $option) = @_;
		if (!exists ${$config}{$option}) {
				print "ERROR, config file is missing required parameter $option\n";
		}
}


return 1;
