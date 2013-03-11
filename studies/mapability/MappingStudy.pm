#!/usr/bin/env perl

package MappingStudy;


sub LaunchJob {
		my ($command) = @_;
		my $pid;
		defined($pid = fork()) or die "unable to fork $!\n";
		if ($pid == 0) {
				# running in the child process
				print "launching: $command\n";
				system($command);
				exit(0);
		}
		else {
				return $pid;
		}
}

sub LaunchJobs {
		my ($commands, $nProc) = @_;
		@pids = ();

		my $curCommand = 0;
		my $numCommands = scalar @{$commands};

		# fill the queue
		while ($curCommand < $numCommands) {
				$p = 0;
				while ($p < $nProc && $curCommand < $numCommands) {
						&LaunchJob(@{$commands}[$curCommand]);
						$p++;
						$curCommand++;
				}
		}
		# replace a running command with a new one
		while ($curCommand < $numCommands && wait() != -1) {
				&LaunchJob(@{$commands}[$curCommand]);
				$curCommand++;
		}
		# wait for the last few jobs to complete
		while (wait() != -1) {};
}

1;

