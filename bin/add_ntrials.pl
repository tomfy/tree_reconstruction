#!/usr/bin/perl -w
use strict;

while (my $l = <>) {
   if ($l =~ /bl: ([0-9.]+)\s+cols: ([0-9]+)\s+subst:\s+\S+\s+TvR:\s+([0-9.]+)\s.*\s(\S+)\s*$/) {
      my ($bl, $ncols, $fc, $rfdist) = ($1, $2, $3, $4);
      my @rfds = split(",", $rfdist);
      my $ntrials = ($fc > 0)? int($rfds[0]/$fc + 0.5) : 1000;
      chomp $l;
      print "$l  $ntrials \n";
   } else {
      print "# line of unexpected format: $l";
   }
}
