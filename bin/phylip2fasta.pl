#!/usr/bin/perl -w
use strict;

my $first_line = <>;

while(<>){
my ($id, $seq) = split(" ", $_);
print ">$id\n", "$seq \n";
}
