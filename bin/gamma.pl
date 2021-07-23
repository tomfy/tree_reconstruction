#!/usr/bin/perl -w
use strict;
use Math::GSL::SF qw( :all );

my $alpha = shift;
my $beta = shift;
my $dx = shift // 0.1;

for(my $i=0; $i<1000; $i++){

my $x = $i*$dx;
my $y = gsl_sf_gamma_inc_P($alpha, $beta*$x);
print "$x  $y \n";

last if(1.0 - $y < 0.001);
}
