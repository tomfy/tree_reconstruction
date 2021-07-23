#!/usr/bin/perl -w
use strict;
use lib '/usr/local/share/perl/5.14.2';
use List::Util qw(min max sum);
use Graphics::GnuplotIF qw(GnuplotIF);
use Math::GSL::SF  qw( :all );


my $fcvcols_file = shift;
my $alpha = shift;
my $beta = shift;

my $persist = 0;


my $plot1 = Graphics::GnuplotIF->new( persist => $persist, style => 'lines lw 2');
$plot1->gnuplot_cmd(" e(n,p) = sqrt(p*(1.0-p)/n) ");

$plot1->gnuplot_cmd("plot [0:10] sin(x) with lines");
$plot1->gnuplot_pause(0);

exit;

open my $fhin, "<", $fcvcols_file;
my @lines = <$fhin>;
my $n_cols = scalar split(" ", $lines[0]);
my @col_data = ();
for ( 1 .. $n_cols ) {
   push @col_data, []; # each element of @col_data is an array ref which will
   # hold the data for a column
}

while (my ($il, $line) = each @lines) {
   my @cols = split(" ", $line);
   while (my ($ic, $value) = each @cols) {
      push @{$col_data[$ic]}, $value;

   }
}

my $n = 160;
my $max_x = max( @{$col_data[5]} );
my $dx = $max_x/$n;
my $x = 0;
open my $fhtmp, ">", "tmpfile";
print STDERR "$alpha $beta \n";
for(0..$n){
   print $fhtmp "$x  ", gsl_sf_gamma_inc_P($alpha, $beta*$x), "\n";
   $x += $dx;
}close $fhtmp;

$plot1->gnuplot_cmd('set key bottom');
my $gnuplot_cmd = " plot [0:*][0:1.02] \'$fcvcols_file\' using 6:10:(e(1000,\$10)) ps 1 with errorbars ";
$gnuplot_cmd .= ",  \'tmpfile\' using 1:2 with line";
$plot1->gnuplot_cmd($gnuplot_cmd);
#$plot1->gnuplot_plot_xy($col_data[5], $col_data[9]);
$plot1->gnuplot_pause(0);
