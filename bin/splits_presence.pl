#!/usr/bin/perl -w
use strict;
use POSIX qw ( ceil );

use warnings;
use List::Util qw( min max sum );
use Devel::Cycle; 

use File::Basename 'dirname';
use Cwd 'abs_path';
my ( $bindir, $libdir );

BEGIN {	    # this has to go in Begin block so happens at compile time
   $bindir =
     dirname( abs_path(__FILE__) ) ; # the directory containing this script
   $libdir = $bindir . '/../lib';
   $libdir = abs_path($libdir);	# collapses the bin/../lib to just lib
}
# print "$libdir \n";

use lib $libdir;

use Getopt::Long;
use CXGN::Phylo::BasicTree;
use Hash::Ordered;
#use TomfyMisc qw ( increment_hash  add_hashes  format_newick_species_info median  store_gg_info );
#use Grouper;

# the input filse should have a newick expression on each line,
# true_tree_file has just one - the 'true' one,
# recon_tree_file has the other, reconstructed trees to be compared to it.

my $true_tree_file = undef;
my $recon_trees_file = undef;
my $scale = 1.0; # branch-length scale factor (multiplies true tree branch lengths)
my $max_trees = 100000;
my $nbins = 30;
my $lower_limit = 0.0;
my $upper_limit = 3.0;
GetOptions(
           'true_tree_file=s' => \$true_tree_file, # e.g. '*.newick'
           'trees_file|recon_trees_file=s' => \$recon_trees_file, # e.g. '*.newick'
           'scale=f' => \$scale, 
           'max_trees=i' => \$max_trees,
           'nbins=i' => \$nbins,
           'lower_limit=f' => \$lower_limit,
           'upper_limit=f' => \$upper_limit,
          );

my $true_tree = undef;

my %stsize_nsplits = ();
if (defined $true_tree_file) {
   open my $fhin_true, "<", $true_tree_file or die "Couldn't open $true_tree_file for reading.\n";;

   while (<$fhin_true>) {
      next if(/^(\s*#.*|\s*)$/); # skip comments, all-whitespace lines.
      my $newick = $_;
      $newick =~ s/\s+//g;
      $newick =~ s/:1[.]0([,\)])/:$scale$1/g;
      print STDERR "# true tree:  $newick \n";
      $true_tree = make_tree($newick);

      #   $copy_of_true_tree = $true_tree->
      my ( $distance, $sym_diff, $branch_score, $TL1, $TL2 ) = $true_tree->RF_distance_x($true_tree, \%stsize_nsplits);
      # %stsize_nsplits = map(($_ => scalar $stsize_nsplits{$_}), keys %stsize_nsplits); 
      $sym_diff /= 2;
      # $sym_diff is topological rf distance; $distance is one with branch-lengths
      #    $distance_histogram[$sym_diff]++;
   #   printf STDERR ( "%6.3f  %3d  %6.3f \n", $distance, $sym_diff, $TL2);
 printf STDERR ( "true tree vs self: %6.3f  %3d  %6.3f \n\n", $distance, $sym_diff, $TL2);
      # my @ssizes = sort { $a <=> $b } keys %stsize_nsplits;
      # for my $stsize (@ssizes) {
      #    print "[$stsize]", " ", $stsize_nsplits{$stsize}, "  ";
      # }
      # print "\n";
   }
   close $fhin_true;
}

open my $fhin_recons, "<", $recon_trees_file;

my $tree_count = 0;
my %stsize_count = ();
my @distance_histogram = (0,0,0,0,0);
while (<$fhin_recons>) {
   next if(/^(\s*#.*|\s*)$/);   # skip comments, all-whitespace lines.
   my $newick = $_;
   $newick =~ s/\s+//g;
   if (! defined $true_tree) {
      die "must specify true tree like: -true <true_tree_filename> \n";
      # $newick =~ s/:1[.]0([,\)])/:$scale$1/g;
      # print "# true tree:  $newick \n";
      # $true_tree = make_tree($newick);
      # #   $copy_of_true_tree = $true_tree->
      # my ( $distance, $sym_diff, $branch_score, $TL1, $TL2 ) = $true_tree->RF_distance_x($true_tree, \%stsize_nsplits);
      # $sym_diff /= 2;
      # # $sym_diff is topological rf distance; $distance is one with branch-lengths
      # $distance_histogram[$sym_diff]++;
      # printf STDERR ( "true tree vs self: %6.3f  %3d  %6.3f \n\n", $distance, $sym_diff, $TL2);

      # my @ssizes = sort { $a <=> $b } keys %stsize_nsplits;
      # for my $stsize (@ssizes) {
      #    my $count = scalar $stsize_nsplits{$stsize};
      #    print "$stsize $count    ";
      # }
      # print "\n";
   } else {
      my $recon_tree = make_tree($newick);
      my ( $distance, $sym_diff, $branch_score, $TL1, $TL2 ) = $true_tree->RF_distance_x($recon_tree, \%stsize_count);
      #   print "n sizes: ", scalar keys %stsize_count, "\n";
      $sym_diff /= 2;
      # $sym_diff is topological rf distance; $distance is one with branch-lengths
      $distance_histogram[$sym_diff]++;
      printf STDERR ( "%6.3f  %3d  %6.3f \n", $distance, $sym_diff, $TL2);
      $tree_count++;
      last if($tree_count >= $max_trees);
   }
}


@distance_histogram = map($_ // 0, @distance_histogram);
print "# $tree_count  ", join(",", @distance_histogram), "   ";
my @ssizes = sort { $a <=> $b } keys %stsize_count;
# print "N sizes: ", scalar keys @ssizes, "\n";
for my $stsize (@ssizes) {
   my $count = scalar @{$stsize_count{$stsize}};
   #   print "YYY $stsize ";
   # print $stsize_nsplits{$stsize}, "\n";
   print "$stsize $count ", $count/(scalar @{$stsize_nsplits{$stsize}})/$tree_count, "   ";
}
print "\n";


#my $nbins = 25;
#my ($lower_limit, $upper_limit) = (0.0, 2.5);
my $width = ($upper_limit - $lower_limit)/$nbins;

my %size_dlhist = ();
for my $size (@ssizes) {

   my $dls = $stsize_count{$size};
   #  print "size: $size  ls:  ", join(",", @$dls), "\n";
   $size_dlhist{$size} = histogram($dls, $nbins, $lower_limit, $upper_limit);
   #  print join(",", @{$size_dlhist{$size}}), "\n";
}
# exit;
print "# underflow: \n ";
print $lower_limit - $width, "  ", $lower_limit - 0.5*$width, "  ";
for my $size (@ssizes) {
   my $absents = scalar @{$stsize_nsplits{$size}}*$tree_count - 
     ( $size_dlhist{$size}->[0] // 0 + sum(@{$size_dlhist{$size}->[1]}) + $size_dlhist{$size}->[2] // 0);
   print "$size  ", $absents + $size_dlhist{$size}->[0] // 0, "    ";
}
print "\n";
for my $ibin (0..$nbins-1) {
   my $bin_lower_edge = $lower_limit + $ibin*$width;
   print $bin_lower_edge, "  ", $bin_lower_edge + 0.5*$width, "  ";
   for my $size (@ssizes) {
      
      print "$size  ", $size_dlhist{$size}->[1]->[$ibin] // 0, "    ";
   }
   print "\n";
}
print "#overflow: \n ";
print $upper_limit, "  ", $upper_limit + 0.5*$width, "  ";
for my $size (@ssizes) {
   print "$size  ", $size_dlhist{$size}->[2] // 0, "    ";
}
print "\n";

# ***********************************************8
sub histogram{
   my $xs = shift;              # array ref
   my $nbins = shift;
   my $lower_limit = shift;
   my $upper_limit = shift;
   my $bw = ($upper_limit - $lower_limit)/$nbins;
   my @h = ();
   my $under = 0;
   my $over = 0;
   for (@$xs) {
      if ($_ > $upper_limit) {
         $over++;
      } elsif ($_ <= $lower_limit) {
         $under++;
      } else {
         my $bin = ceil( ($_ - $lower_limit)/$bw ) - 1;
         #  print STDERR "ABC $_  $lower_limit  $bw $bin \n";
         $h[ $bin ]++;
      }
   }
   return [$under, \@h, $over];
}


sub make_tree{
   my $the_input_newick = shift;
   my $parser = CXGN::Phylo::Parse_newick->new( $the_input_newick, 0 );
   my $tree = $parser->parse( CXGN::Phylo::BasicTree->new() );
   $tree->make_binary();
   $tree->impose_branch_length_minimum(); # using default min branch length
   $tree->show_newick_attribute("species");
   $tree->set_show_standard_species(0);
   $tree->get_root()->recursive_implicit_names();
   $tree->get_root()->recursive_implicit_species();
   return $tree;
}
