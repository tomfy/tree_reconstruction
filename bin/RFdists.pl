#!/usr/bin/perl -w
use strict;
use List::Util qw(min max sum);

my $true_tree_file = shift;
my $recon_trees_file = shift;
my $true_tree = `cat $true_tree_file`;
$true_tree =~ s/\s+//g;
my $comma_count = () = $true_tree =~ /,/g;
my $tree_size = $comma_count + 1;
$true_tree = "'" . $true_tree . "'";

my $true_and_recons_trees = 'true_and_reconstructed_trees';
my $raxml_file_ending = 'raxml_rfds';

system "cat $true_tree_file $recon_trees_file > $true_and_recons_trees";
my $rax_stdout = `raxmlHPC-PTHREADS-SSE3 -T 1 -m PROTCATWAGF -z $true_and_recons_trees -f r -n $raxml_file_ending `;
my $raxml_rfds_string = `cat RAxML_RF-Distances.$raxml_file_ending`;



my ($max_true_v_recon_rfd, $max_recon_v_recon_rfd) = (0, 0);
my @true_v_recon_rfd_histogram = (0) x $tree_size;
my @recon_v_recon_rfd_histogram = (0) x $tree_size;

my @rf_lines = split("\n", $raxml_rfds_string);
for my $aline (@rf_lines){
   next if($aline =~ /^\s*#/);
   my @fields = split(" ", $aline);
#   print STDERR "XXXXXX: ", join(", ", @fields), "\n"; 
#exit;
   my $rfd = $fields[2];
   $rfd = int($rfd/2 + 0.1);
   if ($fields[0] eq '0') {      # true v recon
      $max_true_v_recon_rfd = $rfd if($rfd > $max_true_v_recon_rfd);
      $true_v_recon_rfd_histogram[$rfd]++;
   } else {                     # recon v recon
      $max_recon_v_recon_rfd = $rfd if($rfd > $max_recon_v_recon_rfd);
      $recon_v_recon_rfd_histogram[$rfd]++;
   }
}

my $n_tvr = sum( @true_v_recon_rfd_histogram);
my $n_rvr = sum( @recon_v_recon_rfd_histogram);
# print STDERR "n tvr, n rvr: $n_tvr  $n_rvr \n";q
my ($tvr_mean, $tvr_med, $tvr_stddev, $tvr_stderr) = mean_median_and_stddev_of_index(\@true_v_recon_rfd_histogram);
my ($rvr_mean, $rvr_med, $rvr_stddev, $rvr_stderr) = mean_median_and_stddev_of_index(\@recon_v_recon_rfd_histogram);
# print  $true_v_recon_rfd_histogram[0],  "  $n_tvr  ",  $recon_v_recon_rfd_histogram[0], "  $n_rvr \n";
my $hist_str = join(",", @true_v_recon_rfd_histogram[0..5]);
printf ( #"%s  %s    "
         "TvR:  %6.4f  %6.4f  %4.2f  %6.4f   RvR:  %6.4f  %6.4f  %4.2f  %6.4f   %s\n", 
             # $true_tree_file, $recon_trees_file,
         $true_v_recon_rfd_histogram[0]/$n_tvr, $tvr_mean, $tvr_med, $tvr_stderr,
         ($n_rvr >= 1)? $recon_v_recon_rfd_histogram[0]/$n_rvr : '---', $rvr_mean, $rvr_med, $rvr_stderr,
         $hist_str);

my @files_to_rm = split(" ", `ls *.raxml_rfds`);
#print join(", ", @files_to_rm), "\n";
unlink(@files_to_rm);



sub mean_median_and_stddev_of_index{
   my $array_ref = shift;       # values are weights.
   my $q = shift // 0.5;        # i.e. median by default
   my $array_size = scalar @$array_ref;
   my $n = sum( @$array_ref );
   return (undef, undef, undef) if($n <= 0);
   my $wsum = 0;
   my $wsumsq = 0;
   my $wcume = 0;
   my $w_threshhold = $q*$n;
   my $quant = undef;
   while( my ($i, $w) = each @$array_ref) {
      # my $w = $array_ref->[$i];
      $wcume += $w;
      if (! defined $quant) {
         if ($wcume > $w_threshhold) {
            $quant = $i;
         } elsif ($wcume == $w_threshhold) {
            for (my $j=$i+1; $j < $array_size; $j++) {
               if ($array_ref->[$j] > 0) {
                  $quant = 0.5*($i + $j);
                  last;
               }
            }
         }
      }

      $wsum += $i*$w;
      $wsumsq += $i*$i*$w;

   }
   my $mean = $wsum / $n;
   my $var = $wsumsq/$n - $mean*$mean;
   my $stddev = sqrt($var);
	return ($mean, $quant, $stddev, $stddev/sqrt($n));
}

