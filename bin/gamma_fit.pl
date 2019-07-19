#!/usr/bin/perl -w
use strict;
use lib '/usr/local/share/perl/5.14.2';
use List::Util qw(min max sum);
use Graphics::GnuplotIF qw(GnuplotIF);
use Math::GSL::SF  qw( :all );
use Math::GSL::RNG  qw( :all );
use Math::GSL::Randist  qw( :all );
use Math::GSL::CDF  qw( :all );
use Getopt::Long;

my $fcvcols_file; # input file
my $alpha = 0.4; # initial guess
my $beta = 0.01; # initial guess
my $ytarg_str = '0.5,0.6,0.7,0.8,0.9';
my $mcmc_steps = 1000;
my $burn_in_factor = 0.1; # do a burn-in of this fraction of $mcmc_steps.

GetOptions(
           'input_file=s' => \$fcvcols_file, # 
           'alpha_init=f' => \$alpha,
           'beta_init=f' => \$beta,
           'ytargets=s' => \$ytarg_str, # e.g. '0.8,0.9'
           'mcmc_steps=i' => \$mcmc_steps,
);



my @ytargs = split(/[, ]+/, $ytarg_str);

my $the_rng = Math::GSL::RNG->new();

my %ncols_nc = ();
my %ncols_ntrials = ();
open my $fhin, "<", $fcvcols_file;
my @lines = <$fhin>;
for my $l (@lines) {
   if ($l =~ /bl: ([0-9.]+)\s+cols: ([0-9]+)\s+subst:\s+\S+\s+TvR:\s+([0-9.]+)\s.*\s(\S+)\s*$/) {
      my ($bl, $ncols, $fc, $ntrials) = ($1, $2, $3, $4);
    #  my @rfds = split(",", $rfdist);
    #  my $ntrials = ($fc > 0)? int($rfds[0]/$fc + 0.5) : 1000;
    #    print "XXX: $rfdist  ", $rfds[0], "   ncols $ncols ntrials $ntrials fc: $fc \n";
      $ncols_ntrials{$ncols} = $ntrials;
      my $nc =  int($fc*$ntrials + 0.5);
   #   print "# bl: $bl  ncols: $ncols  fraction correct: $fc \n";
      $ncols_nc{$ncols} = $nc;
   }
}
close $fhin;
my @sncols = sort {$a <=> $b} keys %ncols_nc;

#exit;


# my $Likelihood = 1;
# for my $ncols (@sncols) {
#    my $nc = $ncols_nc{$ncols};
#    print "$ncols  ", $ncols_nc{$ncols}, "\n";
#    my $pr_gamma = gsl_sf_gamma_inc_P($alpha, $beta*$ncols);
#    my $Lnc = gsl_ran_binomial_pdf($nc, $pr_gamma, $ntrials);

#    $Likelihood *= $Lnc;
#    print "$alpha $beta $ncols  $pr_gamma  $ntrials  $nc  ", $nc/$ntrials, " $Lnc  $Likelihood \n";
# }
# burn-in
my ($alpha_pbi, $beta_pbi, $logLmax_pbi, $x, $xx, $xxx, $xxxx) = mcmc_fit(\%ncols_nc, $alpha, $beta, \%ncols_ntrials, int($burn_in_factor*$mcmc_steps + 0.5), \@ytargs);

my ($alpha_opt, $beta_opt, $logLmax, $pa_alpha, $pa_beta, $pa_aprb, $states_str) = mcmc_fit(\%ncols_nc, $alpha_pbi, $beta_pbi, \%ncols_ntrials, $mcmc_steps, \@ytargs);
print $states_str, "\n";
print "# $alpha_opt $beta_opt  $logLmax   $pa_alpha $pa_beta $pa_aprb \n";

for my $ytarg (@ytargs) {
   my ($X, $deltay) = inc_gamma_inv($ytarg, $alpha_opt, $beta_opt);
   print "# $ytarg  $X  $deltay \n";
}


my $npts = scalar @sncols;
for (my $i=0; $i<$npts; $i++) {
   #my $ncols = $sncols[$i];
   my @ncols_sel = @sncols;
   my $ncols_i = splice(@ncols_sel,$i,1); # remove ith element.
   my %ncls_nc = map(($_ => $ncols_nc{$_}), @ncols_sel); 
   my ($alpha_o, $beta_o, $logLmax, $pa_alpha, $pa_beta, $pa_aprb) = mcmc_fit(\%ncls_nc, $alpha_opt, $beta_opt, \%ncols_ntrials, $mcmc_steps, undef);
   my $nc_i = $ncols_nc{$ncols_i};
   my $ntrials = $ncols_ntrials{$ncols_i};
   my $Pr =   gsl_sf_gamma_inc_P($alpha_o, $beta_o*$ncols_i);
   my $cdf_lt = ($nc_i > 0)? gsl_cdf_binomial_P($nc_i-1, $Pr, $ntrials) : 0;
   my $cdf_gt = gsl_cdf_binomial_P($nc_i, $Pr, $ntrials);
   my $pdf = gsl_ran_binomial_pdf($nc_i, $Pr, $ntrials);
   
   printf("%2i %5i  %5i %5i  %6.3f  %7.5f  %8.5g  %7.5f  %7.5f  %7.5f   %7.5f \n", 
          $i, $ncols_i, $nc_i, $ntrials, $alpha_o, $beta_o, $logLmax, $cdf_lt, $cdf_gt, 0.5*($cdf_lt + $cdf_gt), $Pr);
}
# exit;

my @col_data;

my $persist = 0;


my $plot1 = Graphics::GnuplotIF->new( persist => $persist, style => 'lines lw 2');
$plot1->gnuplot_cmd(" e(n,p) = sqrt(p*(1.0-p)/n) ");


while (my ($il, $line) = each @lines) {
   my @cols = split(" ", $line);
   while (my ($ic, $value) = each @cols) {
      push @{$col_data[$ic]}, $value;

   }
}

$alpha = $alpha_opt;
$beta = $beta_opt;
my $n = 160;
my $max_x = max( @{$col_data[5]} );
my $dx = $max_x/$n;
my $xi = 0;
open my $fhtmp, ">", "tmpfile";
print STDERR "# alpha, beta of curve to show with data: $alpha $beta \n";
for (0..$n) {
   print $fhtmp "$xi  ", gsl_sf_gamma_inc_P($alpha, $beta*$xi), "\n";
   $xi += $dx;
}
close $fhtmp;

$plot1->gnuplot_cmd('set grid');
$plot1->gnuplot_cmd('set key bottom');
my $gnuplot_cmd = " plot [0:*][0:1.02] \'$fcvcols_file\' using 6:10:(e(1000,\$10)) ps 1 with errorbars ";
$gnuplot_cmd .= ",  \'tmpfile\' using 1:2 with line";
$plot1->gnuplot_cmd($gnuplot_cmd);
#$plot1->gnuplot_plot_xy($col_data[5], $col_data[9]);
$plot1->gnuplot_pause(0);


sub logL{
   my $x_y = shift;
   my $alpha = shift;
   my $beta = shift;
   my $x_ntrials = shift;
   my $logLikelihood = 0;
   #for my $ncols (@sncols) {
   while (my($ncols, $nc) = each %$x_y) {
      #my $nc = $ncols_nc{$ncols};
      #  print "$ncols  $nc \n";
      my $ntrials = $x_ntrials->{$ncols};
      my $pr_gamma = gsl_sf_gamma_inc_P($alpha, $beta*$ncols);
      # print STDERR "# Z: $alpha $beta $ncols  $pr_gamma \n";
      my $Lnc = gsl_ran_binomial_pdf($nc, $pr_gamma, $ntrials);
      # print STDERR "#  L: $Lnc \n";
      return 'neg_inf' if($Lnc == 0);
      $logLikelihood += log($Lnc);
      #    print "$alpha $beta $ncols  $pr_gamma  $ntrials  $nc  ", $nc/$ntrials, " $Lnc  $logLikelihood \n";
   }
   return $logLikelihood;
}


sub mcmc_step_alpha{
   my $x_y = shift;
   my $alpha = shift;
   my $beta = shift;
   my $dalpha = shift;
   my $x_ntrials = shift;
   my $logL = shift;

   my $rn = gsl_rng_uniform($the_rng->raw());
   my $alpha_prop = $alpha + $dalpha*($rn - 0.5);
   # print STDERR "alpha:  $rn  $dalpha $alpha $alpha_prop   \n ";
   my $logLprop = logL($x_y, $alpha_prop, $beta, $x_ntrials);
   return (0, $alpha, $beta, $logL) if($logLprop eq 'neg_inf');
   #  print STDERR "$logL  $logLprop \n";
   if ($logLprop >= $logL  or  log( gsl_rng_uniform( $the_rng->raw())) < $logLprop - $logL) {
      return (1, $alpha_prop, $beta, $logLprop);
   } else {
      return (0, $alpha, $beta, $logL);
   }
}

sub mcmc_step_beta{
   my $x_y = shift;
   my $alpha = shift;
   my $beta = shift;
   my $dbeta = shift;
   my $x_ntrials = shift;
   my $logL = shift;

   my $rn = gsl_rng_uniform($the_rng->raw());
   my $beta_prop = $beta + $dbeta*($rn - 0.5);
   # print STDERR "beta:  $rn  $beta $beta_prop    \n";
   my $logLprop = logL($x_y, $alpha, $beta_prop, $x_ntrials);
   return (0, $alpha, $beta, $logL) if($logLprop eq 'neg_inf');
   # print STDERR "b lL lLp:  [$logL]  [$logLprop] \n";
   if ($logLprop >= $logL  or  log( gsl_rng_uniform( $the_rng->raw())) < $logLprop - $logL) {
      return (1, $alpha, $beta_prop, $logLprop);
   } else {
      return (0, $alpha, $beta, $logL);
   }
}



sub mcmc_step_aprb{
   my $x_y = shift;
   my $alpha = shift;
   my $beta = shift;
   my $daprb = shift;
   my $boa = shift;             # fixed approx beta/alpha
   my $x_ntrials = shift;
   my $logL = shift;

   my $rn = gsl_rng_uniform($the_rng->raw());
   my $alpha_prop = $alpha + $daprb*($rn - 0.5);
   my $beta_prop = $beta + $boa*$daprb*($rn - 0.5);
   # print STDERR "alpha:  $rn  $dalpha $alpha $alpha_prop   \n ";
   my $logLprop = logL($x_y, $alpha_prop, $beta_prop, $x_ntrials);
   return (0, $alpha, $beta, $logL) if($logLprop eq 'neg_inf'); # reject
   #  print STDERR "$logL  $logLprop \n";
   if ($logLprop >= $logL  or  log( gsl_rng_uniform( $the_rng->raw())) < $logLprop - $logL) {
      return (1, $alpha_prop, $beta_prop, $logLprop); # accept
   } else {
      return (0, $alpha, $beta, $logL); # reject
   }
}

sub simulate_data{
   # simulate numbers with correct topo. at given xs (numbers of align cols)
   my $xs = shift;
   my $alpha = shift;
   my $beta = shift;
   my $x_ntrials = shift;
   my %x_y = ();
   for my $ncols (@$xs) {
      my $pr_gamma = gsl_sf_gamma_inc_P($alpha, $beta*$ncols);
      my $nsuccess = gsl_ran_binomial($the_rng->raw(), $pr_gamma, $x_ntrials->{$ncols});
      $x_y{$ncols} = $nsuccess;
      print STDERR "# $ncols  $nsuccess \n";
   }
   return \%x_y;
}


sub mcmc_fit{
   my $ncols_nc = shift;
   my $alpha = shift;
   my $beta = shift;
   my $x_ntrials = shift;
   my $mcmc_steps = shift;
   my $ytargs = shift;

   my $boa = $beta/$alpha;
   my $logL = logL($ncols_nc, $alpha, $beta, $x_ntrials);
   #  print "# init alpha, beta, logL: $alpha $beta $logL   ";
   my $f = 0.15;
   my ($dalpha, $dbeta, $daprb) = ($alpha*$f, $beta*$f, $alpha*0.6);
   #  print  " dalpha, dbeta, daprb:  $dalpha $dbeta $daprb\n";

   my ($logLmax, $alpha_opt, $beta_opt) = ($logL, $alpha, $beta);
   my ($accept_alpha, $accept_beta, $accept_aprb) = (0, 0, 0);
   my ($nacc_alpha, $nacc_beta, $nacc_aprb) = (0, 0, 0);
   my $states_str = '';
   for (1..$mcmc_steps) {

      if ($_ % 2 == 0) {
         ($accept_beta, $alpha, $beta, $logL) = mcmc_step_beta($ncols_nc, $alpha, $beta, $dbeta, $x_ntrials, $logL);
         $nacc_beta += $accept_beta;
      } else {
         ($accept_alpha, $alpha, $beta, $logL) = mcmc_step_alpha($ncols_nc, $alpha, $beta, $dalpha, $x_ntrials, $logL);
         $nacc_alpha += $accept_alpha;
      }
      ($accept_aprb, $alpha, $beta, $logL) = mcmc_step_aprb($ncols_nc, $alpha, $beta, $daprb, $boa, $x_ntrials, $logL);
      $nacc_aprb += $accept_aprb;

      $states_str .=   "$accept_alpha  $accept_beta  $alpha  $beta  $logL   ";
      if (defined $ytargs) {
         for my $ytarg (@$ytargs) {
            my ($X, $dy) =  inc_gamma_inv($ytarg, $alpha, $beta);
            $states_str .= "$X $ytarg  ";
         }
      }
      $states_str .= "\n";
      if ($logL > $logLmax) {
         $logLmax = $logL;
         $alpha_opt = $alpha;
         $beta_opt = $beta;
      }
   }
   my ($pa_alpha, $pa_beta, $pa_aprb) = (2*$nacc_alpha/$mcmc_steps, 2*$nacc_beta/$mcmc_steps, $nacc_aprb/$mcmc_steps);
   # print "#  $pa_alpha  $pa_beta  $pa_aprb   $alpha_opt  $beta_opt   $logLmax \n";
   return ($alpha_opt, $beta_opt, $logLmax, $pa_alpha, $pa_beta, $pa_aprb, $states_str);
}


sub inc_gamma_inv{
   my $ytarget = shift;
   my $alpha = shift;
   my $beta = shift;

   my $dy = 1e-6;

   my $xlb = 0;
   my $xub = 1000;

   my $ylb = 0;
   my $y = gsl_sf_gamma_inc_P($alpha, $beta*$xub);
   while ($y < $ytarget) { # increase $xub as necessary so it is really upper bound.
      $xub *= 2;
      $y = gsl_sf_gamma_inc_P($alpha, $beta*$xub);
   }
   my $yub = $y;

   while (1) {
      my $x = 0.5*($xlb + $xub);
      $y = gsl_sf_gamma_inc_P($alpha, $beta*$x);
      if ($y > $ytarget + $dy) { # $x is new upper bound
         $xub = $x;
      } elsif ($y < $ytarget - $dy) { # $x is new lower bound
         $xlb = $x;
      } else {                  # close enough
         return ($x, $y - $ytarget);
      }
   }
}
