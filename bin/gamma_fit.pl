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
my $p_ratio = 1.0; 
my $ytarg_str = '0.7,0.8,0.9';
my $mcmc_steps = 1000;
my $burn_in_factor = 0.1; # do a burn-in of this fraction of $mcmc_steps.
my $rng_seed = undef;

GetOptions(
           'input_file=s' => \$fcvcols_file, # 
           'alpha_init=f' => \$alpha,
           'beta_init=f' => \$beta,
           'ytargets=s' => \$ytarg_str, # e.g. '0.8,0.9'
           'mcmc_steps=i' => \$mcmc_steps,
           'seed=i' => \$rng_seed,
);



my @ytargs = split(/[, ]+/, $ytarg_str);

my $the_rng = Math::GSL::RNG->new($gsl_rng_default, $rng_seed);

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
my ($alpha_pbi, $beta_pbi, $p_ratio_pbi, $logLmax_pbi, $x, $xx, $xxx, $xxxx) = mcmc_fit(\%ncols_nc, $alpha, $beta, $p_ratio, \%ncols_ntrials, int($burn_in_factor*$mcmc_steps + 0.5), \@ytargs);

my ($alpha_opt, $beta_opt, $p_ratio_opt, $logLmax, $pa_alpha, $pa_beta, $pa_aprb, $states_str) = mcmc_fit(\%ncols_nc, $alpha_pbi, $beta_pbi, $p_ratio_pbi, \%ncols_ntrials, $mcmc_steps, \@ytargs);
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
   my ($alpha_o, $beta_o, $p_ratio_o, $logLmax, $pa_alpha, $pa_beta, $pa_aprb) = mcmc_fit(\%ncls_nc, $alpha_opt, $beta_opt, $p_ratio_opt, \%ncols_ntrials, $mcmc_steps, undef);
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


my $plot1 = Graphics::GnuplotIF->new( persist => $persist, style => 'lines lw 6');
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
$gnuplot_cmd .= ",  \'tmpfile\' using 1:2 with line lw 2. lc 'sea-green' ";
$plot1->gnuplot_cmd($gnuplot_cmd);
#$plot1->gnuplot_plot_xy($col_data[5], $col_data[9]);
$plot1->gnuplot_pause(0);


sub logL{
   my $x_y = shift;
   my $alpha = shift;
   my $beta = shift;
   my $x_ntrials = shift;
   my $p_ratio = shift // 1.0; # (alpha2-1)/(alpha1-1) # typically < 1
   # piecewise gamma-distributed with 2 pieces (or 1 if $alpha_ratio == 1)
   # i.e. for x <= ($alpha-1)/$beta  pdf = C*exp(-beta*x)*x^(alpha-1) / exp(-(alpha-1))*x^(alpha-1)
   # and for x >= ($alpha-1)/$beta  pdf = C*exp(-beta*p_ratio*x)*x^((alpha-1)*p_ratio) / exp(
   my $logLikelihood = 0;
   #for my $ncols (@sncols) {
#print STDERR "abr: $alpha $beta $p_ratio \n";
   if ($p_ratio == 1) {
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
   }else{
      my $mode = ($alpha-1)/$beta;
      my $pdf1_at_mode = gsl_ran_gamma_pdf($mode, $alpha, $beta);
      my $alpha2 = ($alpha - 1.0)*$p_ratio + 1.0;
      my $beta2 = $beta*$p_ratio;
      my $pdf2_at_mode = gsl_ran_gamma_pdf($mode, $alpha2, $beta2);
      my $relp_ltm = gsl_sf_gamma_inc_P($alpha, $beta*$mode);
      my $relp_gtm = gsl_sf_gamma_inc_Q($alpha2, $beta2*$mode) * $pdf1_at_mode/$pdf2_at_mode;
      print "$relp_ltm  $relp_gtm ", $relp_ltm + $relp_gtm, "\n";
   }
}


sub mcmc_step_alpha{
   my $x_y = shift;
   my $alpha = shift;
   my $beta = shift;
   my $p_ratio = shift;
   my $dalpha = shift;
   my $x_ntrials = shift;
   my $logL = shift;
  

   my $rn = gsl_rng_uniform($the_rng->raw());
   my $alpha_prop = $alpha + $dalpha*($rn - 0.5);
   # print STDERR "alpha:  $rn  $dalpha $alpha $alpha_prop   \n ";
   my $logLprop = logL($x_y, $alpha_prop, $beta, $x_ntrials, $p_ratio);
   return (0, $alpha, $beta, $logL) if($logLprop eq 'neg_inf');
   #  print STDERR "$logL  $logLprop \n";
   if ($logLprop >= $logL  or  log( gsl_rng_uniform( $the_rng->raw())) < $logLprop - $logL) {
      return (1, $alpha_prop, $beta, $p_ratio, $logLprop);
   } else {
      return (0, $alpha, $beta, $p_ratio, $logL);
   }
}

sub mcmc_step_beta{
   my $x_y = shift;
   my $alpha = shift;
   my $beta = shift;
   my $p_ratio = shift;
   my $dbeta = shift;
   my $x_ntrials = shift;
   my $logL = shift;

   my $rn = gsl_rng_uniform($the_rng->raw());
   my $beta_prop = $beta + $dbeta*($rn - 0.5);
   # print STDERR "beta:  $rn  $beta $beta_prop    \n";
   my $logLprop = logL($x_y, $alpha, $beta_prop, $x_ntrials, $p_ratio);
   return (0, $alpha, $beta, $logL) if($logLprop eq 'neg_inf');
   # print STDERR "b lL lLp:  [$logL]  [$logLprop] \n";
   if ($logLprop >= $logL  or  log( gsl_rng_uniform( $the_rng->raw())) < $logLprop - $logL) {
      return (1, $alpha, $beta_prop, $p_ratio, $logLprop);
   } else {
      return (0, $alpha, $beta, $p_ratio, $logL);
   }
}



sub mcmc_step_aprb{
   my $x_y = shift;
   my $alpha = shift;
   my $beta = shift;
my $p_ratio = shift;
   my $daprb = shift;
   my $boa = shift;             # fixed approx beta/alpha
   my $x_ntrials = shift;
   my $logL = shift;

   my $rn = gsl_rng_uniform($the_rng->raw());
   my $alpha_prop = $alpha + $daprb*($rn - 0.5);
   my $beta_prop = $beta + $boa*$daprb*($rn - 0.5);
   # print STDERR "alpha&beta:  $rn  $daprb   $alpha $alpha_prop    $beta $beta_prop  \n ";
   my $logLprop = logL($x_y, $alpha_prop, $beta_prop, $x_ntrials, $p_ratio);
   return (0, $alpha, $beta, $logL) if($logLprop eq 'neg_inf'); # reject
   #  print STDERR "$logL  $logLprop \n";
   if ($logLprop >= $logL  or  log( gsl_rng_uniform( $the_rng->raw())) < $logLprop - $logL) {
      return (1, $alpha_prop, $beta_prop, $p_ratio, $logLprop); # accept
   } else {
      return (0, $alpha, $beta, $p_ratio, $logL); # reject
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
   my $p_ratio = shift;
   my $x_ntrials = shift;
   my $mcmc_steps = shift;
   my $ytargs = shift;

   my $boa = $beta/$alpha;
   my $logL = logL($ncols_nc, $alpha, $beta, $x_ntrials, $p_ratio);
   #  print "# init alpha, beta, logL: $alpha $beta $logL   ";
   my $f = 0.14;
   my ($dalpha, $dbeta, $daprb) = ($alpha*$f, $beta*$f, $alpha*0.6);
   # print  STDERR " dalpha, dbeta, daprb:  $dalpha $dbeta $daprb\n";

   my ($logLmax, $alpha_opt, $beta_opt, $p_ratio_opt) = ($logL, $alpha, $beta, $p_ratio_opt);
   my ($accept_alpha, $accept_beta, $accept_aprb) = (0, 0, 0);
   my ($nacc_alpha, $nacc_beta, $nacc_aprb) = (0, 0, 0);
   my $states_str = '';
   for (1..$mcmc_steps) {

      if ($_ % 2 == 0) {
         ($accept_beta, $alpha, $beta, $p_ratio, $logL) = mcmc_step_beta($ncols_nc, $alpha, $beta, $p_ratio, $dbeta, $x_ntrials, $logL);
         $nacc_beta += $accept_beta;
      } else {
         ($accept_alpha, $alpha, $beta, $p_ratio, $logL) = mcmc_step_alpha($ncols_nc, $alpha, $beta, $p_ratio, $dalpha, $x_ntrials, $logL);
         $nacc_alpha += $accept_alpha;
      }
      ($accept_aprb, $alpha, $beta, $p_ratio, $logL) = mcmc_step_aprb($ncols_nc, $alpha, $beta, $p_ratio, $daprb, $boa, $x_ntrials, $logL);
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
   return ($alpha_opt, $beta_opt, $p_ratio_opt, $logLmax, $pa_alpha, $pa_beta, $pa_aprb, $states_str);
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
