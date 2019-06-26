#!/usr/bin/perl -w
use strict;
use Algorithm::CurveFit;
use Algorithm::CurveFit::Simple qw( fit );
use Data::Dumper;
use Math::GSL::RNG;
use Math::GSL::Randist qw ( :all );

my $RNG = Math::GSL::RNG->new();
my $input_pattern = shift;
my @input_filenames = split(" ", `ls $input_pattern`);
my $threshhold = shift // 0.9; # we will interpolate to find the alignment length s.t. that topo is correct this fraction of time.
my $Ndraw = shift // 100;
my $Nreps = shift // 1000;
my $verbose = shift // 0;

for my $input_file (@input_filenames){
my %ncols_fc = ();
my $branch_length = undef;
my $treename = undef;
my $model = undef;
open my $fhin, "<", "$input_file";
while (<$fhin>) {
   my @cols = split(" ", $_);
   my ($tn, $bl, $ncols, $mod, $fc) = @cols[1,3,5,7,9];
   die "tree name inconsistency: $tn $treename \n" if(defined $treename  and  ($tn ne $treename));
   die "branch length inconsistency: $bl $branch_length \n" if(defined $branch_length  and  ($bl != $branch_length));
   die "model inconsistency: $mod $model \n" if(defined $model  and  ($mod ne $model));
   $treename = $tn;
   $branch_length = $bl;
   $model = $mod;
   my $tree_size = ($treename =~ /t(\d+)u/)? $1 : die "input line has unexpected treename: $treename \n";
   my $stderr = binomial_stderr($Nreps, $fc);
   $ncols_fc{$ncols} = "$fc $stderr";
}
close $fhin;

my @ncols = sort {$a <=> $b} keys %ncols_fc;
print STDERR "#  ", join(", ", @ncols), "\n";

my %ncols_lofc = ();
my %ncols_hifc = ();

for my $icols (@ncols) {

   my ($fc, $se) = split(" ", $ncols_fc{$icols});
   if ($fc < $threshhold) {
      $ncols_lofc{$icols} = "$fc $se";
   } else {
      $ncols_hifc{$icols} = "$fc $se";
   }
}

my $nlo = scalar keys %ncols_lofc;
my $nhi = scalar keys %ncols_hifc;
my @locols = sort {$a <=> $b} keys %ncols_lofc;
my @hicols = sort {$a <=> $b} keys %ncols_hifc;

my @cols_to_use = ();
if ($nlo == 0  or  $nhi == 0) {
   print STDERR "points don't surround threshhold. $nlo  $nhi \n";
   next;
}
if ($nlo == 1) {
   if ($nhi == 1) {
      # linear interpolation:
      @cols_to_use = (@locols, @hicols);
   } else {
      if ($nhi == 2) {          # 3 points total
         @cols_to_use = (@locols, @hicols);
      } else {                  # $nhi >= 3
         @cols_to_use = (@locols, @hicols[0..2]);
      }
   }
} else {                        # $nlo >= 2
   if ($nhi == 1) {
      if ($nlo == 2) {
       #  print STDERR "2,1 case\n";
         @cols_to_use = (@locols, @hicols);
      } else {                  # nlo > 2
       #  print STDERR ">2,1 case \n";
      #   print join("; ", @locols), "\n";
         @cols_to_use = (@locols[-3..-1], @hicols);
      #   print join("; ", @cols_to_use), "\n";
      }
   } else {
      @cols_to_use = (@locols[-2,-1], @hicols[0,1]); # 2 of each - usual case
   }
}


my @fcs = ();
my @ses = ();
for (@cols_to_use) {
   print STDERR "# $_ ", $ncols_fc{$_}, "\n";
   my ($fc, $se) = split(" ", $ncols_fc{$_});
   push @fcs, $fc;
   push @ses, $se;
} 

my $n_points = scalar @cols_to_use;
if ($n_points == 1) {
   die "only 1 point? \n";
} elsif ($n_points == 2) {

} elsif ($n_points == 3) {
   my $formula = 'c + b*x + a*x^2';
   my $variable = 'x';
   for (@cols_to_use) {

   }
} elsif ($n_points == 4) {
   my $solution;
   my @ns = map(int($Nreps*$_+0.1), @fcs);
   my ($max_dev, $avg_dev, $fopt) = fit(xdata => \@cols_to_use, ydata => \@fcs, terms => 2, impl_lang => "coderef");
   print STDERR "$max_dev $avg_dev \n" if($verbose); # it;
   for my $k ($cols_to_use[0]..$cols_to_use[3]){
      print STDERR "$k  ", $fopt->($k), "\n";
   }
   my ($x_soln, $dy, $iterations) = solve_binary($fopt, $threshhold, $cols_to_use[0], $cols_to_use[-1]);
   $solution = $x_soln;
   
 #  my $s = '';
   # for my $i (0..$n_points-1) {
   #    $s .= $cols_to_use[$i] . " ";
   #    $s .= $fopt->($cols_to_use[$i]) . " ";
   #    $s .= $fcs[$i] . " ";
   #    $s .= ($fcs[$i] - $fopt->($cols_to_use[$i]))/$ses[$i] . "  ";
   # }
 #  print "$s \n";

   my $M = $Ndraw;
   my ($xsum, $xsumsq) = (0.0, 0.0);
   for my $i (1..$M){
    my $ppfcs = get_pp_fs(\@ns, $Nreps, $RNG);
      my ($max_dev, $avg_dev, $fopt) = fit(xdata => \@cols_to_use, ydata => $ppfcs, terms => 2, impl_lang => "coderef");
      my ($x_soln, $dy, $iterations) = solve_binary($fopt, $threshhold, $cols_to_use[0], $cols_to_use[-1]);

 #   print STDERR "iterations, time:  ", Algorithm::CurveFit::Simple->$STATS_H->{fit_iter}, "  ", Algorithm::CurveFit::Simple->$STATS_H->{fit_time}, "\n";
      print STDERR "$x_soln $dy $iterations \n" if($verbose);
    $xsum += $x_soln;
    $xsumsq += $x_soln*$x_soln;
   }
   $xsum /= $M;
   $xsumsq /= $M;
   my $sigma = sqrt($xsumsq - $xsum*$xsum);
   my $unc = sqrt(($solution - $xsum)**2 + $sigma**2);
   print "$treename $model $branch_length $threshhold   $solution +- $unc   ";
   print "$xsum +- $sigma $M  $iterations \n";
}

}



# *********************************************************8


sub get_pp_fs{ # for each data point sample posterior distrib.
#my $xs = shift;
my $ns = shift;
my $N = shift;
my $rng = shift;
my @ppfcs = ();
for my $n (@$ns){
push @ppfcs, gsl_ran_beta($rng->raw(), $n+1, $N-$n+1);
}
# print stderr "ns: ", join(", ", @$ns), "\n";
# print stderr "ppns: ", join(", ", @ppfcs), "\n";
return \@ppfcs;
}



  sub binomial_stderr{
     my $n = shift;
     my $p = shift;
     return sqrt($p*(1.0-$p)/$n);
  }

sub solve_binary{
   my $f = shift;
   my $Y = shift;
   my $xL = shift;
   my $xR = shift;
   my $yL = $f->($xL);
   my $yR = $f->($xR);
   my $x;
   my $dy = shift // $Y/1000;
   my $res = 10000000;
   for my $i (1..100) {
      $x = 0.5*($xL + $xR);
      my $ynew = $f->($x);
      print STDERR "$x $ynew $i \n" if($verbose);
      $res = $ynew - $Y;
      return ($x, $res, $i) if(abs($res) <= $dy);
      if ($yL < $Y  and  $yR > $Y) {
         if ($res > 0.0) {
            $xR = $x;
         } else {
            $xL = $x;
         }
      } elsif($yL > $Y  and  $yR < $Y) {
         if ($res > 0.0) {
            $xL = $x;
         } else {
            $xR = $x;
         }
      }
      $yL = $f->($xL);
      $yR = $f->($xR);
    # print "$xL $xR  $yL $yR \n";
   }
  #exit;
   return ($x, $res, 101);

}
