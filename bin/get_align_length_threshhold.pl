#!/usr/bin/perl -w
use strict;
use Algorithm::CurveFit;
use Algorithm::CurveFit::Simple qw( fit );
use Data::Dumper;
use Math::GSL::RNG;
use Math::GSL::Randist qw ( :all );

my $MAXITERATIONS = 100;

my $RNG = Math::GSL::RNG->new();
my $input_pattern = shift; # each file should have several lines (from *_rfds files) with the same bl, various numbers of cols
# typical input line:
# 1000cols/1000cols0.025-0.28_rfds:tree: t16u  bl: 0.07  cols: 1000  subst: HKY1   TvR:  1.0000  0.0000  0.00  0.0000   RvR:  1.0000  0.0000  0.00  0.0000   1000,0,0,0,0,0

# my $branch_length = shift // undef;
my @input_filenames = split(" ", `ls $input_pattern`);
my $threshhold = shift // 0.9; # interpolate to find the alignment length s.t. that topo is correct this fraction of time.
my $Ndraw = shift // 100;
my $Nreps = shift // 1000;
my $verbose = shift // 0;

for my $input_file (@input_filenames) {
   my %ncols_fc = ();
   my $treename = undef;
   my $branch_length = undef;
   my $model = undef;
   print STDERR "Input file: $input_file \n";
   open my $fhin, "<", "$input_file";
   while (<$fhin>) {
      my @cols = split(" ", $_);
      my ($tn, $bl, $ncols, $mod, $fc) = @cols[1,3,5,7,9];
      die "tree name inconsistency: $tn $treename \n" if(defined $treename  and  ($tn ne $treename));
      die "branch length inconsistency: $bl $branch_length \n" 
      # next 
        if(defined $branch_length  and  ($bl != $branch_length));
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
         if ($nhi == 2) {       # 3 points total
            @cols_to_use = (@locols, @hicols);
         } else {               # $nhi >= 3
            @cols_to_use = (@locols, @hicols[0..2]);
         }
      }
   } else {                     # $nlo >= 2
      if ($nhi == 1) {
         if ($nlo == 2) {
            #  print STDERR "2,1 case\n";
            @cols_to_use = (@locols, @hicols);
         } else {               # nlo > 2
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
      # my ($max_dev, $avg_dev, $fopt) = fit(xdata => \@cols_to_use, ydata => \@fcs, terms => 2, impl_lang => "coderef");
      # print STDERR "max, avg devs: $max_dev $avg_dev \n" if($verbose); # it;
      # # for my $k ($cols_to_use[0]..$cols_to_use[3]){
      # #    print STDERR "$k  ", $fopt->($k), "\n";
      # # }
      # my ($x_soln, $dy, $iterations) = solve_binary($fopt, $threshhold, $cols_to_use[0], $cols_to_use[-1]);
      # print STDERR "$x_soln  $dy  $iterations \n";

      get_init_quadratic(\@cols_to_use, \@fcs);
    #  if ($iterations > $MAXITERATIONS) {
         my $x_soln = linear_interp(\@cols_to_use, \@fcs, $threshhold);
    #  }
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
      my $good_count = 0;
      my ($xsum, $xsumsq) = (0.0, 0.0);
      for my $i (1..$M) {
         my $ppfcs = get_pp_fs(\@ns, $Nreps, $RNG);


        #  my ($max_dev, $avg_dev, $fopt) = fit(xdata => \@cols_to_use, ydata => $ppfcs, terms => 2, impl_lang => "coderef", iterations => 1000);
       # #  print STDERR "max, avg devs: $max_dev $avg_dev \n" if($verbose); # it;
       #   my ($x_soln, $dy, $iterations) = solve_binary($fopt, $threshhold, $cols_to_use[0], $cols_to_use[-1]);

         #   print STDERR "iterations, time:  ", Algorithm::CurveFit::Simple->$STATS_H->{fit_iter}, "  ", Algorithm::CurveFit::Simple->$STATS_H->{fit_time}, "\n";
       #  if($iterations > $MAXITERATIONS){
            my $linear_x_soln = linear_interp(\@cols_to_use, $ppfcs, $threshhold);
       #  }
     #    print STDERR join(" ", @cols_to_use), "    ", join(" ", @$ppfcs), "  ", "   $linear_x_soln \n" if($verbose);
         if (defined $linear_x_soln){ #  or  $iterations <= $MAXITERATIONS) {
            $xsum += $linear_x_soln;
            $xsumsq += $linear_x_soln**2;
            $good_count++;
         }
      }
      $xsum /= $good_count;
      $xsumsq /= $good_count;
      my $sigma = sqrt($xsumsq - $xsum*$xsum);
      my $unc = sqrt(($solution - $xsum)**2 + $sigma**2);
      print "$treename $model $branch_length $threshhold   $solution +- $unc   ";
      print "$xsum +- $sigma $good_count \n"; #  $iterations \n";
   }

}



# *********************************************************8

sub get_init_quadratic{
# idea is to get initial quadratic guess to fit
# to 4 points by dropping out each of 4 points
# one at a time, getting 4 solutions and combining ....
my $xs = shift;
my $ys = shift;
my @xs = @{$xs};
my @ys = @{$ys};

if(scalar @xs == 4  and scalar @ys == 4){
   my @x123 = @xs[1..3];
my @y123 = @ys[1..3];
   my ($a0, $b0, $c0) = quad3(\@x123, \@y123);
  my @x012 = @xs[0..2];
my @y012 = @ys[0..2];
   my ($a3, $b3, $c3) = quad3(\@x012, \@y012);

# drop point 1
  my @x023 = @xs[0,2,3];
my @y023 = @ys[0,2,3];
   my ($a1, $b1, $c1) = quad3(\@x023, \@y023);
  my @x013 = @xs[0,1,3];
my @y013 = @ys[0,1,3];
   my ($a2, $b2, $c2) = quad3(\@x013, \@y013);
}else{
die "get_init_quadratic needs 4 data points.\n";
}

}

sub quad3{
   my $xs = shift;
   my $ys = shift;
   my @xs = @{$xs};
   my @ys = @{$ys};
   #  print join(", ", @xs), "\n";
   #  print join(", ", @ys), "\n";

   my $s01 = ($ys[1] - $ys[0]) / ($xs[1] - $xs[0]);
   my $s12 = ($ys[2] - $ys[1]) / ($xs[2] - $xs[1]);
   my $s02 = ($ys[2] - $ys[0]) / ($xs[2] - $xs[0]);
   my $D = ($s01 - $s02) / ($xs[2] - $xs[1]);

   my $a = $ys[0] - $xs[0]*$s02 - $xs[0]*$xs[2]*$D;
   my $b = $s02 + ($xs[0] + $xs[2])*$D;
   my $c = -1.0*$D;

   my $m = 20;
   my $dx = ($xs[2] - $xs[0])/$m;
   for my $i (0..$m) {
      my $x = $xs[0] + $i*$dx;
      print STDERR $x, "  ", $a + $b*$x + $c*$x**2, "\n";
   }
   print STDERR "\n";
}




sub get_pp_fs{         # for each data point sample posterior distrib.
   #my $xs = shift;
   my $ns = shift;
   my $N = shift;
   my $rng = shift;
   my @ppfcs = ();
   for my $n (@$ns) {
      push @ppfcs, gsl_ran_beta($rng->raw(), $n+1, $N-$n+1);
   }
   # print stderr "ns: ", join(", ", @$ns), "\n";
   # print stderr "ppns: ", join(", ", @ppfcs), "\n";
   return \@ppfcs;
}

sub linear_interp{
   my $xs = shift;
   my $ys = shift;
   my $ytarget = shift;
   my $xresult;
   my $n = scalar @$xs;
   for my $i (0..$n-2) {
      my ($y0, $y1) = ( $ys->[$i], $ys->[$i+1] );
      # if (abs($y0 - $y1) < 0.001) {
      #    return 0.5*($xs->[$i] + $xs->[$i+1]);
      # } els
        if (($y0 <= $ytarget and $y1 >= $ytarget)) {
         $xresult = $xs->[$i] + ($ytarget - $y0)*($xs->[$i+1] - $xs->[$i])/($y1 - $y0);
    #     print "$i  ", $xs->[$i], "  ", $xs->[$i+1], "  ", "$y0  $y1    $ytarget    $xresult \n";
         return $xresult;
      }
   }
   return undef;
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
   for my $i (1..$MAXITERATIONS) {
      $x = 0.5*($xL + $xR);
      my $ynew = $f->($x);
      #     print STDERR "$x $ynew $i \n" if($verbose);
      $res = $ynew - $Y;
      return ($x, $res, $i) if(abs($res) <= $dy);
      if ($yL < $Y  and  $yR > $Y) {
         if ($res > 0.0) {
            $xR = $x;
         } else {
            $xL = $x;
         }
      } elsif ($yL > $Y  and  $yR < $Y) {
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
