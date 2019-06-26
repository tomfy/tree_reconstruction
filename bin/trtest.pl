#!/usr/bin/perl -w 
use strict;
use Getopt::Long;
use English;


my $tree_file = undef;
my $branch_length_string = '0.07, 0.01, 0.14, 0.2, 0.28'; # can be list to loop over, e.g.: '0.05, 0.07, 0.1, 0.14, 0.2'
my $n_alignment_columns = 100;
my $substitution_model = 'WAG'; # or 'JC'
my $n_cats = 1;                 # number or rate categories
my $n_replicates = 10;
my $seed = undef;

GetOptions(
           'tree_file=s' => \$tree_file, # e.g. 't8.newick'
           'Ls|branch_lengths=s' => \$branch_length_string,
           'cols|columns=i' => \$n_alignment_columns,
           'substitution_model=s' => \$substitution_model,
           'ncats|categories=i' => \$n_cats,
           'reps|replicates=i' => \$n_replicates,
           'seed=i' => \$seed,
          );

if (defined $seed) {
   srand($seed);
} else {
   $seed = srand();
}
$branch_length_string =~ s/\s+//g;
my @branch_lengths = split(",", $branch_length_string);
#print join("; ", @branch_lengths), "\n";
#exit;
#$substitution_model = 'HKY' if($substitution_model eq 'JC'); # We are not specifying frequencies (so they will be equal) and not
# specifying transition/transversion ratio (so will be 1); as a result HKY reduces to JC

for my $branch_length (@branch_lengths) {
   my $seqgen_substmodel = $substitution_model;
   my $seqgen_bl = $branch_length;
   my $temp_alignment_filename = 'tmp_align_' . $PID;
   if ($substitution_model eq 'JC') {
      $seqgen_substmodel = 'HKY';
   } elsif ($substitution_model eq 'CFN') {
      $seqgen_substmodel = 'HKY';
      $seqgen_bl *= 1.5;
   }
   my $seqgen_options = "-m" .
     "$seqgen_substmodel ". 
       "-l $n_alignment_columns " .
         "-s $seqgen_bl "; 
   #	"-wr ";

   print STDERR "seqgen options: $seqgen_options \n";
   #   my $seqgen_substmodel = $substitution_model;
   #    my $seqgen_bl = $branch_length;
   # if($substitution_model eq 'CFN'){
   #    $seqgen_substmodel = 'HKY';
   #    $seqgen_bl *= 1.5;
   #    }
   my $params_out = "tree: $tree_file  bl: $branch_length  cols: $n_alignment_columns  subst: $seqgen_substmodel$n_cats";

   my $param_str = 
     $tree_file . "_" . $branch_length . "_" . 
       $substitution_model . "_" . $n_cats . "_". 
         $n_alignment_columns;

   my $true_tree_string = `cat $tree_file`;
   $true_tree_string =~ s/\s+//g;
   #my $outgroup = ($true_tree_string =~ /(A0+):/)? $1 : die "No all zero leaf. ?\n";

   my $output_trees_string = "#  $n_replicates $seed   seq-gen $tree_file $seqgen_options \n";
  
  
   for my $irep (1..$n_replicates) {

      my $temp_alignment_string = `seq-gen $tree_file $seqgen_options -z $seed `;
     #  print "AAA: $substitution_model \n", "$temp_alignment_string \n"; #exit;
      
      if ($substitution_model eq 'CFN') { # get binary string:  A, G -> 1, C, T -> 0;
         my @new_lines = ();
         my @lines = split("\n", $temp_alignment_string);
         for my $a_line (@lines) {
        #    print "XXX: $a_line  ";
            if ($a_line =~ /^(A\S+\s+)(\S+)\s*$/) {
               
               my ($id, $s) = ($1, $2);
           #    print "s: $s \n";
               $s =~ s/[AG]/1/g;
               $s =~ s/[CT]/0/g;
             #  print "aa $id  ss: $s \n";
               $a_line = $id . $s;
            }
            push @new_lines, $a_line;
         }
         $temp_alignment_string = join("\n", @new_lines);
      }
      open my $fhout, ">", "$temp_alignment_filename";
      print  $fhout "$temp_alignment_string\n";
      close $fhout;
    #  exit;
         

      $seed += 100;
      my $run_param_str = $param_str . "_$irep";
      my $model_str;
      if($substitution_model eq 'JC'){
         $model_str = ' -m GTRCAT --JC69 -V ';
      }elsif($substitution_model eq 'CFN'){
         $model_str = ' -m BINCAT -V ';
      }elsif($substitution_model eq 'WAG'){
         $model_str = ' -m PROTCATWAGF -V '
      }
      my $rax_stdout = `raxmlHPC-PTHREADS-SSE3 $model_str -n $run_param_str -s $temp_alignment_filename -p $seed -c $n_cats -T 2 ` ; # -o $outgroup `;
      $seed += 100;
      print STDERR $rax_stdout;

      my $raxml_output_treefile = "RAxML_bestTree." . $run_param_str;

      $output_trees_string .= `cat $raxml_output_treefile`;
      system "rm RA*";
   }

   my $output_trees_filename = $param_str . "_" . $n_replicates . ".trees";
   open my $fhout, ">", "$output_trees_filename";
   print $fhout $output_trees_string;
   close $fhout;

   my $x = "$params_out   ";
   $x .= `perl ~/Tree_reconstruction_studies/bin/RFdists.pl $tree_file $output_trees_filename `;
   print "$x";

}

