#!/usr/bin/perl -w
use strict;
use lib '/home/tomfy/Orthologger/lib';
use TomfyMisc qw ( fasta2seqon1line );

my %codon2aa = qw(
                    TCA  S  TCC  S  TCG  S  TCT  S  TTC  F  TTT  F  TTA  L  TTG  L
                    TAC  Y  TAT  Y  TAA  _  TAG  _  TGC  C  TGT  C  TGA  _  TGG  W
                    CTA  L  CTC  L  CTG  L  CTT  L  CCA  P  CCC  P  CCG  P  CCT  P
                    CAC  H  CAT  H  CAA  Q  CAG  Q  CGA  R  CGC  R  CGG  R  CGT  R
                    ATA  I  ATC  I  ATT  I  ATG  M  ACA  T  ACC  T  ACG  T  ACT  T
                    AAC  N  AAT  N  AAA  K  AAG  K  AGC  S  AGT  S  AGA  R  AGG  R
                    GTA  V  GTC  V  GTG  V  GTT  V  GCA  A  GCC  A  GCG  A  GCT  A
                    GAC  D  GAT  D  GAA  E  GAG  E  GGA  G  GGC  G  GGG  G  GGT  G
		    ---  -	
                );

my %aa2codon = ();

while(my($codon, $aa) = each %codon2aa){
   if(! exists $aa2codon{$aa}){
      $aa2codon{$aa} = [$codon];
   }else{
      push @{$aa2codon{$aa}}, $codon;
   }
}

#while(my($aa, $codons) = each %aa2codon){
#   print "$aa ", join(', ', @{$codons}), "\n";
#}
# exit;

my $fasta = '';
while(<>){
$fasta .= $_;
}

my $f = fasta2seqon1line($fasta);

my @lines = split("\n", $f);
for(@lines){
   my $id = '';
   if(/^>(\S+)/){
      $id = $1;
      print ">$id \n";
   }else{
      my $l = $_;
      $l =~ s/\s*$//;
      while($l){
         my $aa = substr($l, 0, 1, '');
         my $codons = $aa2codon{$aa} // ['XXX'];
         my $n_codons = scalar @$codons;
         my $chosen_codon_index = rand($n_codons);
         my $codon = $codons->[$chosen_codon_index];
         print "$codon";
      }print "\n";
   }
}
