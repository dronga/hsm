#!/usr/bin/perl -w 

use strict;


my $infofile = shift or die "info file not specified\n";
open my $infofh, "$infofile" or die "could not open haploview info file $infofile\n";

my @info;
my @rs;
while(<$infofh>) {
#        $_ =~ m/^(\d+)/;
#        push @info , $1;
        chomp $_;
	my @infopos = (split(/\t/, $_));
        push @info , $infopos[1];
	push @rs , $infopos[0];

}
close($infofh);


my $methfile = shift or die "methyl.list file not specified\n";
open my $methfh, "$methfile" or die "could not open haploview info file $methfile\n";

my @cpgsnps;
#only read the first column and store in an array
while(my $methline = <$methfh>) {
	$methline =~ m/^(\d+)/;
	push @cpgsnps , $1;

}


my $blockfile = shift or die "GABIRELblocks file not specified\n";
open my $blockfh, "$blockfile" or die "could not open haplotype block file $blockfile\n";

my $snp = shift or die "SNP ID not specified\n";
my $snpmarker = which($snp, @info) + 1; #haploview marker,  1 indexed

print $snpmarker;

my $psid = shift or die "Process ID (file prefix)  not specified\n";

print "we are looking for $snpmarker...\n";

my $blockcount = 0;
my @haps;    #haplotype strings
my @mhaps;   #haplotype strings reduced to only CpG SNPs 
my @markers; #SNPs per block
my @afs;      #allele frequencies per block

my @hblocks;   #AoA of the strings representing the full haps
my @mhblocks;  #AoA of the strings representing the methylation haps
my @mkrblocks; #AoA of the markers in each block
my @afblocks;  #AoA of the alle frequecies

my @indices; #the indices in of the SNPs in @markers which have CpGs

my $snpblock      = -1; #the block in which the SNP is, to be found later
my $snpblockindex = -1; #where is the marker relative to the hap block?;

my @hapcpgs; #SNPs in CpGs per blocl
my @hcgblocks; #AoA of the CpG markers in each block

while (my $line = <$blockfh>) {
	
	#when arriving at the next haplotype block put everything collected into the container arrays
	if(substr($line,0,5) eq "BLOCK"){
	if($. != 1) { save_block(); }
	#extract the marker numbers and save
	$line =~ m/^BLOCK\ (\d+)\.\ +MARKERS:\ (.*$)/;
	@markers = split(/\ /, $2);

	#which markers are cpgsnps?
	@hapcpgs = intersection(@cpgsnps,@markers);
	

	#and what's the index of these markers?
	@indices = ();
	foreach(@hapcpgs) {
		push @indices, which($_, @markers);
	}
	
	#print "$blockcount: @markers; @indices; @hapcpgs\n";

	print "$info[$markers[0] - 1]  $info[$markers[$#markers] - 1]\n";
	print "$#markers\n";
        #is this marker my snp?	
	#if (grep {$_ eq $snpmarker} @markers)
	if ($info[$markers[0] - 1] < $snp && $info[$markers[$#markers] - 1] > $snp) 
	{
		$snpblock = $blockcount;
		$snpblockindex = which($snpmarker,@markers);
		print "The SNP is marker $snpblockindex in block $blockcount\n";
		print "$blockcount: @markers; @indices; @hapcpgs\n";
		
	}

	#delete everything that was stored from previous haplotype block
	$blockcount++;
	@haps  = ();
	@mhaps = ();
	@afs   = ();
	
	}
	
	#read in the haplotypes
	elsif($line =~ m/^([1-4]+) \(([0-9]+.[0-9]+)/) {
		push @haps, $1;
		push @afs, $2;
		push @mhaps, getchar($1, @indices);
	}
}
	


#put last collection of haps into container

save_block();
close($blockfh);

#my $temp = length(@{@hblocks[$snpblock]}[1]); #$snpblockindex, $snpblockindex+1);
#print "$temp\n";

#write output to file

if($snpblock==-1) {
	print "SNP is in no hap block";
      }
else {
    print ">meth.haps/" . $psid . "." . $snp . ".meth.haps";
  open my $outfh, ">meth.haps/" . $psid . "." . $snp . ".meth.haps" or die "could not open file temp.meth.haps for writing";
  print $outfh "snp.allele meth.hap hap.allele.freq cpg.score\n";

  #get the single alleles of the SNP corresponding to each hap;
  my @mhapalleles =();
  foreach (@{@hblocks[$snpblock]}) {
    my $allele = substr($_, $snpblockindex, 1);
    push @mhapalleles, $allele;
  }

  for my $i (0..$#mhapalleles) {
    my $count = ()= @{@mhblocks[$snpblock]}[$i] =~ /[23]/g;
    print $outfh "@mhapalleles[$i] @{@mhblocks[$snpblock]}[$i] @{@afblocks[$snpblock]}[$i] $count\n";
  }
  close($outfh);

open my  $outfh, ">hap.pos/" . $psid . "." . $snp . ".hap.pos" or die "could not open file temp.hap.pos for writing";
  print $outfh "start end\n";
  print $outfh "@info[@{$mkrblocks[$snpblock]}[0]] @info[@{$mkrblocks[$snpblock]}[$#{$mkrblocks[$snpblock]}]]\n";
  
  close($outfh);

open my  $outfh, ">hap.markers/" . $psid . "." . $snp . ".hap.markers" or die "could not open file temp.hap.markers for writing";
##  print $outfh "start end\n";
##  print $outfh "@{@hcgblocks[$snpblock]}\n";

for my $i (0..$#{@hcgblocks[$snpblock]}) {
    print $outfh "${@hcgblocks[$snpblock]}[$i] @rs[${@hcgblocks[$snpblock]}[$i] - 1]\n";
  }

  close($outfh);

}



sub save_block {
	push @hblocks,   [ @haps ];
	push @mhblocks,  [ @mhaps];
	push @mkrblocks, [ @markers ];
	push @afblocks,  [ @afs ];
	push @hcgblocks,  [ @hapcpgs ];
}

#returns values that are in both arrays (adapted from perldoc)
sub intersection {
    my (@array1, @array2) = @_;
    my @union = my @intersection = my @difference = ();
    my %count = ();
    my $element;
    foreach $element (@array1, @array2) { $count{$element}++ }
    foreach $element (keys %count) {
    push @union, $element;
    push @{ $count{$element} > 1 ? \@intersection : \@difference }, $element;
    }
    return sort @intersection;
}

#find the array index corresponding to a value
sub which {
	my($search, @array) = @_;
	my %index;
	@index{@array} = (0..$#array);
	my $index = $index{$search};
	return $index;
}

#create a new string from single characters at postions in string given an array of indices
sub getchar {
	my($string, @charpos) = @_;
	my @chars = split(//,$string);
	return join "", @chars[@charpos];
}
