#!perl -w
use strict;
use warnings FATAL => 'all';

print "parameters: fasta kmerPropfile k-mer-length countThreshold output alt\n";
print "            k-mer with frequency >= countThreshold will be masked (replace with alt)\n";

my $fasta = $ARGV[0];
my $kmerPropfile = $ARGV[1];
my $kmerLength = $ARGV[2];
my $coutThreshold = $ARGV[3];
my $output = $ARGV[4];
my $alt = $ARGV[5];

# read genome begin
my %seqs;
my $name="";
my $seq="";
open INPUT, "$fasta";
while( my $line=<INPUT>  ){
    if( $line=~/>(\S+)/  ){
        if( length($name) > 0  ){
            $seq=~s/\s//g;
            $seq = uc ($seq);
            $seqs{$name} = $seq;
        }
        $name=$1;
        $seq="";
    }else{
        $seq = $seq . $line;
    }
}
close INPUT;
if( length($name) > 0  ){
    $seq=~s/\s//g;
    $seq = uc ($seq);
    $seqs{$name} = $seq;
}
# read genome end

#masking begin
my %ifMask;
while ( my ($k, $v) = each %seqs ){
    $ifMask{$k} = "0" x length($v);
}
#masking end

#read k-mer begin
my $count="";
my $kmer="";
my %kmers;
open INPUT, "$kmerPropfile";
while( my $line=<INPUT>  ){
    if( $line=~/>(\S+)/  ){
        if( length($count) > 0 && $count >= $coutThreshold  ){
            $kmer=~s/\s//g;
            $kmer = uc ($kmer);
            $kmers{$kmer} = $count;
        }
        $count=$1;
    }else{
        $kmer =  $line;
    }
}
close INPUT;
if( length($count) > 0  && $count >= $coutThreshold ){
    $kmer=~s/\s//g;
    $kmer = uc ($kmer);
    $kmers{$kmer} = $count;

}
#read k-mer end

my $maskMer = '1' x $kmerLength;

#change masking string begin
while( my ($k, $v) = each %seqs ){
    print "$k\n";
    my $v2 = $ifMask{$k};
    for( my $i=0; $i<= length($v)-$kmerLength; ++$i ){
        if (exists $kmers{ substr($v, $i, $kmerLength) }){
            substr($v2, $i, $kmerLength, $maskMer);
        }
    }
    $ifMask{$k}=$v2;
}
#change masking string end

#output result begin
open OUTPUT, ">$output";
while( my ($k, $v) = each %seqs ){
    print OUTPUT ">$k\n";
    my $v2 = $ifMask{$k};
    for( my $i=0; $i<length($v2); ++$i ){
        if ( substr($v2, $i, 1) ==0){
            print OUTPUT substr($v, $i, 1);
        }else{
            print OUTPUT $alt;
        }
    }
    print OUTPUT "\n";
}
close OUTPUT;
#output result end
