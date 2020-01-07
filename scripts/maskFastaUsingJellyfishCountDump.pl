#!perl -w
use strict;
use warnings FATAL => 'all';

print "parameters: fasta kmerPropfile countThreshold output alt\n";
print "            k-mer with frequency >= countThreshold will be masked (replace with alt)\n";

my $fasta = $ARGV[0];
my $kmerPropfile = $ARGV[1];
my $coutThreshold = $ARGV[2];
my $output = $ARGV[3];
my $alt = $ARGV[4];


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

my %ifMask;
while ( my ($k, $v) = each %seqs ){
    $ifMask{$k} = "0" x length($v);
}

my $count="";
my $kmer="";
open INPUT, "$kmerPropfile";
while( my $line=<INPUT>  ){
    if( $line=~/>(\S+)/  ){
        if( length($count) > 0 && $count >= $coutThreshold  ){
            $kmer=~s/\s//g;
            $kmer = uc ($kmer);
            my $maskMer = '1' x length($kmer);
            while( my ($k, $v) = each %seqs ){
                my $v2 = $ifMask{$k};
                while( $v=~/$kmer/g ){
                    substr($v2, $-[0], length($kmer), $maskMer);
                }
                $ifMask{$k}=$v2;
            }
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
    while( my ($k, $v) = each %seqs ){
        my $v2 = $ifMask{$k};
        my $maskMer = '1' x length($kmer);
        while( $v=~/$kmer/g ){
            substr($v2, $-[0], length($kmer), $maskMer);
        }
        $ifMask{$k}=$v2;
    }
}

open OUTPUT, ">$output";
while( my ($k, $v) = each %seqs ){
    print OUTPUT ">$k\n";
    my $v2 = $ifMask{$k};
    for( my $i=0; $i<length($v2); ++$i ){
        if (substr($v2, $i, 1) ==0){
            print OUTPUT substr($v, $i, 1);
        }else{
            print OUTPUT $alt;
        }
    }
    print OUTPUT "\n";
}
close OUTPUT;
