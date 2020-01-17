#!perl -w
use strict;
use warnings FATAL => 'all';


# AN ARRAY
my @seqs;

my $chrName="";
open INPUT, "msa";
while( my $line=<INPUT> ){
	if( (length($line) < 2)  ){ ## empty line
		if( 3 == (scalar @seqs) ){
			my $outputSeq="";
			for( my $i=0; $i < length($seqs[0]); $i = $i + 1 ){
				my $all_identical = 1;
				for my $seq( @seqs ){
					if ( substr($seq, $i, 1) eq substr($seqs[0], $i, 1) ){

					}else{
						$all_identical = 0;
					}
				}
				if( 1== $all_identical ){
					$outputSeq = $outputSeq . substr($seqs[0], $i, 1);					
				}elsif ( length($outputSeq)>0) {
					print "$chrName\t$outputSeq\n";
					$outputSeq = "";
				}
			}
			print "$chrName\t$outputSeq\n";
		}
		@seqs = ();
		$chrName="";
	}elsif( $line=~/^>(\S+)/ ){
		if( length($chrName)<1 ){
            $chrName = $line;
		    $chrName =~s/\s$//g;
	    }
    }else{
		$line=~s/\s//g;
		push(@seqs, $line);
	}
}
close INPUT;
