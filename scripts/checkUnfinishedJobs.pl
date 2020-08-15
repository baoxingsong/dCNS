#!perl -w
use strict;

open INPUT, "$ARGV[0]";
while( my $line=<INPUT> ){
    if( $line=~/\s(\S+)$/ ){
        my $file = $1 . "o";
        if (-e $file) {
            my $good = 0;
            open IN , "$file";
            while( my $lin = <IN> ){
                if( $lin=~/^#done/ ){
                    $good = 1;
                }
            }
            close IN;
            if ( 0 == $good ){
                print "$line";
            }

        }else{
            print "$line";
        }
    }
}
close INPUT;
