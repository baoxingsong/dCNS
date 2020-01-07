#!perl -w
use strict;

open INPUT, "$ARGV[0]";
while( my $line=<INPUT> ){
    if( $line=~/\s(\S+)$/ ){
        my $file = $1 . "queryo";
        if (-e $file) {

        }else{
            print "$line";
        }
    }
}
close INPUT;
