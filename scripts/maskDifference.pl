#perl
use strict;
use warnings;

my %codes;

open INPUT, "$ARGV[0]";
while( my $line=<INPUT> ){
    my ($name, $sorghum, $maize) = split(/\s+/, $line);
    my $code = 0;
    if( $sorghum=~/n/ || $sorghum=~/N/ ){
        $code = $code+1;
    }
    if( $maize=~/n/ || $maize=~/N/ ){
        $code = $code+2;
    }
    if( exists $codes{$code} ){
        $codes{$code} = $codes{$code} + 1;
    }else{
        $codes{$code} = 1;
    }
}
close INPUT;

while( my ($key, $value) = each %codes ){
    print "$key\t$value\n";
}
