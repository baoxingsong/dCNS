#!perl
use strict;
use warnings FATAL => 'all';

opendir(DIR1, $ARGV[0]) || die "Can't open directory $ARGV[0]";
my @dots = readdir(DIR1);
foreach my $file (@dots) {
    $file=~s/\s$//g;
    my $output = 0;
    if ( $file=~/^\d+\.5query.sam$/ ) {
        system("cat $ARGV[0]\/$file");
    }
}
