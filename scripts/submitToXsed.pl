#!perl
use strict;
use warnings FATAL => 'all';

opendir(DIR1, "./") || die "Can't open directory ../";
my @dots1 = readdir(DIR1);
foreach my $file (@dots1){
    if( $file=~/_prime$/ ){
        my $inputfile = $file;
        my $bashfile = $file . ".sh";
        my $outfile = $file . ".and_CNS";
        if(-e $outfile) {
        }else{
            open OUTPUT, ">$bashfile";
            print OUTPUT "#!/usr/bin/bash\n";
            print OUTPUT "../and_CNS mult-cns -i $inputfile -o $outfile\n";
            close OUTPUT;
            system("sbatch  --mem=4GB --ntasks-per-node=1 -t 20:30:00 $bashfile");
        }
    }
}
