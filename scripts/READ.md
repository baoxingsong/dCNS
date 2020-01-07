maskFastaUsingJellyfishCountDump.pl is used to do hard masking given a kmer frequency profile producde by  jellyfish
````
kat_jellyfish count -m 31 -s 100M -t 12 -C Zea_mays.AGPv4.dna.toplevel.fa
kat_jellyfish dump mer_counts.jf > mer_counts_dumps.fa
````
you could mask all the high frequency k-mer in the fasta file using command
````
perl /home/bs674/Dropbox/andropogoneae-conservation/SSW/scripts/maskFastaUsingJellyfishCountDump.pl Zea_mays.AGPv4.dna.toplevel.fa mer_counts_dumps.fa 9 B73_v4_k31_9.fa n
````
All the sequence will be up cased.
All the k-mer with frequency equal to or larger than 9 will be replace with low case `n` and output to B73_v4_k31_9.fa.
