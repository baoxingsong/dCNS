# dCNS
A k-mer free and dynamic programming sequence alignment approach to detect CNS elements.



## Install
### Dependencies
GNU GCC >=6.0 
Cmake >= 3.0
```
cmake CMakeLists.txt
make
```
usage and parameters are available from the command line.

#### About the default parameters:

Majority of the sequence alignment parameters are related with each other. If you want to change a single parameter, other parameters might should be adjusted accordingly.

Here explains how did we optimize the parameters for the sequence alignment of andropogoneae non-coding sequence. If you are working different population, the default parameters might does not work very well.

The match score, mis-match, open gap penalty, extend gap penalty was modified from minimap2 parameters by keeping in mind the high diversity of the andropogoneae non-coding regions. We fit the score distribution into a non-linear model using the non-linear least square model to test do those parameters work for our data.

The k and lambda are used for p-value calculation, they were calculated using a non-linear least square regression approach with smith-waterman scores of random sequence fragments. 

Random sampling was performed between maize and sorghum, sugarcane, miscanthus, setaria. All the pair-wise species gave similar k and lambda estimations. The default k and lambda parameters were determined using maize against sequence extract from all the other 4 species randomly.

The zDrop value allows a large gap and a few of nearby penalties.

The w and xDrop values were set by referring parameters in the example code of Seqan library.

The significant smith-waterman score should be ~54, we used a smith-waterman score 40 as a minimum score of seed. The seed window size 38 was selected to make sure there is only one seed in each window. The step_size value 8 to make sure there is no seed missing for each window sliding. The openGapPenalty2 value -45 basing on the guess that the normal gap are in general  <20bp and TEs are >25bp
 

K-mer masking, we used 20-mer, since 20 is the minimum seed size. (40/2), 40 is the minimum seed score, 2 is the match score.


## If run dCNS on different machine
It is good to recompile the code, if you run it on different machine. dCNS uses hardware instructions (SSE4 and AVX2) to speed up.

## Run dCNS under docker

```
docker1 run -d -i -t ubuntu /bin/bash
docker1 images
docker1 ps
docker1 exec -it 3f97507f8d9a /bin/bash
```

go to docker
```
apt update
apt install build-essential -y
apt install cmake -y
apt install xjobs -y
git clone https://github.com/baoxingsong/dCNS.git
cd /workdir/dCNS
cmake ./
make
```


#k-mer masking

```
kat hist -t 12 /media/bs674/2t/genomeSequence/Setaria_italica/Setaria_italica.Setaria_italica_v2.0.dna.toplevel.fa -m 20 -o setaria.kat.m20.hist
```

Check the distribution using R

```
library(ggplot2)
data = read.table("setaria.kat.m20.hist")
ggplot(data=data, aes(x=V1, y=V2)) + geom_line() + xlim(3, 1000) + ylim(0, 1000000)
ggplot(data=data, aes(x=V1, y=V2)) + geom_line() + xlim(3, 200) + ylim(0, 1000000)
ggplot(data=data, aes(x=V1, y=V2)) + geom_line() + xlim(3, 100) + ylim(0, 1000000)
ggplot(data=data, aes(x=V1, y=V2)) + geom_line() + xlim(3, 50) + ylim(0, 1000000)
```

According to this plot, pickup a k-mer frequency threshold.

Currently, we use the sencondary drivate to pickup the threshold. (It maybe not the best solution)
```
data = read.table("setaria.kat.m20.hist")
d = data.frame(x=diff(diff(data$V2)), y=1:(nrow(data)-2))
which.min(d$x^2+d$y^2)

kat_jellyfish count -m 20 -s 100M -t 12 -C /media/bs674/2t/genomeSequence/sugarcane_tareploid/Sspon.HiC_chr_asm.fasta -o Sspon_k20_count.js
kat_jellyfish dump Sspon_k20_count.js > Sspon_k20_count_dumps.fa
```

mapping the CDS sequence of reference to the query genome sequence
```
python3 longestTranscript.py -g /media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.dna.toplevel.fa -f /media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.34.gff3 -t false -o gene.fa

minimap2 -ax splice -a -uf -C 1 -k 12 -P -t 12 --cs /media/bs674/2t/genomeSequence/Setaria_italica/Setaria_italica.Setaria_italica_v2.0.dna.toplevel.fa ../gene.fa > setaria.sam
```

mask the query genome
```
and_CNS maskGenome -i /media/bs674/2t/genomeSequence/Setaria_italica/Setaria_italica.Setaria_italica_v2.0.dna.toplevel.fa -o masked_setaria_k20_50.fa -b setaria.sam -c ../gene.fa -k setaria_k20_count_dumps.fa -f 33
```

mask the reference genome
```
and_CNS maskGenome -i /media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.dna.toplevel.fa -o masked_B73_v4_k20_46_cds.fa -g /media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.34.gff3 -k b73_m20_mer_counts_dumps.fa -f 46
```

prepare the sequence for alignment
```
python3 /home/bs674/Dropbox/andropogoneae-conservation/SSW/scripts/extractInterGeneticSequence/sequenceUpStreamGeneAndDownStreamV2.py -g /media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.34.gff3 -r /media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.dna.toplevel.fa -c ../gene.fa -q /media/bs674/2t/genomeSequence/Setaria_italica/Setaria_italica.Setaria_italica_v2.0.dna.toplevel.fa -s /media/bs674/1_8t/AndCns/maskGenomeForGenomeAlignment/Setaria_italica/setaria.sam -o and_cns_setaria_maize_V2
```



go the folder and generate commands, the genome file must be the masked genome

```
cd and_cns_setaria_maize_V2
ls | awk '{print("and_CNS cut1Gap -ra masked_B73_v4_k20_46.fa -qa masked_A1025_k20_57.fa -i "$1" -r reference -o "$1".5")}' > command1
```

run all the commands with in file `command1` file. I use GNU parallel to run it.

check if every command finished successfully
```
perl ./scripts/checkUnfinishedJobsV2.pl command1  > command1_missing
```

goto the farther folder and combine those sam outputs
```
cd ../
perl /home/bs674/Dropbox/andropogoneae-conservation/SSW/scripts/combineCnsSamFiles.pl  and_cns_setaria_maize_V2 > 5.sam
```

reformat sam file into bam file
```
cat 5.sam| sort | uniq >5_uniq.sam
samtools view -o 5_.bam -O BAM --reference /media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.dna.toplevel.fa  5_uniq.sam; samtools sort -O BAM  -o 5.bam 5_.bam; samtools index 5.bam
```

# Citation
manuscript under preparation
