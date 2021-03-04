# dCNS
dCNS is a sensitive sequence alignment implementation to detect conserved non-coding sequence (CNS) by using a k-mer match free and dynamic programming sequence alignment strategy.



## Install
### Dependencies
GNU GCC >=6.0 
Cmake >= 3.0
```
git clone https://github.com/baoxingsong/dCNS.git
cd dCNS
cmake CMakeLists.txt
make
```
usage and parameters are available from the command line.

#### About the default parameters:

**The majority of the sequence alignment parameters dependent on each other. If you want to change a single parameter, other parameters might should be adjusted accordingly.**

Here explains how did we optimize the parameters for the sequence alignment of *Andropogoneae* non-coding sequence. If you are working with a different population, the default parameters might do not work very well.

The match score, mis-match, gap open penalty, gap extension penalty were modified from minimap2 parameters by keeping in mind the high diversity of the *Andropogoneae* non-coding regions.

The k and lambda are used for p-value calculation, and they were calculated by fitting the maximum smith-waterman scores of random sequence fragments in to a non-linear least square regression model. The non-linear least square model is available as an [R script in this release](/home/bs674/Dropbox/dCNS/scripts/scriptToestimateKandLambda.r).

Random sampling was performed between maize and sorghum, sugarcane, Miscanthus, Setaria. All the pair-wise species gave similar k and lambda estimations. The default k and lambda parameters were determined using maize as the reference against sequence extract from all the other four species randomly.

The w and xDrop values were set by referring parameters in the example code of [Seqan library](https://www.seqan.de/).

In our dataset, the significant smith-waterman score is set as ~54, and we used a smith-waterman score 40 as a minimum score of seed. The seed window size 38 was selected to make sure there is only one seed in each window. The step_size value 8 to make sure there is no seed missing for each window sliding. The openGapPenalty2 value -45 basing on the guess that the normal gap is in general <20bp and TEs are >25bp

K-mer masking, we used 20-mer, since 20 is the minimum seed size. (40/2), 40 is the minimum seed score, 2 is the match score.

**Another example for setting parameters for Brassicaceae species is available at [here](./docs/parameters.pdf).
## If run dCNS on different machine
if you run it on different machine, it would be computationally efficient to recompile the code for each machine. dCNS uses hardware instructions (SSE4 and AVX2) to speed up.

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


# k-mer masking
We are using [KAT](https://kat.readthedocs.io/en/latest/) to count k-mers.
```
kat hist -t 12 Setaria_italica.Setaria_italica_v2.0.dna.toplevel.fa -m 20 -o setaria.kat.m20.hist
```

Check the distribution using [R](https://www.r-project.org/)

```
library(ggplot2)
data = read.table("setaria.kat.m20.hist")
ggplot(data=data, aes(x=V1, y=V2)) + geom_line() + xlim(3, 1000) + ylim(0, 1000000)
ggplot(data=data, aes(x=V1, y=V2)) + geom_line() + xlim(3, 200) + ylim(0, 1000000)
ggplot(data=data, aes(x=V1, y=V2)) + geom_line() + xlim(3, 100) + ylim(0, 1000000)
ggplot(data=data, aes(x=V1, y=V2)) + geom_line() + xlim(3, 50) + ylim(0, 1000000)
```

According to this plot, pickup a k-mer frequency threshold.

Currently, we use the sencondary drivate to pickup the threshold. (It maybe not the best solution. Other ideas are welcome to test.)
```
data = read.table("setaria.kat.m20.hist")
d = data.frame(x=diff(diff(data$V2)), y=1:(nrow(data)-2))
which.min(d$x^2+d$y^2)

kat_jellyfish count -m 20 -s 100M -t 12 -C Setaria_italica.Setaria_italica_v2.0.dna.toplevel.fa -o setaria_k20_count.js
kat_jellyfish dump setaria_k20_count.js > setaria_k20_count_dumps.fa
```

mapping the CDS sequence of reference to the query genome sequence using [minimap2] (https://github.com/lh3/minimap2)
```
python3 ./scripts/longestTranscript.py -g Zea_mays.AGPv4.dna.toplevel.fa -f Zea_mays.AGPv4.34.gff3 -t false -o gene.fa
minimap2 -ax splice -a -uf -C 1 -k 12 -P -t 12 --cs Setaria_italica.Setaria_italica_v2.0.dna.toplevel.fa gene.fa > setaria.sam
```

mask the query genome
```
dCNS maskGenome -i Setaria_italica.Setaria_italica_v2.0.dna.toplevel.fa -o masked_setaria_k20_50.fa -s setaria.sam -c gene.fa -k setaria_k20_count_dumps.fa -f 33
```

mask the reference genome
```
dCNS maskGenome -i Zea_mays.AGPv4.dna.toplevel.fa -o masked_B73_v4_k20_46_cds.fa -g Zea_mays.AGPv4.34.gff3 -k b73_m20_mer_counts_dumps.fa -f 46
```
#### CNS calling
prepare the sequence for alignment
```
python3 ./scripts/extractInterGeneticSequence/sequenceUpStreamGeneAndDownStreamV2.py -g Zea_mays.AGPv4.34.gff3 -r Zea_mays.AGPv4.dna.toplevel.fa -c gene.fa -q Setaria_italica.Setaria_italica_v2.0.dna.toplevel.fa -s setaria.sam -o dCNS_setaria_maize_V2
```



go the folder and generate commands, the genome file must be the masked genome
```
cd dCNS_setaria_maize_V2
ls | awk '{print("dCNS cut1Gap -ra masked_B73_v4_k20_46.fa -qa masked_A1025_k20_57.fa -i "$1" -r reference -o "$1".5")}' > command1
```

run all the commands with in file `command1` file. I use [GNU parallel](https://www.gnu.org/software/parallel/) to run it.

check if every command finished successfully
```
perl ../scripts/checkUnfinishedJobsV2.pl command1  > command1_missing
```

goto the farther folder and combine those sam outputs
```
cd ../
perl ./scripts/combineCnsSamFiles.pl  dCNS_setaria_maize_V2 > 5.sam
```

## OUTPUT
The output file is in sam format, and it works with majority functions implemented in tools compatible with sam format.
The 5th column is the sequence alignment score. There would be some information lost when converting the sam files into bam files.
The 6th column is always start with regex `[0-9]+H` , which tells the coordinate where query sequence alignment starts from, the value is 1 based coordinate.

Those SAM files could be reformatted into bam using [samtools](http://www.htslib.org/)
```
cat 5.sam| sort | uniq | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}'  | samtools view -O BAM --reference masked_B73_v4_k20_46_cds.fa | samtools sort > 5.bam
samtools index 5.bam
```

To reformat bam file into bed format, we excluded those reference genome masked regions. We did this by using [bedtools](https://github.com/arq5x/bedtools2/) and [SeqKit](https://github.com/shenwei356/seqkit).
```
seqkit locate -F --only-positive-strand --bed -m 0 -p n masked_B73_v4_k20_46_cds.fa > ns.bed
bedtools merge -i ns.bed > ns_megered.bed
bamToBed -i 5.bam | bedtools sort -i | bedtools merge > 5.bed 
bedtools subtract -a 5.bed -b ns_megered.bed > 5_nons.bed
```
* please NOTE. Each line in the bed file could NOT be interrupted as a CNS
## Multiple sequence alignment
If you are interested in multiple sequence alignment, Firstly, you should perform pair-wise sequence alignment for each species against reference species.
Give each sam file a unique name. And use this command to generate multiple sequence alignment.
`dCNS multCns -i Zea_mays.AGPv4.dna.toplevel.fa -o msa.fasta -s sorghum.sam setaria.sam miscanthus.sam sugarcane.sam 1013.sam 1025.sam`

# Funding
This work is funded by NSF #1822330

# Citation
The dCNS manuscript is under review.
