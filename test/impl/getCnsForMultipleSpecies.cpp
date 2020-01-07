//
// Created by Baoxing song on 2019-01-03.
//

#include "include/gtest/gtest.h"
#include "../../impl/impl.h"
#include <string>
#include <iostream>
#include <sstream>

TEST(getCnsForMultipleSpecies, c1){
//TODO compare the output coordinate with the real sequence
// check the ...Rference sequence carefully

std::map<std::string, std::string> sequences;
    std::vector<std::string> seqNames;
    //std::string filePath = "/home/bs674/Desktop/SORBI_3001G063600_SORBI_3001G063800";

    std::string filePath = "/media/bs674/2t/pairConservedGeneNeighbhood/sequenceUpStreamAndDownStream/SORBI_3001G457700_5_prime";
    readFastaFile( filePath, sequences, seqNames );
    int8_t ** seqs = new int8_t * [seqNames.size()];
    int8_t ** seq_revs = new int8_t * [seqNames.size()];
    SequenceCharToUInt8 sequenceCharToUInt8;
    std::vector<int32_t> lengths;
    std::vector<std::string> seqs_string;
    for( int i=0; i <seqNames.size(); ++i ){
        seqs[i] = sequenceCharToUInt8.seq_to_int8(sequences[seqNames[i]]);
        seq_revs[i] = sequenceCharToUInt8.rev_comp(seqs[i], sequences[seqNames[i]].length());
        lengths.push_back(sequences[seqNames[i]].length());
        seqs_string.push_back(sequences[seqNames[i]]);
    }

/*
    uint32_t windowsSize = 12;
    uint32_t mini_cns_size = 8;
    int32_t mini_cns_score = 6;
    int matrix_boundary_distance=3;
    int _open_gap_penalty = -1;
    int _extend_gap_penalty = -2;
    */


    int matchingScore = 2;
    int mismatchingPenalty = -3;
    int openGapPenalty = -4;
    int extendGapPenalty = -2;

    int32_t windowsSize = 10;
    int32_t mini_cns_seed_size = 7;
    int32_t mini_cns_score = 9;
    int32_t step_size = 2;

    int matrix_boundary_distance=3;

    int minimumNumberOfSpecies = 2; // the number of non-reference species
    std::string output = "test";

    Scorei m(matchingScore, mismatchingPenalty);

    double lambda = 0.382291;
    double kValue = 0.006662;
    int w = 10;
    int xDrop = 20;


    int32_t mini_cns_size = 20;
    double outputWithMinimumLengthPercentage = 0.8;
    getCnsForMultipleSpecies ( seqs, seq_revs, lengths, windowsSize, mini_cns_seed_size, mini_cns_score,
                               matrix_boundary_distance, openGapPenalty, extendGapPenalty, matchingScore, mismatchingPenalty, false, output,
                               sequences, seqNames, step_size, seqs_string, minimumNumberOfSpecies, m, mini_cns_size,
                               outputWithMinimumLengthPercentage, lambda, kValue, w, xDrop);

    ASSERT_EQ(0, 0);
}



TEST(getCnsForMultipleSpecies, c3){ //todo test it
//TODO compare the output coordinate with the real sequence
// check the ...Rference sequence carefully

    std::map<std::string, std::string> sequences;
    std::vector<std::string> seqNames;
    //std::string filePath = "/home/bs674/Desktop/SORBI_3001G063600_SORBI_3001G063800";

    std::string filePath = "/media/bs674/2t/pairConservedGeneNeighbhood/sequenceUpStreamAndDownStream/SORBI_3001G457700_5_prime";
    readFastaFile( filePath, sequences, seqNames );
    int8_t ** seqs = new int8_t * [seqNames.size()];
    int8_t ** seq_revs = new int8_t * [seqNames.size()];
    SequenceCharToUInt8 sequenceCharToUInt8;
    std::vector<int32_t> lengths;
    std::vector<std::string> seqs_string;
    for( int i=0; i <seqNames.size(); ++i ){
        seqs[i] = sequenceCharToUInt8.seq_to_int8(sequences[seqNames[i]]);
        seq_revs[i] = sequenceCharToUInt8.rev_comp(seqs[i], sequences[seqNames[i]].length());
        lengths.push_back(sequences[seqNames[i]].length());
        seqs_string.push_back(sequences[seqNames[i]]);
    }

/*
    uint32_t windowsSize = 12;
    uint32_t mini_cns_size = 8;
    int32_t mini_cns_score = 6;
    int matrix_boundary_distance=3;
    int _open_gap_penalty = -1;
    int _extend_gap_penalty = -2;
    */


    int matchingScore = 2;
    int mismatchingPenalty = -3;
    int openGapPenalty = -4;
    int extendGapPenalty = -2;

    int32_t windowsSize = 10;
    int32_t mini_cns_seed_size = 7;
    int32_t mini_cns_score = 9;
    int32_t step_size = 2;

    int matrix_boundary_distance=3;

    int minimumNumberOfSpecies = 2; // the number of non-reference species
    std::string output = "test";

    Scorei m(matchingScore, mismatchingPenalty);

    double lambda = 0.382291;
    double kValue = 0.006662;

    int w = 10;
    int xDrop = 20;


    int32_t mini_cns_size = 20;
    double outputWithMinimumLengthPercentage = 0.8;
    getCnsForMultipleSpecies ( seqs, seq_revs, lengths, windowsSize, mini_cns_seed_size, mini_cns_score,
                               matrix_boundary_distance, openGapPenalty, extendGapPenalty, matchingScore, mismatchingPenalty, false, output,
                               sequences, seqNames, step_size, seqs_string, minimumNumberOfSpecies, m, mini_cns_size,
                               outputWithMinimumLengthPercentage, lambda, kValue, w, xDrop);

    ASSERT_EQ(0, 0);
}

TEST(getCnsForMultipleSpecies, c2) {
    std::string str = "SORBI_3009G024600_setaria_-_III_2725295_2825294";
    str = str.substr (18,5);
//    str.insert(0, "-");
//    str.insert(str.length()-1, "-");
    std::cout << str << std::endl;
}

//
//cd wtdbg
//$Software/wtdbg2 -t 64 -i $workdir/fastq/g.fastq -fo wtdbg_g -L 5000
//$Software/wtpoa-cns -t 64 -i wtdbg_g.ctg.lay -fo wtdbg_g.ctg.lay.fa
//        minimap2 -t 64 -x map-ont -a wtdbg_g.ctg.lay.fa $workdir/fastq/g.fastq | samtools view -Sb - > wtdbg_g.ctg.lay.map.bam
//        samtools sort wtdbg_g.ctg.lay.map.bam > wtdbg_g.ctg.lay.map.srt.bam
//        samtools view wtdbg_g.ctg.lay.map.srt.bam | $Software/wtpoa-cns -t 64 -d wtdbg_g.ctg.lay.fa -i - -fo wtdbg_g_final.ctg.lay.fa
