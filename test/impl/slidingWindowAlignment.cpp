//
// Created by Baoxing Song on 2019-07-09.
//


#include "include/gtest/gtest.h"
#include "../../impl/impl.h"
#include <string>
#include <iostream>
#include <sstream>

TEST(slidingWindowAlignment, c1){

    std::map<std::string, std::string> sequences;
    std::vector<std::string> seqNames;

    std::string filePath = "/Users/bs674/SORBI_3009G024600_5_prime";
    std::string output = "/Users/bs674/SORBI_3009G024600_5_prime_slidingwindow_2";
    readFastaFile( filePath, sequences, seqNames );
    int8_t ** seqs = new int8_t * [seqNames.size()];
    SequenceCharToUInt8 sequenceCharToUInt8;
    std::vector<int32_t> lengths;
    for( int i=0; i <seqNames.size(); ++i ){
        seqs[i] = sequenceCharToUInt8.seq_to_int8(sequences[seqNames[i]]);
        lengths.push_back(sequences[seqNames[i]].length());
    }

    int matchingScore = 2;
    int mismatchingPenalty = -3;
    int openGapPenalty = -4;
    int extendGapPenalty = -2;

    int32_t windowsSize = 150;
    int32_t step_size = 50;
    Scorei m(matchingScore, mismatchingPenalty);
    slidingWindowAlignment ( seqs, lengths, seqNames, windowsSize, openGapPenalty, extendGapPenalty, matchingScore,
            mismatchingPenalty, output, step_size, m);

    ASSERT_EQ(0, 0);
}
