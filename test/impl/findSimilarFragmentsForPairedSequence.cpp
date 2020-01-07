//
// Created by Baoxing song on 2019-01-03.
//


#include "include/gtest/gtest.h"
#include "../../impl/impl.h"
#include <string>
#include <iostream>
#include <sstream>

#include "../../InputParser.h"
#include "../../model/model.h"

#include <limits>
#include <iomanip>
#include <cmath>
#include <algorithm>
/*

TEST(findSimilarFragmentsForPairedSequence, c2){ // just to make sure that every line has been analysed

//    int matchingScore = 4;
//    int mismatchingPenalty = -3;
//    int openGapPenalty = -6;
//    int extendGapPenalty = -4;

    int matchingScore = 2;
    int mismatchingPenalty = -3;
    int openGapPenalty = -4;
    int extendGapPenalty = -2;
    double lambda = 0.382291;
    double kValue = 0.006662;

    uint32_t seed_window_size = 10;   //40s
    uint32_t mini_cns_size = 7;
    int32_t mini_cns_score = 9;
    uint32_t step_size = 2;

//    uint32_t seed_window_size = 10; //todo test the time costing
//    uint32_t mini_cns_size = 22;
//    int32_t mini_cns_score = 44;
//    uint32_t step_size = 5;
//

    int matrix_boundary_distance = 0;
    bool onlySyntenic = false;
    double pvalue = 0.1;

    std::string input = "/Users/bs674/SORBI_3009G024600_5_prime";
    std::string reference = "SORBI_3009G024600_sorghum";
    std::string query = "SORBI_3009G024600_maize_V4_+_8_135910122_136010121";
    std::string output = "/Users/bs674/SORBI_3009G024600_5_prime_rev_o";

    Scorei m(matchingScore, mismatchingPenalty);

    std::map<std::string, std::string> sequences;
    std::vector<std::string> seqNames;
    readFastaFile( input, sequences, seqNames );

    SequenceCharToUInt8 sequenceCharToUInt8;
    std::string seq1_string = sequences[reference];
    std::string seq2_string = sequences[query];

    //seq2_string = getReverseComplementary(seq2_string);

    int8_t * seq1 = sequenceCharToUInt8.seq_to_int8(seq1_string);
    int8_t * seq2 = sequenceCharToUInt8.seq_to_int8(seq2_string);
    int8_t * seq1_rev_com = sequenceCharToUInt8.rev_comp(seq1, seq1_string.length());
    int8_t * seq2_rev_com = sequenceCharToUInt8.rev_comp(seq2, seq2_string.length());

    uint32_t length1 = seq1_string.length();
    uint32_t length2 = seq2_string.length();


    int lDiag = 0;
    int uDiag = 20;
    int xDrop = 20;

    std::vector<PairedSimilarFragment> pairedSimilarFragments0 =
            findSimilarFragmentsForPairedSequence (seq1, seq1_rev_com, seq2, seq2_rev_com, length1, length2,
                                                   seed_window_size, mini_cns_size, mini_cns_score,
                                                   matrix_boundary_distance, openGapPenalty, extendGapPenalty, matchingScore,
                                                   mismatchingPenalty, m, step_size, seq1_string, seq2_string, pvalue,
                                                   lambda, kValue, lDiag, uDiag, xDrop);

    if ( onlySyntenic ){
        std::vector<PairedSimilarFragment> pairedSimilarFragments = syntenic(pairedSimilarFragments0);
        pairedSimilarFragments0 = pairedSimilarFragments;
    }

    std::ofstream ofile;
    ofile.open(output);
    for( int i=0; i<pairedSimilarFragments0.size(); ++i ){
        std::cout << "i " << i << " start1 " << pairedSimilarFragments0[i].getStart1() << " end1 " << pairedSimilarFragments0[i].getEnd1()
                  << " start2 " << pairedSimilarFragments0[i].getStart2() << " end2 " << pairedSimilarFragments0[i].getEnd2() << " score " << pairedSimilarFragments0[i].getScore() << " cigar ";
        ofile << pairedSimilarFragments0[i].getStart1() << " " << pairedSimilarFragments0[i].getEnd1()
              << " " << pairedSimilarFragments0[i].getStart2() << " " << pairedSimilarFragments0[i].getEnd2() << " " << pairedSimilarFragments0[i].getScore() << " " << pairedSimilarFragments0[i].getPValue() << std::endl;

        for( int j=0; j<pairedSimilarFragments0[i].getCigar().size(); ++j ){
            uint32_t cigarLength = pairedSimilarFragments0[i].getCigar()[j]>>4;
            uint32_t cigarType = pairedSimilarFragments0[i].getCigar()[j]&0xf;
            std::cout << cigarLength << "MID"[cigarType];
        }
        std::cout.precision(10);
        std::cout << " evalue: " << pairedSimilarFragments0[i].getEValue() << " pvalue: " << pairedSimilarFragments0[i].getPValue() << std::endl;

        //        std::cout << getSubSequence(seq1_string, pairedSimilarFragments0[i].getStart1(), pairedSimilarFragments0[i].getEnd1()) << std::endl;
        //        std::cout << getSubSequence(seq2_string, pairedSimilarFragments0[i].getStart2(), pairedSimilarFragments0[i].getEnd2()) << std::endl;
        std::cout << pairedSimilarFragments0[i].getAlignment1() << std::endl;
        std::cout << pairedSimilarFragments0[i].getAlignment2() << std::endl;
    }
    ofile.close();
}

*/
/*

TEST(findSimilarFragmentsForPairedSequence, c22){ // just to make sure that every line has been analysed

//    int matchingScore = 4;
//    int mismatchingPenalty = -3;
//    int openGapPenalty = -6;
//    int extendGapPenalty = -4;

//    int matchingScore = 2;
//    int mismatchingPenalty = -3;
//    int openGapPenalty = -4;
//    int extendGapPenalty = -2;
//    double lambda = 0.382291;
//    double kValue = 0.006662;


    int matchingScore = 2;
    int mismatchingPenalty = -3;
    int openGapPenalty = -3;
    int extendGapPenalty = -1;
    double lambda = 0.0574;
    double kValue = 0.00001807;


    uint32_t seed_window_size = 10;   //40s
    uint32_t mini_cns_size = 7;
    int32_t mini_cns_score = 9;
    uint32_t step_size = 2;

//    uint32_t seed_window_size = 10; //todo test the time costing
//    uint32_t mini_cns_size = 22;
//    int32_t mini_cns_score = 44;
//    uint32_t step_size = 5;
//

    int matrix_boundary_distance = 0;
    bool onlySyntenic = false;
    double pvalue = 0.1;

    std::string input = "/home/bs674/SORBI_3001G121600_5_prime";
    std::string reference = "sorghum@SORBI_3001G121600";
    std::string query = "maize_V4@SORBI_3001G121600_+_1_270453646_270553645";
    std::string output = "/home/bs674/SORBI_3001G121600_5_prime_o";





    Scorei m(matchingScore, mismatchingPenalty);

    std::map<std::string, std::string> sequences;
    std::vector<std::string> seqNames;
    readFastaFile( input, sequences, seqNames );

    SequenceCharToUInt8 sequenceCharToUInt8;
    std::string seq1_string = sequences[reference];
    std::string seq2_string = sequences[query];

    //seq2_string = getReverseComplementary(seq2_string);

    int8_t * seq1 = sequenceCharToUInt8.seq_to_int8(seq1_string);
    int8_t * seq2 = sequenceCharToUInt8.seq_to_int8(seq2_string);
    int8_t * seq1_rev_com = sequenceCharToUInt8.rev_comp(seq1, seq1_string.length());
    int8_t * seq2_rev_com = sequenceCharToUInt8.rev_comp(seq2, seq2_string.length());

    uint32_t length1 = seq1_string.length();
    uint32_t length2 = seq2_string.length();


    int lDiag = 0;
    int uDiag = 20;
    int xDrop = 20;

    std::vector<PairedSimilarFragment> pairedSimilarFragments0 =
            findSimilarFragmentsForPairedSequence (seq1, seq1_rev_com, seq2, seq2_rev_com, length1, length2,
                                                   seed_window_size, mini_cns_size, mini_cns_score,
                                                   matrix_boundary_distance, openGapPenalty, extendGapPenalty, matchingScore,
                                                   mismatchingPenalty, m, step_size, seq1_string, seq2_string, pvalue,
                                                   lambda, kValue, lDiag, uDiag, xDrop);

    if ( onlySyntenic ){
        std::vector<PairedSimilarFragment> pairedSimilarFragments = syntenic(pairedSimilarFragments0);
        pairedSimilarFragments0 = pairedSimilarFragments;
    }

    std::ofstream ofile;
    ofile.open(output);
    for( int i=0; i<pairedSimilarFragments0.size(); ++i ){
        std::cout << "i " << i << " start1 " << pairedSimilarFragments0[i].getStart1() << " end1 " << pairedSimilarFragments0[i].getEnd1()
                  << " start2 " << pairedSimilarFragments0[i].getStart2() << " end2 " << pairedSimilarFragments0[i].getEnd2() << " score " << pairedSimilarFragments0[i].getScore() << " cigar ";
        ofile << pairedSimilarFragments0[i].getStart1() << " " << pairedSimilarFragments0[i].getEnd1()
              << " " << pairedSimilarFragments0[i].getStart2() << " " << pairedSimilarFragments0[i].getEnd2() << " " << pairedSimilarFragments0[i].getScore() << " " << pairedSimilarFragments0[i].getPValue() << std::endl;

        for( int j=0; j<pairedSimilarFragments0[i].getCigar().size(); ++j ){
            uint32_t cigarLength = pairedSimilarFragments0[i].getCigar()[j]>>4;
            uint32_t cigarType = pairedSimilarFragments0[i].getCigar()[j]&0xf;
            std::cout << cigarLength << "MID"[cigarType];
        }
        std::cout.precision(10);
        std::cout << " evalue: " << pairedSimilarFragments0[i].getEValue() << " pvalue: " << pairedSimilarFragments0[i].getPValue() << std::endl;

        //        std::cout << getSubSequence(seq1_string, pairedSimilarFragments0[i].getStart1(), pairedSimilarFragments0[i].getEnd1()) << std::endl;
        //        std::cout << getSubSequence(seq2_string, pairedSimilarFragments0[i].getStart2(), pairedSimilarFragments0[i].getEnd2()) << std::endl;
        std::cout << pairedSimilarFragments0[i].getAlignment1() << std::endl;
        std::cout << pairedSimilarFragments0[i].getAlignment2() << std::endl;
    }
    ofile.close();
}
*/

// most general one
TEST(findSimilarFragmentsForPairedSequence, pairCnsXExtend){ // just to make sure that every line has been analysed

    int32_t seed_window_size = 38;
    int32_t mini_cns_score = 40;
    int32_t step_size = 8;
    int matrix_boundary_distance = 0;

    int matchingScore = 2;
    int mismatchingPenalty = -3;
    int openGapPenalty1 = -4;
    int extendGapPenalty1 = -2;

    int openGapPenalty2 = -45;
    int extendGapPenalty2 = 0;

    int zDrop = 50; // this one should be slightly larger than openGapPenalty2
    int32_t bandwidth = 20000; // using very large value, 12000 there are some strange result, should be checked very carefully

    double lambda = 0.382291;
    double kValue = 0.006662;

    bool onlySyntenic = false;
    double pvalues = 0.1;

    std::string input = "/Users/bs674/SORBI_3009G024600_5_prime";
    std::string reference = "SORBI_3009G024600_maize_V4_+_8_135910122_136010121";
    std::string query = "SORBI_3009G024600_sorghum";
    std::string output = "/Users/bs674/SORBI_3009G024600_5_prime_many_to_many_seed_xtend.sam";
    std::string output2 = "/Users/bs674/SORBI_3009G024600_5_prime_many_to_many_seed_xtend";

    std::string refChr = "8";
    int32_t refStart = 135910122;
    int queStrand = 0; //1 is positive strand, 0 is negative strand

    std::string queChr="9";
    int32_t queStart=2183621;
    int refStrand = 1; //1 is positive strand, 0 is negative strand

    Scorei m(matchingScore, mismatchingPenalty);

    std::map<std::string, std::string> sequences;
    std::vector<std::string> seqNames;
    readFastaFile( input, sequences, seqNames );

    SequenceCharToUInt8 sequenceCharToUInt8;
    std::string seq1_string = sequences[reference];
    std::string seq2_string = sequences[query];

    int8_t * seq1 = sequenceCharToUInt8.seq_to_int8(seq1_string);
    int8_t * seq2 = sequenceCharToUInt8.seq_to_int8(seq2_string);
    int8_t * seq1_rev_com = sequenceCharToUInt8.rev_comp(seq1, seq1_string.length());
    int8_t * seq2_rev_com = sequenceCharToUInt8.rev_comp(seq2, seq2_string.length());

    int32_t length1 = seq1_string.length();
    int32_t length2 = seq2_string.length();

    int w = 10;
    int xDrop = 20;

    std::vector<PairedSimilarFragment> pairedSimilarFragments0 =
           findSimilarFragmentsForPairedSequence ( seq1, seq1_rev_com,
           seq2, seq2_rev_com,  length1, length2, seed_window_size,
           mini_cns_score, matrix_boundary_distance, openGapPenalty1, extendGapPenalty1, matchingScore,
           mismatchingPenalty, m, step_size, seq1_string, seq2_string, pvalues,
           lambda, kValue, w, xDrop);


    if ( onlySyntenic ){
        std::vector<PairedSimilarFragment> pairedSimilarFragments = syntenic(pairedSimilarFragments0);
        pairedSimilarFragments0 = pairedSimilarFragments;
    }

    std::ofstream ofile;
    ofile.open(output);
    std::ofstream ofile2;
    ofile2.open(output2);
    ofile << "@SQ\tSN:8\tLN:181122637" << std::endl;
    for( int i=0; i<pairedSimilarFragments0.size(); ++i ){
        if( refStrand == 1 && refStrand==queStrand ){
            int32_t quePos = queStart + pairedSimilarFragments0[i].getStart2() - 1;
            int32_t refPos = refStart + pairedSimilarFragments0[i].getStart1() - 1;
            ofile << queChr << "\t0\t" << refChr << "\t" << refPos << "\t" << pairedSimilarFragments0[i].getScore() << "\t" << quePos << "H";
            for( int j=0; j<pairedSimilarFragments0[i].getCigar().size(); ++j ){
                uint32_t cigarLength = pairedSimilarFragments0[i].getCigar()[j]>>4;
                uint32_t cigarType = pairedSimilarFragments0[i].getCigar()[j]&0xf;
                ofile << cigarLength << "MID"[cigarType];
            }
            ofile<< "*\t0\t0\t" << getSubSequence(seq2_string, pairedSimilarFragments0[i].getStart2(), pairedSimilarFragments0[i].getEnd2()) << "\t*" << std::endl;
        }else if ( refStrand == 0 && 1==queStrand){
            int32_t quePos = queStart + pairedSimilarFragments0[i].getStart2() - 1;
            int32_t refPos = refStart + seq1_string.size() -1 - pairedSimilarFragments0[i].getEnd1()+1;
            ofile << queChr << "\t16\t" << refChr << "\t" << refPos << "\t" << pairedSimilarFragments0[i].getScore() << "\t" << quePos << "H";
            for( int j=pairedSimilarFragments0[i].getCigar().size()-1; j>=0; --j ){
                uint32_t cigarLength = pairedSimilarFragments0[i].getCigar()[j]>>4;
                uint32_t cigarType = pairedSimilarFragments0[i].getCigar()[j]&0xf;
                ofile << cigarLength << "MID"[cigarType];
            }
            ofile<< "\t*\t0\t0\t" << getSubSequence(seq2_string, pairedSimilarFragments0[i].getStart2(), pairedSimilarFragments0[i].getEnd2(), 0) << "\t*" << std::endl;
        }else if( refStrand == 1 && 0==queStrand ){
            int32_t quePos = queStart + pairedSimilarFragments0[i].getStart2() - 1;
            int32_t refPos = refStart + pairedSimilarFragments0[i].getStart1() - 1;
            ofile << queChr << "\t0\t" << refChr << "\t" << refPos << "\t" << pairedSimilarFragments0[i].getScore() << "\t" << quePos << "H";
            for( int j=0; j<pairedSimilarFragments0[i].getCigar().size(); ++j ){
                uint32_t cigarLength = pairedSimilarFragments0[i].getCigar()[j]>>4;
                uint32_t cigarType = pairedSimilarFragments0[i].getCigar()[j]&0xf;
                ofile << cigarLength << "MID"[cigarType];
            }
            ofile<< "\t*\t0\t0\t" << getSubSequence(seq2_string, pairedSimilarFragments0[i].getStart2(), pairedSimilarFragments0[i].getEnd2()) << "\t*" << std::endl;
        }

        std::cout << "i " << i << " start1 " << pairedSimilarFragments0[i].getStart1() << " end1 " << pairedSimilarFragments0[i].getEnd1()
                  << " start2 " << pairedSimilarFragments0[i].getStart2() << " end2 " << pairedSimilarFragments0[i].getEnd2() << " score " << pairedSimilarFragments0[i].getScore() << " cigar ";
        ofile2 << pairedSimilarFragments0[i].getStart1() << " " << pairedSimilarFragments0[i].getEnd1()
              << " " << pairedSimilarFragments0[i].getStart2() << " " << pairedSimilarFragments0[i].getEnd2() << " " << pairedSimilarFragments0[i].getScore() << " " << pairedSimilarFragments0[i].getPValue() << std::endl;

        for( int j=0; j<pairedSimilarFragments0[i].getCigar().size(); ++j ){
            uint32_t cigarLength = pairedSimilarFragments0[i].getCigar()[j]>>4;
            uint32_t cigarType = pairedSimilarFragments0[i].getCigar()[j]&0xf;
            std::cout << cigarLength << "MID"[cigarType];
        }
        std::cout.precision(10);
        std::cout << " evalue: " << pairedSimilarFragments0[i].getEValue() << " pvalue: " << pairedSimilarFragments0[i].getPValue() << std::endl;

        //        std::cout << getSubSequence(seq1_string, pairedSimilarFragments0[i].getStart1(), pairedSimilarFragments0[i].getEnd1()) << std::endl;
        //        std::cout << getSubSequence(seq2_string, pairedSimilarFragments0[i].getStart2(), pairedSimilarFragments0[i].getEnd2()) << std::endl;
    }
    ofile.close();
    ofile2.close();
}



//general + two gap penalties

TEST(findSimilarFragmentsForPairedSequence, c3){

//
//    uint32_t seed_window_size = 10;   //40s
//    uint32_t seed_max_size = 7;
//    int32_t mini_cns_score = 9;
//    uint32_t step_size = 2;

//
//
//    uint32_t seed_window_size = 35;   //40s
//    uint32_t mini_cns_size = 1000;
//    int32_t mini_cns_score = 25;
//    uint32_t step_size = 7;
//
//

//    uint32_t seed_window_size = 75;   //40s   // those scores are not reasonable
//    uint32_t mini_cns_size = 1000;
//    int32_t mini_cns_score = 50;
//    uint32_t step_size = 10;
//
//    uint32_t seed_window_size = 20;   // there would not two seeds in the same windows, since 25/2=12 > 20/2
//    int32_t mini_cns_score = 25;  // the overlap size of two windows is 20-5, 15 is enough for a seed
//    uint32_t step_size = 5;

//    uint32_t seed_window_size = 40;
//    int32_t mini_cns_score = 50;
//    uint32_t step_size = 8;

    int32_t seed_window_size = 38;
    int32_t mini_cns_score = 40;
    int32_t step_size = 8;
    int matrix_boundary_distance = 0;


//    int matchingScore = 4;
//    int mismatchingPenalty = -3;
//    int openGapPenalty = -6;
//    int extendGapPenalty = -4;

/** those parameters could generate very good result, but very slow
    uint32_t seed_window_size = 9;   //40s
    int32_t mini_cns_score = 9;   //(9-2)*2-5 ==9   a 7 bp with one mis-match. 4*2=8<9 there could not be two seeds in the same window.
        //  9*2-2-4-2-4 = 6 <9 , there could never be a insertation + a deletion in the seed
    uint32_t step_size = 2;
**/


    int matchingScore = 2;
    int mismatchingPenalty = -3;
    int openGapPenalty1 = -4;
    int extendGapPenalty1 = -2;

    int openGapPenalty2 = -45;
    int extendGapPenalty2 = 0;

    int zDrop = 50; // this one should be slightly larger than openGapPenalty2
    int32_t bandwidth = 20000; // using very large value, 12000 there are some strange result, should be checked very carefully


    double lambda = 0.382291;
    double kValue = 0.006662;


    bool onlySyntenic = false;
    double pvalues = 0.1;

    std::string input = "/Users/bs674/SORBI_3009G024600_5_prime";
    std::string reference = "SORBI_3009G024600_maize_V4_+_8_135910122_136010121";
    std::string query = "SORBI_3009G024600_sorghum";
    std::string output = "/Users/bs674/SORBI_3009G024600_5_prime_many_to_many_2gaps.sam";
    std::string output2 = "/Users/bs674/SORBI_3009G024600_5_prime_many_to_many_2gaps_o";


    std::string refChr = "8";
    int32_t refStart = 135910122;
    int queStrand = 0; //1 is positive strand, 0 is negative strand

    std::string queChr="9";
    int32_t queStart=2183621;
    int refStrand = 1; //1 is positive strand, 0 is negative strand



    Scorei m(matchingScore, mismatchingPenalty);

    std::map<std::string, std::string> sequences;
    std::vector<std::string> seqNames;
    readFastaFile( input, sequences, seqNames );

    SequenceCharToUInt8 sequenceCharToUInt8;

//    std::string seq1_string ="CGCCGCCGCCGCGTCCTCCTCCTCCTCCTCCTCCAGGCTCTTCACGGCGTGGCTGGTGGCGTCGTGGTACGCCTCCAACATCGGGGTGCTCCTCCTCAACAAGTACCTCCTCTCCGTCTACGGCTTCCGCTTCCCCATCCTGCTCACCGCCTGCCACATGACCGCCTGCACGCTCCTCTCCGCGCTCGTCCACCACCACCACCACCACCGCTCCTCCTCCCGCAGCCGCGGCTCCAGGTCCAGGGCGCAGCTCGCGCGCGTCGCCGTCCTCGGCGCCGTCTTCTGCGCCTCCGTCGTCGCGGGGAACGTCTCGCTCCGCCACCTCCCGGTCTCCTTCAACCAGGCCGTCGGCGCCACCACGCCTTTCTTCACCGCGCTCCTCGCCTACGCCGTTGCTGGTCGCCGCGAGGCCTTCGCCACCTACGCCGCGCTCGTCCCCGTCGTCGCCGGCGTCGTCATCGCCACCGGGGTGAGTCAGTACTCAGTAGCGCCGCGACGAGCGATTCGAG";
//    std::string seq2_string ="CGCCGTCGCCGTCGCCGTCGCCGTCCTCGTCCAGGCTCTACACGGCGTGGCTGGTGGCGTCGTGGTACGCGTCCAACATCGGCGTGCTGCTGCTCAACAAGTACCTGCTCTCCGTGTACGGCTTCCGCTTCCCGCTGCTGCTGACCGCGTGCCACATGTCCGCCTGCGCCGTCCTCTCCACGCTCGCCCAGCACGCCTCGCCCCGCCCCCGCAGCAGCTCCTCCCCGCGCTCCCACCGCCAGCTCGCTCGCGTCGCGCTCCTCGGCGCCGTCTTCTGCGCGTCCGTCGTCGCCGGGAACGTCTCGCTCCGCCACCTCCCCGTCTCCTTCAACCAGGCCGTCGGCGCCACCACGCCTTTCTTCACCGCGCTCCTCGCCTACGCCGTCGCCGCCCGCCGCGAGGCCTGCGCCACCTACGCCGCGCTCGTCCCCGTCGTCGCCGGCGTCGCCATCGCCACCGGGGTGAGGACGCTAGGCTAGTAGCAACGCGGCGGATGGGGATTAGAG";
//
    std::string seq1_string = sequences[reference];
    std::string seq2_string = sequences[query];

    int8_t * seq1 = sequenceCharToUInt8.seq_to_int8(seq1_string);
    int8_t * seq2 = sequenceCharToUInt8.seq_to_int8(seq2_string);
    int8_t * seq1_rev_com = sequenceCharToUInt8.rev_comp(seq1, seq1_string.length());
    int8_t * seq2_rev_com = sequenceCharToUInt8.rev_comp(seq2, seq2_string.length());

    int32_t length1 = seq1_string.length();
    int32_t length2 = seq2_string.length();

    int w = 10;
    int xDrop = 20;

    std::vector<PairedSimilarFragment> pairedSimilarFragments0 =
            findSimilarFragmentsForPairedSequence (seq1, seq1_rev_com, seq2, seq2_rev_com, length1, length2,
                                                   seed_window_size, mini_cns_score,
                                                   matrix_boundary_distance, openGapPenalty1, extendGapPenalty1,
                                                   openGapPenalty2, extendGapPenalty2, matchingScore,
                                                   mismatchingPenalty, m, step_size, seq1_string, seq2_string, pvalues,
                                                   lambda, kValue, zDrop, bandwidth, w, xDrop);

//    std::vector<PairedSimilarFragment> pairedSimilarFragments0 =
//            findSimilarFragmentsForPairedSequence (seq1, seq1_rev_com, seq2, seq2_rev_com, length1, length2,
//                                                   seed_window_size, mini_cns_size, mini_cns_score,
//                                                   matrix_boundary_distance, openGapPenalty1, extendGapPenalty1,
//                                                   openGapPenalty2, extendGapPenalty2, matchingScore,
//                                                   mismatchingPenalty, m, step_size, seq1_string, seq2_string, pvalue,
//                                                   lambda, kValue, zDrop, bandwidth);

    if ( onlySyntenic ){
        std::vector<PairedSimilarFragment> pairedSimilarFragments = syntenic(pairedSimilarFragments0);
        pairedSimilarFragments0 = pairedSimilarFragments;
    }

    std::ofstream ofile;
    ofile.open(output);
    std::ofstream ofile2;
    ofile2.open(output2);
    ofile << "@SQ\tSN:8\tLN:181122637" << std::endl;
    for( int i=0; i<pairedSimilarFragments0.size(); ++i ){

        if( refStrand == 1 && refStrand==queStrand ){
            int32_t quePos = queStart + pairedSimilarFragments0[i].getStart2() - 1;
            int32_t refPos = refStart + pairedSimilarFragments0[i].getStart1() - 1;
            ofile << queChr << "\t0\t" << refChr << "\t" << refPos << "\t" << pairedSimilarFragments0[i].getScore() << "\t" << quePos << "H";
            for( int j=0; j<pairedSimilarFragments0[i].getCigar().size(); ++j ){
                uint32_t cigarLength = pairedSimilarFragments0[i].getCigar()[j]>>4;
                uint32_t cigarType = pairedSimilarFragments0[i].getCigar()[j]&0xf;
                ofile << cigarLength << "MID"[cigarType];
            }
            ofile<< "*\t0\t0\t" << getSubSequence(seq2_string, pairedSimilarFragments0[i].getStart2(), pairedSimilarFragments0[i].getEnd2()) << "\t*" << std::endl;
        }else if ( refStrand == 0 && 1==queStrand){
            int32_t quePos = queStart + pairedSimilarFragments0[i].getStart2() - 1;
            int32_t refPos = refStart + seq1_string.size() -1 - pairedSimilarFragments0[i].getEnd1()+1;
            ofile << queChr << "\t16\t" << refChr << "\t" << refPos << "\t" << pairedSimilarFragments0[i].getScore() << "\t" << quePos << "H";
            for( int j=pairedSimilarFragments0[i].getCigar().size()-1; j>=0; --j ){
                uint32_t cigarLength = pairedSimilarFragments0[i].getCigar()[j]>>4;
                uint32_t cigarType = pairedSimilarFragments0[i].getCigar()[j]&0xf;
                ofile << cigarLength << "MID"[cigarType];
            }
            ofile<< "\t*\t0\t0\t" << getSubSequence(seq2_string, pairedSimilarFragments0[i].getStart2(), pairedSimilarFragments0[i].getEnd2(), 0) << "\t*" << std::endl;
        }else if( refStrand == 1 && 0==queStrand ){
            int32_t quePos = queStart + pairedSimilarFragments0[i].getStart2() - 1;
            int32_t refPos = refStart + pairedSimilarFragments0[i].getStart1() - 1;
            ofile << queChr << "\t0\t" << refChr << "\t" << refPos << "\t" << pairedSimilarFragments0[i].getScore() << "\t" << quePos << "H";
            for( int j=0; j<pairedSimilarFragments0[i].getCigar().size(); ++j ){
                uint32_t cigarLength = pairedSimilarFragments0[i].getCigar()[j]>>4;
                uint32_t cigarType = pairedSimilarFragments0[i].getCigar()[j]&0xf;
                ofile << cigarLength << "MID"[cigarType];
            }
            ofile<< "\t*\t0\t0\t" << getSubSequence(seq2_string, pairedSimilarFragments0[i].getStart2(), pairedSimilarFragments0[i].getEnd2()) << "\t*" << std::endl;
        }


        std::cout << "i " << i << " start1 " << pairedSimilarFragments0[i].getStart1() << " end1 " << pairedSimilarFragments0[i].getEnd1()
                  << " start2 " << pairedSimilarFragments0[i].getStart2() << " end2 " << pairedSimilarFragments0[i].getEnd2() << " score " << pairedSimilarFragments0[i].getScore() << " cigar ";
        ofile2 << pairedSimilarFragments0[i].getStart1() << " " << pairedSimilarFragments0[i].getEnd1()
              << " " << pairedSimilarFragments0[i].getStart2() << " " << pairedSimilarFragments0[i].getEnd2() << " " << pairedSimilarFragments0[i].getScore() << " " << pairedSimilarFragments0[i].getPValue() << std::endl;

        for( int j=0; j<pairedSimilarFragments0[i].getCigar().size(); ++j ){
            uint32_t cigarLength = pairedSimilarFragments0[i].getCigar()[j]>>4;
            uint32_t cigarType = pairedSimilarFragments0[i].getCigar()[j]&0xf;
            std::cout << cigarLength << "MID"[cigarType];
        }
        std::cout.precision(10);
        std::cout << " evalue: " << pairedSimilarFragments0[i].getEValue() << " pvalue: " << pairedSimilarFragments0[i].getPValue() << std::endl;

        std::cout << getSubSequence(seq1_string, pairedSimilarFragments0[i].getStart1(), pairedSimilarFragments0[i].getEnd1()) << std::endl;
        std::cout << getSubSequence(seq2_string, pairedSimilarFragments0[i].getStart2(), pairedSimilarFragments0[i].getEnd2()) << std::endl;
    }
    ofile.close();
    ofile2.close();
}





//weighted + 2 gap penalties
TEST(findSimilarFragmentsForPairedSequence, weightedApproach){

    int32_t seed_window_size = 38;
    int32_t mini_cns_score = 40;
    int32_t step_size = 8;
    int matrix_boundary_distance = 0;

    int matchingScore = 2;
    int mismatchingPenalty = -3;
    int openGapPenalty1 = -4;
    int extendGapPenalty1 = -2;

    int openGapPenalty2 = -45;
    int extendGapPenalty2 = 0;

    int zDrop = 50; // this one should be slightly larger than openGapPenalty2
    int32_t bandwidth = 20000; // using very large value, 12000 there are some strange result, should be checked very carefully

    double lambda = 0.382291;
    double kValue = 0.006662;

    double pvalues = 0.1;

    Score score("/home/bs674/scoreMatrix");

    std::map<std::string, std::string> sequences0;
    std::vector<std::string> seqNames0;
    //readFastaFile( "/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.dna.toplevel.fa", sequences0, seqNames0);
    readFastaFile( "/home/bs674/Zea_mays.B73_RefGen_v4.dna.toplevel.fa", sequences0, seqNames0);


    std::map<std::string, int32_t>  chrSize;
    for( std::map<std::string, std::string>::iterator it = sequences0.begin(); it!=sequences0.end(); ++it ){
        chrSize[it->first] = it->second.size();
    }
    std::map<std::string, int16_t *> weight;
    std::map<std::string, int16_t *> weight_rev;
    //readGffFileWithEveryThing("/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv3.31.gff3", chrSize, weight, weight_rev);
    readGffFileWithEveryThing("/home/bs674/Zea_mays.B73_RefGen_v4.42.gff3", chrSize, weight, weight_rev);


    std::string input = "/home/bs674/SORBI_3009G024600_5_prime";
    std::string reference = "SORBI_3009G024600_maize_V4_+_8_135910122_136010121";
    std::string query = "SORBI_3009G024600_sorghum";
    std::string output = "/home/bs674/SORBI_3009G024600_5_prime_weighted_2gap.sam";
    std::string output2 = "/home/bs674/SORBI_3009G024600_5_prime_weighted_2gap";

    std::string refChr = "8";
    int32_t refStart = 135910122;
    int queStrand = 0; //1 is positive strand, 0 is negative strand

    std::string queChr="9";
    int32_t queStart=2183621;
    int refStrand = 1; //1 is positive strand, 0 is negative strand

    Scorei m(matchingScore, mismatchingPenalty);

    std::map<std::string, std::string> sequences;
    std::vector<std::string> seqNames;
    readFastaFile( input, sequences, seqNames );

    SequenceCharToUInt8 sequenceCharToUInt8;

    std::string seq1_string = sequences[reference];
    std::string seq2_string = sequences[query];

    int8_t * seq1 = sequenceCharToUInt8.seq_to_int8(seq1_string);
    int8_t * seq2 = sequenceCharToUInt8.seq_to_int8(seq2_string);
    int8_t * seq1_rev_com = sequenceCharToUInt8.rev_comp(seq1, seq1_string.length());
    int8_t * seq2_rev_com = sequenceCharToUInt8.rev_comp(seq2, seq2_string.length());

    int32_t length1 = seq1_string.length();
    int32_t length2 = seq2_string.length();


    int w = 10;
    int xDrop = 20;

    std::vector<PairedSimilarFragment> pairedSimilarFragments0 = findSimilarFragmentsForPairedSequence_wighted ( seq1, seq1_rev_com,
            seq2, seq2_rev_com, length1, length2, seed_window_size,
            mini_cns_score, matrix_boundary_distance,
            openGapPenalty1, extendGapPenalty1, openGapPenalty2,
            extendGapPenalty2, matchingScore, mismatchingPenalty,
            m, step_size, seq1_string, seq2_string,
            pvalues, lambda, kValue, zDrop,
            bandwidth, w, xDrop, score,
            weight[refChr]+135910121, weight_rev[refChr] + (chrSize[refChr]-136010121));

    std::ofstream ofile;
    ofile.open(output);
    std::ofstream ofile2;
    ofile2.open(output2);
    ofile << "@SQ\tSN:8\tLN:181122637" << std::endl;
    for( int i=0; i<pairedSimilarFragments0.size(); ++i ){

        if( refStrand == 1 && 0==queStrand ){
            int32_t quePos = queStart + pairedSimilarFragments0[i].getStart2() - 1;
            int32_t refPos = refStart + pairedSimilarFragments0[i].getStart1() - 1;
            ofile << queChr << "\t0\t" << refChr << "\t" << refPos << "\t" << pairedSimilarFragments0[i].getScore() << "\t" << quePos << "H";
            for( int j=0; j<pairedSimilarFragments0[i].getCigar().size(); ++j ){
                uint32_t cigarLength = pairedSimilarFragments0[i].getCigar()[j]>>4;
                uint32_t cigarType = pairedSimilarFragments0[i].getCigar()[j]&0xf;
                ofile << cigarLength << "MID"[cigarType];
            }
            ofile<< "\t*\t0\t0\t" << getSubSequence(seq2_string, pairedSimilarFragments0[i].getStart2(), pairedSimilarFragments0[i].getEnd2()) << "\t*" << std::endl;
        }


        std::cout << "i " << i << " start1 " << pairedSimilarFragments0[i].getStart1() << " end1 " << pairedSimilarFragments0[i].getEnd1()
                  << " start2 " << pairedSimilarFragments0[i].getStart2() << " end2 " << pairedSimilarFragments0[i].getEnd2() << " score " << pairedSimilarFragments0[i].getScore() << " cigar ";
        ofile2 << pairedSimilarFragments0[i].getStart1() << " " << pairedSimilarFragments0[i].getEnd1()
               << " " << pairedSimilarFragments0[i].getStart2() << " " << pairedSimilarFragments0[i].getEnd2() << " " << pairedSimilarFragments0[i].getScore() << " " << pairedSimilarFragments0[i].getPValue() << std::endl;

        for( int j=0; j<pairedSimilarFragments0[i].getCigar().size(); ++j ){
            uint32_t cigarLength = pairedSimilarFragments0[i].getCigar()[j]>>4;
            uint32_t cigarType = pairedSimilarFragments0[i].getCigar()[j]&0xf;
            std::cout << cigarLength << "MID"[cigarType];
        }
        std::cout.precision(10);
        std::cout << " evalue: " << pairedSimilarFragments0[i].getEValue() << " pvalue: " << pairedSimilarFragments0[i].getPValue() << std::endl;

        std::cout << getSubSequence(seq1_string, pairedSimilarFragments0[i].getStart1(), pairedSimilarFragments0[i].getEnd1()) << std::endl;
        std::cout << getSubSequence(seq2_string, pairedSimilarFragments0[i].getStart2(), pairedSimilarFragments0[i].getEnd2()) << std::endl;
    }
    ofile.close();
    ofile2.close();
}





//weighted + 1 gap penalty
TEST(findSimilarFragmentsForPairedSequence, weightedApproach1gap){

    int32_t seed_window_size = 38;
    int32_t mini_cns_score = 40;
    int32_t step_size = 8;
    int matrix_boundary_distance = 0;

    int matchingScore = 2;
    int mismatchingPenalty = -3;
    int openGapPenalty1 = -4;
    int extendGapPenalty1 = -2;

    int openGapPenalty2 = -45;
    int extendGapPenalty2 = 0;

    int zDrop = 50; // this one should be slightly larger than openGapPenalty2
    int32_t bandwidth = 20000; // using very large value, 12000 there are some strange result, should be checked very carefully

    double lambda = 0.382291;
    double kValue = 0.006662;

    double pvalues = 0.1;

    Score score("/Users/bs674/scoreMatrix");

    std::map<std::string, std::string> sequences0;
    std::vector<std::string> seqNames0;
    //readFastaFile( "/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.dna.toplevel.fa", sequences0, seqNames0);
    readFastaFile( "/Users/bs674/Zea_mays.B73_RefGen_v4.dna.toplevel.fa", sequences0, seqNames0);


    std::map<std::string, int32_t>  chrSize;
    for( std::map<std::string, std::string>::iterator it = sequences0.begin(); it!=sequences0.end(); ++it ){
        chrSize[it->first] = it->second.size();
    }
    std::map<std::string, int16_t *> weight;
    std::map<std::string, int16_t *> weight_rev;
    //readGffFileWithEveryThing("/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv3.31.gff3", chrSize, weight, weight_rev);
    readGffFileWithEveryThing("/Users/bs674/Zea_mays.B73_RefGen_v4.42.gff3", chrSize, weight, weight_rev);


    std::string input = "/Users/bs674/SORBI_3009G024600_5_prime";
    std::string reference = "SORBI_3009G024600_maize_V4_+_8_135910122_136010121";
    std::string query = "SORBI_3009G024600_sorghum";
    std::string output = "/Users/bs674/SORBI_3009G024600_5_prime_weighted_1gap.sam";
    std::string output2 = "/Users/bs674/SORBI_3009G024600_5_prime_weighted_1gap";

    std::string refChr = "8";
    int32_t refStart = 135910122;
    int queStrand = 0; //1 is positive strand, 0 is negative strand

    std::string queChr="9";
    int32_t queStart=2183621;
    int refStrand = 1; //1 is positive strand, 0 is negative strand

    Scorei m(matchingScore, mismatchingPenalty);

    std::map<std::string, std::string> sequences;
    std::vector<std::string> seqNames;
    readFastaFile( input, sequences, seqNames );

    SequenceCharToUInt8 sequenceCharToUInt8;

    std::string seq1_string = sequences[reference];
    std::string seq2_string = sequences[query];

    int8_t * seq1 = sequenceCharToUInt8.seq_to_int8(seq1_string);
    int8_t * seq2 = sequenceCharToUInt8.seq_to_int8(seq2_string);
    int8_t * seq1_rev_com = sequenceCharToUInt8.rev_comp(seq1, seq1_string.length());
    int8_t * seq2_rev_com = sequenceCharToUInt8.rev_comp(seq2, seq2_string.length());

    int32_t length1 = seq1_string.length();
    int32_t length2 = seq2_string.length();


    int w = 10;
    int xDrop = 20;

    std::vector<PairedSimilarFragment> pairedSimilarFragments0 =
            findSimilarFragmentsForPairedSequence_wighted_1gap ( seq1, seq1_rev_com,
                        seq2, seq2_rev_com, length1, length2, seed_window_size,
                        mini_cns_score, matrix_boundary_distance,
                        openGapPenalty1, extendGapPenalty1, openGapPenalty2,
                        extendGapPenalty2, matchingScore, mismatchingPenalty,
                        m, step_size, seq1_string, seq2_string,
                        pvalues, lambda, kValue, zDrop,
                        bandwidth, w, xDrop, score,
                        weight[refChr]+135910121, weight_rev[refChr] + (chrSize[refChr]-136010121));

    std::ofstream ofile;
    ofile.open(output);
    std::ofstream ofile2;
    ofile2.open(output2);
    ofile << "@SQ\tSN:8\tLN:181122637" << std::endl;
    for( int i=0; i<pairedSimilarFragments0.size(); ++i ){

        if( refStrand == 1 && 0==queStrand ){
            int32_t quePos = queStart + pairedSimilarFragments0[i].getStart2() - 1;
            int32_t refPos = refStart + pairedSimilarFragments0[i].getStart1() - 1;
            ofile << queChr << "\t0\t" << refChr << "\t" << refPos << "\t" << pairedSimilarFragments0[i].getScore() << "\t" << quePos << "H";
            for( int j=0; j<pairedSimilarFragments0[i].getCigar().size(); ++j ){
                uint32_t cigarLength = pairedSimilarFragments0[i].getCigar()[j]>>4;
                uint32_t cigarType = pairedSimilarFragments0[i].getCigar()[j]&0xf;
                ofile << cigarLength << "MID"[cigarType];
            }
            ofile<< "\t*\t0\t0\t" << getSubSequence(seq2_string, pairedSimilarFragments0[i].getStart2(), pairedSimilarFragments0[i].getEnd2()) << "\t*" << std::endl;
        }


        std::cout << "i " << i << " start1 " << pairedSimilarFragments0[i].getStart1() << " end1 " << pairedSimilarFragments0[i].getEnd1()
                  << " start2 " << pairedSimilarFragments0[i].getStart2() << " end2 " << pairedSimilarFragments0[i].getEnd2() << " score " << pairedSimilarFragments0[i].getScore() << " cigar ";
        ofile2 << pairedSimilarFragments0[i].getStart1() << " " << pairedSimilarFragments0[i].getEnd1()
               << " " << pairedSimilarFragments0[i].getStart2() << " " << pairedSimilarFragments0[i].getEnd2() << " " << pairedSimilarFragments0[i].getScore() << " " << pairedSimilarFragments0[i].getPValue() << std::endl;

        for( int j=0; j<pairedSimilarFragments0[i].getCigar().size(); ++j ){
            uint32_t cigarLength = pairedSimilarFragments0[i].getCigar()[j]>>4;
            uint32_t cigarType = pairedSimilarFragments0[i].getCigar()[j]&0xf;
            std::cout << cigarLength << "MID"[cigarType];
        }
        std::cout.precision(10);
        std::cout << " evalue: " << pairedSimilarFragments0[i].getEValue() << " pvalue: " << pairedSimilarFragments0[i].getPValue() << std::endl;

        std::cout << getSubSequence(seq1_string, pairedSimilarFragments0[i].getStart1(), pairedSimilarFragments0[i].getEnd1()) << std::endl;
        std::cout << getSubSequence(seq2_string, pairedSimilarFragments0[i].getStart2(), pairedSimilarFragments0[i].getEnd2()) << std::endl;
    }
    ofile.close();
    ofile2.close();
}




// cut and link and 1 gap
TEST(findSimilarFragmentsForPairedSequence, cutAndLinkdApproach){ // just to make sure that every line has been analysed

    int32_t seed_window_size = 38;
    int32_t mini_cns_score = 40;
    int32_t step_size = 8;
    int matrix_boundary_distance = 0;

    int matchingScore = 2;
    int mismatchingPenalty = -3;
    int openGapPenalty1 = -4;
    int extendGapPenalty1 = -2;

    int openGapPenalty2 = -45;
    int extendGapPenalty2 = 0;

    int zDrop = 50; // this one should be slightly larger than openGapPenalty2
    int32_t bandwidth = 20000; // using very large value, 12000 there are some strange result, should be checked very carefully

    double lambda = 0.382291;
    double kValue = 0.006662;

    double pvalues = 0.1;

    std::string output = "/home/bs674/SORBI_3009G024600_5_prime_cutAndLinkdApproach_2gaps.sam";
    std::string output2 = "/home/bs674/SORBI_3009G024600_5_prime_cutAndLinkdApproach_2gaps.o";

    std::string refChr = "8";
    int32_t refStart = 135910122;
                     //135917208
    int queStrand = 0; //1 is positive strand, 0 is negative strand

    std::string queChr="9";
    int32_t queStart=2183621;
                   //2266095
    int refStrand = 1; //1 is positive strand, 0 is negative strand

    Scorei m(matchingScore, mismatchingPenalty);

    std::map<std::string, std::string> sequences0;
    std::vector<std::string> seqNames0;
    readFastaFile( "/media/bs674/1_8t/maskGenomeForGenomeAlignment/B73_v4_k31_9.fa", sequences0, seqNames0);

    std::map<std::string, std::string> sequences1;
    std::vector<std::string> seqNames1;
    readFastaFile( "/media/bs674/1_8t/maskGenomeForGenomeAlignment/sorghum_v4_k31_5.fa", sequences1, seqNames1);

    SequenceCharToUInt8 sequenceCharToUInt8;
    std::string seq1_0_string = getSubSequence( sequences0["8"],refStart, refStart+99999, 1);
    std::string seq2_0_string = getSubSequence( sequences1["9"],queStart, queStart+99999, 0);
    std::string seq1_string=seq1_0_string;
    std::string seq2_string=seq2_0_string;

    seq1_string.erase(std::remove(seq1_string.begin(), seq1_string.end(), 'n'), seq1_string.end());
    seq2_string.erase(std::remove(seq2_string.begin(), seq2_string.end(), 'n'), seq2_string.end());

    int8_t * seq1 = sequenceCharToUInt8.seq_to_int8(seq1_string);
    int8_t * seq2 = sequenceCharToUInt8.seq_to_int8(seq2_string);
    int8_t * seq1_rev_com = sequenceCharToUInt8.rev_comp(seq1, seq1_string.length());
    int8_t * seq2_rev_com = sequenceCharToUInt8.rev_comp(seq2, seq2_string.length());

    int32_t length1 = seq1_string.length();
    int32_t length2 = seq2_string.length();

    int w = 10;
    int xDrop = 20;

    std::vector<PairedSimilarFragment> pairedSimilarFragments0 =
            findSimilarFragmentsForPairedSequence (seq1, seq1_rev_com, seq2, seq2_rev_com, length1, length2,
                                                   seed_window_size, mini_cns_score,
                                                   matrix_boundary_distance, openGapPenalty1, extendGapPenalty1,
                                                   openGapPenalty2, extendGapPenalty2, matchingScore,
                                                   mismatchingPenalty, m, step_size, seq1_string, seq2_string, pvalues,
                                                   lambda, kValue, zDrop, bandwidth, w, xDrop);

    std::ofstream ofile;
    ofile.open(output);
    std::ofstream ofile2;
    ofile2.open(output2);


    ofile << "@SQ\tSN:8\tLN:181122637" << std::endl;
    for( int i=0; i<pairedSimilarFragments0.size(); ++i ){
        if( refStrand==1 && 0==queStrand ){

            int32_t quePos = queStart + pairedSimilarFragments0[i].getStart2() - 1;
            int32_t refPos = refStart + pairedSimilarFragments0[i].getStart1() - 1;
            std::cout << "line 651 quePos:" << quePos<< std::endl;
            std::cout << "line 652 refPos:" << refPos<< std::endl;

//            uint32_t start = pairedSimilarFragments0[i].getStart2();
//            uint32_t end = pairedSimilarFragments0[i].getEnd2();
            uint32_t start = 1;
            uint32_t end = 1;
            uint32_t number_of_chr=0;
            uint32_t number_of_n=0;
            int started=0;

            for ( int j=0; j<seq2_0_string.size(); ++j ){
                if ( seq2_0_string[j] == 'n' ){
                    ++number_of_n;
                }else{
                    ++number_of_chr;
                }
                if ( number_of_chr == pairedSimilarFragments0[i].getStart2() && started==0){
                    start += j;
                    quePos += j;
                    started=1;
                }
                if ( number_of_chr == pairedSimilarFragments0[i].getEnd2() ){
                    end += j;
                    break;
                }
            }

            number_of_chr=0;
            number_of_n=0;
            uint32_t ref_start = 1;
            for ( int j=0; j<seq1_0_string.size(); ++j ){
                if ( seq1_0_string[j] != 'n' ){
                    ++number_of_chr;
                }else{
                    ++number_of_n;
                }
                if ( number_of_chr == pairedSimilarFragments0[i].getStart1() ){
                    ref_start += j;
                    refPos += number_of_n;
                    break;
                }
            }
            std::cout << "number_of_chr:" << number_of_chr << " pairedSimilarFragments0[i].getStart1(): " << pairedSimilarFragments0[i].getStart1() << std::endl;
            ofile << queChr << "\t0\t" << refChr << "\t" << refPos << "\t" << pairedSimilarFragments0[i].getScore() << "\t" << quePos << "H";
            std::cout << "line 696 quePos:" << quePos<< std::endl;
            std::cout << "line 697 refPos:" << refPos<< std::endl;
            std::cout << "line 698 start:" << start<< std::endl;
            std::cout << "line 699 end:" << end<< std::endl;

            int que_i = start-1; // the index of reference
            int ref_i = ref_start-1; // the index of query


            std::vector<uint32_t> newCigar;
            for( int j=0; j<pairedSimilarFragments0[i].getCigar().size(); ++j ){
                uint32_t cigarLength = pairedSimilarFragments0[i].getCigar()[j]>>4;
                uint32_t cigarType = pairedSimilarFragments0[i].getCigar()[j]&0xf;
                uint32_t op = 0;
                uint32_t length = 1;
                for ( uint32_t h=0; h<cigarLength; ++h ){
                    op = cigarType;
                    if( op==0 ){
                        while( seq1_0_string[ref_i]=='n' ){
                            op=3;
                            if( newCigar.empty() || op != (newCigar[newCigar.size() - 1]&0xf) ){
                                newCigar.push_back(length << 4 | op);
                            }else{
                                newCigar[newCigar.size()-1] += length<<4;
                            }
                            ++ref_i;
                        }
                        ++ref_i;
                        while( seq2_0_string[que_i]=='n' ){
                            op=6;
                            if( newCigar.empty() || op != (newCigar[newCigar.size() - 1]&0xf) ){
                                newCigar.push_back(length << 4 | op);
                            }else{
                                newCigar[newCigar.size()-1] += length<<4;
                            }
                            ++que_i;
                        }
                        ++que_i;
                    }else if(op==1) {
                        while( seq2_0_string[que_i]=='n' ){
                            op=6;
                            if( newCigar.empty() || op != (newCigar[newCigar.size() - 1]&0xf) ){
                                newCigar.push_back(length << 4 | op);
                            }else{
                                newCigar[newCigar.size()-1] += length<<4;
                            }
                            ++que_i;
                        }
                        ++que_i;
                    }else{
                        while( seq1_0_string[ref_i]=='n' ){
                            op=3;
                            if( newCigar.empty() || op != (newCigar[newCigar.size() - 1]&0xf) ){
                                newCigar.push_back(length << 4 | op);
                            }else{
                                newCigar[newCigar.size()-1] += length<<4;
                            }
                            ++ref_i;
                        }
                        ++ref_i;
                    }
                    op = cigarType;
                    if( newCigar.empty() || op != (newCigar[newCigar.size() - 1]&0xf) ){
                        newCigar.push_back(length << 4 | op);
                    }else{
                        newCigar[newCigar.size()-1] += length<<4;
                    }
                }
            }

            for( int j=0; j<newCigar.size(); ++j ){
                uint32_t cigarLength = newCigar[j]>>4;
                uint32_t cigarType = newCigar[j]&0xf;
                ofile << cigarLength << "MIDDSHI=XB"[cigarType];
            }
            ofile<< "\t*\t0\t0\t" << getSubSequence(seq2_0_string, start, end) << "\t*" << std::endl;

        }


        std::cout << "i " << i << " start1 " << pairedSimilarFragments0[i].getStart1() << " end1 " << pairedSimilarFragments0[i].getEnd1()
                  << " start2 " << pairedSimilarFragments0[i].getStart2() << " end2 " << pairedSimilarFragments0[i].getEnd2() << " score " << pairedSimilarFragments0[i].getScore() << " cigar ";
        ofile2 << pairedSimilarFragments0[i].getStart1() << " " << pairedSimilarFragments0[i].getEnd1()
               << " " << pairedSimilarFragments0[i].getStart2() << " " << pairedSimilarFragments0[i].getEnd2() << " " << pairedSimilarFragments0[i].getScore() << " " << pairedSimilarFragments0[i].getPValue() << std::endl;

        for( int j=0; j<pairedSimilarFragments0[i].getCigar().size(); ++j ){
            uint32_t cigarLength = pairedSimilarFragments0[i].getCigar()[j]>>4;
            uint32_t cigarType = pairedSimilarFragments0[i].getCigar()[j]&0xf;
            std::cout << cigarLength << "MID"[cigarType];
        }

        std::cout.precision(10);
        std::cout << " evalue: " << pairedSimilarFragments0[i].getEValue() << " pvalue: " << pairedSimilarFragments0[i].getPValue() << std::endl;

        std::cout << getSubSequence(seq1_string, pairedSimilarFragments0[i].getStart1(), pairedSimilarFragments0[i].getEnd1()) << std::endl;
        std::cout << getSubSequence(seq2_string, pairedSimilarFragments0[i].getStart2(), pairedSimilarFragments0[i].getEnd2()) << std::endl;
    }
    ofile.close();
    ofile2.close();
}




// cut and link and 2 gaps
TEST(findSimilarFragmentsForPairedSequence, cutAndLinkdApproach2gap){ // just to make sure that every line has been analysed

    int32_t seed_window_size = 38;
    int32_t mini_cns_score = 40;
    int32_t step_size = 8;
    int matrix_boundary_distance = 0;

    int matchingScore = 2;
    int mismatchingPenalty = -3;
    int openGapPenalty1 = -4;
    int extendGapPenalty1 = -2;

    int openGapPenalty2 = -45;
    int extendGapPenalty2 = 0;

    int zDrop = 50; // this one should be slightly larger than openGapPenalty2
    int32_t bandwidth = 20000; // using very large value, 12000 there are some strange result, should be checked very carefully

    double lambda = 0.382291;
    double kValue = 0.006662;

    double pvalues = 0.1;

    std::string output = "/home/bs674/SORBI_3009G024600_5_prime_cutAndLinkdApproach_2gap.sam";
    std::string output2 = "/home/bs674/SORBI_3009G024600_5_prime_cutAndLinkdApproach_2gap.o";

    std::string refChr = "8";
    int32_t refStart = 135910122;
    //135917208
    int queStrand = 0; //1 is positive strand, 0 is negative strand

    std::string queChr="9";
    int32_t queStart=2183621;
    //2266095
    int refStrand = 1; //1 is positive strand, 0 is negative strand

    Scorei m(matchingScore, mismatchingPenalty);

    std::map<std::string, std::string> sequences0;
    std::vector<std::string> seqNames0;
    readFastaFile( "/media/bs674/1_8t/maskGenomeForGenomeAlignment/B73_v4_k31_9.fa", sequences0, seqNames0);

    std::map<std::string, std::string> sequences1;
    std::vector<std::string> seqNames1;
    readFastaFile( "/media/bs674/1_8t/maskGenomeForGenomeAlignment/sorghum_v4_k31_5.fa", sequences1, seqNames1);

    SequenceCharToUInt8 sequenceCharToUInt8;
    std::string seq1_0_string = getSubSequence( sequences0["8"],refStart, refStart+99999, 1);
    std::string seq2_0_string = getSubSequence( sequences1["9"],queStart, queStart+99999, 0);
    std::string seq1_string=seq1_0_string;
    std::string seq2_string=seq2_0_string;

    seq1_string.erase(std::remove(seq1_string.begin(), seq1_string.end(), 'n'), seq1_string.end());
    seq2_string.erase(std::remove(seq2_string.begin(), seq2_string.end(), 'n'), seq2_string.end());

    int8_t * seq1 = sequenceCharToUInt8.seq_to_int8(seq1_string);
    int8_t * seq2 = sequenceCharToUInt8.seq_to_int8(seq2_string);
    int8_t * seq1_rev_com = sequenceCharToUInt8.rev_comp(seq1, seq1_string.length());
    int8_t * seq2_rev_com = sequenceCharToUInt8.rev_comp(seq2, seq2_string.length());

    int32_t length1 = seq1_string.length();
    int32_t length2 = seq2_string.length();

    int w = 10;
    int xDrop = 20;

    std::vector<PairedSimilarFragment> pairedSimilarFragments0 =
            findSimilarFragmentsForPairedSequence (seq1, seq1_rev_com, seq2, seq2_rev_com, length1, length2,
                                                   seed_window_size, mini_cns_score,
                                                   matrix_boundary_distance, openGapPenalty1, extendGapPenalty1,
                                                   openGapPenalty2, extendGapPenalty2, matchingScore,
                                                   mismatchingPenalty, m, step_size, seq1_string, seq2_string, pvalues,
                                                   lambda, kValue, zDrop, bandwidth, w, xDrop);

    std::ofstream ofile;
    ofile.open(output);
    std::ofstream ofile2;
    ofile2.open(output2);

    ofile << "@SQ\tSN:8\tLN:181122637" << std::endl;
    for( int i=0; i<pairedSimilarFragments0.size(); ++i ){
        if( refStrand==1 && 0==queStrand ){

            int32_t quePos = queStart + pairedSimilarFragments0[i].getStart2() - 1;
            int32_t refPos = refStart + pairedSimilarFragments0[i].getStart1() - 1;

//            uint32_t start = pairedSimilarFragments0[i].getStart2();
//            uint32_t end = pairedSimilarFragments0[i].getEnd2();
            int32_t start = 1;
            int32_t end = 1;
            int32_t number_of_chr=0;
            int32_t number_of_n=0;
            int started=0;

            for ( int j=0; j<seq2_0_string.size(); ++j ){
                if ( seq2_0_string[j] == 'n' ){
                    ++number_of_n;
                }else{
                    ++number_of_chr;
                }
                if ( number_of_chr == pairedSimilarFragments0[i].getStart2() && started==0){
                    start += j;
                    quePos += j;
                    started=1;
                }
                if ( number_of_chr == pairedSimilarFragments0[i].getEnd2() ){
                    end += j;
                    break;
                }
            }

            number_of_chr=0;
            number_of_n=0;
            int32_t ref_start = 1;
            for ( int j=0; j<seq1_0_string.size(); ++j ){
                if ( seq1_0_string[j] != 'n' ){
                    ++number_of_chr;
                }else{
                    ++number_of_n;
                }
                if ( number_of_chr == pairedSimilarFragments0[i].getStart1() ){
                    ref_start += j;
                    refPos += number_of_n;
                    break;
                }
            }
            std::cout << "number_of_chr:" << number_of_chr << " pairedSimilarFragments0[i].getStart1(): " << pairedSimilarFragments0[i].getStart1() << std::endl;
            ofile << queChr << "\t0\t" << refChr << "\t" << refPos << "\t" << pairedSimilarFragments0[i].getScore() << "\t" << quePos << "H";

            int que_i = start-1; // the index of reference
            int ref_i = ref_start-1; // the index of query

            std::vector<uint32_t> newCigar;
            for( int j=0; j<pairedSimilarFragments0[i].getCigar().size(); ++j ){
                uint32_t cigarLength = pairedSimilarFragments0[i].getCigar()[j]>>4;
                uint32_t cigarType = pairedSimilarFragments0[i].getCigar()[j]&0xf;
                uint32_t op = 0;
                uint32_t length = 1;
                for ( int32_t h=0; h<cigarLength; ++h ){
                    op = cigarType;
                    if( op==0 ){
                        while( seq1_0_string[ref_i]=='n' ){
                            op=3;
                            if( newCigar.empty() || op != (newCigar[newCigar.size() - 1]&0xf) ){
                                newCigar.push_back(length << 4 | op);
                            }else{
                                newCigar[newCigar.size()-1] += length<<4;
                            }
                            ++ref_i;
                        }
                        ++ref_i;
                        while( seq2_0_string[que_i]=='n' ){
                            op=6;
                            if( newCigar.empty() || op != (newCigar[newCigar.size() - 1]&0xf) ){
                                newCigar.push_back(length << 4 | op);
                            }else{
                                newCigar[newCigar.size()-1] += length<<4;
                            }
                            ++que_i;
                        }
                        ++que_i;
                    }else if(op==1) {
                        while( seq2_0_string[que_i]=='n' ){
                            op=6;
                            if( newCigar.empty() || op != (newCigar[newCigar.size() - 1]&0xf) ){
                                newCigar.push_back(length << 4 | op);
                            }else{
                                newCigar[newCigar.size()-1] += length<<4;
                            }
                            ++que_i;
                        }
                        ++que_i;
                    }else{
                        while( seq1_0_string[ref_i]=='n' ){
                            op=3;
                            if( newCigar.empty() || op != (newCigar[newCigar.size() - 1]&0xf) ){
                                newCigar.push_back(length << 4 | op);
                            }else{
                                newCigar[newCigar.size()-1] += length<<4;
                            }
                            ++ref_i;
                        }
                        ++ref_i;
                    }
                    op = cigarType;
                    if( newCigar.empty() || op != (newCigar[newCigar.size() - 1]&0xf) ){
                        newCigar.push_back(length << 4 | op);
                    }else{
                        newCigar[newCigar.size()-1] += length<<4;
                    }
                }
            }

            for( int j=0; j<newCigar.size(); ++j ){
                uint32_t cigarLength = newCigar[j]>>4;
                uint32_t cigarType = newCigar[j]&0xf;
                ofile << cigarLength << "MIDDSHI=XB"[cigarType];
            }
            ofile<< "\t*\t0\t0\t" << getSubSequence(seq2_0_string, start, end) << "\t*" << std::endl;

        }


        std::cout << "i " << i << " start1 " << pairedSimilarFragments0[i].getStart1() << " end1 " << pairedSimilarFragments0[i].getEnd1()
                  << " start2 " << pairedSimilarFragments0[i].getStart2() << " end2 " << pairedSimilarFragments0[i].getEnd2() << " score " << pairedSimilarFragments0[i].getScore() << " cigar ";
        ofile2 << pairedSimilarFragments0[i].getStart1() << " " << pairedSimilarFragments0[i].getEnd1()
               << " " << pairedSimilarFragments0[i].getStart2() << " " << pairedSimilarFragments0[i].getEnd2() << " " << pairedSimilarFragments0[i].getScore() << " " << pairedSimilarFragments0[i].getPValue() << std::endl;

        for( int j=0; j<pairedSimilarFragments0[i].getCigar().size(); ++j ){
            uint32_t cigarLength = pairedSimilarFragments0[i].getCigar()[j]>>4;
            uint32_t cigarType = pairedSimilarFragments0[i].getCigar()[j]&0xf;
            std::cout << cigarLength << "MIDDSHI=XB"[cigarType];
        }

        std::cout.precision(10);
        std::cout << " evalue: " << pairedSimilarFragments0[i].getEValue() << " pvalue: " << pairedSimilarFragments0[i].getPValue() << std::endl;

        std::cout << getSubSequence(seq1_string, pairedSimilarFragments0[i].getStart1(), pairedSimilarFragments0[i].getEnd1()) << std::endl;
        std::cout << getSubSequence(seq2_string, pairedSimilarFragments0[i].getStart2(), pairedSimilarFragments0[i].getEnd2()) << std::endl;
    }
    ofile.close();
    ofile2.close();
}






// cut and link and 1 gaps
TEST(findSimilarFragmentsForPairedSequence, cutAndLinkdApproach1gap){ // just to make sure that every line has been analysed

    int32_t seed_window_size = 38;
    int32_t mini_cns_score = 40;
    int32_t step_size = 8;
    int matrix_boundary_distance = 0;

    int matchingScore = 2;
    int mismatchingPenalty = -3;
    int openGapPenalty1 = -4;
    int extendGapPenalty1 = -2;

    int openGapPenalty2 = -45;
    int extendGapPenalty2 = 0;

    int zDrop = 50; // this one should be slightly larger than openGapPenalty2
    int32_t bandwidth = 20000; // using very large value, 12000 there are some strange result, should be checked very carefully

    double lambda = 0.382291;
    double kValue = 0.006662;

    double pvalues = 0.1;

    std::string output = "/home/bs674/SORBI_3009G024600_5_prime_cutAndLinkdApproach_1gap.sam";
    std::string output2 = "/home/bs674/SORBI_3009G024600_5_prime_cutAndLinkdApproach_1gap.o";

    std::string refChr = "8";
    int32_t refStart = 135910122;
    //135917208
    int queStrand = 0; //1 is positive strand, 0 is negative strand

    std::string queChr="9";
    int32_t queStart=2183621;
    //2266095
    int refStrand = 1; //1 is positive strand, 0 is negative strand

    Scorei m(matchingScore, mismatchingPenalty);

    std::map<std::string, std::string> sequences0;
    std::vector<std::string> seqNames0;
    readFastaFile( "/media/bs674/1_8t/maskGenomeForGenomeAlignment/B73_v4_k31_9.fa", sequences0, seqNames0);

    std::map<std::string, std::string> sequences1;
    std::vector<std::string> seqNames1;
    readFastaFile( "/media/bs674/1_8t/maskGenomeForGenomeAlignment/sorghum_v4_k31_5.fa", sequences1, seqNames1);

    SequenceCharToUInt8 sequenceCharToUInt8;
    std::string seq1_0_string = getSubSequence( sequences0["8"],refStart, refStart+99999, 1);
    std::string seq2_0_string = getSubSequence( sequences1["9"],queStart, queStart+99999, 0);
    std::string seq1_string=seq1_0_string;
    std::string seq2_string=seq2_0_string;

    seq1_string.erase(std::remove(seq1_string.begin(), seq1_string.end(), 'n'), seq1_string.end());
    seq2_string.erase(std::remove(seq2_string.begin(), seq2_string.end(), 'n'), seq2_string.end());

    int8_t * seq1 = sequenceCharToUInt8.seq_to_int8(seq1_string);
    int8_t * seq2 = sequenceCharToUInt8.seq_to_int8(seq2_string);
    int8_t * seq1_rev_com = sequenceCharToUInt8.rev_comp(seq1, seq1_string.length());
    int8_t * seq2_rev_com = sequenceCharToUInt8.rev_comp(seq2, seq2_string.length());

    int32_t length1 = seq1_string.length();
    int32_t length2 = seq2_string.length();

    int w = 10;
    int xDrop = 20;

    std::vector<PairedSimilarFragment> pairedSimilarFragments0 =
            findSimilarFragmentsForPairedSequence ( seq1, seq1_rev_com,
                                                    seq2, seq2_rev_com,  length1, length2, seed_window_size,
                                                    mini_cns_score, matrix_boundary_distance, openGapPenalty1, extendGapPenalty1, matchingScore,
                                                    mismatchingPenalty, m, step_size, seq1_string, seq2_string, pvalues,
                                                    lambda, kValue, w, xDrop);

    std::ofstream ofile;
    ofile.open(output);
    std::ofstream ofile2;
    ofile2.open(output2);



    ofile << "@SQ\tSN:8\tLN:181122637" << std::endl;
    for( int i=0; i<pairedSimilarFragments0.size(); ++i ){
        if( refStrand==1 && 0==queStrand ){

            int32_t quePos = queStart + pairedSimilarFragments0[i].getStart2() - 1;
            int32_t refPos = refStart + pairedSimilarFragments0[i].getStart1() - 1;

//            uint32_t start = pairedSimilarFragments0[i].getStart2();
//            uint32_t end = pairedSimilarFragments0[i].getEnd2();
            uint32_t start = 1;
            uint32_t end = 1;
            uint32_t number_of_chr=0;
            uint32_t number_of_n=0;
            int started=0;

            for ( int j=0; j<seq2_0_string.size(); ++j ){
                if ( seq2_0_string[j] == 'n' ){
                    ++number_of_n;
                }else{
                    ++number_of_chr;
                }
                if ( number_of_chr == pairedSimilarFragments0[i].getStart2() && started==0){
                    start += j;
                    quePos += j;
                    started=1;
                }
                if ( number_of_chr == pairedSimilarFragments0[i].getEnd2() ){
                    end += j;
                    break;
                }
            }

            number_of_chr=0;
            number_of_n=0;
            uint32_t ref_start = 1;
            for ( int j=0; j<seq1_0_string.size(); ++j ){
                if ( seq1_0_string[j] != 'n' ){
                    ++number_of_chr;
                }else{
                    ++number_of_n;
                }
                if ( number_of_chr == pairedSimilarFragments0[i].getStart1() ){
                    ref_start += j;
                    refPos += number_of_n;
                    break;
                }
            }
            std::cout << "number_of_chr:" << number_of_chr << " pairedSimilarFragments0[i].getStart1(): " << pairedSimilarFragments0[i].getStart1() << std::endl;
            ofile << queChr << "\t0\t" << refChr << "\t" << refPos << "\t" << pairedSimilarFragments0[i].getScore() << "\t" << quePos << "H";

            int que_i = start-1; // the index of reference
            int ref_i = ref_start-1; // the index of query

            std::vector<uint32_t> newCigar;
            for( int j=0; j<pairedSimilarFragments0[i].getCigar().size(); ++j ){
                uint32_t cigarLength = pairedSimilarFragments0[i].getCigar()[j]>>4;
                uint32_t cigarType = pairedSimilarFragments0[i].getCigar()[j]&0xf;
                uint32_t op = 0;
                uint32_t length = 1;
                for ( uint32_t h=0; h<cigarLength; ++h ){
                    op = cigarType;
                    if( op==0 ){
                        while( seq1_0_string[ref_i]=='n' ){
                            op=3;
                            if( newCigar.empty() || op != (newCigar[newCigar.size() - 1]&0xf) ){
                                newCigar.push_back(length << 4 | op);
                            }else{
                                newCigar[newCigar.size()-1] += length<<4;
                            }
                            ++ref_i;
                        }
                        ++ref_i;
                        while( seq2_0_string[que_i]=='n' ){
                            op=6;
                            if( newCigar.empty() || op != (newCigar[newCigar.size() - 1]&0xf) ){
                                newCigar.push_back(length << 4 | op);
                            }else{
                                newCigar[newCigar.size()-1] += length<<4;
                            }
                            ++que_i;
                        }
                        ++que_i;
                    }else if(op==1) {
                        while( seq2_0_string[que_i]=='n' ){
                            op=6;
                            if( newCigar.empty() || op != (newCigar[newCigar.size() - 1]&0xf) ){
                                newCigar.push_back(length << 4 | op);
                            }else{
                                newCigar[newCigar.size()-1] += length<<4;
                            }
                            ++que_i;
                        }
                        ++que_i;
                    }else{
                        while( seq1_0_string[ref_i]=='n' ){
                            op=3;
                            if( newCigar.empty() || op != (newCigar[newCigar.size() - 1]&0xf) ){
                                newCigar.push_back(length << 4 | op);
                            }else{
                                newCigar[newCigar.size()-1] += length<<4;
                            }
                            ++ref_i;
                        }
                        ++ref_i;
                    }
                    op = cigarType;
                    if( newCigar.empty() || op != (newCigar[newCigar.size() - 1]&0xf) ){
                        newCigar.push_back(length << 4 | op);
                    }else{
                        newCigar[newCigar.size()-1] += length<<4;
                    }
                }
            }

            for( int j=0; j<newCigar.size(); ++j ){
                uint32_t cigarLength = newCigar[j]>>4;
                uint32_t cigarType = newCigar[j]&0xf;
                ofile << cigarLength << "MIDDSHI=XB"[cigarType];
            }
            ofile<< "\t*\t0\t0\t" << getSubSequence(seq2_0_string, start, end) << "\t*" << std::endl;

        }


        std::cout << "i " << i << " start1 " << pairedSimilarFragments0[i].getStart1() << " end1 " << pairedSimilarFragments0[i].getEnd1()
                  << " start2 " << pairedSimilarFragments0[i].getStart2() << " end2 " << pairedSimilarFragments0[i].getEnd2() << " score " << pairedSimilarFragments0[i].getScore() << " cigar ";
        ofile2 << pairedSimilarFragments0[i].getStart1() << " " << pairedSimilarFragments0[i].getEnd1()
               << " " << pairedSimilarFragments0[i].getStart2() << " " << pairedSimilarFragments0[i].getEnd2() << " " << pairedSimilarFragments0[i].getScore() << " " << pairedSimilarFragments0[i].getPValue() << std::endl;

        for( int j=0; j<pairedSimilarFragments0[i].getCigar().size(); ++j ){
            uint32_t cigarLength = pairedSimilarFragments0[i].getCigar()[j]>>4;
            uint32_t cigarType = pairedSimilarFragments0[i].getCigar()[j]&0xf;
            std::cout << cigarLength << "MIDDSHI=XB"[cigarType];
        }

        std::cout.precision(10);
        std::cout << " evalue: " << pairedSimilarFragments0[i].getEValue() << " pvalue: " << pairedSimilarFragments0[i].getPValue() << std::endl;

        std::cout << getSubSequence(seq1_string, pairedSimilarFragments0[i].getStart1(), pairedSimilarFragments0[i].getEnd1()) << std::endl;
        std::cout << getSubSequence(seq2_string, pairedSimilarFragments0[i].getStart2(), pairedSimilarFragments0[i].getEnd2()) << std::endl;
    }
    ofile.close();
    ofile2.close();
}







// cut and link and 1 gaps
TEST(findSimilarFragmentsForPairedSequence, cutAndLinkdApproach1gap2222){ // just to make sure that every line has been analysed

    int32_t seed_window_size = 38;
    int32_t mini_cns_score = 40;
    int32_t step_size = 8;
    int matrix_boundary_distance = 0;

    int matchingScore = 2;
    int mismatchingPenalty = -3;
    int openGapPenalty1 = -4;
    int extendGapPenalty1 = -2;

    int openGapPenalty2 = -45;
    int extendGapPenalty2 = 0;

    int zDrop = 50; // this one should be slightly larger than openGapPenalty2
    int32_t bandwidth = 20000; // using very large value, 12000 there are some strange result, should be checked very carefully

    double lambda = 0.382291;
    double kValue = 0.006662;

    double pvalues = 0.1;


    std::string refChr = "8";
    int32_t refStart = 1;
    //135917208
    int queStrand = 0; //1 is positive strand, 0 is negative strand

    std::string queChr="3";
    int32_t queStart=13728310;
    //2266095
    int refStrand = 0; //1 is positive strand, 0 is negative strand

    Scorei m(matchingScore, mismatchingPenalty);

    std::map<std::string, std::string> sequences0;
    std::vector<std::string> seqNames0;
    readFastaFile( "/media/bs674/1_8t/maskGenomeForGenomeAlignment/B73_v4_k31_9.fa", sequences0, seqNames0);

    std::map<std::string, std::string> sequences1;
    std::vector<std::string> seqNames1;
    readFastaFile( "/media/bs674/1_8t/maskGenomeForGenomeAlignment/sorghum_v4_k31_5.fa", sequences1, seqNames1);

    SequenceCharToUInt8 sequenceCharToUInt8;
    std::string seq1_0_string = getSubSequence( sequences0["8"],refStart, 77950, 1);
    std::string seq2_0_string = getSubSequence( sequences1["3"],queStart, 13828309, 0);
    std::string seq1_string=seq1_0_string;
    std::string seq2_string=seq2_0_string;

    seq1_string.erase(std::remove(seq1_string.begin(), seq1_string.end(), 'n'), seq1_string.end());
    seq2_string.erase(std::remove(seq2_string.begin(), seq2_string.end(), 'n'), seq2_string.end());

    int8_t * seq1 = sequenceCharToUInt8.seq_to_int8(seq1_string);
    int8_t * seq2 = sequenceCharToUInt8.seq_to_int8(seq2_string);
    int8_t * seq1_rev_com = sequenceCharToUInt8.rev_comp(seq1, seq1_string.length());
    int8_t * seq2_rev_com = sequenceCharToUInt8.rev_comp(seq2, seq2_string.length());

    int32_t length1 = seq1_string.length();
    int32_t length2 = seq2_string.length();

    int w = 10;
    int xDrop = 20;

    std::vector<PairedSimilarFragment> pairedSimilarFragments0 =
            findSimilarFragmentsForPairedSequence ( seq1, seq1_rev_com,
                                                    seq2, seq2_rev_com,  length1, length2, seed_window_size,
                                                    mini_cns_score, matrix_boundary_distance, openGapPenalty1, extendGapPenalty1, matchingScore,
                                                    mismatchingPenalty, m, step_size, seq1_string, seq2_string, pvalues,
                                                    lambda, kValue, w, xDrop);

    for( int i=0; i<pairedSimilarFragments0.size(); ++i ){


        std::cout << "i " << i << " start1 " << pairedSimilarFragments0[i].getStart1() << " end1 " << pairedSimilarFragments0[i].getEnd1()
                  << " start2 " << pairedSimilarFragments0[i].getStart2() << " end2 " << pairedSimilarFragments0[i].getEnd2() << " score " << pairedSimilarFragments0[i].getScore() << " cigar ";

        for( int j=0; j<pairedSimilarFragments0[i].getCigar().size(); ++j ){
            uint32_t cigarLength = pairedSimilarFragments0[i].getCigar()[j]>>4;
            uint32_t cigarType = pairedSimilarFragments0[i].getCigar()[j]&0xf;
            std::cout << cigarLength << "MIDDSHI=XB"[cigarType];
        }

        std::cout.precision(10);
        std::cout << " evalue: " << pairedSimilarFragments0[i].getEValue() << " pvalue: " << pairedSimilarFragments0[i].getPValue() << std::endl;

        std::cout << getSubSequence(seq1_string, pairedSimilarFragments0[i].getStart1(), pairedSimilarFragments0[i].getEnd1()) << std::endl;
        std::cout << getSubSequence(seq2_string, pairedSimilarFragments0[i].getStart2(), pairedSimilarFragments0[i].getEnd2()) << std::endl;
    }
}






TEST(calculateK, c1){ // just to make sure that every line has been analysed

    std::string seq1 = "ACACTATTTTTTAAGGGTCTAATTTTATCTTATATAAGGGAAAATTTCGATAATCCTATGAGGCTGTTTAGGGCAGTTTGAAGGGATGATCTCTTGCATTTAAGGTAAATTCTTGAAGGGATGAGAATAATTTTGATGAGAGGTTTTGCAAAACCTGAATCTGGTGGATCCTTGATTACTTTGTTTTTTTTTATGTTGTACTCTCTTATATTTGCTCACTCCCCGTGCCAATCTTGGTCATCCCATCATTTCCATTCTACTCGTGACATATCTCCTTCTACCCATAATTCAACTATCAGTCTACAACACAACCACCTTTTTCAAGGATAGAACTTACCAAATGCTCAGAGAACTGAGTTGTGGTCTCATTGTTATGACACTATGATCAGTTTATAATTGTCAAGGACTGCATTTTTGAAATGTCAAAACCATGGAGTCAAAAATATAGGGCCAGACTAGATTCCCACCCCTAGCCCGTGCAACTTCTCGCTCTCCAGCTAGAGGGCAGAAGGTGGGATTGCCTCTTTGTGTTTTGTTTTGCTCTTGTGGACTTGAATTTTGTAAGAGTATGTGTATGTTTTGCTCTTGTGGACTTGAATTTTGTAAGAGTACGTGTATGCGGTTCCATTTAGAAGTGTAGGACATCTAAAAACCCAATGATTAGGGGAATTTTGAAAAAAAAAACTCAAATACCTAGAATCCATGACTAAATGGTAAAGATCAAAATGATTAAAGATTTTTATTGTTCCCGCACCTTTGAATAAAAGATTACACGGACAATTTTATTTTCACATCTGTATTTTAGGATTAGGATTCTATGTGCTATCCAAAAATTTATTTCTAACTTCTACACATTCTAGATTTCCTATAGGATTTCAATTGCACAACATTCTAATTTTTTTTTTCAAATGTCGTAACACCGTAACACGATGAACTGTCACCCGGGAGAGGGGGTGGGGTTGGCTTTTGTTCAGGGATAACGCTACTGAAACTGATACCCCGGTGATAACGGTGCTTTGTCATGATATCGAACAATGCTACGAAAACTGACACCTCTATCCTGGTGATAACGGTAATTTGTCATGATATCGGTGGATGAAATGGACATTAGACTTCACCACATCATCGATCTGTGGGATGAGATCACCGCCTGTTAGCGATTGATTATTTAACACAAAAACTCGGGCACAAAGCAGACTCCTTTGAAATTCCTCTTTTCCCTACACATTTTCAACATTAAAACTCTCCAAAAACAAGAACCGTTAAGGCCTTGTTTAGTTAAAAATTTTTTGCAAAACAGGTACTGTAGCACTTTCGTTTGTATTTGACAAATATTATCCAATCATACACTAACTAGGCTCAAAAGATTCGTCTCGTCAATTTCGACCAAATTGTGTAATTAGTTTTTTTTCGTCTATATTTAATATTTCATGCATGCGTCTAAAGATTCGATGTGACGGAGAATGTAAAAAATTTTGCAAAATTTTCTGTGAAGTAAACAAGGCCTAAATGGAATTCCTTTTCTTTTGGTAACTTCTTGGGATTCCGCTTCCTTATTTTTAGAAACTTGTACAGCTCTGAAAGCCAGGGAGAAAAAACAGCTGTTCTGAGCAGATGTAGGGCTGGACGACCAACACATATATAGAAACAAGAGGCTTTTGCGCAGGAAATGCCGGATGCATTTCTTCAGGGCCGGCCAATGCTTTTTGTGCCCACAGAAGACGCCATCTCTACTCCAGTCCAGCTAGTACACACACTGGCCATGCATGTGGACCGTGGATTCCTCAACTTCCTCACTTACTTGGGAATTGGGACAATGTGTGATTATGACGGTGATGGCCGGTGGCTGTTTCAAACCGTAAAATAAAAACCTTGCAACAACTTGGATTAGTTATTCAAGTTGCTTCACCTGGATGGAAGATGCAAAGTCAGCATGCTATCATATGGTGTGGACTGACAGGAGGATTAGAATACTAGTGTCCAGGTAAATTTCAGGTTCATTCTGTCTTGACTAGTTGGCGCCAAACTATTCATGGAGGAAGGAGTAGGGAAATTCCTCAGAGGTATTGAGTGCAGTGCACACATTATGGGCGTCAAAAAAGATGACCTAGCACATGATACATCCACGTCAGAGAGTGGGCCCGTTGATCCATGGACACTGATGAAGGCGAAAATATCCTTCATCTTTTCATCTTATGAACAAGTGTCAGTGTCAGTGTGTCACACAGGTGAGAGGCACATGACCTTCAGGTTGTAAATAAAATACTGACAGCCAACAAAAAAGAACTATAACATGCACTTTGAGTTTTGTAATTTACAGTGGAGACAACATAGTTGATGAAAGTTTTGATGGAAAAAACAACCATATATCATGTCATGTAGATCTTTCATAGACTGACAGGACTTATGCCAGTAGATGCTCTGAATCTGGATGATACTCCGAGTAGATTGAAGTCAAATTTATCATAAAATAAAACGAAATCCAACATGCATCCTCCCTGTATATGAACCGTGTATATGGATCTGAGCGGCGTGCCGACATGCACGGAGACGCTCCTACGCCATCCTCTCTGTATGACGTGTCACGCCAGATTGACACGGCTCTGGATTGTGGATAAATTAGGCCAACCAATCTCAGGCTCGTTCGTTCCTCCTTCCACCTGGATTCTCCCAGGAGATGAAATCCACGATAAATTATCAGAAGAAGATTCCAGGTTATAATTAATTCAGGAGTCCAATCCGGATTAACCGAACATGCCCCGGCTTTTATGATCCAACTATCAATACGCTCATATTTTGAAAATAATAGGAAAGAAGCACACTTGAAGAGCGCGGTACAACAGAGAGTTGTATGCTGCGTTCGGGAAGGATGAATCACTCCCGAAAAGGAGTCTATTGGAGAGAAGCACAGAGGTCAACCCGCTTCAAATATGGAACATGAATTCTGGCAATGCAACGGAGTCCAAAAACCCAAGGTTTCCTCCGCAAGGTTCGTCCACGAAGGGTGAGTCAGGGTCTAAGATCAGGCCGAAAGGCGCGAGTTCAGGGCAAACCCTAGCATTGTTTGGCGCCAAAGGAAACTAGAAATTAGTTGATTTAATTGAAGGAAGTGAGCCCTTAGTACAGACAGACTGATTGGGTGCTAAGATCCAAAGTCGAGAGGGAAGCAACTCTAGTAGAGCGAGACGCGAAGAGTGACAAGTGGTCCGGCCGGGTCAAGTGCTAGAGCTTGTGGTAAGCACTCCACGTGGGAGAGTGTGACTTGTGGGTCGATGAATTATCGGAGGAGGATCGTTTAACTGTAGCAAGAGCACGAAAAATTGAGCGCCGAAGGCGTCCTTGGGGTGATCTCGTAGTTTCTATGGGGTGGAGACAATGGGGTCGGTCCATGGATTTTCCTTCCTTTTGCCACATTTCGCTTCCTTTTGCCACATTTCGCTCAAAGGGTTGAAGGGAGATAGTGCATCAAGCTATTCGCAAGGGCCAACTTTCCTCTTCCCCAGGGATCCCAGATGAGGGAAGCCTAGGAGAGCCGCCAGCTTCCGCGATATACGACGATTCGCTCGTTATATGAGGCGCACATGAAACTTTTATTACCTCTCTTCATATCATGTCAGTATATTCAATATTCATAAAACTCTCATGAAAACTTGGACTTATGAATTAACTGCTTTACTTATTGGAATATAAATGTTGACTATTTTCTATAAATTTTGTTAAATTTAAAATTGTTTGGCTCCTCGTGAAATTTCCATTCTCTTATTCATGTAACTCCAACTCCATTTTGTGCTGAACAATTGGAATCAAATAAAATGGATTGAAGTGCCTGAGGAAATTGCAATTCTATGGTCATTTGTCTCATTGAGCAATAATGATAGGCACTGGCTACTTGCCTACTTGCATTGTGAGCTGCAGTGTCATCATGTCGTTCATGGGCACAAGTTTAGCTACTACCAGTACCAGTGCATGATGCACTAGTCGTCATCTTCCTTTTATTTCAATGAATAAAAAAACCATCGTCTCGATCTTACAGGCTTGCCTAGAGCATGTGAATCGCGCATAGATAAGACTGACAGGTACATGGCATCCTTCAGGGTTAAAATATGAGAAAACCCAAAAGATGGAAATAACGGGTCCTACATCATTTAAATTTAATTAGAAAGGGAGACAACAAAGTTCATGAAGAAAACAATAGGACTTTCTGATACTTATAAATAAATTGCTGAAACTTACTTGTCATACCGAGCTAAAGCGTTGATGGATACAATATATATATGTCTTGTTTATATCCAAAAAGTTTTTGGATTTTGACAATATAGCACTTTCGTTTTTATTTGACAAACATTGTCCAATCATGGAGTAACTAGATCTAAAAGATTCGTCTCGCGATTTACAGGCAAACTAGGCTATTAGTTTTTGTTTTCATCTATATTTAATGCTCCATGCATGTGCCACAAGATTTGATGTGACCAGGAATCTTGAATTTTTTTTTGATTTTTAGGGTGAACTAAACAAGGCCTACTAGTTACATGCTCTTATCAATTTCATATGAAACATATAATCTCCCATCTAAAAAAATCTCATTTCTTGAGGAGTCTTTAAATATAACCAAATATATAATAAATATCACTAGATCAACATTAGAATACATTTTTATAACAAATACATTTGGAGATACAAATATTGATATTATATATTATAAATCTAGTCAAATTTAAAAAAATAATTCTAAAAAAAAGAGATGTGTATTCTTTTTGTGGCAGAGGAGTACATATATTTAAAGAGACAACGAATTGAAGGTAGCCTAATTGGTAAGCTTCCTTATAAATTATTTCGATCGGTGATAGTCAAGATGTGTCTTGGTGAACATTGAGGTGAGTGGCTATGGTGGTAGCTTGGGGGCGGGGGGGGGGGGGGGGGAGGGGAGGGGAGGGGAGGGAGAAATGCAAGTCGCAATTTGTATGACAGAATATATACAATGGGCAGGTCCCAAGATGAAGTATGCATGGTATCATGAAAACCAACCTTTATTTGGGATATCTGTTAGCAAAGAAAGCCGAACTACAGTCTTCCATTTATTAATTTAAACATTGCTGCAGAGATTAGAGATTTCACTAAAATTGAGAGAGAATATCTTAAGTCAATCTCAGGATCACCCGACTAAGGCCTTGTTTAGTTCACCCTAAAAACCAAAAACTTTTCAAGATTCTCCATCATATCGAATCTTGCGGCACATGCATGAAGCACTAAATATAAACGAAAACAAAAACTAATTGCACAGTTTACCTGTAAATCGCGAGATGAATCTTTTGAATTCCATGATTGGACAATATTTGTCAAATAAAAACGAAAGTGTCACACGAAACTAGTTCTTAGGTGAGCACTGCTAAATAAAAAGAATCTTAAGTAGGTGACTACTACTGCCCACAAGTCCACAACCACAAGCAGTGACAGTCTGACAGATTATGTATTAGTGCTTCACGCAGTAGTGTGTGTTATGATCACCAACTGAATAGTATAATCAGAGAGTAGATAGCAGTGTACGTACTCAATGAACTAGCTACTATAGTGCCTAGCTAGTTGTGGCAGCTGATATATAGGACATACTAGGAGTCACCAACAATGCATGATTGCAAGTGACAGGAGGATTGCTACTAATTACTAATGATAATCCAAATAAATCCAAACGACGACTAACATGCGGGTGGCTGTCAATTCATTACTCCCTAGCTCGATCGGCCGTGTGTGTTGCTGTCAGCCGTCTAACCAACCCTGGCGACCAACGGAGGTTGCGTGCGTGGTGTGAGGTGTCGCCGGTAGTAGACACGACACGGCGAGGCATCGATCACCGGCCGGCCGGCGACGTACGTAGCAGCACGTGTCTGCCACCCGACAGCGGCACGTAGGATTATTTTGACTACTGGGATCTGGATGAAGGCAATTCTGGACTGGCATATATATTCCGGTCACGCATACGTGTATAACAAAGGAAAGCAGTAGCTAGCGCTTATCCGTTCCCTGCGTGTGCATGGTTGCAGGCTGCTATATCGATCCGTTCCAACACCAATGTTTTAAATCATCACGAACTGTAACATTTAGCAGCCTCCTCAAACAGTTCATCTTTTGTTAAATCAGACTATAACGGTTAAATTGTACTCCATATTAAAGAAAATAACTTACTTCTTTCATCCTTAAAAGAATGTAAGGGCATGTTTAGCATAGCTTCGGCATCAATTTTTTTGCAGAGAAATCATAGCGAACATTTTAATGGAGTGGGAGTTGTTTTTCCAGTTTCAAGTCGTCGTAGCTGCAAAGAGTGGTTGTTCCTCTTCGAGTACTCTTTCCTAAATCTCATTTCTCGACAAGTCCAACAATTATAAGTTTGACTCTAGGTTTATAGAATAAATATATGGTACTAATAACATTTATATCTTTCAATAAATTAGTTATAAAATATATTTCTTTTTTATATTTAATCTATAGTGATATTTATTATGCATCATAAATATTATTTCCTCCATCCTAAATTATAAGACATTCGTTTTTTTACTCTAAATTTGACCACTCGTCTTATTTTTGTACAAATATAGCAAAATTTAAGTAATTTTTGAAGAACTTTTATTAATAAATCAAGTCACGATAAAATAAGTGATATTTTGTATAATTTTTTTCTTAATAAGACGAGTGATTAAACTTGGTGTCAAAAGTCAAACATTTTATAGTTTGCGATGGATTGAGTAGTATTTTTTTTATTGTATATATCTGGTTAAATATATATTATTTGACTCTTAGAGAAGTGAGACTTCTTTAAAAAAACCGGAGGGAGTAGATCTTTTGTTCAAGCAGCTATATTATTTGGAGATACCACATGAGCTTTAAAATCTCAGCTCCGACATAGGATTGACCTTGAATAGATCGCATGGCATCGGTACAACCTAGTAGGAGTACGTATTCATGATTCCAAGTGCTCGAAACAAGCGAGAGAGCAGTAATGAATCTGCAGTCAAACAAGGAATCACACTCGCGCACGTACAACTAACCAACTTCGTACGTACAAAAAAAATACTCCTTTCATCCAAAAAGTGAAAATAAAAAACCGATGAAAAAAAAACAGACGGACAAAAGGGAGAAAATGGGGAATGGAAATAAAAAACAGATGGAAAAAGGATCGGAAAAAAAAAGGAAAATGAAAATAAAAAACAGATTGGAAAAAGGGAGAAAATTTTGAAATGAAAGTAAAAAACAGATGGAAAAAAGGAAAATGAAAATAAAAGGCATGGAGAAAAGGGAGTAAATGGGAAATGAAAGATAAAAACTGAAGGAAAATGAAAAATTAAAAACAGATGAAAAAAGGAAATGATAAAAAAGCGCTAGAAGGAGGGCTTGAACCTCCGACCTTGTGGTTAACAGCCACACGCTCTAACCAACTGAGCTATTCCAGCTGCTGTTCTTGTAGAGTATACTTTAATATTTGTCAGAAGCATGAAAGTGTCTACATGTTGAAAGGTTTGCATTCAAATATTTTTGTGAACTTGTCGAGTGAAACCATGGTTTTGAGGACACTTTGGTACATCTTATAAGGTGCCTATGGCCAGAATAAGTTGTAACGAGAAGCTATACTGATCATATAGCACAACAGAGAGTCTCAGAAATTAGAAGTTTTACATAGATTCTCGTAAGCAGACCTTGTTCACGTTAGGCACATATTTCTATAGAATTAACCAACGGAGAAAACAAATATGACTATTGAATGGCAACAGTTGCACCCAGCCATGGTACTCTTGCTAGAAGCTATATACTGATCATAGAACAACCACTCAATCAGTCACCTTTCTTCAACATTGTTAAGGAACACAGTGCTCGGTGATGCAAGGCTAGCTCGTACATGATTACTCTGCAGAAAACTCTATTAAAAACTTTAACTGTTGCTATGCAGCCTTAAAAGAAAATCAGCCTCAGTACTACGGGCATGCCATGAACATGATACAAAACATACAAGACAACAAGCATATTCACGCTCTATAGTTACACCTCATGTCATCATTTTTGTGAGAGAACATCATGTCAATGCCATCATCTGTTTGGTGCAGCCGAGGCAGATTTTTATCAGTCATGATGAATGCTAGATGATGATGTGGAAACGCAGAACTTGCAGCTTCGCTGAATCAGGTGCTTTCGTTCTGCACTTGGAACCATGCGAGCCCCGGTTGACCACACTGCCTCGACTCTGGGATAGGAACTATTGGCCCCAAAAATATGCGCCATAGGTCGAGTATTATCGTGACGAGCATCGCTGCGCCCAGAGAGTCAACAGCTACCCTGGATAGGTTCCACTTCTCATCTTTTTCATCAATCCTGAACCACAAGTAACAGAGATACACTAAGTTGCTGAAAAGAGTTATAATAAGTCCTCATTTTCTGTATAAAACAGCACCAAATCAAATCATGGTCATGACATGTCTGCCAACATCAATCACGGGAAGTACTATTGGTCTTCAGACAAGTGTTATTTGTCTCAATAATCAGAGTGCCATTTATTCCATTATTTCATCAGGTGATGAGTGGGAACAGGGAAGATATAAGCTTGCTACATAAGATTTTTTACTTAATCTCAGACTGAAGGAAAAAATAAGTTACCAAGTTTACCACAAAAGCAAAAAATATTCACGTAGACAAAATACGAAACAACAAAATTCACAAAGCTCAATGAGCCAGACACTCCTACAAACCCAACTGAAATGCTGAGCAAGTTCATTCTATTTCGAGATTTGAGACTGGAGAGATGTTATACCTGCTCAATATATATATCTATTTAGTAAATAGATCTTTTGTGTAAAATTAAGACTGGTAAACATAGTCTGATTTAGTCATTGAAGGCAAACAGACAGACGGATAAAGCACAAAAAGAAATAGTAGAAGCATCGCTGTTCTAAAGGAGTTCATAGGAATGAAGAGCAAGTGAGTTCAAGTACAAGAAACCAAAATACTAGAATATCATTATGGAACTAGTTGATCAGTGGTCACATCTAAGCTCAGATTCATTTTACCTTGAAAACATTGGGAAGCTGACGATGAAGTATATTGCGTAGAACAAGGAACCTACTTTGTACATTATATCCCGGTCAACAAATTCATAGTATGGGAACTGCCAAAGTCCACACAGAATATCAACCATAAAGATGTGAAAAAGAGAGAGAAAGTAATGGACAGAGATACATATATATACAGACATCAGGAATGGAACACTTAGCTTCTTAGATAGCATTGATTAAATATGAAAAGCAAGACCAGCTGTATTTCTTATTGAATGCCTCTACTATTCCCAGAAATAAAATGTATTGAACAGGCAACTCATACTTGCTATAAATTATGGTCAGATTCGAATATTTATCCTAGATTTGGTATATATAGTTAATCCAAAAGAATCCTGACTCTTTCATCGTCGACTGAATATTAAATATTTGCAGACATCGTGCAGGTTACTACCTTAAAGCTGTTGTTTCTCACTTGACTAATAGGAAGGATGTCATGCTGAAATAAATTAAGAGATGGGTGCTCTGAATTATTTAAACTCTATGACTATAGTTGCATGATCTTATGATTACATTAGTACTACAATTTCAAGTTTTCTTGGAGTGCGATAATGGTCAAATGAAACTAACACATTGTTCTCTCAAAAAGAACTGGAACTAGAATATTATTCAAAATATGCAACTGATGGTACACGCACATGACTTACTGAGATAACGAGTCCTAAAAACTATGGATTATCTAGAAGTTTTTATCTGCTATTGTCAGCGACATGTCTAAATGTCCTGAAATGTCTACCAGGTCCAGCATGTCTAGAAGTTTTTATCTGCTATTGTCAGCGACATGTCTAGAAGTTTTTATCAACACTAATTTGCAAATGCCCAATGAACTAAACAAGGGATATGCAACATTACTACTACAAATTGGGGTGAGCATTAAACTATTTAGTGTAATATAAGACATTAATTTTTCAAACTATCTATGCAAAAGGAAAACTGATCATACTAGAACTAGAAAGTGGAAGTTAAATATAAGCTAGAGCATACATTTGCAATGGCTAGGGTCTCCAAGTACGCTATGAAGTATGACAGTGCTAAAATCCATGCAGCTTCAAATAACCAGCGAATAGACTGTGGTAAGTGAGCTGCAGAGTGACGTAATCTACGGAGTGTCATATTTGATGCCATATGATAGAACAGGAAACAAGCATGGGTCAGGAGAAACGTCGTATGGGGTACCTGAAACATGAAAGCAATTAATAAAATCATAGATAACACGGCAACTGTGATCAATCATTAACTCAGAAAAAAATGGGGACAAGGTTCTTGCACACATTTCTTCTTCTTAAAACAAATGCTGAACATACATGCCAGCTCAACTATCTAGTTTCTATTGTCAATTCTTCTGAACGTGCTAATACAGGTGAACCGTGCAAGCTTGGATTCAAATTCTGAACCAGAGTGGTTCCACATACCAGTTTATCTGTCATTTACTGAAAGTAAGGAACTAAGGAAACAGGAACATGTCATCCCAAACCTTTTTTTAAAAAAAATATGTTAATTAAACTAATAATCCAGTCACTCCAAACTGTGATGTGATCCAAGTAACAGATTTCCCCTGTTCCACAACCATCCTGCAACTTCCATGTCAATAAAATAATAACAATAACTTCTCTATTGAATTTTATCCAAGTAGACCTTAAGAACAGCACCAACAAGAAGTTGTTGGCTGTACATAGTGATCTATTACTATGACATCACATCAAATTTCAGTACATTGAATTACCCAGAATTTATTTATGACAGAGAAAATTCAGTTTTAACTGCAGTGACGTAGTATTGGAAACATGAAACCTCTAGTTGTACTATCACATGCATGTCGGACCACATGACAACAAGTTGCAGTACTTGTAACTGTGCGTTTGCATAACAAAGATGCCACTGAAAATAAGGGAGAACATCTGCAAGAACAACTAGAAGCCATCAAGTCTTACATTATTCATCCTCCATGATGGAAAAGTATATGACGCGCCAAGAACTGTGAAGAAATAATGTGTCCAGAAGTAGTTACCAACATAGCTAAAAATTATAATCCAAATATTAGCCTGCAATGGCAAGATGTGTTTTGGATTTAGCAAATGGTTTTGAAGAGATAAACTAAACCAACAAATCTCACTAGTAAAAGAAAATCCCTAAACAATTGAGAAGTACAGGTCTAGTTTTTTCTTGTATTGGTAGCTTGCTACAGAATCATTAGAGTAGCTGCACACTGAATTCAAATCATAATGATTTCAGAGGGAAACTGGTAATTTTTCCCCAAATTCGCCACATAGTTTTTCTGCACTCACTAGGTGAAATTAAAACTCATTGAAGTTAACTAGCTGGATATTTTAGGATCTGATGTACAGTAACAAATATAATATACACAATGTTGAGTCATTCCATGAAAATAATTATTTAATATTTATAAACCTAAATAGTCACGAATGATTCACAGATGAAAATAATTACCTTGACCCAGTATCGATCTTTGAAATTTCTAACACTATCTGCCTGAAAAGTGAAGTTACTGTTAAATTAGTTTCTGAACATATAGATTCACTTTAGGTATTTAAGAAGAAAAGTTCACATTAAATAATATGAAGAAAATGACTGCATTTTCAATTAAAAAATGAACAGCAATATCATGAAGGGAACCAAAGAAAGTCAAAAACCAATTTATTTTAATCTCCAAAATGATAGCCCTGGAAGCAAAGTACTGGTACAACAGAACCAAGGACACTATTCAGGGTAGCAAGGAGGACATATTGAACATGACCTGGTATGAATTGTACAAATCACCCATAAGGTCGCAAATCACAAAAATGAACCTGTTTTGATACATGATTAAAAATAACAGGTGCTATCAGAGGTGAAAAGAAGAGCACGTGATCGGCAGCATCAACACAGGCATGCCACAGGGCCGCTGACTGACAAATGCAATGTCTCAGGACAGAGGGAAAAAAAAACAGACACCAAAGAGGGCCCATAGTGGAGTGGGCTTAAATCTAATGATCTTCCAGCAACAGATCATCCTGTATCTTGCCCCCTCCCCGTACCAGTGTGCTCTGCATTATTATTATTATTATTTACAACAAGTCTTCTCTGCATTTTAATTGTCTGCCAACTACTTGTTATTCCTTTCCATCTCACATTCTCACCTAATGGTAGTTGGTTGGTATGTTATTATGAAAATGTTACCACGAGCATTGTCCCTTCACTGGTGTGGCTTACAGACAAATCACATGTGCACCTTATACTCTGAAGAATTAAGGAAAATGATACATATCCTGCGACCTTGCTTGAACAGTAAAAGATTAAGGCATACATATGTTGGCGAGATGAGCAGCCAGAGCCAGCTGAGGGAACAGAACAGAACAGAACAGGATGGCGGTATTTATTATTATTATTATGGATTACCTTCCCTACGAGGAATAGGGGGATGAGGAACGCCGGCACGGTGGACACCAAGCCAAGGATCAGGTACTCCAGCTCCGTGAACCTCTGCACATGGCACAAGAAGCCAAGAAGAATGGGTAGACAGATAAAAATGCATTAAGCAGTTTATTCACCCGCACCGATCTCGTGGTGGTGGCCACATCGTCCAAGATCTCCAATGGAGAACCCAATCGATGGACCGGCCAACAACAAACCAACCGAGAATCATGCGTCTGGCTCTAGATGAACACAACACAAGATGGATGAGACGAGACGGAGGAGCACCTCGTAGAGCTTGAAGGGGACGACGACACCGAGGCAGAGAGTGAGCCAGAAAGGCGTGTAGAGCAGGAAGAAGGCCTCCCCCCAGCGCTTGCTCCCGTCCGCCGCCAGCCACGCGCTCCTCCTCCTCCCCGCGCCCCCGCCGCCGCCGCCCTTGGGCTTCGTCGCCCCCGGCCGCCGCGCAGCTGCGAGAGGGAGGAGAGGAGAACCCGGTCAGATACGCGGAAAGAAGCGCGGGCGCGGGCGCGGGCGCGGGACGCACCTGCCATTGGCGGGCTGGTTCGGCTTGCACCTCCCTCCCGCAGGACGCCGGCCGAGGGATGGAATCGAATCGGATCGGAGCCGAAGAAGAGGCCAACGAGAGTGGAAGTGTAGTGTGTGGGCCGGCGGTGGTTGGCTGGCGCGGCGCGGACGGAGGAGCCGAGGTGCTTCTGAATTCTGATAGAGGGGAGGGAGGGGAAAAAAAGTGCGGAGGGAGGCATCAGATAGGCCTTCCCGCGCTTGCGCTCCATTGCCGTGCTTGCATTGCATGGGCCGACCCGTTCAACTTTGCGAGTAACCATTGTAACTGGGCTGGGCTTTGTTTTGTCCCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCCCAGACTATCAGAGATGGTGGATAGCCAGCTTCATAGATCCAAATAAGTACTTCTAATGGTGTGTAAGGAATTGACTCTAATGTCAAAGGTTGCCAAGTAAGCATATAAATATTTTATAATATGAATGTTTGGGATTTTTTTGATAGACCTACGACGCAGAGACTTATATACATGGAGTACACATTTATCCCCATAAACTCATGCACGTACACCCTACCCCTATGAGGATCTACAAAGGGACTGAATTGGATTGACAAATCTCAAGATTGGTGAGGGCCAAGCAAAGAAAAAAGAAAAGAATTAAAGGATATTATAATCATTTTACAATTTTAGCACATAGGAGAAACTATTTTACCAAACAAATTTCTAAAACAGCTTCGACTGTATTAGTAGAGCTGCTCCACTAACTGAACCTAGAAAAAACCGCTTCACCAATGAAGAAGAGTTGTGCCAATCGGAACCTAATAAATTATATATAAGAAATAACTATTTCTTCTGTTCTATAAAGAAGGGCATGTTTAGCATGGCTCCAATAGTGACCTCTGTAGTGGAGCCATTCCAAACATTATAAAAGGAGTTGTTTTCTCAAAAATGGTAGAGAAGTCAAAGTTGTTTTGTATGTTAAAGCTGGAGATGCAAAAGATAGCTCCTCCTTAAGTACTCCTTCCTAAATCTCGCTTCTTGTGGAATCCGACAATTATAACTTTAATTAGGTTTATAGAATATGGTATCACCTAATTCAAAGCAATGCACGAGAATGGAGGACAAGGATAGGTTGTCACCGAGCATATTGGGGAGAACAACTCTCTACTTGGATATGCTCACATTTTTTATGTTATTGTTGTTTTTTTTCTCTCACTAAAGTATCAATCCATTTTTCAAATCATGTGCCGAAGGCATGTAAACTTGGATATGATCGAATGGAGCACAAATATTGTTTATATAGAGCTTGTTCGCTTGGCGGATTAGCTATGGCTAAAAATATTATTCGTTGATTTATTGTGAGAGAAAGACACTGTTGAATGACTGATAGATTGGACAGATAAACTCAAGTAAATAGGCTCATAGTTACTTCATCTAAGTTTGGAATTTGAAAAATAGAAGAAATATGAAATAGGGAGAAATGGAAAAGGGAGAAAATGATTAGGTGGAATGGAACTTGTCAATGAAAAATAAAAGAAATATGAAAGGAGCAATTTGTAAACACAAAAGAGACAAATGAAAAAAAAACTTGTGAATAAAAATAAAAGAAAGAGGAAAAGGAGCAATTGGAAAATAAAAGATGAAAAATAAAAGAAAGAGGAAAAGGATCAGTTTGAAAATAAAGAAGGCAAATGGAAAAAAAACAAACTTGTGAATGAAAAATAAAAGGAAGAGAAAAAGGAGCCATTGGAAAATAAAGAAGTCATAAAAATATGACGCTAGAAGGAGGGCTTGAACCTCCGACCTTGTGGTTAACAGCCACACGCTCTAACCAACTGAGCTATTCCAGCTGTTGATCTTGTGGTCTACAGTTTATTATTTGTCATACTAAGTGCCTACGTCTTCCAAAAGCATGCATATTAAAACCTTTGGTTTTAAAGACGCTCCGGTACAGCTTATAAGGTGCCTATGGCCAAAATAAGTTGGAACTAGAAGCTAAAGTGATCGTACAGCACCACACAGAGTTTCTGACATTAGAACAGTTTTACATACCCCTCCATCCTATAAATAAATACATTTCTTTGACTTCTATGATTAATGTTTGATCATTCGTCTTATTCAAAAAAATAATAAATATTATTTTTTATAACTTATTTTATTGTTAGATATATTTTTATAATGATTTATTTATTTTATTATTTATAAAAAACTTTAAATATTTAATTATGAGACAGAGAGAGTAATTCTCAACAAACAGTCCTTGATTACGTTTTTTTTTTGAAAAGGGAATATATTAAAATATCTATGCACGGAGACACTTCCATGAATTCACCAAAGGAGAAACAAATATAACTATTGAATGACAACAGTTGCGCCCAGCCATGGTACTCTTATAGTTTCCATCAAGCCATCAAGGGTATATCATTTCCATTTATCAAAAAGAAGTGCTGCTATGATTTCTCTCAAGGGTCCAATTCATTCGCCAATTATCCACTGCTCTTTTCAGTGTGGCCCTGGCCTCGAAGATCTCCGCGTCTGCGAATGTGTGCTCCACAATTCCTATTGGACCCGGCTGGAGGACGTGGATTACAGTCTCAACATCTTGCTCCACCTGCCTGAGAGAGAGAGAGCATCAAGAAAGGAAATGAAGTGTCCATCATATATATACTGCACTGCACGGGTTAGTGAAAATTCAAAGGTAAAACTAAAAGATAATACCTGATAATCTGGGGCAGTCTTTGAGGCAACAGTGGCAGGAAAGACATCACGGTAATCTTTCTCTTAACTCAGATTTATGTAACAAAAGTGACCGAGTATATTAAAATATAGAAGATTTTATTCAATTCTGTTGAGCAAGCATTGAAATATTAAATATAGATTATTTATAAAACTGAAAACATAGAGAATAATTTGCAAAATAAATCTTTTGAGTCTAATTAGTCTATAATTAAACACTAATTATTAAATAAAACAAAAATACTACAATATCTGTTAAACTTTAACACCTCTAACAAACACGCCATAAGATCCAAACACGTCAGAAGATGGAAGGATAAAGGGGGTAGATACTGACCTGAGCAATCTCAACAAGATTTGACCAAAATAGCTAGCTAAATTCAGCTGCACCAGCAGGCAACAGTGGTGCAAAGTCACTTTCCTTGCCATTAGTACTAATAATTATCTAGTGGTCTGATTAAACTAAGCACACACAATGCTTGGCCTGGCACGCGTCCTAGATGCGAGGCCAGGTGTCGCACCACTGAAGCTCAGCGCTGCCACGTCTTCTCTCTGCTGGCAGCTTTCTTTGTTTCTTTCCCTGGGACTGGGATCAGCATCTGCTTGCTACTACGCTTCCTCATCTGCTGGCACTGGCGTCTCTCCTCCCCTGGTACCTGAGCTGAGCTGAACCGGGGCAAAATTTTGCACTGAGATTGAGCTACTACCTCCTCCAGCGCCTCTTCTTGCGCACCTAGAATAGAAAGAGAACTACTGTAGAAATATACAGCTACGAGAACTAGAAGGCGGCGAGCACATCGGCGGCAATGGGGTCGGAGACGTTCCTGGAGATCCTGCTGGCCATCCTGCTGCCGCCGGTCGGCGTCTTCCTCCGCTACGGCATCGGGGTAACTCAGCTTTGCATTTTTCTGATGCATCATCATTTGTCAGATGTAGGAGTAATTGCTTGCAGTTTACATATCTAGTTTGTGTCTGTGACTGCAGGTTGAGTTCTGGATCTGCCTGCTGCTCACCATCCTGGGATACATCCCGGGCATTATCTATGCCGTCTACGTCCTCGTCGCATGATTATGATTATTATATATAATTATGATGTGTGATAATGACCGGCCGGGAGCTAGATTGGTTCGAGACCACGGGTTTCTGAATTCTGATCCATCCATTCCATACATGTAGCCCATGCTGTATGTTCTTGTATCTTCTAAATTAGCTTAGCTAGCGGTTGTTTCTCTTGTGTCCTGATCCATTTTGTTGTGCCCATGGACAAGGATGGGAGAAAGGAACTGTAGCTTCGATGCTTCGTCTATGCAAGTGTGCTAGATCTAGTCTAGAGTTTCTACTCTACTAATTTCGAGTCCAGACATGTTTCTGAAATCTGCATGCACTGCGATGTTTGTTTAATATAATATATGTGGTGATGCTGATGCATATGAATGGATCGTACCGCCCATGCTGAACCCTCTACGGCCTCGTTTAGTTCTAAAATTTTTTAAGTTTTTCATTACGTTGAATCTTACACATGGAGCATTAAATATAGATAAAAAAATAATTAATTGTCCGTTTTGTCTTGTGAGACAAATCTTTTGAGTCTAGTTAGTCTATGATTGTACAATAATTCTTAAATACATACGAAAGTGCTAATTTAGTCAAACTTAAATTTTTTCACCAACTAAACAAGGTCTAAGAAAACATCTTTAACAATCCTCCAATATGACCTCTTATCTTTATAAAATTGTGATATTGTGATATTGGACTGTGATATGAGGGGAGATGGTCACTTTTCTTCTTCGACAGTCTCTTAGGGCATGTTTGGAACAGATGAATAAATTACGCAGGAATAGGAGAAAATAGAGAAATAGGACAAAAATTGAGATGCAAAACAGAGAATTGGAAAACACATGAGTTTTTTTAAGAATAGGTGTTTGGAACGCAGGAATATGACTAGTCTATGAATTTTAGAGTTAGCTCTAAACTCAATTGTTTAGATCTAGTCGTAAATATTTTTTTACATTTACAAGGCTATAGATAGTTTTTTCTTTTATACAAAGTACAAAGTCAACCCCACATGTTATCGTTTATTTTTCTACTACAGTACTAGTGTACTACAGTGCACTTTTTTTTTTCCTTGACTAACAGCCTCAGCGTCAAAGTGTATGGCAGCTCATCAGCCTAGGGACTAACACCAGCTCGTAAGCCCTCATCATCGACAGCAATTGTAAAGGAAAAGAAAACATGGGACTAGAGTGGATTCTTTTTTCCTGTGAAATCGACATCTGGCTACAGGAAATCAAAGGAAAGGAAATCTCGTTTTCTCTGTTCCAAACACATACTTCAACTTCTTTTTTCCGTTTCTTGATTCCTACGATTTTTCTTTGGATTTCCTCCAATCCAGACAGGCCCTTAATATGGTTCTGCCACAAAAAAATCTACCCATCTCCTCTCAAAGGATGGGGAGGGGGAGTTGCATCTCCTCCAATAGCTTGGGTCGCAATGTTTCTTCTTTCTCTTTTTTCACCCCGGTTGTTGGTTGGATCGGATCGGAGCTCCTATGGGCTATGGAGGACCTCATCCGACCTAGCTCGCCGGATCCGCACTAGATGCAAGAGGAAGGAGGCGATCGCCGTTGGGTCTTGGCGCTGCCAAGCTCGCCGCCCATGGTTGTGCTCAAGCTTTTCTTATCCATGGCTGCGCTCGAGAGCTCCACAACGAGCCTACGTCGGCCTTCGGTGCCTGCAGATCTCGTGGAGAGAGAATAGAGATGCAAACCTCGGTTGGGGACTAGATCGTGGGCGCTATGGCTGAGGTGGAACGGGTGGAGGAGAAAAGGTATGGTGGCTGCGCTCGAGAGCTCCACAACGAGCCTACGTCGGCCTTCGGTGCCTGTAGATCTCGTGGAGAGAGAATAGAGATGCAGACCTCGGTTGGGGACTAGATCGTGGGCGCTATGGCTGAGGTGGAACGGGTGGAGGAGAAAAGGTATGGTGGGAGGAGACAGGGAGCGGCCTTCGACTCCAGTGGCCATACCAATGGTGAGTGAGGGGACTGTTGGGGTTTAGTCCCACATCGTGTAGCGATAGTGGGGGAGCACAACAGATAAGGTGGGGGAAGTTCTCACCTGACGTCTTTTGGATTGAATTGGGCCCAAGACCTTATATGTTATCAGCTCTAATGGACTTGAACATGTGGAGCGTACCTGTTCTAACAGGGACATGTGGTCCTCATGGGAAAAGATGAAAGAAAAAAAAGAAAAGAGAAAAAAAATGAAGAGATTATTGGTCTAATAAACTCAGAAAATGATATTAGGTCATCTCAATACCATTTTATTGGGCTATTGCTGCTCAATTCTAGTTGTGCCCTTGTATGGAGTCATTTTTTTTAACAACTCGTATGGAGTTAGTTGGGCACACCGCACACGAAAACGTTTATGGCCTGCGAAAATGTAGAGATCGCTATGGTTCACGTGACCGAATGCCAAGACAACGAACTCGGCGATCCCAATCCGATGAATTTTTCGCTCAACAAATCTCTAGCGATATTCAGTCCAATGAATTTTGAATTCCTCTCGCTCAACAAAACTTCAGCGAGATTCATGGAGGTCCCTCAACAAATCTCCAGCTAGATTCATGGAGGTAGTTTCAAGGTTTGGGTTGAATTCATGAAAAGATTTCAGATTGCCGGCGGTGACATTGCTCTGAACTCCGAATCTCCGATATATCCAAATCGTCAGCCATAAACGCACCCATTTATTTTGAGGCAACTCGCGTATATGATGTATAGCTTCTGAAGGGAAAATAAAAAACGCTTCAAAAACATGGAAAGAACCAAGGAACATCTGGGATACTGGCTGATTCAGATCATAGTAAGGTACAGAATAGGCAGCCAAAAACAGCTGCTCCTCAACACAAGGTTCATTCAAGTACCAGCCAAATGTTATAAGACGGATCATGCTTGATTCCTACCAGCTACTCTTGTCCAGCTGCGTTCATCAGCTTACCACTGCTCATCCTGCTGCTTGTTGCACCTCGCACTACTACCACTAACCATACCCCATCCACCTCCTCTCACTTGGGGCTCTTCTTCCATAGCCCACTCATCCTCCTCAGCACGCCGCCACCACTGCTGTTCCGCCGCTTGGGCGACTCCTCGTCCACCTCATCCGAGAACGGCGAGCTCATCAGCTTCCCACCAATGGAGTTGTTGTACTCAGGCTCCAGCGACCCGGTCCTCTCCACCGCCCTCCCACCATTGTTCTTCTCCCCACCCCCTCCACCCAGCACGGCGGCGGCAGCCTCGGCGGCCTTGCGCCACTGGTCGGACTGCACGCGCAGCCGCCGTAACTCCGCCTCCATCTCCCCGGTGGCCACCTGCGCCGCGTCCAGCTGCTCCGAGGCTCGTGCCGCACGGCGGCTGCTCTTGTCGGCTTCTTCGGTCACCAGGCCCAGCCGCATCCGGACGTCCTGCTCTGCGGCCTTGGCCAGCTCCAGCTCCGCCACCGCCGCCTCGTACTTGTGCTGCAGCTCCGCCTCCACCCTGCCTGCCTGCGTCTTCAGGGTCTCATTCTCCTCCGCCAGGCTCTGTAAAGCGTTCTCCTTGTCCATCAGGGACGCTTTCAGCTCCGCGACGTCCGTGATCGACTTCATCAGCTTTGCCTCCAGCTCCGGCTGCTGCATGACGTCTGCCCTGTGAGGATCCTGCTGCGCTTTGCCCGCTGCACAGGCTCCTGCTGCCATTAGCTCCATGAGTGCATCGTTCTTGCTCTTCAGCTCGAGCTTAAAGTCACACAGTTTGCGGGCGTACTCAGCCTTCACACTCTCCAGTGTCTCGTACAAGCCTTTATTCTCGATAGTCATCCGAGTTTGCTGCTCCTGGTACCTGATCTCAGCAACCTCAAGAGCAGTTCGCAACTGCTCAATCTCTGGGTCAGAACTTTCACAGCCCTCATTGGCATCGGCGTGGCTAGTGGTGATCAGCACCTTCTTCAGAGATTCAGGCTCAGGGTTGCCGAAGCTGACATCAGTACTGGCCATGTCCTGTGCCTTCTTCAGATCCTCCTGGAGCGACGCTATGAGAGTCTTGGAGTCATTGAGCTCCATGTCCTTTGACCTCAGGCATTCTTGCATGCGGACACCCTCAGCTACAAGTGAGTCAATGGTAGCCTTTCCTGTCTCCAGCTGCTGCTTGGTCTCATTTGCCATGGCATTTGCCTCAGCTGCGACCTTGTCACTCTCACTGACATTGACCTTCAGGCCCTCGATGGTTGCTAGGCGGATCTCCATCTCCTGCTTTAGGCTCTCAATCTCAGCTTCCACATACTCGCATTGTTTCGCACGGGCTGCGTCAGACTCTGCCGAGGCCTCCAGCTGCTGCTTCAGCCTCTGAATCTCGGACATGGCTGAGCTGAGCGCAGCTGCATCTACTGATTTCTGCTTCTGCACTGCCTCAAGTTCTGACTCCCAGGCACGGTCACGTTCCTGTGAGATCTTGCGCAGTTCCTGGAGACGGGATTCCTCTGCAGCTGAGGACTCTTCAAGCTGGCGCTGCAAGTCCTCCAGCTTCGAGGTTGCAGCCTGCTCTTGCATTTTTGCTTCCTCAGCTTCCTGCTGCACATGGCGCCTTCGTGCCTCTGATGAGCTAAGTTGCTCCTTGGCCTTCTTCAGTTCATCTTGCAGCTGGCTAACCTTGGACTCCAATTCAGATAGCCTACTTGGGCGCTTCTTCTGAAACAGAGCGTAAATTATGTGACTAAATGATTCAATATGACGCAAACCATGTTACAGACAATTATGTAGTTTCTATGATGAACAATATGCAAATGTTGAAATAATCATAACAAAATAGAGTTGTTTGCTCAAAACAGAGTAGTTTCCAATCATGAAGCTAACAAAAGGAAAATTGAATGTTAAAAAATAGGACAAGGAGTGAAGATAATGTGCTGGCTGCAAGATATACCGCTACCTCAGTAATTGGGCTTCGTGGAGACCGCCGCTCAGTGACCTTTGGGCTCCTCTCCGTAGGTGTTCTTGTTGGGGTAACTCCTGTAGAATCAGTTTCATTCCCTCCAGTTTTGGCCACACGGGATGCCCGCGGTGTCCTTGGGGAAGTCCTCTGCGGGGCATCTGAAGACCCGCTCCTGCAGAGCTCAAATGTAAACACCTAGCATCAAGTCGTGACAAAAGGAGCACTAACAATCTAAGACAGGCTAAAAATATTTTACAAAACATTGATACAGAAATCGGATATTTTATAAATTGTGATTGCCAACATGTATGTACTAATCAAGTTAGCTGCATACTATACAAGGTATGGCTACAGAGAGATATGTACATAACAAAACATCCATGATCCTAGCATCTAACAGTTAAATTTGGAAATCCTTGTTTTTCTGTTATGTAAGTGGCGACTTCGAGTACAAGGAAATGTCTGGTACCAAAGCCCAACTCCAAGTAATCGCATTAACTCATTGAGCCTCAAACCTCAATGAGTCCACAGAAGAACTAAATACTGATAAGATAAACATACAAAGGTGTTGCAACCATAGTCCACCCAAGGAGATGCGAATTCAGGAAATAAACATATGCAGATACACAGAACATGGCAAGATGACATCTTTTTATTTTATGAACTGGAAATACTGGAACCCAACTGTTTATCGACTCTTAGAAAGACAGTCAATCTGATGCCTCTTCTTTGAAGAGGTTCAAGCCAGAACAATTCCATTATCTAAAAAGACACAGATTATCTTTCAGGATTCACAAAGATGCAATGACAACATATTCAGTTACTTGCTGTAGGGGTGTGAACCTCTCCAGTTTCCATAAAAAAACAAATATTCAGTTACTCTAGAAGAACTTACAATGTTGTGAAGAGTGACTACAATCTAGAAGCATGTAATTTTGTGTGAAAACAGGACACGGTTCATTTATTAACAAAATGAGCTTCCCAAAGAACTAATACAACAGCTGTATGACACTAAGACACAGACAAACAAAAGGATTGAAATCCTACATATCCATATTTCAGTAGTACTTGCATAGATCACAGCTCGAGAAACGATCCCCAAGGACGGATGTAACATAAATGCAATGCAGTGCTAGTCACTACATAAATCTAGCTTGGTTTGGTAGGGCAGTAATGTCAAGCAAAATGGTACATACGTACACGGAGTAAGTCCCTACCTTCACATACAAATATGGACCTAGCCTACATCGGTAGCATTTCATCCTTAAACTATACTACTAGGGTAAAGAGCAATAAGGATCTACTAAGGCAATAGCAGATTTCAAATTCCAAAAATATTTTCATGCTTCCAGAATTGGTCCAGCAGATCATAGGTTGAGGGGTGATAAGCATTTGTACAGTTGTTTGTGTAGTGTACAAGAACCAATGTAACATAAATGCAACCAGTACCAGTCAGTACACATTCAGCACATCAAAAATGGTACATACATCCTAGCAGCTTCAGATACAAATTATGGGCAGCATGTCACCCTAAATTTTGTAAATCATTTTGATATCTCCAGACATGGTGCTGCAGTTGCAAAATGCATGTTGATTCAATTAACAACAGGAACAAAACAGAAATGTACTACAGGTATTTTGACTTTTTGTATTGCCCTAAACAAAATGACCTGCATTGGTCCCTCGGCCACAGTAGAATAGCTTTTTTTTTCTTTTTTTTGAAAGGGACAGTGCAAGAGTTAGATCCATTAAATACGAACCAAGTATTAATCCTTGGATGCCCACTGTCTACTGCGGTCGTAAACAGGGCCACCCAGAAATGCACATGCATGATGCATGTGCAGAACCCAGTGCAGCAGTTCTTGGGCAAACACAAAAGAAGAGAAGGAAATACAACCCAAGAAAGTCCCAGTGCATGTATGTAGAGAGAGAAGAGTGTGAGAGAGCGCTAAAGCTAACCTTGTTTTTGAAGTCTGCATGATGAGATTCCTCACAGAGTCCTTATGCACTCACCAAACTGCTCGGGCAAACAAGAAAGTAGTGATAAAGTCAGTACAAAATTGCACAACTTTCTTTATGTACAAACAAAAGCATGCCTAGTTAATGTTGTACTTAGCATCTGTAAGTATGTAATAAGACACAAAAGCTGAAAGAGCTAAATGGCATCTAAACTTCCCCCCCCCCCCTTTTTTTTTTCAAAAAAGGCATCTATAACTGAATGAACCTGAAATAATGGGCAGCCACTGCCTTTTGATTAGAAAAAAAAAATGGAAACTTCGCTACTCTGCCTTCAGGATTTGTGAAAAATAATAGTACTCGCTATCTCTGCCTGCCTTCACATCAAATTTAAAAGTATGAGGTGTCGTTTCCTGCCCATGCCCTATACCCTGACCCTGAGGCAGAGGCTACAACAGTGGCACCAGCAAGGACCCCCAATCCAAGGACGCCTGTGGAGGGAATCCAGCAAAGAGCCAGGAACCTACTTCCACGACAAACCTAAAAAAAACTATTGGGAACATACCAAACTGCCCTCGGTTCAACCACCCCAGCACACCTCCACGGCGACGGCCAAAACGGGGATTCCCTCCGGTGCTGCTGGTGCAGCACTGGAGAAATAGCAATGGCAGCCTTTTTCCTACAAGGAGTGGTGGAGGAGCAATAGCATGGGCAGCCAGCCTCATGGCTCAAAGTAAGTACTGTAACCAAAACACCTAGGAGAGGGTGAAAAATTCAAACGGGGAAGACGCAGTAAAACTCGCACTGTTTACTTTACTGTCTCCAAAGAAGAAGAGGGGAAAACTTCCACGGTTTTACATGGAAGGAAGGATGTGGTAGATCCTGAAGCGTGCCTACCTGCTCAGCTCCCATGGTTGTGTGGCTGAGGAGCGGGGAAACGCGTCCATCGATGGACAAGAAGAACCGGTGGAAAACAGCCAGAGAAGTTCAAGTGCAGATGCAGGTAATGTGATTGATGCAGCATCAGAGAGAGAGACACACACACAGAGACACCCAGAGAGACACACAGAGGGAGGTTTCTTTAAAAGACGAAGCTTTTGAAGAAGAACAGCTGCAAGCAATAGGCAGGAACGCCTGCAATGAGGGCTGGTAACGCATCGTTTGCGGTTGGGATTGAAACTTAAAAGCGTCGTTAGTTTGTGATTAGTGTACGAGCAGCGCTAAGAGAGGCGCGTTATACTGTACTTATTACTAGTATTCTGTTTCCCGCATGCTCTGTAAATGACGCATTGGCTGCGGATTAAAAGCAATAATAAACGCCCGTACCAGCAGTACTACGGCAGTGTCTTAGTTTAGCAGCATGGCTATGGTTTCAGTTGCCCAAACGCAACCAAACCAAATGGCAGAGTGTGGGTAGCAAGAGCAGAGAGGACGAATGGCAGATGAAACTGAGCTGAGCAGTAAGACAAAAGAAAACCACCCCAAAACTTGAATCTGAAACCAAAGCGCACCCCCAAAAAAAAAACATGACAAGCATCCAAACCAAACCATTCACAGATAAACAGGGCAAGAAAGAGAGAGAAAAAATGGGGTAGAGGAGAAACGGAAGCAAGCATGGTGAACAACGGGGACAAGAACAAACAGGAAGAACAGCTTGCATTTTTTTTTCCCGCGAGGAACGAGCTAGATTCAAGGATGCTGCTGAATGCTCACCTCTGCCTCTGCTTCTTCTCCCTCGTCGCTCGCTCGCTCTTCTTGGAAATCTTGGCGCCGGGCTCCTCTCTCTCTCTCTCTCTCTCTCTAACTCCTGCTGGAGCAGCAGGTGGTGTGGGTACAGTGAGAGTCTCTCCTCCCAGTTCGAAGAGAGCTCTGAAAGGGAAGGGAGTGAGAGGTGAGTTCGTCCCGGGCTCCCTCGCTTTCGCTTGCTTCTATGGCGGGGTCCGGGGTGGCCGGGTGGGGGGGGATGGGCCCGCCGCCGGGGCGTCGTCGATGCTGCAGTGCAGTGCGGCCGCTGCGCCTGTGCGCCTCCCTGTGCCTGTGCGTGCGTACACTACTCCGGCTCGGCTCGGGCTCCCCGCACGAGGCCAGCGACGGGCACGTGTCCCGGGGCGTGCGGGGTCCTGTGGGGCCCACGAGGCAGCGAGACACGGGGGGCGGGCGCCTGGCACACCGCCGCTTCCGCTTCCGCCTTCCTTGGGCCCTGGTCTCGAGGGCTTTCCATTCCTTTCCCGCCGCACAGAACCTGGCACGTCACCGACAGGTGGGGACGGATCGGAGTCCGGGCCCACATGTCGGTGGGTTGCGATGGAAGAGCGCGGGGTGCGGCGCCCCGGCCCGGTCCCGACAAGGACTGGCCAGTCAAACGAACGGGGATCTGTCGCCACTCCTGCCCCCTAGTAGTACGCGCAACCATACGGAACGGGGAAAAAAAAAGAAAAGCCTTTTCTTTTTTTTGTCTTTTCTTGGCCAAATACTAGTGCTAGGTTTGTGTTGATGGCTAGTAATAAGTACTCTTCCATGCTCTATTTGAAAACTTTTTTTTTTGCTTCTGAAATCTATGCATTTTTCTTGTGGAAATGAACTGCTTCGTTAGGAATGTCCATTATCGAATAGACCAAAAAATAAAGGTTGAGCCGTGTGTTTGTTTTGTTACAATTAGAAGAGCCAAGACCTTTGGTTGCCCTCGACCAACCTAGGCCTTTGTTTATTTCATGGAAAGAAAAATTTATAGTATATTTATAACTATAGTAATTTGTATTTTGAAAGATTTGTCTCGCAAAATATATATATTAGATAATTAGTTATTTGTTAATCTATATTTAATGTTCTACGTATATGTCTAAATATTTGATGTGACGGAGAGAAAAAATAGATAAAAAATAAATAGAGCCTTCACTGCCCCGTAGGTACGCATGGACTAGTTCTAATCAATCATGCTCCTTACAGTACAGAAGAAGATGGCACCATGACAAGTTGACAACAAGAAATCACTAGCCTACGTGGCGAGTATTTACGCATACCTCGTTGATTCTCAAACGATGGAATCGTGTATGCAAGGGCCTACAGTCAATCCATCACTCATGTGTAGTCCCTCCCTACTGACTAAACTATGATTTATTCTAGCTTTGTTATAAATCAAACTTATCTAACTTTGACTATGCTTTTAAACAGAAATGTATTATCAAAATATATTTTGTAGAGGATTAAATGAAACTAAATTGGAATAGTCAAATATCATTGTACAAATTGGAATATAGTCCAATATCATTGTAGTGCTGAGCTCATAGTTTACTCATGAATGAATAAACGAAAGATTGTCATGGACTAAATTAAATTACATATATCACTAAGTAATGAATCGTGCAGTCTCTCTTAGTTAGGAACTCATAAGAGGACAAACATGTGTCAATTGAGTGGTAGCTAAGTTTCTATTTCTCTCTAGCTCATAACTTCTTTCTCCCTCTTCTTCACCGTCTCTTCTCTTGTGAATGCTCTCTACCTCTTTCCCTTCCTCCTCCCACTCTCCTCTCTCCACATGACATCAGTAGGGAGGAGGATTAAGCCCCTCTCCCCCTTAGACCATCCTTTTTTGAAATATGCAATGCACAACTCAAAGACGTTTTGTTGTTTACATGTAAGCCCAAGTGCCATGTACACTCCATCTTAAAATTAAGGAGACAGAGCATTCAAATTAGAATGAGAGATATGAAGTCAGACCAAGGTTTAATTCCATAGCCCGACCTGACCAGTGCTAAAGTCAGAACGCCTATCCTCGGGCAAACATAGAATCCAAGGTGGTGACTCCAGCGGCTGGTCAATGCATGAGCGTCATTGTCCACATCATCATCATCATCTTTTGGGGGCAATGATATGGTTACTCTCTAGTTACGTGTCTTCTTACGACTATAGATGTCGAGCCTTTCAATGTTGTCGGCATTGTTGATTCATTGGTGTTCCCTTTGGCATGCCTAGGTGTGGCATCTTGGCTGACATGCATCTCACTTCACCCACCATGGGTCTGTAATAGCATCTCCAAGAGTTTTTTCCTTAAACTTGCTCTTTAAATTCTTATGCGCGGGGAGTCTTTTGAAAACTTGTATGCACTCTAGAGAGTCAACCCTCTCTTTCTTTTTCTCCATTTTTCAACTAGCGAGAAATCCAGTAGAGAATGGAAATATTTAGGAATCTACGTTTAAAAAAAGCTACTAGAAGGTATTATTTCACCAAAATATCTTATATATATATATATATATATATATATATAAGTGCTATTATACACCCTAGTGTAGAATACTAATCTACACCACAAGTCAAAACTGAGGAGAAAAAAAAACTATTCTACATCCTAGGTGTAGAATACTATTCTACACCACAAGTCAAAACCTAGGGTGTAGAATAGCGCTGCCATATTATATATATTCATGTTTTTAACAATAATAGCTACTAAAAATACAAACCATCTTGAAACCATCGCAAAATAGGATATGTAATGTTAGGTAATTTGAACAGTTTTGAGAGCTCAGACTATCTCCAACAATCGTCACCCAAAATACAAGACTCATTTGTCCTTTGGGTAGCGCTACAGGTAAAAGGTTCCATACCTATTTTTAGTCTTCTCCAACAACAAGACCTAAAAGACAACACTCTCTGCAAATGGGTCTCGAGGAGAGAGGATACCCAAATTTGGGTTATGCCTCCCCTGATACCCAAAATAGGTTTTCTGTATGGATACTCTGTTGGAGACTATAGGTATTGTGTTGGAAACCTATTTTGGGTTTGGGTTCCCAAATGCATAGGCTATTGGAGACAGCCTGAGAGCAAGATTTGGGGAGAAGAATCAGACTAGGACTATCCTTTGGTTGAAATGGAATTTTTCCATGGTAGAGGCAGGGCACTCTCTACACGGGACGGAATTGAGGCCCACGGTATTGACATGACGGAATTGAGGCCCACGGTATGGCCGGTCCATGGGCTACACGGGACAGCCGTGCCACGTTGTAGCCCAGCCCTGTACAGCCGTCGTCACATCATCATCGTATAGTATGTACGAAACCAAAGACGTGTACAGTGTACTGTAGCTTCATATATACGGTAGAAGTCGCAAGTTAACTACACGACGACGAGAGGCGACGATGTTTGATCGAGTCTCATCTCAACGTATACCTCATCCAGTAAAAACAAAAGGATGCAGAAACCACCGGCCGTTTCGTGATCGTCCTCGCCATGGGGACGCCGGGGCAGGCAGGAACTCTACTCTGGACCTCTGTTAGCAGTGCCGCCGGCCAGAACCCAAACAGGCAAAAGATGCATGCAATGCAGCCATGGAGGCACGAAAAAAGCCAAGATCTCCGAGTGTTGAGTTTGTCCAAGTGAATTGTTACGAGAGAGGCCGCTGAATTGTACGTACTGAGCATTCAGGAAGCATCGATCACCTGCTGTAACTGTAATGTAATGTAAAGGCTGTAATGTAACCTCTTCTTCATCTCGCTGAAGCTCGATTTGATGCTGCCATGGTCATGGTAAAAATAAAATAGTGCCGCTGTTGGAATGTGGGTGCTTTAACGATGACCACCAGTTATCCCAAAAAAAAAAAAACCCCCGATGACCACCGGCTCTTTGCTCTGCCCTGCTCTGCCATGCATGCATGATTCTTGCTCCATCCATCATTCTGCATTTTTGCTTCTTTTTGGGTTGACAACAGTCCCTCTGCTATTTGCAATCAGTTACGGTCTTAATAACAACATTCCAATTTGCCTTTTAAAAAAAAACATTCCAAATCCATCGCATGACGAAAGATTTCCCCACAGGCAAAAATATCTCCGAATCTATAGCCGAACACGCCCAAATTAACCTGTACCCTCCAGCGGTCCAGCCTCCGCTTTGCGCCAACTATTACTCTATTAGCTGCCTACTTACCTTACCTCTATGAAATTCCACTTTCATCAACAGCACGACAACACAATATACTTATATGTAAACATGCATGGATAGAATTCACAAGCAGAACCATGTTTCGGAATTAGGAACATGTGCACTCAACCGTTTCTAGCGGGTTCCTTTTTTAGTAGCTGATGAAACATGTTTCTTATGTAGTATGTAGTTATAATGTACAGCAGATATACTGTGTATATGTCTGTATTCCACAATCTGATGCAGCACATGAACAAGGCCCCTGTTTGGTTTACCTTTCCTAAATTTTAGCCAGCTAAACTTTAGTCACTTTAATAGCTAAAGTTCTAAATACATTGACTAAAAAAGAAGCTAGCTACTAAAATTTAGCAAGGGAAAACAGGGTCTAAGTGTTGTTTTGGTCTGGGTTGGACTATGTTAATCAGGTTTCCTCAGCGTTGTTGAACGACAAGTTGATCAGGGGAGAAAACGGAATCGGTTTGGATAAAAGCAGATGCTCCTGCAATCTGCACAGTATCAGTAACCGCCTGTGTGGAACGGAACAACAAAAGGCTTGGAGGGAGGACGCAGGCAGATTCCTCTCGTGCAATTGTGAGTTTCTTGCTGGAGCCGGACAGAATGCAACTTTTTTTTTAAAAAAAAAATAACAATGGACAGAATGCAACTCGTTCTGCGCCAAAACAGACACACACCCCACACACAAACACAGAGGAGATGCAGCCAAATCCAAACCGAAAAAGGGGACTGCAGAGAGATAGAAAGAGAGACATCAAAGGCTCTTTGGAAATGAAAAAGGAAAAGCAACCTGAGTTTAATGGTGTGCACCAAGCCACCAACCAACTTTCAAGTTTGGTGCACAACATTTGATCCCAAATCCACAGCCACACCAATTGCATTGTCCTTTTGAAAAATCTTCTTGCATTGCTTGTCACAATTTCACAAAACACAAGCTAGCATCATTCCAATATATGTATGATTGAATGGCAGCAACCACAAGACCACACCTACATAGTACTACTTCTTTAAAAGCCTCTTTTTTTCAAAACAAAAAGGAAAAAAATGTAAAAAAAATGCAGAAACAAAAAGAAGGAAGAAAGTAGTCTCTGCAGTCTGCACATCTTCACTGCAGACGTGAGAGACAGAGCTCATCCTCATGCAAGAATCTCTTTTTTTTTTGCACTAACCCCCCTTGCCTACCCACAATCACACACTCCACTCCACCACCCACCTTCTCCTGAGTCCTGATCTCCCTCCTCTCTTCTCCTTCCTCCCCTGCCGCGCCTCCCTCCCAAACCCCCAACCAAAACCCCACCACCTCCTCTTCTCCCTCTCTCTCGCGGGGTCTCGCTGGCTTCCGCCGCCACCCGTCGCGGCGTTCTCTCGCCCCGGACGCCTTCCCGCGCCGCACCAGTCCGCCCGCTCCTCGGTAAGCGCCTCGATCGCCGCCTCCCTCCGCCGTCCGTACCGTACGTTCCCCTCTCCCTCTCGCTCTCCATAAGAGGCGGCCGCGGAGACGGAGCTCGCGCCATTCCCGTACTCCCTCGTGCTGATCGGTTCGGTTCGTTTCGTTCTTTCTTTCTTTCTTCCCAGGCGTGTCGCGCTACATGGCGTCGTCGGCGCCGGCGGAGGAGGGCGACGTCGCGGCGACGCTTCTCGGCAAGGACGCTTCTTTCGGGTCCAGAAGCACGGTGCTGCGCCCTGGGAGGCTGCGCGTTATGCACCCGCACGTGGCGGAGCTCTTGCGGAGCCCGCGCCGGCACGCGCGTCCGGCGGCGAAGCCCGCGGCGACGCGGCCGACGCAGACGCTGTGCACCGAGCGCGCGAGGTACGCGTGCGCGTTCGAGGAAGACGTCGGCGGGGGCGTCGCGGCACCCGGCCGCCTCGTGTGGGGCAAGGTCAGGGACCACCCGTGGTGGCCCGCCCAGGTGTTCGACGCGGCGGACGCGTCGGCGGACGCGCGGGCGCTCCGCAGGCCGCGCGGCGCCGTCCTCGTCGCCTACTTCTGGGACAAGACGTTCGCGTGGAACGACGCCGCGGCGCTCCTCCCGTTCCGCGCAGGGTTCCCGGGCCTCGCGGCCATGGCCCCCGTCGCTGCCGCCGTGGACGCCGCGCTCGCCGAGGTCGCGCGACGCGTGGCCGCCGGCCTCTCGTGCTGCTGCGGGGGTGGTGCCAGCGCCAACTACAGGCAGGTGATCGACAACGCCGGCGTCCGTGACGGCGCGTACGGGGCGCCCGTGGACGCGGCCTTCGCCAGGGGCGCGCTCCAGGCCGAGGCGCTCGTCGGGTACCTCTCCGCGCTCGCCACCAAGCCGCGCGCCGGGGCCGACAGGGTGGACCTCACCGTTGCCGCGGCGCAGATCGAGGCGCTCGGCCGGTGGAGGCGCTCGACGAGGGGCCTTCCCGAGTACACGGTCGTCCACGGCATCGATGGTGTCGTCACTGCGACAGCCAAGCGGAGGAGGTCGTCGACCGCAGGAGGCGGCTCTGCGAAACGGAGGGTGACCAGCAGGAGCAGCCGCAGCGCGACGAAAGGGAACTCAGCCCGTGACACCGGCGACTACGAAGCTCTGGAACAGGAGGATCTCCCGCTGCCCACACCGGCCCAGCAGATGTTCACCAAGATGGGGAAGCTCATGAGCCGTGCCGCGCAGCAGATGTCGCTCTCGCCGGTGATCCTCAACAGGGCCAACGGCGGCAACTCCTGTCCTCCTCCACCTCCGGCAGTGCCAGCGCCAGACATGGCAAGGTGTGCCATAGCTGCAGATGATGAGCAGTTTCCTCCAGTGAGCAAGAACAATGGCGATCACGAAACTGGGCTGGTGCTCAACTTCAGCAGCGCAAGCGCTGTGCCTTCGGCGCGTCACCTGACCATGATCTTCAGCCGGTTCGGACCTGTCAAGGAAGTCCGGGCAGAGAACTCTACTGCTCTGGTGATCTTCAAGAACAGCGCGCATGCGGACGAGGCGTTCTCAGGCACCGCCAAGATCGGCTCCATCAGGGCATCCCTGGTCAGCTTCCGTATCACCTCGTTGCTTCCAGTTCCAGCTGCTGCACCAGTTGATGATCCTCCACAGCCACAGAGCATGTCATTGGACACCTCCTCACCTGTTGCAGCTCTTCAGTAGCTACTTTGATCGGTAAGCGTATTTACTATCAATGGTAATTTCTTCTTTGATTGGTTCGTGTATTTACTATCATGGTAATTTCTTCTTATTTTGTGTTATATGAGAGTGTAGCTTCATTAGCTGTTCATGTGCTCTTGTTTTTCTTTTTCTTTTTTTAACTTTGGCTAGTAATTACAAATCATGTTTTTCTTTCAGGTGATGGAGACGATGAGGCAGTGAAGAGGTTGCTCGTTCAACCATGCTTCGAACTTTTGGGTTGCCAGGGGATGGTAGTTCGTGAGATGAGATTTTTTGGATTGGAACCTTAGGAAGTAGATCCTTTTTGGCTGTCGTCTACACGTTGACATTTAAATTCAGATGCTGTTGATGATGTAAAAACTTCTGAGTAGTGATGTTTCATGCTCTTATATATAATATGATACATCTGGTTTGGTCTCTAAATGGATTTTTGATGCCACCAGCAGATCATTCTGTTACTAGATTAGCAGAATGTCTTGTCAAGATGCAATCTTTATTTCCGCGTCCAAAATCTAGCAGAATTTGAACTTGTTATCATGCAGGTTTCATACACCAGAGAGTGAGTGCTATCCATCACAAATCCATGTTACTAGACAACAAATATTGTGGCTGAGATTCCGAAGAAATTAGCTTGCATATTCTTTGACAATCAATCATTAAATAGCAGCGACGAAGGGCGTTCGACGCAGTAGTTATATACATACAGTACATACGAGGGGAGAAAACTTGAGCTACGTCTAGAGCAACTAACCTGATCCTCGACCAGAAAAATAGATAAATCTGAGTTCAAACAAATAGCTAAATCAGAATATTGAGAGGTGTAATGGCTATTCACATCAAGTCTAGATCTGTCGTGTGATTCTTTCTAATAAATCACCAATAATCAGTTGCTATCTACTCCTTCTATAATGTACACTGTTCTTCCATCATCCTGTGCACTTCTGTCTTGTTATCATGATAAGTGGAAGGAGCCTAGACCCTCTGTAATGGCATATTTGAGCTTCGTTCTCATTATCTCCTGCGATCAAGAAACAGATACCAATGATAAATAGGTCTCAGTTTTTCATCTTCACGAAAACGAAACAGATTACTAAACTGATATAATTTTTGTGGTACCTTCGAGGAGTAAGGAGGGAGCTTGATAAAGTGCCGGCAGGTATTGACACTTGGCAATTCGTCGTCTACATTGCCATCGCATTGCTGTTGAGCAAGATGCATGGCAGATGATAAGTACTCATCAAAAGAATAGGTCTATGCCTAGTAGTGGTATAATTAGCTATGCTAACAAACCTTGCGTACGACAGTGAGCTTGGGCTCAAGCGAAGCTAGGCCGCCAAGCGGGAGTTGAGGGGCTCCAGTACTGAACTGTATGAAAGCCCGCTGCTCCTCTCTTCCAAACTCTCGCAATATCTCCAGGAACTGCAAGTTGTTAAGTTGAGCTCTTAATTAGTTGATTAGAGTTCAAGGAATACAATGTCTGATAGGCATACTACTTACAATAATAATCGGTTGACTGCTCATATCGTAGCCATGCTCAAACTCCATGTGATCCTCGAGGTTCTTTAACTGCAAAAGATAAACACAATTACAACAATGAAAGGATCCTGTTTGCACAAGACTTTGGGAGTACAGAAAGAAACGTACAGCCCAAGCATCCTGTTCACCACATAGTATGCATTCCATTTCCTCCTCAGTAAACATTTTGAGGGCCTTAAGAGCAAAAACCTGCAAGTGGTAGACCTGAATTAGTATGTCTGCATAGAACGAGAGAAACTCATTTAGCTCTATTTCATGAGAAACAGACCTCATTAATTCCAGACTTGAAAGCCTCGATCTGCTTGGCAATCCCGCTCTTCAGTGTTGCATCAGCAACCAAAGATACATATTCACCTAAACTTTCAAGTGTTACCGGCTTCTCTGAGCCTCCAGGAACAAGTTCATATTCAGGATTTCCAGGAAGAGTGAAGTCAAGACACAAATCCTCCAATTTTACGTTTTTATAGGTCAAATCAACCATAGGGTTAGATGCTCTGGAAGAGGTTTCCAAAAAGTTCTTCCTACTAACAAGTGCTTGAAACTCTATTACTATCTTGCCCAGCTCAGGATCAAACAATGGGATGTCATATATATCGAGCTCCTAAAACAAAAAATACAGTGATTAGAACAACGGAGAAGGAATTCAAAATCAGAAATTTGGAATCAGAAGGGTCATAGCATGGGTGATTGCTGAAAGAATTACCTGCTCAAGCATGACCTTGTAAAATGCTTTAGAAAGTGGAATATCCAGGATCCTTCCATCCAGAACTGCTCTCGCTACAAGATTCCCTAAAAGCTTGAATTTCTGTAGCATATTAGTAAAATCAACTCCTTGTGATGAAGTGCCTGATGGCGACCACGGCTTTGGAAATAATCCAAAAGGAGCATGAATGAATCCATGTTCACCACTGTCTCCTCTCCACATTCCAAGCCCGCCCCTTTGAAGTTCATGAGTAACTGTAGTGTAGAATTCAAATGTTGGACCTCGCCCAGTCCCAACCTCTCCTTCGAACTCTACGTCAATGATTCTACTACTTGAGCCATGACTGGTCATAACTGATACAGCACCTTCAAGGATAGCACTTCGGGCTACCCTGTACTTTTTACTTTGTGGAGGACTCTTAACCTGGTCCAAGTGACCATTAACTTGGTTCAAAGTATTGTGAACCTGTTCTGGAGAACTGTTAACCTGGTCCGGAATGAATGAGCGATGCACAGTCAGGCAGAAGTACTTCCACCGTGTATCGAAGGACAACAGGAATGGGCAGGTCTCAACCAGATAGACACACCATGAGGGTATGAGGCCATCCTCAAAGAAGGTATCTTGCATTTGCACCTCTAGCTTATTTACGAGAAGGCTACTGATGAAGTGATGTCGTGGAACTGGATATATTGACACCTCGAGATCATTGAGGTCTTGTAGAGTTCCTTCAGCAAATTTGTTTATCTGTTCATCCATCAACAGTTGATATGAAAAACGGTTCAGCCCTTCAAGAACTTTCAGCATGAACAAGAGATTGTAGGACAGATCAGACACATCCAGATCACCAGGAAGTTTGCCATGAAGTATGGCAGTGAAAAAAGGGTCCTTTAACCATGCGCGTTGAAGGTTTTCTTGCACGTGGGATAACCGAGTATAGTAGGAACTCTGAGAAGAGATCTCCTTGCTTTTCTTTCTTTTCCTATAGGTTATGTTGTGTTCTTTTTCCCAGAAAGATGGGTCGATCAAGAGATCAGATTGTCCCTTGTTCATCAGACGGAGTATCGATTCGAAAAATGTAACAGATGGTTGGAGTGTCACGCCATTGTACAAGAATCTCAGTCTTGAGGATCTGTTCCCATCATCTAAAAGATACGCGGGTGTTAATTTATTTATGTTGCTATCTTATTATAATGAACGTCAAGCAAGACTAACATAACAAGGAAAGCACGAGGGAATGAGACATGAAAGAGAGATGGGTGTATACCATTTTTGGACTCTACTAATTTTTTGCTCCCATTTGCTTCTATCTCTTGATTGGAATTCTGAAATCATAGAAACAAGAGGGAACTATTAAATATTGCATACAAGAGTTTACAATACCGAAAAAATAAAAGCAAAAACATGCAAGATTTATCCTAAATGGTTATGATTATATGCATGTATTCCCTCCATTCCAAATACTATGTATCTAGACATATTGTATATCTCAGTGCATAGCAAAAGATATGTATCTAGAGAAGTCATATAATTTGGAACAGAGGAAGTAAACTGTATTGGTATGATAGGATAATTTCAAGATTAGAAGGGAATAATTCAAGTTCAAGAGAATTAGCCCATTTTTTGCATTAACATTGAAGAAATTATATTCTTGAGCAACTGTGTATTGATGAATACCTTGGATGCAGCTTCCTGATCAGTTCTTCTAAAGACCCTAGGAAATAGAATTGGCTCAGTGGTATCTGGGGTTGAGAAAAGGTCGACACTAAGAACACCATTGTAGTTCCGTAATTCTCTCTCCTTCCGCGACCTCCGAAATTTTAATTCCAATAACGTTGGTTCCTGGGCCTCAGGATATCTCAGCGGAATCATCATGCTCTCACGTGAAATCTGTTCATCAGATAGCATTACAGGAAAACTGTCATAGCACATGTGCAAAGTGTCTAGCAGCTTCTCCACCAAAATCCCAAGTGGATTTGCTGAACTTTTATTGGAGAGGGTGAGAGCCAAATTAGCAAATTGCTGCAGCCGACTCTGCACCTTACTAAGCTGCTCAAGTACAGGTCCATCATTGAGGTCTGCATTGCAGTACGCCCCATTTGAGAGATAAACTGCAAAACATTTGATAGATCCACTCTGCACAAACTCAAAAGTGGAAGTAACTGGCAGTTCATCTGACAATAATCTCCTTGAAAGATCAGACAACTGCTTGCAGCTATCAGGATTTTCAGCAGGGTGTGTCAAGGCATATACATTCAACCTAGCAAAAAAGTCCCTGACACTTCTGAGAGTAAATCCAAACCTATGTGAAGTTTTTTTGCTCCGTTTTACCGAGAAGAAGTTTTTCTTTATCTCTTCTGCTAGCTCCATGACAGCATCATTCTCAATGCTGCATGCCTCAACAGTCGAGGAAGATTCCAGGTCAAAGCAAAGGCAACTTTCTTTCATGCTGTTGTTTCTTTCTGACTGGTGACTGGTGTTTTTTGATTGTGATACTATATACTCAATTGCATGCTTTACACCCTCCTTGGTAAAGGTCTCCAAGGAAAAGAGCTGGTCTTTTTCCAAAAGGGTCCTTGAAATCTTGAGTGTTTCAAATAAGATGTGGTGGTTTTTCCGGGCCAACAAACAATTAAGAAAGCTGATACAAAGACGAGAGATAAACATGAACATCATGCATGGTGACAAGGCAGCACATTTCAAAGCTCAAAATCTTACTTAAGAAAAAAAAGTAAATTTTGTAGAAATTTTACCTTGAGAGATTTACAGTCCTCTGTAACTCCACTAGGAAATTAGGTGTGATTAATTGAACAATCTTGCCAATGAGAACAACACAGCTGTGGCAAATTGATGATAGTGCAGCAGATTTTGCAACCTACATTGGACAAAAGAAGTGCACATGTAAAGGGAAATTCTGGATTTCTGCAGGCACCCAGAGGTCCTCCTGTAGAGGACCTGAAAAAATACCTCGTACCTGTACTATAAGGGCGACGATGCTAGCAAGCTGGTTCATGTATCGACTTTGTCTCATAATGAGTTTCCTCTTTGCAGTGACTAGTTTAGCATGTTGATCAGAGGTTTTAAGAGGTGGCATAAGCTGGTAAATGAGCTCCACAAGTATCTGGACCTGCAGTTGGGAATTGAATAATGTCAACTGCTAAGAGATTCACCAAGATGTTTGCAGGCTATAGCATGCATATAAATTGCATCCACATATATATATCTACTATCTTGTGATATCACAGAAGTCTACTGGTGGCTTGAGTGTAAAAGGGCTTCCCTTTGTGCACTCTCTGAAAGTAAATGCAACTATACTCTGATTTACTTTTTTTTTACAAATATGTGGCACATGATGTTAACAATATTCCCTCATGGAACTTGGCACGTTAGGGGAACCTCAATGGGGAGTTAAGGCACCTCAATCAAGGAACCTCAATGGGGGAGTTAAAGCACCTCAATAATAGGTGCATGACCACAAAATTGCACTAGGTAATAAAAAGAGAAGTAGAGTACCTTATGACTATCGTGGTGCAAATAGCTGTAATAAGTTATCATCTGTTCAAGCAATTCACAAAAATTTAACTCAAAAAGAGACTTCACAGCCTTTGCGGAGACAGAAGCAAGACCTTTAAGCAGCCCCAGAATGCCCTGCCAAACATGTTCGATTGAGTAAAGCAAAAAGATAACAACGAATTCAAGGAAGGATGATGAGACATACAGTTAAAGTGGCATCATTGAGGCTCTTCCACCCCTCATCAACCATCAAAGTCATTGTCACCTGGACCACATTCGATTCACAGAGCTTACTCATGTGTTTAGTTCCGGAAGCACCAGCAGCAGCCAATGCCAAGCACGAAATAGTGGATTCGAGAATCTAGCATCCACAAGGATCAACGGCCCATCAATATATCTTCCCATATAAAAGAAACAAACAAAGAAACCATCAATATATAATTATTACCATCCTATCAGAATACTGCAGAAGGTTGCACAGTGCAGGCACAGCCTCCATGGCCTTAGCAGCATTCACTTCATCATAATCGGTGAAGACATTCCAGACGATCTGGAGCGCCACTTTCTGTAGGCATATCCACAAACAGTTGGTGAGCACAGGACATACGAGACAGAGAACGAACAATGCACAAGCAGATCAAGCCACAAAACATTCCTGCTTAATGAAATACATGCCACGCATGTTCTCGAAAAAAAAAGCCACAAAACATTCGCTCATGTGCCATTACCTGCTTGTTCGTCGAAAAGAAGTCGAAGAACTGCAGCACAGCAGCTGCCACCCCCCTTCTGAGGCACTCGTCTGGACATTCCAAAGATATGGCGTCCAGAGCTCGCAAGCACTGTCAAAGTGATAGCAAAGAGAGATGGAGATGTCTTTTAGCACCAACTTGTCAAAATAAGAACAGCAAGCGAAGGTAGAAGGACACTTATTACAACACAAGGAGTTGAGTTTAATCAAATATATAAACGGACTTATAAACAAATATATGCACCTGCACTTGTGCAAAAGAGAATAAAGGAACCCAAACTTAGTACATACCATCCCGTTTTTCATCAGGTCGCACCCAAGACCCAACTGCAGTTTTATGCCCCCTGGATACAGTAGCAACTTTACTATGGCTTTAGTGAACTGGAGTGTGTGCAAACCACCTACCCAACCACCGCCCACAATTTTAATTTAAAACTTGACCACTGCACAACGACTACTAGATTATATGATCGCGGAACAAATCCTGTCAATATTAATTGCAGTCCACGTATGAACAAAGAGAGGCTCAGAGCTGGTTAGCCGGCGTGACTACACGGCACATAGTACGTAGATGGAAACTAACCAGATTTTATTTGCAATAACAAATTAACAATCCGCCCGGCAACAGATTAATCACAGCCCGCGACACTCAGTATTAATGGATCGGTACGACTAATCTAGACGAGCAACCTTAATTAGGAGCAGCCTTATTACGATCGAGACAAGCAGACAGGGTTTAAGAAAAAACAGTGGTGGGAGGTTTATTAGGGTCTTTAGGTTTATGCAGCTTTATTGGATTAATGGGGTTTTAGTACGACGGCACGGCACGTCACTAGGAGGGGGAAAAGAAAAAAAAGGAGAGAGACGAACAAGAGTGCGCAGCAGAAAAAAAAAAGCAACCAACCTCCTCGGCGAGTTCGATGCAGTCGACGGCGAGGAGCCTGTCGCGGAGCGCGTCGACGGCGCCGTGTCGCGCGAATCGCGGGGCCCACTGTGGCGCGGCCTCGCAGGCCTCGGCGATGGCGCGCGCGGCGAGCAGCGGCACGTCGCTGCTGCCGCAGGCGCCGGACCCTGAGCCGGAGCCGGAGCCGCCGGCGAGGAGGGCGGGGAGGCGCGCGGCGAGGCCCGCGTGCGGGATGGCGAGGATGAAGTCCGGCCCCGACACGGCGAGCACGTCGCAGAGCGCGGCCAGCCCCGCCGCCAACTCCCCCGCGCCCGCGGACTCGTCGGCGATCGCCGCCACCGCGTACTGCAGCACCGCGTCGTCGTCGCCCGCGTCCATGCCGAGACGGACGCGCGCGCGCTCCCCTTCCCGCGCTCGTGGGGGTTCCGGACTACGCTAGTAAACCTAAGTAGCGGGCGGCACGGCAGCAGCAGCAGCAGCCGCAGCAGGGCGCAGGTGCCTACCCACCGATGGCGGCGGCTGGAGCTAGGGCGAGCGTGGCCGATCGGTCGGTGGGCGTGGGGCGTTGCGATCGTTCGTTCGCGGGCGCGCGCGGACTCGTGCCGTTAGTTGGTTGGTGGTAGTGCTGGTAGGGTACGGCAGAGCGGCGAGGTGCGAGAGCCAAGGCTTATATTGCGGCGGCGAGGGGGACGGGCGGCGGACCCGGACTCCGATGCGGCGCATGCACGGGACCCGACCCCGAGAGGCCGGACAGGCACGGCCCGGCCCGGCACGGGAGTCCGTGACGCCCCCGCCTCCGCGGCCTGCGGGGGTGGCTCGTCTCGTCGCACGGATCCCGTATCCGACTGGGACTCGAGCGCGGCCGAGTCCGTGCGGGCACGTGCCACGTCAGCACGTAGGTGGGGGCCGCGCGCACAGGCAGCCTCCCCTCCCCCCACGCCCACGGGAACGTCGACAAAGAGACGGAGAGGGGGACACGCGCGGTACCCCACGAGGCCAGAGGCTTTTTGCCTGTGTGCCGATTTGCCGGGCGCTTGCAGCGCCGCAGGAACAAGGAAGTCGACCGAGGCGGCCGCATGGATTCGCCGTCGGAGGAGGTCGCCACCACGGGCCGTGGGCGCGCGCGTGGGGCGTAAAGCGGTCACGGCGGCGGACGTTGCCATTCAGGGTTTCAGTTCCGATGGATGGTTGCAATTGGTTGGATTTTTTAAAATCTAACGGCCCAAAAATGTGTTCTAGTCTAATCGCACAAAAGTTCTGCACCATGAGTTCCAATCCTATGATATAACGTAGCAATTCTAAGAAAGTACGTTACACGAAAGTAAATGTGCAATGGAAAAAAAAAAGAGAAGAGACGTGATGATGTTTTGTTGAGGTATTGAAAAGTTTGTATTATCCACTAGTTCTTATTAGCGCGCCAATGAAAAGGTACCAACCACCTTAGACCTGCACAAAGGTTGATTCTTCTCCGCTTCAGAAAAACAGATCACTGTACAACAACAACTTTTAGGTCTCTCGGAGTTGCACCAACAGGCAAGAGATCACCAAAATGCTCCTCAATACCACCAAACCATGCCAGTTCGAGCCCATAGGTCACCGAGGGCTCGGGCTGCATGGCCATATATATGCCTACACACGGGAGATATCTTTGGTACCACACACAAGAGATGTCGCCGCCCTCGATGACCTATGGGCTCAAACTAGCATGGCCCGAGGGAACAAACCATGCCTGGGCCACACCCCAGCCACGTGGCCCGACATGGCATGGCTCGCAAAATAGCCTCGGCCTCGTGTTGGCCCGTCAGATTGATATCCCTGCACATATTTGGTCTGGTTCATCACCAGATTTCATGTTTTGTTTAGTTGGGGAAAATTTTTGGTTTTGGCTACTGTAGCACTTTCGTTTGTAGTAGTAATTATTATCCAAACATGGACTAACTAGACTTAAAAGATTCGTCTCGTAAATTACAGACAAATTGTGTAATTAATTATTTTTATCTATATTTAGTACTCTGTGCATGCGTACAAAGATTCAATGTGATGGAGAATCTTGATTTTTTTTTGAAACTAAACAAGGCCTCATATATATACAACATTTTCACTCCAAAACACTAACTTCTCCCCTATCCAAGCAACGACAACAGCTTTTTAATTTGTGATGGTTTACATGTCCAGATGTTCAATATATATTTTCTAGATCGTATTTGGCATGGGTTAACAATGTTGAATAGAGACAAAACCCACATCAAAGTTGTAGTACTTAACGTGATCTATAACTGTGTAAGGGCACTCACAATGCAAGACTCTATCACAGAGTCCAAGACAATTAATTACATATTATTTATGATATTTTGCTGATGTGGCAGCATATTTATTGAAGAAAGAGGTAGAAAAAATAAGACTCCAAGTCTTATTTAGACTCTAAGTCCACATTGTTCGAGGTAATAAATAACTTTAGACTCTATGATAGAGTCTGCATTGTGAGTGCCCTAATTAACTTTTTTTTCATTTAAGATGATTACTGTTTGGAATAAAATATGTGTTGGTTCGCTTAAATATCTAACGTTTTGCAGCCAACGGTTTGAAGAGTAGAGCAGCTCGGCAAAAGGAATCAAGGCGGAGGAAAAAAAAAACATAGGTATTGCGAGATAGTATTTATTTTTCCTTTTTTCTTCTAGGAGTTGATGAGAGTGGTTTGAAACTGAAATGAATTATTTGACTTATACAGACGTGAGACGTGCAATCCGCGAACTTGCGTTCCGTCAAACTCGGGTATATACGTGAAAGCAGGGCGAACGAGTGAATAAATTAGTGAGGTAGCATTTGTCGTGCCGCTGCGGTCGTGGGATGTCGTCTGTCTATACCAACATTATTTTATACAATAAAGAAAGAATTTACATTCGAAAAGAAAAAAATAATAAAAAAATGTACTCCTAAATTGAGTTGTACGTACGTACTGCATTCAAAAGAATACACTATATATTGAATTTCCCGGTGTTATCGTGTCGTTCTTCCACCTCTCGAGGTTGACGGAGTCCGTACGGGAGGAGGAGCAGGTGCGGTGTGTGTGTCGCCGTCGTCTATCTACGCTCTCGACGATGTTCCCGTGCATGTCGCTGCACGGCGCGTTGCTGTTGACTCGCGGTGGTCCTATCCCGGGAACGGGAAGATGGCCTATGCTTTCATTTCTATATTTTACGTTGCGGATGTGAATATTCCTGCTTTCACGTTATAAAATAAGGGCACTCACAATACAGTCTAAAGTTATTTATTACCTCGAACAATGTGGACTTAGAGCCTAAATAAGACTTGGAGCCTTATTTTTTCTACATCTTTCTTCAATAAATATGCTACCATATCAGCAAAATACCATAAATAATACGTAATTAATTGTCTTGGACTCTATGATAGAGTCTTGCATTGTGAGTGCCCTAATAATCACTTTTTGTTGGGTTAGTTTTCCGAGAACAGGGATGAGAACCTTGACTTGCTAAGCTAGGCTACTATTAGTCCAATGGCAAGGATTCTATTGGCTGGCTGCAGCCCCTGTGCGATGGCGCAAATGTTCTTCGCGGACGGCGAGCGGGTGAGTACGCGCGCCGTACGTCACGCTGTGACCTTGGTTAACCAGGCGAGGCGAACGGTGGATTGCTGTCAACGGATCGGAGAGCAGGTAGCCCGTGTCACGGGTGACGACGGGCTCGTGCGTGCCGTCCCCGATGCCATGCCGGACGCCGTTTTGGCGTGTCGTCTCGTCTCGTCTCTCGTGTCGTCTGTGTCCCGCCACCCCGCGCAGGCCGCCGCCGCCGCCAGCAGGCCCCGCAGCAGTACAGTTGAACTCACACGCCACGTGTCCCCCTCGCCCCCTTCACGCCCAACTGCAGCTGCAGGCTGCAGCCCCGCCCCCTGACACGCCATTAGCGCCTCCGCCCCGCCTACTTTAGCAGGCCATCACATTCCATTACCAGCAACCCCTCCCAAATCCCGGCAGCCCGGCTCCCTCCTCCGGTCCTCCCAATCTCCAGCCGACCAAAAGCCTCCTCCGCCTCCGCCTCCGCCTCCGCACTTCCGTTCCGATATCGATGGGGCAGCAGGAGGCCCAGACTCACCACCCCGACGAGGCGGCGCTCCCGACAACCACCGCCGCCGCCGCGTCCTCCTCCTCCTCCTCCTCCAGGCTCTTCACGGCGTGGCTGGTGGCGTCGTGGTACGCCTCCAACATCGGGGTGCTCCTCCTCAACAAGTACCTCCTCTCCGTCTACGGCTTCCGCTTCCCCATCCTGCTCACCGCCTGCCACATGACCGCCTGCACGCTCCTCTCCGCGCTCGTCCACCACCACCACCACCACCGCTCCTCCTCCCGCAGCCGCGGCTCCAGGTCCAGGGCGCAGCTCGCGCGCGTCGCCGTCCTCGGCGCCGTCTTCTGCGCCTCCGTCGTCGCGGGGAACGTCTCGCTCCGCCACCTCCCGGTCTCCTTCAACCAGGCCGTCGGCGCCACCACGCCTTTCTTCACCGCGCTCCTCGCCTACGCCGTTGCTGGTCGCCGCGAGGCCTTCGCCACCTACGCCGCGCTCGTCCCCGTCGTCGCCGGCGTCGTCATCGCCACCGGGGTGAGTCAGTACTCAGTAGCGCCGCGACGAGCGATTCGAGCAGAACCTGAACCTTTCGCCGCGAGAGATTGGAATTGGAAAATTTGGAATGCCTGACACTCCCAACCCAAATACACAATGCAGGGCGAGCCGAGCTTCCACCTCTTCGGATTCATCATGTGCGTCGCCGCCACGGCCGGGCGCGCGCTCAAGTCCGTGCTGCAGGGGATCCTGCTGTCGTCCGAGGAGTGAGTAAACTTTCCCAGCCCCATTGTTTTGTGTTTCTGGCCGCGCCGCCCGTCGGGGTCGGAGCAGACCAGAGGCTTCCGTTCCGTGCCGCTAGCTGCTGCAACATTCTGACACGGCCGCTCCCCGCGCAGGGAGAAGATGGACTCCATGGACCTCCTCCGCTACATGGCGCCGGTGGCCGTCCTGCTGCTGGTGCCGGCGACGCTGGCCATGGAGCGCGACGCGTTCGGGGTGGTGGCGGACCTCGCCCGGGTGGACCCTAGCTTCCTCTGGATCCTGCTCTGCAACTCCTGCCTCGCCTACTTCGTCAACCTCACCAACTTCCTCGTCACCAAGCACACCAGCGCGCTCACGCTCCAGGTGAATTTTTTTATTATTCAAGCTCACCGCATGCTCCGATGGCTGCACTGCAAGCCGATGCTGCCGACCTCCTCCTGCTCCTCCTCCTCCTGCTGCTGCCATTTGTGTCATTGTGTGGTGCTGATGCCCTCCTCCTTCGCTCCTCTGTCCTTGCAGGTTCTTGGGAACGCGAAAGGCGCCGTCGCGGTGGTCGTCTCCATCCTCATCTTCAGGAACCCCGTGACGGTCGTGGGGATGCTGGGCTACGGTGTCACGGTCGCCGGCGTCGTCTTGTACGGCGAGGCCAAGAAGAGGAGCAAGTGAAACTTGGTGTGGTGGTGACTCAACAAGCTCCTGTACCTTAACCTGTACATTGGTAACCGGCGGTGGCCTGAGAGTTTTGTTCCGATGGAATGGAACGCCGGCGATGATAGACAAGCAAGTAGGTTTGTTCCTTGTTAACCTGCTCCAGTGTCGCAGTCTGAGCTCAAACTGTAGGGAATTATATTAAACTTGCATACCGTCCTATAAGAAGCGTAGTTAACTTGGTCATTGGCACATGCGGTGAGACGTACTGATCGATCCTTTTTTGTAACAGTGTTTGTTTATGTACATTACTATTACTACAACAAATGTACGACAAAGCGAACTGGATGTATGTTCATAGAGCTTATTGGCACTAGCCATAACTAGAGAGGAGGAATGCAGTTTCTTTGCTTGTGCTCTGCTGCTGTGCACCTTCAGATTCGCTCTAAGATACGTGCACAGCACGTGCCCTTTGGGGCAATTTTAACTTCCCTGAGACGAGGCAGTGAGGCACCATGGTCCATGGTACTGGCTGCAAGGTTTACAGTGGATCCGTTGCTTTACTGATTGCTTTTACGATTGCAAACGGCACGGGCTGCGCCATCAGCCTGCTGGGCCGGACCGTGCCGCGAGCAAAAATCAGCCCGCACAGCAGCCCATCTATGCGAAGCGGGTTTAGGGAGCGGACACATGAAGCCCAACGGCGCCGCTCCAACCCGTTAGACTCGCGAGCGACTCATCCCAGACACGAACACGATGACGCGAGGGCTGCCGCCGCCGCCGCCGCCGCCGTTGCAACCTCGATGGGCGGCTCTTCCCATGCCGCAGACATCGCCACGGGCGAATCCGCCATCGGCCCCCGCCACGGCCCCCCTACACCCCGGCGCTGTTGGGACCTTCGCCGGCAATGAGCAGCCAGCAGGTATGCCCAGCCCCTTGCCTGTAAATTTATCTCCTTGATTTGCAATGGTCGTCCGATGATTGTTATTTGGCTCTTGTCCTTTCTAGCAATAAATTGATGGTTGGGGCATTGTCGCTGTGCAGTTTTGATTTAGATTTGATTAAAGGTGGTCAGAAAAAAGATCGGCTTCTTTTTGCTATAAATTGTCTTCAGTTCCGTGAGGATCTTATGTCCCCTTCATCAATCGTGTAGGATATTCCTGAATTCCTTCATCTAAAACTGTATTTCTATAATTTATTGTGATATTCCGATAGCAATATTAGTGTTGAGGATATATGATTGTTCAGGACAGGGAATCCAGAATTTTGTGGTTCCAGGAAAGCATTAGCTGCTCCTTTGTAGTCCTGAGCATAGAGGTATAACTACTGTGAGATACAAGTGTAGAGTATATTAAGTTGCGAGGCTAAATTTGACGCAGCTGCTACACTGCCTAGTATATTGTATACCTTGAGCAAGCTGAGAGCTGCACCCTCACAGTACCACGCCTTTGGCCACCGAGGACGCATGGCTTTGCATTGTTGAGTATCTGACAGGGCACACTTTCCATCTCCCAGCCGCAACCAGCATAAGCTTCTGTTGGAAAACACGGTAGCATCAAGAAGATTTTTGAATAATGCCTACATTAGTAATCAGAAATGTGCAGGGCCTTAGATTATAAATGAGTATAGAAGCAGACTGTCTAGTGATACTGTAGAAGCACTTGTCTGTCTTCAAGACTGGCTTCGAGGTAGATTACATATTCTTTTGTGTCCTAACATGTTTACTATAATATAGTTGAATAACTTTTGTGTCCTAACATGTCTACTATAATTTGTTTTGATTGTTTTAGATTGTAGTACCCATGATGACATAGCTGGATGTATTATGGGAGATGAAGATGATGATGAGATCTAACTCAGGTATAAATTTGCAGATGAATTCGCTTTATGAGATGATCTTGCTTTTCTCAGTAAATGATTTGCTTTTCTCAAATCTCAATACTATGTTACAACAGTCAGTCCATATGATTGGTTGCTCTCCTTGTTTGCTGCCAATTTCATTGTCTAGATCTAAAAGACCAGTTGGTACCTTGCAGTATTTGGACCTTTATTTTCAATAAGTACCGTCTGGCACATCATTTATATGTTGTAGAATGGAAGCAGTCTTTGTCCCTTTGTCATGCTATCATGATGAAAAGTGAATATTTCTTGGTAGTTTGGAAGAATGCTTCTTCCTCAATCTTATTACTCAACATGCGGCAGGAAACATCAAGCTTTTGCAAGGCTACAATAGTGAAAAGTGACACGAGTTGAAGTAATCATTGATAGGCAATATGTTTGAGTTTAATACCATTGATGAATCTATTACAAACTTGTGATGATATGGATGGCAGTCACCTTTTTTATTTTTTTAGTTTGTTGTATGAACTTGTGATGACACTGCAAGGTTGCAAATGCCTCTTATTTATATTTCCCTTGTCTTATTGTACCACTTTCTAAAACATTTGGAACTATACTGCAGTCAGCAGCAGCATGTAATTAAGTAGTAAGTAGCAGTTTGCTTATTGCCTTAGAATTATTCAGAGATTTCCCTTTTAGATCAAGACTTAAGAGCACTGAATTATATTGCAACTAGCAACAACTCAATTTTTCTTATTTTTCTTTGATTTCTTGTTCAGAACTTGTAATTTAATGCAGTATGCTTTTGTTCAGCAGCAGCAACCAGCTTTGACTTATTGTTCAGAGCTTACGGCGAGATAGTGACCGCAACAATTAGCTTATGATTTTCTTAACAATTAGCTTATGATTTTCTTAGATGCTAACTTGCATGTTTTTATGATAACTTATAAATGGTGAAAAACTGTAAAAGAAAATTCAATCATATATTTCTAATGATACATAATATTTTATGTAATTTACATATATTTTTTAGCTAACATATAGATCACTAAAAAATTTATAAAGTAATAAATAACAAAATAAACACACACATACATAACTTACCACGATTCTTGTTTCATGCAAGGCTGCAAGGCAGCGCGGGCTTGCGGCCGCCGCAACAACGCCGCACGCGCACGCGGATGGAAACGGGCCGCCGTATAAGCCCGTGTACGATGATGCGGCCTTTGCGGACGCGGGACGGGCTGCCCGTCGCCGCGCCGTCTGCAAGTTTAATTGCTTTCCACCACAGTTTATGCAATGCATTCCACGCTAACTAGCGGTTCAGTTTGGCCATGACGTGCTTCTTGGGCCGGGCCTTCTCCGTTTGACGGGCCGTCTTCCTTGCCTTTTCATTGAACAATTTCTTTATACGAGCAAATGGAACAAATTGTATAATGACAAGGAAAATATTGACAACTCAAGTGCTATATCGTACTACGAAATATTACTATCCATGCACGTCCAACAATATACTTACTACTACGACAATATATACAGCCACAAGGTCGCTACGTACCAAAATTCAACTCTGTTTACTCATGACTAGCGTGTTGCTCAGAAGATTATCTTGCACACTCACTTTTCGTTCACCTTTGCTCCTGTCACATCTTTGCGCTTAATGCAACTCCCTGAAGGTTCTGCATTCTTTTGGATCCACGGATGCTTCATGATGTCTTCAAGACAAAGCCTCTTGCTTGAATCCTTAACTAGCAGCTGCACACACAATCACACAAAAAAAAGGCAACAGTCCAGCCCACAGCCATTTAAATTCGTTAATATCCTATATAATAACCCACTACTTTTAGACGTCACCTTGGAAATGAGATCCTTAGCCTCTGAAGACACACGAGGAGTTGCAGGGAATGCCCAATCCACCTTAGCTATCCTGCACAAAATTCATGTGGCAGTGGCACACAAAATGTAAGGAATGAGTACCAGGGAGTAAAGAAAAACAGACGCTTCAGCCTTCAGCCCCATACCTTCTCAAGGTATCATCCTGTTCGTCAGCTTCGAAAGGAGGTGAACCATACAAGAACTCATAGCACAGGATTCCCAGGGTCCAGTTGTCCACGGCATGGTCATGAGCTTTCTTCTCTACCATCTCTGGTGCCAGATAGTCGATGGTGCCGCAGAGTGTGTGTCGTTTAGCATCAGAACGAGCTGCCCATCCAAAATCTGCAATTTTAAGTCGGCCCTGGAGAAAAAAAGAAGCTAGCTAGTTCCTCAAATCAACGCCTAGCATAATAGCATTTGATGCATCAGTCAAGAGTAAAATGCCCATCACAACAACAGCAGTTCTTCAGAAATACTGCAGCATTGAAGAAGCATTTGTAGAACCATTTCTCCACATATGTAAGGGAATAGCAAAAACACGTGTTTCCAAATTCATTCATAGAGGCATGGTATGTATCCTATCCAGGTGTGACTTGTATGCAAGATGAATTATCCCTACCCTTATCGAAATGTAACAACAAATATCATCTGCTTCGTACATGTATACAGTGCCTCGAGTTTCTTCTTTTTTTATGATCAATGGTAGTAGAGAATCCCAACCTGATTATCATTCATTCCTGAGAGTTTGTTCTAAGCTCAGGTTCGGCAGTTTTACAAGAAGAAAAAAAGAAACAGGGAAAGATCAAACACATGGGCTAAGTATTTTAAGTAAACGCCGGGCACCATGCCTCAATGGTAGCTAGATCATCTGTGCTGTGCAGCAATCTAACTTGATTTTCATATTTTCATGTTGTCTATGAATTGACCTTCCAAGACAACGTTTCCGTTTTGACACTTACAAAGAGCCAAGTCTTACTGTCAGCATTGTTTGCATAAACACTGAAATACTTCTATAATGTACACTTAATAAGGGCAAAACATTAGAAAAAATTTCCAACCTATGAGATCAGCTGTAGATTGGATTCTTGAGCAAACTTAAAGTAATGGTTCTATTAGATCCTATTGTATGCAGAGACAAGACCATCAGCCTGCCTATTAAAATTATAACAAAATGGAGCAAAAACTACAGGCACTGCCATTGAGCAAAGTATTTTACAGTAATTATTATGTTAGAATCCATTTGTGTAGAATAAGACCATCCGGTAAACAAAACTACAACAACACGAAAAACAGTCAAAGTACAAAGACAGGCACAGATAAATTTGCTTGCAACAGCCAACAATAAATAATATCAGAATAAATGCATGAATCAAAAATGGTGCTAGAAATAATAACAGCCGCAACAGCTTCATAAAGAAGAACCAACCTCGAGATCGAGTAGCAAATTCTCAGGCTTGATGTCCCTATGTATGATTCCCTTCTTGTGGCAGTATGCCAGTGCCCCAGCAAGGCTCGCAACATACTAAAACACATAAGATTCGTTGGCATCACAATCTAACAAAAGAATGGCAAATCAGAAATTCAGAATGAAAGCTGCATATGTATATCATGGTCTAGAGCGTGGTATGGCAATGCCAGCACATGTTCTTGGCCATTATGGCAAACACCAAAAGTTGAGACAGACAGTGTCATGCTAACATAGATGCCATGAATTCACACATTGTTTCGTGTATATATGTATGGGTTTAGGGGAATAAAGTCATAAACTGCAGAGTGCACTGTATGTGCGGCACATGGGTACCAGTCCATACAAAAAATTGGTTCAGCTTAACTTGTAAAATTGGATTGATCCTTTGCCTTGTATTTCTTGGAAATGCATTGGAGTTTCTAAAGAGTACTTAATACTGCAATAAGATACTGCAAAGATCCATAAAATACAAATATCATTCTGGAAGCTTATAAAATGGCACTAAATTTTTTCAAATATCAAATTATTTTCAACGTCTTTTTAATGAGCTGTAACTTTAACTTGTATCCATTCCACAGAGTGAACGTCTTTTTAATGATATCTGATTCAAGTAAAATCAAGTTGGGCCACTACTCAAGACTTCCTTGCTACCATTCATGACTTGAGTGAAAATAATAGCAGCCTAGTTCTCTGAAACTTAGTTTGTACAGGTGTTTCTCTCCCTTTCCCTTTCGCTTTCCCTCTCTGTTTCTCTCCCTTTCCCTTTCCCTCTCTCTCTCCTCTGTCTTTTTCATCCTGTTTTGTTTTAATCAAGTAATAGTAGAGGGCTCACTCCTCACTGTATTTCCCCTAAAAAGTATCACACAATTGCAGTCACATAAGACTTTGCCTCGACTATTATGGTTCAGATTAACATAACTTCAGTCAACAATAGTAGGAGGAACGGCAGTCCAAGCGCAGACACAATCACACAACCGGGAAAATAACAGTGAAGCAAAAAGGGAAACGGGGAGGGGGTGGGGGGTTCTTACAGTGGCGGCCGTGCGCTCGTCGAAGCGGCCAGCGGCGCGAAGGACCTTGTAGAGCTCGCCACGGGCCGCATACTCGAGGACGAGGACGACTTTCTCTTCGTCGTGGAACCAGGTGAAGAGGCGGAGCACGTTGGGGTGGTCGAGGTCGCGCTGGATCTCGACCTCCCGCCGCAGCTGCGCGTGGAAGCGGTACTTCTCCAGCTTCGCCTTGAATATCACCTTCAGCGCCACCACGTACCCGCTCTGGATCATCATCGCCACGCGCAGCACAGATTAACTCCCATTGGAGGCGAAACATTAGGTTAATCGAGAGAATCGGAGAGTAGAGTTGATCGGGCCGACCTGTTTCTCGCGGGCGAGGTACACCTTGCCGAACTTGCCCTCACCGATGTACCTGCCGATCTCGAAGTCGGACATGCTCCACTCCTCCTCGCGCGCCGCCATCGTAGGCAGGGGTGGCGAATGGTGGAGCGAGGCGCGAGCCGAAGAAGACGGCGACGGGTTCCGTCGAGGAGGAGCCCGAACCGGAGGAGGGGGTGGAACGGGTTTGGAATTTGAAATGGAGAGCGGCGGCCGGCGAATTGAGTTGGCGATCGTGAAGACGGGAACTGCCGGAAGCCGCAGCCCGCAAGGGGGAGGAGAGAGGGGCAGGCTTCCTATGTGGGCTTCGCAGTCTGGCCCACTGTACCATCTCGGCCCACTAGCTTTTACTATCTCGGTCCTCGACGTACGCCTAGACCCCACGCGCAGACTGGCCCAGTGATAACAGACCTTCTGTAGAAGAGAACAGAACTGGAGAAAAGAGCAGGCTGAGCAACACCCGGACGCAATTTGCCGGTGCCATGACGTGCTTCCTGGGCCGAGCCTTCTTTGTTTGACGGGCCGGGCCAGGTCAACTGTTCTACGTAGCGCCATATAAACCAAGAAAAAGGTAATTTGGACTTTCATAATAACTTAAGAGGTCAATTATTGACTTCAAAAAACCTATGAAATATATTTGAAAAATAATCCAAATATTTTTCTATAAAATGTACATTTTATCATCACACTATAAAAGAGCGAGTTTTTTTAGAAAACTCTCACTATCTAATAACCGTTCGATCTAAATATCGGGCTATATGAGCCCTCCGTCATCCATGATCTTATATGATATAAGATGTAAGACTCCTCCTAAATCAAAATTCCTCCTAGAGTCCAAAAACACATGGATATAAGACTCTTTCTAAATTAAGATTTCTCCTAGACTCCAAAACCAAAAAACAATAGAGTTATAAATACCTCAGCTTTTGTTCCCGAGTAATATTAGTTCTTTATCCATATTTGTGTATCTTTTCTAATTCTTTTAGGCCTTGTGGTGGCACAAAGAGTCCGGCGGCCATCTTTTTTTACTGTCGGCGGCCACTTTGGATACGTCGTGATGGTAAAGGGAGGAAGGCCTTTTTGGCAACTTAAATTTTTTTAAAGAAAAAATTAATCTTTTATCACTTCCTTCCTACGACCGGAAAAAAAATATATGCAATCACCACACTATAAAAGAGCGAGTTACTTTAGCGAGCTCTTCTTCGTCACTTATGATCTTATATGATGTAAGATATAAGACTACTTATAAACCAAATTTCCTCCAAAACTCCAAAAATAGAAAATAATAAAGTTACAAATACCTCAGTTTTTGGTTCCTAGATAATCTTGCTTCCTCAGTATTTGAGTTCGACCCCAACCCACATGTATAGCTTTTCTAATTTTTTACATCCAATTCATTAGCACCGTAACATGCTTTGCACGAGAATAGTAAATAGATTAATTTTAAAAAAAATATGGATGCACCGCCACTTGCCGCGGATTGTTTGGAAAAAAATAATACATATCAAGATGCATGCATGCAACTAATAAATAGATTAATTGGAAAAAGGAAAATATGAATGCATGCATGCAACTATAATAAATAAACTAATTTAGAAAAATAGAATAATACGAAAGCACATCAATACGAATACATGCGAGTGTATGTATTACAATTAAAAAGAAAATAATAAAATTTGAGCATACATGATTTTATTTATTATGACAAAAAAGGAAGAAAAAGATAGCAAACATCTTAAATTACTTATAGTAGAATGAACAACAACTCTAACCTGGAGTGAAAATAAAATCAACGAACGTCCAACAATACAAACACGTGTGTGTCGCACGTGCATATTTGCTAGTAAAAAAAATACTCAAATCTAACTAACATTAGAAATTGTATAGAAAATTTTTTATAACTTAGAAAAATATAAAACCAATTCTTTTTTCTAAAAACATGATGTGTATTTTTCTTATAACCTAAAAAATGAAATTAGCTTCATGTTCGCATATTTACTCATAACCTCAGTGTATCTCAGGTACCTAGCAGTTACCCATGGGTATAATTAAATATTATTATTTTTGGTTTCAAATTACAAAATGTTTTTACTTTCTATTGTACCCGGACCCATCGGATATGAAACGTGCAGCAGCTACGCGCACATGTCTAATTGCCATCCCTAGTCTTAAAAGACAAGTCACCAACTCCTGGGGACATGGAAGTTCTTTTTCTATAAATTCTACGTCAATTCATGTATAATACAATTTGTGTGTTGCAGGTAAAGAAAAAATCCGCTTGGGCACTGTGGCAGTAGCACACAAACGAAACATTTCGAATCCCAATTCCCCACACCTCAAAAATTATATAATAAAATATTAAAATTTATAATACAAATAAGTATCATTAGACTCTTTATTAACTATATTTTTTTATAAAACTATTAGCTATATTTTTATAAATCTTTATGTTTTTTTCCATAATTTTTAATTAAACTGAAGTTGGAATGATTTGCAATTTGGATGTGTTGCAATTGCAAAAAAGGCAAAAATCCTCTTGGGCAGCCGTGGACTGTGGACTGTGGTAGACTGGTAGTAGCACAAGCACAGGAAGCGCGCAAAGAAAACATTTCGAATCCCAATTCCCCCCAAACAGCAAAACCCCCTCCTCGTCGCCGTCGTGCTCCGCCGCCTCGAGCCCCTCTCTCGGTGATCGCGATGGCGAACATGCCGGTCTCGGAGATGCGGCTTCCTCCCCACCTCGCCCACCTCCTCGCCGCGCGCCGCCTCGACACCGCCAAGGTCAGCCTCCCGACTCCATCTCCTCCATTGCTGCCTCCACCAAAACCCCACGTCCCCTAACCGCGCGCGCCTCGGGCCTGCGCTCTGTCCTTGTCCAGGACGTGCTGTCGCTGCCGGAGGTGGAGCTCATGGCCATCCTCGACGCCGGCCTCCCCACCGCACGCGCCGCCGTCGCCCTCGTCAGCGAGGCCGCCTGCCCGCCCTGCCAGACGGTACCCTCTCGCCCACCTGTGCTCCGCTTGCGCCCCTGCCTCGGATTTAGTGGGGTTTCGACTCAGTGGTCGATTTGGGGTTGTAGTCTGTGACTCTGGATTAGGGTTTACACTGCATTGACTTGCCCTGGCGCTGCGCAGGCGCTCGCGCTTCTGGAGGAGCGCGTCAGGTTAGGAGGCGGCGGCCGGCTGGCCACCACGCTCTGCGGGTTGGACGAGGCATTGGGCGGAGGAATCCCCATGGGAAAGCTCACTGAGGTCGTTGGTCCCTCGGGGATCGGCAAAACGCAGGTGAGCTTCTGGTTTCAGAGATTCACGCTTAATTTGCGCACCTTTTAGGGTTCTGTGGGTGCACGGATTCACCGTGTACGCGTGTGGGCGCCAGGTGTTTGTGAATATGTACCTACATGCTATGCGTCTCCATGAACTTGGTATACAATACAGAACCAGAAGGAAATGTACCTTGTGTATATGCTTCATGTGACGCTGCATCATGCGCCTGCAAGTGTTATGATTCTGGTTCCTAACTCCCTATGATCTTTGGCATCTTGATTGAGTCCACCGATCAATTGACATTTTGTTGAGCCAAATGGTTATCGATCCGGAACATGTTGTGTGCTAACTTTATTGTACATTATCTCATTCTGTTGCTTGTCTTGTATGTTTTTTAAATTTGTAGTTATGGTTTTGCACACTGCAGTTCTGCCTGAAGCTTGCGTTGTTAGCAGCATTACCGGAATATTATGGAGGTTTAGATGGTAGAGTTGTGTATATTGACACAGAATTCAAGTTCTCTTCACGGAGGTAGCTGAAGTTAATGACTACTGTTTACTACTCCCAAGTATATTAAATATTTTTTTCTTTAGTCGAACTTCTCATGTTGTCATTGATATGTTGCAGGATGATTGAGATTGGTCAGAAAAGCTTTCCTCAAATATTTAGACAAGAAGGCTTAGCACAAAAGGTTCTTATGCTCCACATTCTTATACTTGGGAACGTCACCAGGAAACAATGCTATGGAACCTTGCTAATATTTCTTTCAGTTCTTGTAGTACTGTTCTCACAATCCTGAAATGTCTTGCCTGTTGTTTCGTATCCTTTGTTTATTTTTGCAGATGGCTGGGAGGATCCTAGTAATTCGACCGACAACTTTAGCTGATTTCACCAAGAGGTATAGCTGTTCATGAAATCAGCGAACTTCATTACCCATGCTCTTTCCTTTATTTGCAAAATATACAGACTGATATTTTGAGTATTTTAGTTCTTATAGGATAGGCGCATGAAGGACCGAAAGTTCCTTAAAGTATCTTTGAGATTCTAATACCATGCTGTTGAAAATTTGGAATGCTTATATTCTGCTCATCCTCATGCCGTCATAGTTTAGATCTTGTGCTTTAGGTTTATGCTCCTCACTAAACTCGACATGTTGTGCTAACACTGCCCCCCCTCTTTTACTGTCTTTCTACAGTTTGGAAGAGATGAAGGTGACTCTTCTTCAGCATGATGTGAAGTTACTCATTGTTGATAGCATGGCTGCTCTTATGTCTCTGTAAGCTATGGACATGGAATTCCTTAGTCTAGTTTGCAAATAACTTGCTGCATATTTTAAAAGATTGAAGTGTATTCTTGTATACACTGTTCTACATATATAAGTATCTTCTTTTGCTGGCAACCAGATCATCAAAATATCCAATTTCAGATTGAAGTATGTTTGACTATATAATAGCATAATGTCCAAACTGCTATTTTTGTACCAAAACAAGAGTCACAGCCGCTATGTTAATTGACTAGTTTGTGTGCCATCTCACTGCCTGGTACACCTCCAATGCTGGCTATGTTTAGGAGCATTTCACATGTGTAGGCTTTAGGTTTTACATAAGCTGAGTTTATAGTTCCTGAATAGTTGTGCTATCAGTGGATGGAAAACTAATCCTCTATTACATCCTGTCTGTAACAGGGAGAACGAGAAGGTTACAGCAGGTTTCAGCCAACACCCTTTAAGATGGACCCTTTCTTTTCTTAAGTGTGTCTCCCAGTCTCTTAATGCTCACAAAACTTTCTTTATATTGTTTTTGTTTACCAATATTCATCCTGCTGAAATACCTTTTGTTTCCTCTCTTTCTTCTTTCCTTGTGTAGGTCTATAGCAGAGTTCTCAAGAATTCCAGTCGTGGTTACAAACCAAGTACGCAGCCAAAGTAATGACGATGGTTACCATTTTTCCTTTGAAGGTTTGGACTTCTTCCTATAACCGCTATGTTTTAGACAATTTGTTGAACTTGACTAAAGTATCTGTCAAAATCATCTGCAGTGGACAAAAAGGATGGTAATAAATGTGCTGAAAAGTTTGATTCTCATCTTGTTGCTGCACTAGGGATTCAGTGGGCTCATGCCATAACTGTCCGTCTAGTCTTTGAATCTCATTCAGGTTTGTCCTGATATTTTTCTAAATACAGAAGGTGTTTATGTTCTCTTTTCATCTGCCATTTATGGATCGCAAGTTTGCAACACTCATACCTGGATTTTTGAATTCATATCACATCAGTACCTATTATTAAGATCATATGAATGCTTTGGTGTTTTTTTTTCTACTATACTTCAATTTCATCTGGTTTGATAATTGATATGCTGCCTAGTTACGATTGCATCACGTTAGTTTTATCTTAGAACCTGATTTCAGTTTTGATGTTTGCTTGAAGCACGCAGAAGCACCACTAAACACAAGGAGCTTTTAAAATTGCATTGCCTTTTTTTTATTTTATTTTATCTTTCAACATAGCACTACATCAGTGAGTAACATTTAGTGCTGCATCATAATATGGAAAATTTGACATCACTTGGTAAAAAATGGTAGGATAGAATAACGATGAAAGAAACTTACATCTATGTTTAGGAGTCTTGCTCTTTTTTTTATGAACATAAGTGGCATGAGCTCTGCCTTTCAATTGAAATAGAAAAGAACTTCTAAATTATACTGATACAGGATGTGGCCTAAAAAAAGCCCCAACACTCAAGGCACTCTAGGAGTAACCTATTGGCACTACCAGTACCAGACCATTTTTTTTTTGAACGCAAACCAGTACCAGACCATGGGCCAACACAAACATGAAAACTAGAAATTAATTTGGATACAAATCCAAAAGCCCTCTTACATTGGTCCATGTCATCCTTCACTCCCCAGGCCACTTGCTGGGCCGTTTGAAAGTGGTTTTTGAATATGTGTTGGTTTCACTCTTTCCAAATGTTCCACATAATATAGATTGCCATACCATGGAACTTGTGCCGCTGATCTTCTGGGAATCCTCTCGCTGTTTCCTCCCACCAAGAAATGATGGTTATGCATTCACTTGATTGTGAGGGTAGCACCACATTGCAATGTTCCTAAGTTGCAACTTGTAGCCAAACCGCTTGAGCAAAGGGGCAGAGTACGGACAGATGGTGACCTGTTTCTAATGGTCCGTGACATAAGACACATGCATAAGATACACACGTCCTGAAGAGGACAGCCTTGCAAGTTTAAATTGTCCACAGTGAGGATCTTGTCTTGGGTCTGTATCCAAGTGAAAAGTGTGCACTTATTTACTGCGTGTGCCTTCCATATTAGATTGCTATTGTACTTGATGTAGGATCCTTTGAATTGTATGTGATATGCTGAGCCTGTGGAGTACTGACTATCTGAGGTCCACCTCCAGCTGATTGTATCTTGCATTTCAGGCTGTAAGCCCACATCTTGTAAACGAATCCATAAAGCAACAAATTCTTGGACTTGGACTGTGGTAGTTATCTTGTTGCGGAAACCATCGATCCAAGCATCGTTGCTCATCTATTGCTCAATGGTTTTGTTTTTTCTCATGACAAGCTAAAAAATGTGTGGCATGAGGTTCCTAGGCGCTTCTCCATCAAGGAGGCTTGCTCTGACAGTAAAACCTTTATATGAGTTGCTGAAGATGTGCATTTTCGCTAACTAAAATACTCCCTGTTTGCTATCATATTTCTTCATCACTTACTATTTTACTAACCAAATAGCAATTTCACATTTATCCATCGTATCATTTTAGACAAAATAAATTTAAATTGCATTCGCCATAGGCCACAGGTTCATCAAGGTGGCAAAATCACCTATGTCTCCAGCAGTAGCATTTCCATTCGTTGTTGAGTCATCAGGCATTACCTTACTAAGTGACGAAAGCATTGATGTGACAGGTCCTGAGATTACCTCAATTCGTTGTCAAGGTACAATCTATGTTTGTTCTGCAGATAATGGGCTTAAGCGTGCAGTCTTCTTGTTATTTGTCCTTCTAACACAGATGTTTCTCCAGGTCAAAATGTTCTGTCTCGATAACCATTTAAATGCTACCATATCTACTGAAAAGGTGTAATCCTACGCTTTTTGTTTTCTACTGGACAATATGGATAGTATTAATTTAATTTCAGTGCACTCTCCTCAAAGTTCAGGAGTAGTCTTGCTTCTCTATCTTGTTTTCTCTACATAGTGTACCATGTAAAAGATATAAGAGCGTTTTCAAATTGCTGTGATTATAAAAGGGGTTAAATGCTTTGCAATTGTTTGATTGGGATTATTAAGGGGATGAAAATACTTTTGACCAGAAACTGTGGGGTGAGTTCCCTACAGGAGTTTTTTTTTCTTTAAATAAAAGAAAAAAAAATGTACCAAGGGTGTTTCATGGGCTGAGGTCCATTAAAAGGGGAAAGGACCAACATGATATAATTCCAGTGCTGTATAATTATGGTACCAACATGATGAACTCCAGTAATGATCTGGTTCTTTATTTTGCACTATGGGATTGGAATAAGGTTCACCTGCGTGGGTCAAATTGGTTGATGTCAAATCTGGTTATTGTGGGTGGTATGATTTCTTTTCTGTTCAGTTAGTGGGAACGGTATTCTTGCTGCTGTCTACTTTATTCTGCATAAAAAAACGAAAATTATGTCAATCTTTTCCATTATGTTTCTGTGAAACTTCTCCCATGCAACTGTCCAATCAATTTGAAAAGATCAATTCTTTGTAAAGAAACTTATTTCTTTTTGACAATTTTGCAGGTCTCTATCCTTCTTTGATAAGTTATGTTCTAGTGTAATGATTGGGTGTTTTTCCACGAGGATCGAGGAAGGTTGAAATAGCCGGTATTCAGTTGATGTTTCTTTCCAAGTATCCTTTTGTGTGGTAACCTAGGCACCTAGCAGGAACAAGATGATATCATTCAGGGAGAAAATCAGACATGGCCTGTTTTATGAAACATTGTTTAGATTGGTCCTTCATAACAAGGTAAATACTTGTGTACCAAAAACGTTATGAATTTGTCTGCAACCATTAGCTCAACGATTTCAGAAGCAGATCTTTGTTTTCTGGAAGGGGAAATATCTATAGCACCTGAATTCGAGAACTGTCATAGGAAGATAAAGCTAAACGCTGTATGAAGAATGCTAGTAAACTATGAGGTCATGCTCAATGGTTAGAAGAAAAACGGCAGTCATATACAAGAGAAAAATCGCCTTGATGCAGTCTTGTGTAACTCGTGTCAAAGGTTAAAGCAGTAGTTTTTTACTGAATTGGCTGACTTCACTTTTGTATGTACTCTGAACTTGTAGGTTTCGCCCTACCTGTACCTGATGGCGATGGGATCTGCTGTTTTCGTCTGTCGCGGAATCAATACACGGTATAACAAATTTACTCATGAGCTTTGCAGAAAATTTAATTGTCTGTGTTTTCTTATGCATGACTCTGCGGCCTGATCCATTACTCATTAGATAATGTATAAAATCTTATAGAATCCCTGATTGGACATTCAATATCATATACACAGAACCATGACATTTTTAGTGGAGATGATGTCTGACCTGTGGGTTTAATTTGGTCTGGAGTAGTGGACTTTCTAAGAATGGATATCAATCCTCCACGTCTTGCGATAGCATTTTTTATACCTCGAAAACCTCGTATGGTGTGTGCGCCACATCATGCCTGCTTTCCGTTTCCACCCACAGAAAAGAACAGGAGCGACAGCCCACAGTTGTCACCTTGTGGGCTTCTTTTTCCTTCTAAAAAATCTCCAGCAATCCCCAAATCATATGAATCCTTTGAATCCTTCCAGGCTTTTGTTTACTTCAATCTAGGACTTGCCCTTCTCGCCTATTCGCCCTTCCTGGCTGGAATAGGCACCCTCCTCCTGGAAAGAAAGAAAGAAGCTGCTCACTGGGGGAAACTTTTGCGCACATCATCACACCACCACTAACACACGCTTAGCATATGCTGGAGCCTGGGGGTGGTGTTTTACCGAAGGGACAGTGCTTGCTTTTGGGAGGTTAAAGAAACGGAACTAGGTTCTAGGGGAACCAATAGAGGCGAAGCAACTGAACTGTCTTGTTAGTAGACAGGAAGAATCACTTTACGGGTTCTGCTGATCGCTTGGCAATTATACATAAATTTAAATAAATAAATAAACTTAAAGCCCAATGCTTCCTTTCGTTCGATTTTGACGGAACGCTGTCGACTTTTGGGAGTGGTGGTTTATCCTGATTGGACGAAAGGGGTATGTCTCTAGAAGCAACATAATTGATCATGAGACTCATTATTGATGATAGGTCATGGGGAGCACCAGGATAACACACACTGTGCAATCGCGACGACTGCAAAGACCTAACAACAACTTATTCCCCCACATGCCCTGAAAAAGCAGGGAAATATTCACCATATATATGCGCTAGTGCGTTTCCTCATTTGGCCATATAGGGGGTATTTTAACAATAGGCTTGGCCACTGCTATTTTGAGTGAGGGGACATAATAAGCTTTACCTTGGGGCATCTGAAGAAACATGAGATGCATGAAGACTGGATTCTGTCTGGAAACAGCGACTTTAATTGTTCCTGGGGCAGCTGATTCTCGTAATCGTGTTGTGCCCCGTCTGGTTGGAAATTAACACTAAGGGCTAGCTGGATGTAGATCGAGTTAGGTATTTGACCACGCTGTTATGTCTGATGACGCAGTTTTAACTCATTTTAGTAACTAGCTGTGGGTCCCCGGCCTGAACCACTTGGCCGGGTTGATTGCATGCTCTAGCTACTAGCTAGTGTCGAGTCAGCTAACCAAAGTTACCTTTAATTTCGGTTGGTCGTATTCCGTTGCAGACCTAATCATCATCATAGTACTAACATGGCAGTAGCAGCATGCAGTATGGCAGTAGTACACATTAGTATAAGGCCTTGTTTAGATCCAAAAAGTTTTTGGATTTTGACACCATAGCATTTTCATTTTTATTTGACAAACATTGTTCAATCATAGACTAACTAAACTTAAAAGATTTGTCTCGCGATTTATAGACAAACTATGCAATTAGTTTTTATTTTTATCTATATTTAATGCTCCATACATGTGCCGCAACATTCGATTTGACGGGAATCTTGAAAAGTTTTTGGTTTTTGGGTGAACTAAGTAAGGCCTAATTAAAAAAAAGGCCCAAAAAAAAGGTGTTACCAGGTTTCTATCCTAGGAAGGAGTCAAACATGGGCGTCGCGTCTGGTAGTTTGACAGGCTGGGATTGGCTAGTAGTTGCTGAGTCCAATGTTTGTATGGATACCAGTCGTGCCTAGCTAGAGTAGACATGCGTGCCTGAGAGCGAAACGTGAACCAATTGGATTGGATGGTTCGCTAGAATAAGAGCAGGACCATCCCCATCCGTAGGTAATAAACTCCTTTTTCCACAAGAGTAGGTGTTACAAGTATTACACCTCTCCCTCTCCCTCTCCCATGAAATGATTCATCAGTTGCATCTCGATCCCCTTCCCTTGCCTTATCCGCGTGGTCCTTTCCCTGGCCTTATCCGTGGTCCGTGCGGCCCCAACCCTGATAGCCCCAGCAGCCTCCAGTCCAGGCGCCGGCCTTGCTCACCCCACTCGCTACGTCCCAAAATGCAAGCAACCCGGCCAGCGCAATGCAATACCGTGTGTGTGATTGCGTTGCCTGCGGCGCAAGTTGCGTGGCCCGCCGGTGGATGTTGCCAGTGATTTGTCTTTAAACAACATGAGCTGACGCTTTGCACTGCATGGCATTCAAGCTAGCATTCGTTTGTTTGGGATGATAATGATGTATAATAATAAAGTAGCAATAGTATTTAAAAGTACCAGTAGCATGTCGCTACGGGGCGTCGTTTCCAGAGCGTCCGTCTTACTCTATGGAGCTGATCAATCAGTTCAATTCTGGTTCTGGAGGTAGTGGTAGATAATAAAGCGCGCGCCCAGAATTAGCCATAATAATCATTTGTTTAGCGGATGGATAGATGGGTGGATGGATGGACACTGGGTGAGGTGAGGTGACCGGCAAACCATAATGTCCACAAACAAAACATTAGTTTTCTGATTGGTTGAGCAAGTGGCCCGTGTAGGGTATAGAGTATCATTTGATTTGGAGTCATCTGTTTTTTAATTTATTCGGAAAGCTAGCGACGGACAGTTGTGATTTGATTGATGTGTTGTGTGTGGGCTACTATTACCAAACTGTTTCTTCCAGGGTATTCTGATGGTTTCCAATCCAATAAATATCATAATACTAGGAGATTAGCATAACATAGTATAGCTAGGCTAGGCTCCCTAGTAGTATATGTAGTGGTTGGTTGCCTAGCGTTGCCCGGCCACAGTCAGTACAACTTGCGCTTTGTCTATGCATCCTACTGTGTCGGCACTAGCGTTCTTCTCTTCCTTCCGAGCTGGACCCTTTGCATTGGCGATTGCCGGTCGTGGTTTCTTGTTTTGTTTGGTTCACTTTGGCTCTAAAAACTCAGGGATGGATTGGGGATTTGTCCCGCCCGAATGCAGTAATAGTGGTTTGCCTTAATTTGTTTTTTTTTTCAAGTTGGAGGCTTTTGCCCTATGATGACGACGACGACGACGACTGTCTGGAATAGAAAAAAAAAGAACTCTTTTTGGAAAACGTCGTCCATGGTAGCCAGGTGGTGTGTGTGTGCGCGCGCGCATGCCGCGTGTTTGTATGTACTGTACGTGCGTGCCCCGTGACTGGAGCGGCATGCATCCACCAAAGGATGGGAACCGCCGGATCTGAATGTCCGGTGCCGGTACCCCCTCTCCCTCCTCTGCCACTACTGCTACTACTACTAATCATCATGTAGTGCTCGTGTGGATTATGGATCAGGATCACCCTTTCACCCTTCCCTCTGCCTGCCTCGACTAGAACTCTCGAATCCTCCTCGTGTATCTTGTCTGTCGTTAATGGTACCAGCCTGTCTGTCCAAAAGTCACGCTTTTCCATCAACTTTTGCTGTGCTTTTTTTTTTCTCGTACCCGTACTAGTAAAAAGGTCGTGGTGTCTGCCGGCCTAGGTTGTAGACCGCCGCAATTGACGACTGCAGGTAACGCTCTAATCCATGGCGGCAAATAATCTCATCACTCCCATCAGCGGCGTCGTCTCGGTCGTGATCGTGACAGGAGGCGGGATTCGTCGTTCCTCGTCGGTGGGTGGACAAGACGATCCCCATTCTGGGAAGTATACAAAGGTAGTACTACTACTACGTACTCGATCCATTCACAGGAGGCCACGAGGGATCCTGTAGCGCTGTCGAGGCATCATTTCATTGGACAATTGGAGATCCCATGGTAGTAGCATGCATGGTACACTTGCGTCTGCACTCTACAGATTGCGATCGAGGCGAACTACGCGCATGGCCTTGTTCAGATGGTAAAAATCATTGGTTTTAGATACCGTAGCATTTTTATATTTGGTAATTATTGTCTAATTATAGATTTACTAGATTTAAAAGATTTGTCTTGCAAATTATAAAAATAGTTATTTTTTAATTTATATTTAGTGTTCTATACATGTGTGCAAGATTTGATATGATGATAAATTTAGGGAACCAGTAACGCCTGTCACTGGCGGTGGACCACGCGAGCGCCCGCCGCTGGCCCGTGCCTGCCTGCCTGCCTGCCTGCCTCCTGGCTCATCCCATGGGCCGGGCCGGCCGCATCCGCACGCGGCACGCGGTTGTGAAAACCTGAAAGCGATCGAATCTCGCTCTCCTATCGACTGAGTATTCTTCTTCCTCCTCTTTTTGTCTTTTTTTTTATTACGACCGGTTTGTACGTACGTGTACGTATAAGAAAAAAAATGATGGGGGTGCTCCTAGGTTTTCTAAGATGTGGTTGCTACTTTGCTTGCTCCTCACCCGCCGTGGGTATCTTGTTTTTCATGCCGCGCACCACAGATAAGGGCCCCTTTTGGATCGCGGCATCCTTTCTCGTTTCATGCATTTTTCCTCTGAAAATAAACTGATTCCTGTAAAAATTAAAATTCGTGTGTTCTAATTTCTAAAGGGGTCTGGGTGTGATTGGTTGTTTACTCTAACTTTGAGCTAGACTAGTTATATTATATTTTTAGAAGGAAAAGGTCAGCAAAAGCTTTCCTGAAATTATATTAGTAGATGAAAATAAAACAAAGTTATATACAAGGAAGCAGCGGACCCTTTAGGAAGCCAGATAGTCATTATATTATATTACTTGTATCTTAGTATTAGTCTGCCAGCGAAACGACGTTGATTGTTTGCCTGCAAGACAGAGATATGTTGTATATGTGAGATAAAAAAAAAACAATACAAGGTTATGGCCATGATCTTTTTTTATGAATTATTACCATGATTTTATTTTTTTTATTATTAAGTGTAGATGTATCTCAAAAGTGCTTATTTTATCTTGCGCATTTTGTGTCTTAAAAGTCCTTTCAAATATACTAATATATCAATTGATATATACTAATCGTGGCGGTACACTTGTGCATGCGGACGCCGGCTCCCATCTGTGCCTGCCTCTCTGCCGCAGTAGCATAATAATATGTACTGGCACAACGTTTGTTTATTGGCAACTTCTCCGTGTAACTTGGCGGCCATAAATAAATAAATAAGTAAATAAACCTCAGTTTCCTAGAAGATTGGTGTATGCGGCGCGCGATCTTTGTATGTGACCACACTCAGAGATGCATACGCAGATGCTAACCCGCCACCACCCGCATCAACTCCACACTGGCTATGTAGTTGTTTACTCCCTTTATTTTTGAAAGGTTCGACTCCTAGATTTCATACCATTTATTTTTTTTTTAAATTTGACCAACTTTATATAAAAGATCATCAATATCTATGATATAAAATGAGTATCTTATAAAAATATATTCTATGATAAATCTAATGGTATTAATTTGGTACTATAAATCTTAGTGTGTTTTTTAGAAATTTTATCAAAATTTTAAATGTTTGACTTAGAATAACTCTAGAAGTTGAAGCTTTTAGGGAACAGAAGAAGTATTACTAAGCCTCTGTCCAAAGTATAATACAATTCTAAGTTTGAACCAAATTTGAACTAAACAAGGCCTAAACTATATAAGCCATTCTTTCTTGTAATCTTGTTAGATACCTAGCTTTTATCGTATACTCCTAGATATCATTTTAGACAAGGTTTGAGTCAAACATTAGGAATATAAATTTCTCAATAGTTTTTAAATTGTTGAGTTTGCAAATATAAAATTATATAAATAGATGTGTCTTGAAAATACTTTTATAAAAATATACATATATTATTTTTTCTAAATATTTTTATAAAAAATAAGAAATTAAAGTTATTTTTTGAAGAACATATCGCTGTGTTAAACGACGTCTAAATCTCGTCGACAGAGTAGAAAGAAAAACACGGCTGAATGATTTACAAATCAATCTGGCTGATAAGATAGGGTCCACACACCACAACAGTGAGCGACGTTGAGATGGGCATATGATGGGATGGGACGAGAGGGAAATCATTAAAAAATTCACAAGTACGGTGGGCCCCGTCCCGCAGTGGCGTCGCTGTCTGGACTCTGCACCGAGCCTACTCCTCGCCCGCCACCGGCTGCAGCCGGTCAAAGGCGGGGGCCGCTTTTCCGGACTCGGAAGGTCCGACGAGGATCCGGCACCCGATCTCCCTCCATCCCAGGTTACCCGATCGGGAAATTCAATGCAAATAGTACTACTCCAGTAGGAAGCTTGCGGGTTTCGGGTGTACGATTACGAACGAACCCGCTATTGTGCTTAGTTCACTGCCTTGTTTAGATGCGAAAAGATTTTGAATTTTGCTACTGTAGCATTTTCGTTTGTTTGTGGTAAATATTGTCTAATCATGGACTAACTAGGATCAAAAGATTCGTCTCGCGATTTACAGTTAAACTGTGCAATTAGTTTTTATTTTCGTCTATATTTAATGCTTCATGCATGTGCCGAAAGATTCGATGTGATAGGAAATCTTGAAAACTTTTTGGTTTTTGGGTGAACTAAACAAGGCCTAAAAAACCAAGAACTTTTTAAAGATTTCTCATTACATCAAATCTTACAACACATACATCTAGCATTAAATATAGATAAAAACAAAAACTAATTACACAGTTTACCTGTAAATAGTGAGACAAATCTTTTAAACTTAGTTACTCTATAATTAGACAATATTTATCAAATAAAAAAGAAAGTGCTACAATATCAAATTCTAAAAAGTTTTTGATGGGTAGACGAAAAGTTTTTGGGTTTCATTTGTATTTAATAATTATTGTCTAATTATGAACTAACTAGAGTTAAAAAATTTATCTACAATTAATTATTTTTTCTATTTATATTTAATACGAACCTTCAAAAGTTTTTAGATTTTAGTTGGAACTAAACAAGGTCTGATTTGATTAGAGGATTGGATCCGGATCCCCGTGGGCCGCAGCGGAGCTGTCTGGTCACTGGTGGTGGTCTTTCGTGTCGGCGGTTTCCGCTGTCAGCTTGCTCGCTCGCTTGTCTTTTGTTTTGCTGTAGACAGATTCATGTTTATCTTCGGTACATGCCTTGGCGTCCACACAAAACACACGCACAGTAGGCTCACCTCACCGTACCTTTCAAGCACCATCGTCACCATGCTACACCTACATACAATGACAAATAGATTAAAACTAGGACCTTGTTTAGTTCATCTAAAAACCAAAAAAAATTCAAGATTCTCCGTCACATCAAATCTTGCGGCATATGCATGGAGCATCAAATAGATGAAAATAAAAACTAATTACACAGTTTGCCTGTAAATCGCAAGACAAATCTTTTAAACCTTAGACAATGTTTGTTAAATAAAAACGAAAGTGTTAAAGTATCGAAATCCAAAAAAAAATAGATTTAAACAAAGCCTACGTAATATGTACCACATGGGTATTATTCATGGTTGAAACTATATACATACATGCATGCATGTCTCATCAATTTAGGTACTAGCTTGGATAAAGAGATATAGTCGTCGATAAGCTCTCCAATTAATCTCGTGAACAGTACACTTCATCACCAGTTCGTACCGTAATAAGGATCGGCGATGGACCGATTGACTTTGCTATAGTAGTGTGTATGTGTGTATTGATAAGAAAAAGAACAAAAAAAACCCGTATGAATGCCGATTATATAACACTGATTTGTTTATTTGCAACTATACGTGTTGGCTGAAAGGACCGGCCGAAGGACGGACACGTCCGCAGAGGCTTTAAATTAATTATTTCCGTTATTTTTATGAAACAACATCTATCGAAGAAAATATACGTACTGTCACGGCAGTTTCAGAACCGAAACAATCATGTCAATGTTGGAAGAAAATGTTGCCAATGTACTAGTAGCTAGATGATACGGGGGAAGTGGACGCAACGTAACTTTATTTTCACAGCAAAAAAATGCGAGGAGGGAATAAAGATAATGAAATTTACCGGTATTTGAGGTATCAATCCTAATACACTCCACTGTAGTCTCTCAGATGGTTTTTCCAAGAAGCATTTAAGGTGTTGTTTCGGTGAAAAACTAGTGCCGAATCGGTCTTGCCGATCCTATAGTCGGCCTCACCATTTAGACCTACCTCTACTATATATCGATCTTAGCAATATGCGGTCTATGCAATTGTATATTTGTTAATGTGTAACTAAATTAACATGATCAACGACTCTAGATCTAGGGTTGTGAAGTGAAGGTCTGCCGCCATCCGGCATCTTTTTTTTCAGTTTTCTATTAGTATTTTTATCAAATAAACAGTAAAGATAGTTACTATTTTAATTTTTTCCTCTGAATATGCTATCAAACCTCGAATATATATATATATATGCTTCGCTGAATGATATTTCCAGTTTCGAAACCTTATTGTAACCACGAAATAAACAATCCTTTCTTTACAAAGCACCTTTATTTCATTATCTCATCGATATCTATGGACAGCTATGAGATGATATCGATCGGTCGGCCAGCCGGCCGCTCCATTATTGATGGTTATGTCATGCATGCATATATACGGATGTCCCATGCATATATACGACCACATAATAATAAACAGTATAGTATATTGTTTAATTTCCGCCTCGTGCGTTTCTCTTCCTCTATCGACAGGTACACCTTAATTAATAAATTAAATATCATCATGTGGATATAATTCCTGCTATTAAGGATTTAAACAATTTGAATTTGACTAAAGTTATATACAAAAGTACTATTAGTAACTATTCAAAATAAATATTGACTATGGAACATATTTTTCGTAGTGAATTTTATTGAAGATATGAATGGTGACACTAATTATTACTTCCTCCGTCCCAAATTATAAGTCATTCCAAGAATCTTGGAGAGTTAAAGTTTTTCAAATTTGACCAAATTTATATAACAAAATAATAACATTTTTGGTACCAACCAAGTATCATTAGATTCTTTGTTAGTTATATTTTCATAGTGTACCTATTTTATGACATAAATCTTTATATTTCTCTCTATATTTTTGGTCAAACTTGAAAATGTTTTGACTCTCCAAGATTCTTGGAATAACTTATAATTTGGAATGGAGGGAGTACAAATTTGATCAAATTTTAAAACTATGCTGACACACAAGTTTTATAGTTGCATTTTTTTTTGTGGTTTGAGGGAGCATATATGTAACAATAATATATGGTTTTGTGTCTCTTTACGTTTTTCTCAAGGTATACGGATGCATGCATATCCGGTTCCTGTTGAATTCGATTTTCTTGCTGCTTTTTGCAAAGATTCCGTAGCAGTGACACTGTAGCTTTGCCCACCGCGAAAGGCACTGTTACTTCTGGTACCTACTACTAGCTACCGGACAGAATGGACGGTAGCACGTGTCTTCACTGACCACGGATACCCCGGTTCAGCGAAAAAAAAAATCTCGGAACTATTTTTGAATCCTCAGGGCAGTCCCAATGGGGAAACTAAATGTAGTTTCTATTCCTATTAATTATGTTGACAACTCAGCATTTTGATGATGTGGCACAAGAGTTAATAAGGAAAGAGACAGTAAAACGAAGAAATCGTTTCTACTCCACGAAACCACCTCTACCCTCGCAACGAGGCGGATAGACACGACCCATTCTCCGTTGGGAAGAAACCACCGGTGTCGAGGGATCCCGTACGCCGCCGCGCTGCCTGCTTTGCGCTCGCTCGCTCGATCCCGCCCGGCTTCGCTTCCCTCCCTCTCCCGCGCCGTGATGGGAATTGACGACGGGCTCGAGCTTCACGTTGCCGTAGGAGGAGTGTGGCAAGAAGTGGTGGTGGTGGTGCGGCGGCGGGGGCTTGCTGTCGACGAGGACGCCAAACTGCTCAGAGATGACCGGACCGCAGCCGGTGGTGGTGAACGAGGGTGTGGTGGAGGTGGAGGTGAGCACAGACACAGAGGGCGGCGAGGCCGAGGGGTTGCCGCCGCCGCCTGCGGTGTCGGCGTCTGGGAGGTTGACGATGGTGATGTCATGGATGCTGGAGCGGCGCTTGTCCTTGCCGCTGGAGTTGAGGCGGATGAAGTACTTCTGTGCATGGCTGGGCACCTGCGTCGGCGTCCAGCTCGTCACGAAGTTGCGCGAGATGTTCCGCCAGTCCCCGCGCCCGTACTTCTTCAGACCCATCAAGAACAGCCTGCGGTACGGTAGGCAACCAGAGAAAGATTGCAGATTAGCCATGAAATTTTCAGCAATCAAGGGCGCAAGTTCCATTTGTTGTGCGAGAGAGTGAGAGGATCAATTGGGGTCAGGGCGGCTGCTTACTTGTGCTCCTCCTCCGTCCAGGGGATGCATGTTTATGCTGGACTGCAAGTCATTGTTTATTGTCGTCTTTCTAAATTATGTTATTCAAATATTGTCTTGTTTCCTATCATGTGATGAAGCCGCTGATGACTACATGGCTGATTGCCACAGCTAATCTAATCATCTTTTGTGTTGTATTTAATTTATCTTTTGTATTAAATATCTTATAGAAACCATGGGTTGGGAGAAACAGTTTCTATAAGACATTAATTGCTCTTTTTCTCTCTTGGAATCTACAACTAGAAATTTCCATTGGGACTGCCCTAAGATGTACTTTGTTATATATAACAAATGCATTTTTAGATATAAATGCTGGTATTGGTTTTTAAAACTAACCACTTCACTGTATCTTTTCTTTTCGGTCTCTAACGATTTCCCTATACATTTTGGTGTATTCCAGAGAGCACCTTGCTCTCTATTTTTGATTAACGAAAAATCTGAAATAAAAGATTACTACTCCTCCGTTCCTGAATAGAATTAAACTAATTATATCAAAAGATTTGTATTACAAAAGTAGTTTACATTAGATTGAGCTTAGAAATACATTCTAATGATCAAACTTTTGTTGCCATAATTTTTTTTCTCGAACACGCAAGAGAGTTGCGTATCATTGTATTAATAAAAGAAAAGAGTCACAATAGACCCAAACCCACACACACTACATGAACACTGAGGAATTATCGGTAAAACAAAAGGCTAACCCGACATCTATAGCAGTCCTTCAATCTGCTCCTCCAGCAATGAGCAGGGAGATCCCCTTGGCTCCTGCCATGCTCCACCAATGGGCCTCTTCCCTAGCCTAAGCCAGGGTATTTACACTTAGGATTGCACCATTGAATACACAATCATTTTCCATGCCGCCAATGCTCCACGCTCCCAAAATGACTAGAGAATTCAGACCTTTTCTCACATCTCCACTAGCAGCGTTATCCACTTTTCTCCACCAATCATCAAAGGACTGATGAGAAGCTTGTGGAGCAAGTGATGCAAAACCAAACTGCAAAAGGAGGGTGAACCAAAGATCTCTTGACAATACATAGCCCACGACCAGGTGCTGAATATTTTCCTCCTATTGATCCCAGAGAAGGCATTGTTCTAGATGAGGGAGCCACTACTTGCCAATCTTTCCAACATCTGTTATGAGCAACGAACCACATAAAAACTCTGCATTTGGGGGATGCTCAGGATTTCCATATTCTCTCATAAGGTTCAAAGGAGATTGATTCATGAAACAGGGTGCACAACAAAAATTATAGGTTAAAATATGTTTTTTTCAAGACCTCACCAAGTTAAATCGCGTCTTATATTTGAGAACAGAAGGAGTATATTTGAGTATCTAATTAAAGAAACTGTTAAAGAATATTTTTTTTTACCAAAATCTCTTTGCATATACATCAATTAAGAGGAACTTGAGAAGGCTGTTGGAGTTGCTTTTTTTATTCAACTATCCTTGTCCAGTCGTTCATGACATTATGTACCGTGATAGCTGATTTTGATTAAAGTTTGATTTTTGTTTTTGATTTTTTTGTTGCCGTGTAATATGTTGGGTGGTAGTTGGCGCGTCGAAATCGAGGCGTGCACGGACGAGCAAGGACGACGGCACGTGAAACACACGCAGTACTGCTGCGCGGAAAAGGGAGGGGGCCACGTACATGCCTGTTTCCATTTCCAGCGGGTCGTCCGTCGCCATGCATTTAAGGCTTTGTTTATAACAATAACAAAAATGTTTAGAGTTTCGTTTCTGCTTCCTCCATTCTAAATTACAAGTTATTTTAAAATATTAGAGAATCAAATCATCTCAAACTTGACTAAATTTATATAATTAAATAATAATATTTTTGATACGAACTAAATATCATTAGATTATTTATTAATTATACTTTTATAGTATACCAAGTCGATGTAATAAATTTTTGTAATTCTCTTTATAATTTTAATCAAACTTAAAATGTTTTGAATCTATAAGATTTTTAGAATGACTAATAATTTGAAATAGTGGGAGTATTTATTATTATTGTTTAGCTATAGAATTTTTTTTTGATAGTCGTTCAATCATAAGAAGACGTTGATTAATGTTCTGGTAGCACACATCGCTAGACCCACGGTGGTCAGACCTAGATTTCATTGTTCCCTAGGTGGATCGCCTGAATTTGTATGCCTTGGTGCATGCTCTCCTTATTTGTGGCAAGATATGGCGTCATGCTAGCATGAGTCGATGCCAATAGTGTGCCATGCAGTGGGAGGCCTCTCTAAGAAGGCCATACCTCCTTCTGTGGGTGGTCTTGGATGTTAAGGTTGTCATACACTACTTTCATTCTATTTTTTCACCTATCTTCATCAAGGGCTCATCAGTGTGGGTGAATGGAACTTCACGGTTAAGCGTGCTCTTTCTCGCGTGTGGCAGATGTCCTTGAATCCTATGTATGATGCATCGTGTCCTTCTCCTAAGGGTGACCCATGGAGGGCTCCATCATTGAGGTGTGACTGCGTCTAGGTGAATCGCACAATGGTTGTGCCGCCAAAGACTCTATGTAGTAGGTTGCAGGCCCGGTGGCACAACGACATCAACAGAATGTGGTTTTTCCCCTTTGGTTGTTGGTGTGAGGGAGCTATCAGAAAGCATTGCCAAATGTTATCTGTGTTGATGACTATGATACCCTTTGGCGCCATTCACCTCCTTGGAAGTGTCACCAAGATACTCTGTCTAGATCTGCCATAACCACAAGAGACCTTCGATGGAGACCTAGTTTTAATCCTACCAAGCAGGCAATGTACGCTTGTGTAACACCCGTGGCCGCCGACGATCGCCGGCGACGTGAGAGACAGAGTGGATCGCCGAAGTCGGCCTGGCGAGCAGAGGAACCTCTGACCTCGCGGTTACTGTATTCAGTAGATAAAAATGTCGTCCCCTCCTTACAGAGAGCGCTGGTCCTTTAAGGCTTACAGACTGGCCACTCAGGCCCACAAGTCATAGTCTATAATATTACATGCATAATGGAAATACGTCTTCTCTTCTCAACGACGCCTTTCTCCTGTGTCCTCATCCTCACGCCGTGACATCTCCCCCGGCTGAAGAATCTGCTTGTCCCCAAGCAGGTGCAGATGGAAAACGTCGCTTGACCACATGATAATCTTCCCAGGTTGCTGACTCCTTGGGTAGGCCTGTCCAAGTTAGCAAGACTTGTGGAATGGCACTGTTACCTTTCTTTACTAGACGGCGATCAAGTATCTGATTAGGTACTGCTTCAGTAGCTTCGATGTCTGACACAGTGGGTAAGGTTGAGAAGACTGGAGTGTGGTCTTTGTGGAAAGCTTTCAGTTGGGAAATGTGAAATACTGGATGTATATTGCTGTCCGCTGGCAACTTCAGACGATATGCTACTGTGCCCACCTTCTCCAGAATTTCATATGGGCCAAAATACTTGTGAGAGAGTTTGGGAAATGGCCGATTTGCTACTGATGATTGAGTGTAAGGCTGAAGGCGTAACAACACTGAATCTCCCACTGAAAAACTTTTGTCAGTTCTGTTGCGGTCAGCAAATAGTTTCATTCTGTTTTGAGCTTGTTCCAGTTTGTTCTTGAGATGCTGGAGGTGAAGTTCTCTGTTATCGATCAGCTCAGTCACTGATGGTGGTGCATTGTCAGTATCAGTGGGTGGTGTACCCAAATCTGGATCATAGCCATATAATGCCTTGAAAGGTGAGCAACCAATTGAGCTGTGGTATGAAGAATTGTACCACAGCTGGGCCAGGGGTAGCCAGGAGTGCCAGGATTTAGGGGAGTCTTGAACTGAACACCTGAGATACATTTCCAAACATTGATTGACTCTTTCAGTTTGGCCATCAGTTTGTGGGTGATATGCAGTGCTGAGGTTCAGGGTGACTTTGTAAAGCTTGAACAATTCTTTCTAGAATGAGCTCACAAAGATTGTATCCCTGTCTGACACCAGTGAGTGAGGCATACCGTGGAGTTTGACAACTGTGTCCAAGAAAAGCTTGGCTATGGAGGCCGCTGTGTAGGGGTGTTTCACTGGAATAAAGTGAGCAAATTTTGTCAGTCTATCCACTACCACCATGATGACAGAGAACCCATCTGACTTAGGTAAGCCTTCGATAAAGTCCATGGAGAGGTCTCTCCAGATGCCAGATGGTATGGGTAATGGTTGTAAAAGGCCAGCTGGGTGAGTATGTGAGTGCTTGGCACGTTGACATGTGTCACACTGCTTGATGTAAGACTCAACATCAGACCTCATGCCTTTCCAAATAAAGTGGCGTTTGAGGCGGTGATAGGTGGCATTGACCCCGGAATGACCCCCAAGAGCACTGGAGTGGCAAGCAGCTATCAATTTGGTTTGGAGGGCAGAATTGTTACCAATCCAGAGAAGATCATGGAGCTTGATCAAACCTTGCTGGAGACTGTATCCATTGGCATCAGGACTCTTGACAGACAGTTGTGTGATTAAGTGCTGAGCATGTGGATCTGTTGTATAAGAGTTTAGGACTTCTTGCATCCACTCTGGTTTGACCACTGAGATTGTATGCAGGGCCAACATTGTTGCAACTCGAGACAGGGCATCTGCTGCCATGTTCTCTTTGCCTTGTCTGTATTGAATTTTGAATTGGAGTCCCATCAATCTGGCCATTGCCTTCTTCTGCATCTCGGAATGAAGGTGCTGCTCTGACAAGTATGCAAGAGATCTGTGGTCAGTGAGAATAATGAATTGTTGTCTTTGCAGGTAGTGTCTCCATTTTTCTATGGCCATAATTAAGGCTAAGAACTCCTTCTCATAAATGGAAAGCTGCTGGTGTTTACTGCCTAGAGCTTTGCTGAGATATGCTACTGGTCTTCCATGCTGCATTAACACAGCCCCAATACCATCAGCACAAGCATCAGTCTCCACTGTGAATTGGGCTTGGAAATCAGGCAGAGCCAGTACTGGAGTGGTTGTCATTGCTGTTTTCAGGTTATGGAAAGCTATTTGGGCTTGTTCATTCCAAGCAAATTGTTTTTGTCTGAGAAGTGATGTAAGGGGTTTGGCAATCACTCCATAATTTTTGACAAACTTCCTGTAATAGCCTGTTAGTCCCAAGAATGCCCTGAGCTCAGTGAAAGATTGTGGAACTGGCCAATGTGTCATAGCTGCTGTTTTGGCTGGATCTGTTGCAACTCCTTTAGAAGAAATGATGTGGCCTAAATATTCCAAACTGTGCTGTGCAAAGGAGCATTTGGATGCCTTTAGGTAGAGTTTATGGGCTCTCAATGTTTCAAACACCAGCTGCAAGTGTTGTTGATGATCTGCCAAGGATTTACTATAAATGAGGATATCATCCAAGAATACCAGCACAAACTGCCTCAAATATGGTTGTAACACTTGATTCATAATGCATTGAAAGGTTGCTGGTGCATTTGTCAGTCCAAATGGCATTACTTTGAACTGGTAATGTCCTTGGTGAGTTTTGAATGCTGTTTTGGGTTCATCAGCTGGCAGCATTCTGACTTGGTGATAACCAGACCTCATATCCAATTTAGTGAAGATTTTGGAACCCTGCAACTCATCTAGTATTTCCTCAATAATGGGCATGGGAAATTTGTTCTTGATAGTAATAGCATTGAGTTTCCTGTAATCAATGCAGAACCGCCAACTGCCATCTTTCTTTTTGACCAGGAGAACTGGAGAAGCAAAGGGACTATGGCTATGTGTGATTAGTCCATGTTCCAAGAGTTCTTTAACTTGTTTTTCGATCTCAGTCTTGTGTTGAGGTGAATAGTGGTAGGGCCGAGAGTTGAAAGGTATTGCATCAGGTACAAGGGGAATGGCATGGTCATGGCTTCTTTGAGGTGGAAGTGTGTGGGGTTCCTGAAAGATGTCGACATATTGGTGCAGGAGCAGCTGTATGTTTTCCTCTGAGACTCTAGGTACTGTCTGGGGAATTGCAGAGGACTCTGTGTTTGTATGCAGCAACACATAGGCCCAAATGTCATTTCCCTTTGTTGACTTGTATAGGTGCTGAGCTGACATCCTGGTGACCTGAGGTGGGGGTAGTTTAATGCCTTGTAAGGTGACTGGGACACCCTTGTGCTGAAAGGATATAGTTTTGTTTACCCAATCACAAGCCATTGGACTGTGCTGTTCCAGCCAGTCATAACCCAGAATGGCATCATAGGGCAGCAAATCCAGGACAATCATGTTTGCAGTGAAAGTATGTCCCTGAATGTACCATGTGAGGTCCTTCACTTGGGTTGCTGCTGTCATCCATTCACCATTTGCCAGTTTGACTCTTTGTTTGGGAATAGGAGTGGGTATGAGTCCTGCCATCTGCACAAAATGGGCACTCACAAAGCTGTGAGAACTGCCTGTGTCCAGTAAGGTGAGCATAGTTTTGTTCTTGACCGTTGTCTTCAACTTGATACAGTCTTTAGATTCAGTGCCAGCCATAGCATTTAAGGACAATTGACAGAAGTCTTCATTTAACAACTCTTCTATTGCTATTTCATTCAGGAGGTCTTCACTGAGTTCTTTATCTAAATCATTGGCTACCAGAGCATTGAGATGGGGCTTGGTTCTCTTAGGACAGACCTCAGCATGTCCAGGATCAAACTTCTCCCCACAAGAGTAACACAGGTTGTTGGCCTTCCTATAATCTCTTAGCTGTTTGTCTCTCCACAATGTGCCATATGTTGGATTAGGTCTAGCATCTGGTCTCTGGTATTGTTGCCTAGGCTGGTTAGGTTGCTTCTGGTATTTGGCTTTCCCCCTTTCCACTTGCCTCTGTTGAATTTTTGCTATTGTTGCTGCTTTCTTGACAGTAGTTGGTGTGTGGGGCTCCACCACTGCTCTAATATCTTCTTTGAGGCCATTGATATAAGCAGTTGCAAAGAACTGCTCATCCACTTGGGATCCTTCCATTGAGACCTCAAATTGAAGAGCTTGGAACTGTGTTGTATACTCCTCCACTGTACTGGTTTGTTTGAGGTTTAATAAATCTGTCACAGCTGATCTATAGTCATCAGACCCAAATTCCAGATGAACAGCCTCACAGAATGCTTTCCATGAAATCTTCGGACGATGCTTCTTGTAGGCTTGCCACCACTTGGCTGCATTGTTTTGGAGATGCATTGCTGCAGCTGTAACCCATAACTTCTCAGGAATGGAATAGATGCTGAAATAATTCTCACAATTGTCAAACCAGATGGTTGGGTGAACACCATCAAAGGATGGAAACTGCATCTTGGGTAAGGAGTGGTGTGGCAGGTCACCCTTCCCATGGTTGTCAACTCGGTGCTTGTGTTGGTGAGAGCGAGCTGGCTTGTATGGAGCCTTGCCCTTGCCAAACATGTTATGAAAGTCTTCATCTTCAAGCACTTCTGAACTTGATTCCTCACTGACAAGTGATTCTTCTTGATCAAATTGCTTCAAAGTTAGCCGGGCTACGGCCTGTCCATTAGCCTTGACCTGTTGGGAGATGAAATCTTGCTGTGCAGTGAAATTGTCCACCTTCTGATTGTTTGCCATAACCTGAGTCTTCAAGTCCTGTTGTGCCAGATTGATATCATTCATCCTGTCAAACAGAAGGTCCACACTTTCCAGCATCTGATCCCAACGAGACTTGTCATCGTCGCTTTTTCTCGACATCTCCTCCAACAAAATTTGAGTCTGCTTGGAGGGCTTGACAGATGCCGTAGTGGCACTCGGTGGAATGCCTGACGACTCCACGGTGAAGCAACAGGCTTGGTCGAAACCTCCCCTTTTGCGTACTTCACCAAAACCCAACTCCTGGAAGTAAGGCGGGATCGAGGGATTTTACACAGGAATAAGTATCCTGCAACCCAAATCCACGGATACTTCCTCGCCAACGCCTCCGCACCACGCTGGAACCGCTAAACTCCGGAATTTTACACGCGAACAAGTGCCGCGCAACCCTTGTCGCTGTCCAGGTGTGCTCAGTGAAAGAAGTAAGCTAAAAGTGGAAGCTTACAGTGGTTCACGCCGCCGAACCCTCAGTGTCGATGTAGTCCCCCAGTCGAAACAGTGGAAGCCGTAGTCGAAGTCGTGGACGAAGATCCCGCGGCTGCCCCCGTCGCGTCGCCGGAGTGCGTCCGCGTCGCAGGTGATGGAGCCGTGCGTCGGAGTGCCCTCTTCCGTCCTTCACGCGCCAGATCCTGTCGATCTAGCCCGCTTCCCGACACGCCGCCGCCGTGCGCTGCTGCCCGGTGTTGAAACGGGCTCCCGTCACGCCGCCGCCGTGCGTATTTGTGGCTGGTGAATCGGATCGATCCAGCAGGTTGGAAGGTGAGTGTAGGGTTAGGGTGGGTTGCCGAGGAAAGACCGGCTTTGATACCACCTGTAACACCCGTGGCCGCCGACGATCGCCGGCGACGTGAGAGACAGAGTGGATCGCCGAAGTCGGCCTGGCGAGCAGAGGAACCTCTGACCTCGCGGTTACTGTATTCAGTAGATAAAAATGTCGTCCCCTCCTTACAGAGAGCGCTGGTCCTTTAAGGCTTACAGACTGGCCACTCAGGCCCACAAGTCATAGTCTATAATATTACATGCATAATGGAAATACGTCTTCTCTTCTCAACGACGCCTTTCTCCTGTGTCCTCATCCTCACGCCGTGACATCTTGCCATTGCAACGTACCTATAGGAGAGAGGTATGACCTTCTTACTGGTGTTTTCTTGGTTGATGTGACTTACCCCTCGAGATTTCGGTGCCTTGAAGCAAAAATGATTTGTGCTACGAAACTTTGTTGATGTCAAACATTGCTAGATCGAGTCCTTGTAGTATAGATGTTATTTGTGGTTATGCCACACTATTGCCCTTGTGCTGTGACTCTTGTTTTGTTATGTGCTCATGCCACATTTTTGCTCGTATGTTGAGCATGTTGTTCCTTTATTATTCAACCCCTGTTTAGCCCAAGGTCATGCTATACAAAGTGATACTTACCGGCCTCCATTACTTGTCAGTCCAACATAAGTGGAAAATCATTCAAGTGGGAATTCTCTTCCCCAGTGAGCTCTAGAAAATTTGCATTATATCTTTTATAGTGGCCGATGTCATTAGTTTAGATGAAGATTTTTGGGTTTGGGGTCACCTTAGGTAGAGATGCATATTTTTTAGCGAATAAAATACATCTACGACTATTATAATATTGAGAATGATGATTAGAGTAATCAACCAAATTAAATGTGAGAAGCTTTTTTAGAAATTATCATTTCTCTTTGTCTGTTTCATCAAAGAGACTCTTAATTAATTTGGTAATAGTAATGTAATATACGTACTCCCTCCGTTTCAAATTTTAAACCGTTTTCCCTTTTCTTGGTACATAGTTTTCTTATGCACCTAATTATATGTAAACAATACATCTACATACATAGTAAAAACAATAAATGTAGAAGAAACAAAATGACTCACAATTTAGAATGAAGAGAGTAGTATACTCCCTCTATAAACCGTTTTTATAGGAAATAAAATACATCTAACGACTATTATAATATTTCATGTCCATAAAACTATAATAAAATCACATTGAGACTGGCTCTGTGAAACTTGCAAAACAGAAATCTAACAACTTCGACGATCAGGCGTTTAACTCCAAATGTCTTGATAACTTGTTCACTGTCTTCGAAATCCAGACCCAGTTTGATTTATTAATTGCACACAACTTTCCAAACTGAACGGAGTGCGCATGCAGCTGATGACATTTCTATTCGGACACACAAAAAGATCGTACTATTGTATGCATTGTCTTTTGGAATTATATAAAAAGCATGGATCTTCATATGCAAGCGAGTCTCTGTTGAGTCATTTATTGTTTGTTTACACTGTCCTAGTCATCCGTCACAGTCCGTTGGGCCCACGGCCGGATGGTGACCTGTCGATTCACCATTGGCGCTTGCTCGGGACAGCGTCGTTTAATGGCGGCTTCCTGCTGTAGGTGACCACCACGCAGGTAGATAAATCGTTCAAGACATTAACACAGTTTTTAAAGTTTAATTTTTATTCATTATTTTTTATAAAAATATTTATTATAAAATGGTATATATATATTTTTCTAAAAATATTTTGAAGATAAATTTATTCATATGATTTTTATATTTTTAAACTCAATAACTTAAAAGTTATTCATGATTTATATTCTCAATGTTTAACCCAAGTCTTGTCCAAAACAACATATGTAATACATATAGTCTCTCCGTCAAAAAAAAAAAACATTGCTAGAGAGTTCAAAAATTATCCCATAAAGTTAAAATTTGTCCCATAAAAATTGTCAGTGAAAGCGCTGCCACACCATTCTACACTAAAAGAAAATTAAAACATGCTTTTACCCAAAAAAATAATTAGAACATGCTACATAAAGAGAGAATCATAAAGATTCATCCCGCCATGCAGGGGATGAAGTTACAACATGCCACCGGCACCGGGGTCTACAGATCTAGATAAATTTTAGAAATCTTAGTAGAAAATTATTATGCTCTCTCCGTCTCATAAATAAATATATTTCTTTGACGATATGATTCATATTTGATCATTCGTCTTATTTAATTTTTTTTACAAATAATAAAATAAATGAATCATTATAAAAATATATCTAACAATAAAATAAGTCATAAAAAATAAAAATATTTATTAATTTTTTTTAAATAAAACGAACAGTCAAATATGAATTCTAGAAGTAAAAGAAATACGTTTATTGTGGGACGGAGGGAGTATAGTGCTAATAGTTCCTAATGACTCCAAGACCGTGGAGGCCATCTGGTCTAGATTAATTAGCCCTTGCCGCCATGACGAACAAGCGGCTGAGCTTAAGCTGCAGCTAATGAAACATCCACCTCACCTTTGCGCGAGCTAAGACATTGAAAATGAACTAGATGGCCCTAGGAGTGCTAACCAATCACTACCGTGCACACATGTCACACTTATGATTTGACGAACTATGGTGCCGTGAGCCCGTGACGGCGTGAACCAAATTAAACAACCCCTAAAGTAAAAAGGTGAACCAGATCTAGGGTAAATAACTAATTAATTAGTACTCAATCCATCCCAAATTGTAAGTTGTTTGACTTTTTTAATATCAAATTTGACCACTCGTCTTATTTAAAGTTTTGTACAAAATATCACTTTTTTTTTGTTGTGTATTGGTTTATTAATAAAAGTTCTTCAAAAACGACTTAAATTTGACTATATTTGTACAAATTTTTTGAATAAGACAAGTGGTCAAACTTGGGGTAAAAAAGTCAAACGACTTATAATTTAGGATGGAGGTAGTACATACGTAGTAGTACTCACTATTTTCTCTTGCTGATAACAATATGCGAATGGACATGTTTTCTGTGGTGGAACCTATATCGATCAACTTTAAATCAAGTTATTCTTTTCTTTATCAAAGTTAAGATTTATTTTTCACGTAGCTCATGCGACATATATCTATTGGCAACGAGTCACTCATTTTGACTTCATCAAGCTAAAATTTTGACTATATATTCTCGCAAAAAAAATTGATTATATATGATGTGTTTCGTGAAACAAGAAGAGAATGAACAACCGCCCTACCCCTGCGGATCGGAAAGAAATTTAAATGAAGGGGATGCTGTCTGTAATAAAAGTCATAAAGAGGGAACGGGAGATGCGTTGGTGATGAGCACATCTGCCTTCCTTTTTGTCACGTTTGCTTAGGTTCACAAATTTTTAGGGCCTGTTTAGATGGAAGATGAAAAATTTTTAGCTGTCGATATGTCGGAAGGATGTCGGGAGGGATTTTTAGAAATTAATAAAAAACAAATTACATAGCTCGTCAGGAAACTGCAAGACAAATCTATTAAGCATAATTAATCTGTCATTAGCATATGTGGGTTACTGTAGCACTTAAGGCTAATCATGGACTAACTTGGCTTAAAAGATTCGTCTCGTAATTTTCAACCAAACTGTGTAATTAGTTTATTTTTTTTATCTACATTTAATGTTCTATGCATGTGTCCAAAGATTCGATAAGATGGATAAAAATTTTTTGGATGGGAAACTAAACAGTTCTTTATTTTGAGCCTGGGATCACAGAAAGATGCCTGAAGTGGAAGGTGCCCATGAAGAAAAGCCACCGGTGTCGCTGTTGGTTCCCTGGGCCCACCCATCGCCATCGGCGCTAGCTCTTATCTGCAAGTGGGCATATGTAAAAGTTACATGCTCCCTCCGTTCCAAAATAAGTGTCGTTTTCGCTTCCCGAGAAATAACTTTAACTAAATATATAGTAAAAAATATTAATATTTATAGTACATAATTGGTATCATTGGAAAGATCTTTAAGTCTAATTTTTTAACAAATTTATTTAGAGATATAAATATTGCATGTATTTTCTACAAATCGAGTCAAATTTGTGGCACAAAAACCTAAAACGACACTCTTTTCGGAACGGAGGGAGTATGTAGTATGCTGTCATAATATATATAGTATAGACAAAGCTAAAGTGGTAGTATATGCACACTACAACCTTTTGTTTCGAAAAGAATACACCTCTAGCCTAACTCTTGATTAAGCTATCTCACAAGTTTCAACAAAGGACGAAAGTAGAGAGTTGCAATGTACAAAAATGCTCAGGATGTACTGAATATAAGTATCATTAGGTTAACTGTGTTATATTCTTTTCATAGTAAACATTGACGTTATTTTATATAAATCTAGTTGAACTTGAAATTAGTTAACTTCTCGAAATACGAGATATGCAAATTTAAGATGAGGGAGTAGGAGTATATGCTAATCCTCCATCACAATGTCTTGTGTTTTTTTAATTTTTAGATTGGTGTGTATGTATATAAAAAAGGTTAAAGTCCTATATCCAATATATAAGTTGCGTTATTGGTTATATATAAATGCAACTATTTAATTGTAGACCAGGCTATATATAGTACTCCAAGCTTTTAGTAATAAAGAACAAATATATAGTGCTTAGATAATATCATACTACTATATATAGCTGTTTATAAAGAGCATCTATCTTGGTGTACTACCTAGTAGTTAGGTAATATGTGCCCCTAAGTTAAGCACACACTTGTTTTTGATCATTACATTGTACAACATAGACATCCACAACACGCACATACATTTACGAAGTACTGTACCGGTAGATATTGAGAATCACGGAGTATTGTACCGGTAGATCTTGAGAATTACGAAGTCACCATAGACTTTGTCTATCACTAAAAAATATAGTGCCGTTAAATTCTAAAATAAATCAAAGAAAAATACAAGCATCCATACTAAGTCAATGACTTAAATTTAGATGAGCAAGTTCCATCATAAGAAATCTAACCAACTCATATACAATGAGTTCAGAGCTAAGGTCTTGTTTAGTTCCTTTCGTTTTTATTGGACAAACATTATCCAATTATGGAGTAACTAAGCTTAATAGATTCATCTCGTGATTTACAGGTAAACGGTGTAATTAGTTTTTATTTTCGTCTATATTTAATACTCCATGCATGTGGCGTAAGATTCGATGTGACGGAAAATCTTGAAAAGATTTTGGTTTTTGGGTGAACTGAGCATCTCCAACAACTTGGTAAAATTCACTTGTCAAATCTTGAGTTATAGCAAGTTGCTAATTTGTTTAGCAAAAGAAAAAAGGATGTGTCTCCAATAAGTTGCTATTCTCACTTGGTAAAAATAGAGGATGACAGATGGGACCAGCTCACTTGGTAAAAGAATATGTGTCTCCAACAAGTTGCGCTCGGGATAGCAACTTGGGATAGCAACTTGACAAATCTCGGCAGTAGAAACCTCTGGATAGCAACTTGGGAAATTTAGCAAGTTGAGATAAGGAGTTGTTGGAGGAGGATTTTTATGCTTTTAGCAGATTTGGTAGGATAGCATGTTGAAATACCAAGTTGTTGGAGATGCTCTAAACAAGGCCTAACTCAATTTTTAACGCGGCGGGGGGTGGGTGGCATGTGTACAATGAGTACATACTTCTCCCTCCATATTGGATTTAGAAGATGTTTTAGACACAATTTGAGTCAAACATTAGGAATATAAATTTTAATAATTTTTAAATTATTGAGTTTACAAATATAAAAATTATATAAATAAATTTGTCTTGAAAAATCTTATACATATATTATTTTTTAAAATATTTTTATAAAAATAAAAGATCAAAGTTTGTTTTAAAGACTTTGTCGCTGTCCGAAATAAACTTATGTACTAGCAAATCAAACGGAGCGACCTGGGCCCAGCGGGTTGAAGCCCGAAAATAGTAACGTGGGCACTGCCGATATATGGCTGCTGCCCTGCACTGCCAGTCAGCCGGTACGGTGCCAGGGCAGGATGGACAGGACCACTTGTTCGTTCCTCCTCAATAATTCGAGGCCTGGGGCCCACTGCTACTCTTGACCCACCAACAGTGAGCTGTGCTTATCGCTGACAGTAAAAAAATATGACAAGAAGAATAATGGAATTGGTAAGTTCCTGACATTAAGGCCTTGTTTACTTCAAAAAAAAAACAAAAAATTTTCAAGATTCTCCGTCAAATCGAATGTTTAGACGTATGCATGGAGTATTAAATATAGACGAAAATAAAAATTAATTGCACAGTTTGGTCGAAATTGACAAGACGAATCTTTTGAGTCTAGTTAGTCCATGGTTGGACAATAATTACTACAAACAAACGAAAATGTTACAGTGTCACGAAATATTTCCCATCAGGAAGTAAACACGGCCTAATTGGAGCCCGGCAAAGGAAAGGAGAAGGCTTTCCTGTTATTCCGAGACCGGTCGGACTCCAACTTGGTCTTGGAAGGCGTTCAAGAGTCAACCAGGCGAAGAATGTGTGTCGAATCAGCCTTTGCTTGTCAGTTTTGTTGATTTACAACAAGTGCGTGTTAGCCCTTATTGGCGTATAACCCAAATTAGCGTTAACAGCTCTACTTGACCTTACATTGGGTCTTTCTTCTTTACTGGGAGTAGAAGTACCCATATCATCCTATCACTGTAAGTCTGCATCCCAAACACCGAAGGATAGATGATCGATCCTCGGTCCTCAACAGTAAGCTCTGCTCATGGCAACCTCTCTTTCTTCTTCTCTCTTTTGTGTAGCTAAAGATCTTGGCATATACACTCTCGTATGAAATTGAGGATGTGTTCTTTTTATTCATCGTTATAATTCGAGAGTTTCTTTCACCACCACATTTGTCTACTGAGTATGACCTTGTTTTCTGAGTCCGCTTGGATGGAAGAATGGAGGGTGCCTATTATCTCTTCTTTTCTTTTCTTTTTTTAGTAGTGGTGGTTGTTTGTATTGCTACCTTTTTTTTGTTGTAGATGCTACCCGTTATTATTATTGTAGCTTCATCCCTCATTAACATTGGTGAATCTAGCTACTAGCTAAGAACATCTTTAGAGTCATCAGTATCAACTTCATTGCCTTATATAAAGTTTTTACACTTTAGTAGTATTTTGTTAAAAGTTTTATTAATTTCATCGGTGTGGACGAACCCCCCCGATAAACTCCCTAAATTTAATTTGCCCCTACTCACTAACCAAGCGTATTCATAGCATGCTGCCTTTTCTCGAGATAACTAGTTGCATATCATCCTTTTCATCATCGTTTCGAAGGCTGTCCCGGCTTCTTCCGCCATCAATTTTTAGACGCAACATACAGTTGGTTGAGAGATGAACTATCTGGAGCTATACCAAATCCAATATACATGTATGTATGTGTGTATGTGCATATCCCATTATATACGATGCATGACTACTATGTACTTAGCAAATGAGGAAGACTAACTACGTGTGCGCATGCGTATTTGGCCCCACCCGACCGAGAAATTTCAACCATATTATTCTTCGTACGAGCAGAAGAGAAAGATTATTATTGCGGTACACAATGCAAGTCATATGTACGCGTACTAACTCGTATTTATGTCATCGTACACCTACAACTTTATTATATTTAATTAGCAAATATATAGAAACCACTAGTTGCAAGTATATATATTTCTAGAAATATCTAGCTAAAGATTTAGTACTCCCTCTGGTTCAGTATATAAAGTACTAGGTACAACAATTGTATTATAGGTAAGCAAATGCCAATATGGAAGAGAGAAGAGAAAAGAGAGAGGAGAGCTCACCTATTATACAAACAAGAGCTTGGCATAGATGCCTATATATTTGTGGATCCAAACATCAAGTGAATGCGAGGAGGTGTAGAAAAGAGAAGTTGTGAAATGGGTAAAAATAATAATCACTACTACAAAAAACTTAATCATTGCAGCCACCTTCACCACTGCATGTTCGTATTACGAACTGTTGGTGAAAATGGAGGTGGTTATCACCGCTGGCTAGTAACACGAACCGATAGTGATAACCTATCACTGTCGGTGGTGATGGTCGGCCCATTTTCACCGCGGGTTGGTGTTACGAACCTATGATGAAATTCACTTTCATTAACATTTTCAAATATAGTCGGATGAAAACAAAGTTTATAAGTTGATGATTTTTCATTTGAAATCATTTATAGTCCCAGAATTCCGTTTGAAGTTCTCATACTTTAAAATTTGTTTCTTTTGAATTGTACAAACAACATTAGATGAAAAAATTACAAGAAACAAAGTTATAAAGTTGATAATTTTTTATTTGAAATTATTTATGGTCTTAAATTTGTATTTGAAGTTTTAGCATTTTAAAATTCAACTTTTTCAAATGACCTAAGATGAAAAATCAACCAAAACCAAAGTTGTAAATCTCAATGAGAGCAACAACTTTATTGTTGATGACTTTTTATTTGAAATCATGAAGTCCTAAAAGTATAGCTCAAGTTCTCATAATTCAAAATTTAAATTTGACAAATGACCTCGGATGACGAAACCATCGAAACCAAAGTTGTAGATCTCGATGAGAAGTACAACTTTGTTGTTAATGGTTTTTTCATTTGAAGTCATTTGGGGTTCTAAGTGAGAAATGAATACTAGGGGAATCAATATGCTCATGTGACTTAGAGTGGTTACAGTAGCTGCTCTTGAGGACACACATCTCTGGTTCGAATCCGCCCGCCATGTAGCATGGGTATTTTTGTGCAAAAAATCGTGTGACTTGTGCAGTCGGATGATGTGGCTTTCTTAGGATTAAAACAAAATTTTAGCTACTTATTAACTAGAGATTTTAGAATTTCTAGAAATCCTAACCAGCGGTGATGACCATTCTTAACACCAATGACTTTCATCACCGAGGGGCCAAAACTAGTGTTGATGAGGTGTTTGGAACCGGCTTAATCAAGATTATGTAGTAGTGAATAATAAAATTTTACATAGAAATAAGAGACAAATATTATATAGCTTGGCTCTTATAATAACTTGACTTATAACCAAACACTTGCTCTACTATTTTGGGCTTGCTCTAAGACATCTTCTTACATGATATGATCTCTTAAGTTACACTTTGACCTTTTTAAAGCAAGTCGGCAAGGGAGCCCCTTACCTATTTTTTTCTATTTATCAAAAGGGGAGAAAGTTATATTGTTTATGATTATGAACTTATCATTACTATGTAGTATATTTGATTACAAATTTGAAAGTTGATATAGCTAAACTATTGGTCAAAGATTATAAAACTTAATGTTGTTGCTAGCCAAATTTAAAGTGATTTTGGATGGGTCCACACACAAACCAAACTCGTTGTTTTGTGGATGGACAGGAAAGTCCCACCAATAATAAGGTGATGGGTGACTCCATTCTTAATCTCCACATATCCCCTGTTTAATCATCCTTCCCGAAAGAAATAGGCCACATCATGCTGAATTTAATAATTGGTAGGAGTACAAAATTGTGCGGGCACTGATGAGGCACACATGACACATGCAATATACAATAGATTTCTCGGACACGCCTTGATTTCATGATGTTCCAGTTCCAACAACCATGCATGCACTCCCTGGGCCTGCTGCATGATGCTGCCGTGTTTCCTTGTTCATCGTCTTTGTTTCCTTGTTCATCGTCACATCACAGACTAAACGAGAGATCGTCGTCGTCGATTATAAATATATATGATGATGATGGAGAGGGAAAATATATTCCTTTTTCGTAATATACAATATATTGACATCGATCGGTCTAGGAGAGGGGAAATATATCTCTCTGTCGGCTCATATAATGCATACATATGAGTGCAAGTGGACTCGATGGATGGGTCTGCCTTGTTTGGATGTTGTCGAATTCACTTCAATCCATGTGTGTTGGTGTGGAATTTAGTTTAAATTCCACTCTAATCCACTTCAACACATGTGGATTGATGTGAATCCGACTACATCCAAACAAAACCTAAACTGTAATGAATGTAGTATGCCTCAGCTAGAACTTCTTGTGGTTTATAGCGTTTACAAGCTCTAACCTATCCAAACACATATCAGTTTTTACACTAGCTTTTGTAATCTACTCGTAAATATAATGCCCAGATTTGGCATAGCAGGCTCTTTATATATAGCTCTTAACTTTACCAACTTATGATGAGAGGGAAAAAAATTCTAACCCAAACTTAAAATGAGTCTTGTATTTTGAGGAATAATGGTCCGCTCGGCACCCATGACAAAGTTATACTACATGTTGGTAAAATAATAACTAGCTAGATTGTCACAAGTCACGAGCCAGGTTTTTATATAATAATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATCCGGCCGGCCAGCGAGGATGATCTTTGTGATAGGCAAAGGGCAACCGGGCGGGGATGGATGGATCTTGGCGATAAATGTCTGCGGGCCAATCAAGGCGTACGAGAGTGGTGGCGTGCAGGTGCAGGTGCAGCGCCGGCCGGCCGCGCCCACTTGTAGTACGTGCACGTGCTACGTACTGTACGCCGCGCGCCCGCCTCACCACTTCACATACTCCGCTGTTCATAACTTCTTGCCAACTTGCACGTACGGCTGCACGTACACGTACAAAAGGATCAACTGCATGCCACACGTACTAACCAAAGGATGATCTTTTTATTTATTTTTTTGTGTCTGCTTTTGTGCTTTCACGTGACAGTCTTCAAGACCGCAACGAAAGACACAGCCGCTAGCTAGCCCATGGTAAACAAAAGGGTACCTACACCGATAAATTCTATTTTTTAAAAAACAAAGAATCAAAGATTTAATTAGTCATTAGTTTGGTACTCCTTCCGCTATGTAATATAGTGCTTTGAAATTTGTTATGAAACATAAGAGCTTCGAGAGAGAGAAAAAAAATATAATGTTTACTTAATAGAAATACTTGATTATATATATGGGAACTAAGTCCGCCGTTAATCATTGGTAAGAGAATAGGATTTCCATAAACTAGGAAGTTATTACAAGTGATTGATGCATGCCTCTATTAATTACACTCTTCTTTGTTTTGTTTGAGAATAACTAGAATGTGCTATGTTTCAGGATGGAGGGAGTGTAGACTCTATAATAATAATAAGTGGCGGAGATAGGAGCAACTTTTACGAATGGCTAATGAAAAAAGTCTACTTAGCCCCAAGTATTGACCCATGTCTAATCGAGGCCCCCAACTACAAAACCACATAACCTACACTCTAAAGTTTGTAAAACAGGGCAACGAACCCCAGGGCGGTTTCGCTGACGTGGCGTCTAGCGTGTCACTGGTATTGACACGGTGGCGCGGTGGATGTGGCTAGCTAGTGGGTCCCACACGTCAGGCATGCTATCCCCTCAAGAACTCCGTCGCCTCCCCTTCGCGTCTTCGCGTGTCACCTCGCCCGCCTCTTCGCGTCTCGCGACGCCGACGCCACGCTCGTCCTCACACGCGCAACACCGCAAGGCTCGATGGCCCCGTCGCCCGCGCGGACGCCGTGAGCCTCTGCATCCGTATCACCACGCTCACGCTGCCCGCAGAGGCCGCCGGGCTCACGGGTGGTCACGAGAGCGCCGACGACCTTCCCAACCGCATTCATCGGCACGGGGTGCTTCGCGCTCTCGGTGCAGCCTGCTCTGCTGCTCCACAGCCCACCACAGAATGACATGGCGTCCGACACAGAGTCATTCCTCGTTCCCGGCCTCCCCGACTTGGTGCGACTCACCACATCGAGGGTCGCCGAGGCGATGCTCTCGAGGCGGACTCACGTGAGTTCTTGAACCGCATGTTCGACGCCGAGCACGCCACGACTGGGTGGGTCGTCAACTTGTTCGCTGACCTCAAGCAAAGGTACATCGAACACAACGAGAAGGACACGGGGAAGCAGGTGTTCACCGTCGGGCCGGTCTGCCTTGTCAACGACGACGACAACGACACCCTGGAGCGTGGCCACACCGGGGAGGCTGATACGACCGCCGAGGCCATGCGCGTGCTGAGGTGGCTGGACACGAAACCAGCTTGGTCTGTGGTCTACGTCTGCTTCAGTAGCCTCACCCGGTTCCCACGGGATCAGGTTGTGGAGCTCGGTATGGGTCTCGCCGACTCCGGCGCGAACTTCGTGTGGGTGGTCAGGGACAAGAACGCGACGCCGCCGCTCCCGGACATCGACGACACGGTGCCGGGCCACGGGCTTGTGGTCGGGGGGTGGGCCTCGCAGGTGGCTGTGCTGCGGTGGGCGCGTTTGTGACACACTGCGGGTGGGGCACAGTGACCGAGGCAGCGGCGGCGGGCGTCCCGATGGTGGCGTGGCCGGTGTTTGCAGAGTAGTTCTACAACGAGGTGCTGGTGGTGGGGCTCACGGGCACGAGCATCTCCATGGGCGCAGAGAGGGGGGTATGTGTGGGGAGGCGAGGCGCTGGGCGGGGTGGTGGTGGGCAGGGAGGCAGTGGGCTGACATGTGAGACCCACCGGCTAGCCACATCCACCACGCCACCGTGTCAATACCAATGACAGGCTGGACGCTACGTCGGTAAAACCGCCCTCCCAAACCGTCCAGGGGTTCGTTGCCCTGTTTTACAAACTTTAGAGTGTAGGTTATCGAGTTTTGTAATTGGGGAGTCTAGGTTAGACACAGGTCAATACTTCAGGGGTAAGTAGACTTTTTCTAATGGCTAATTTTTCTTATGCTGATCAATGTGTACACAAATGTCGGAGCAATTAGCCACTTACAAATTTCATGCTAACTAATACCTAAAAATGTAACAAAAAAAGCAAGTAGACTATTAAGTATCCAAAAAATTTCATTTAAATAGTTACTGTATAGGCTATATAGCAAATAGTTCTAAAGATTAGTATGAAAATTACAAAAAAATAAAAAGTAATTAAACTTTGATTTTTAGTATTTCTTAGACTGCGGCAAGAGGAGGCCAAGGCCCAAGTTCGTCTGGCCTTAGCTCCGTGTTTTCGTGGGGATTTTTTGTCCCCACTGTTCAAAACTATTATGTGTTAGAAACAGTATAGATAAGTTATATCCCTACAAACTCATTATAATTAATAAAACATTCATCATATTTCTCCATCTATACTATATATGTCATATTAGCTTATTTATTGTCGATAAAATTCTAATGCAACCCCACTGAAACTGGTCTAGTAGCCATAGAATATCCACTCCCTTGCCATGAAAATTAACCCCCTGCTGCTATTACCAAAAGTAATATAACATAACTCACTAAATCAATATCTAAATTACCTCCCTAGTTTTTACTGTTTTCTTAGTAAATAGTAAACGCTAGTGCTAACTTTTTTTTTTGAAACAAGTGCTAACTATTCATATTGTCTTGGTTTAGCACTCTTGGCCCATCCTAGTTGGCACTTCAAGCTTTACCGTAGATCTGAGCATCAGCATACACTTAAACGAGTAAACGACCGTTTATGCCATAGTAACGAAGAGTGCATTGTATGACTTAGGACTTGTTTAGTTGATGACCAAAATTTACCATGCCAAACAAAAATACAAACGAAAAAATTTGAATGAAAGACAAATCACACAAGCTTTTAACAATGGTAGCAACTAATAAAATCAAAATCTCAAGATGTTACGACTTACGAGGCACAAGACTCCTTCCATTTCAATTTTTTTATCATGTATTTAACATACGCCATATCTAAATATATAATAAAAATAATATATCTAGAAAAGTCAAAACTCTTATAATTTAGGCCTTGTTTAGTTCAACCTGAAAACTAAAAAGTTTTTAAGATTCTCCGTCACATCGAATTTTGCGGCACATGCATGAAGTATTAAATATAGACGAAAACAAAAACTAATTACACAGTTTGTCTATAAATCGCGAGATGAATCTTTTGATCCTAATTAGTCTATGATTGGACAATATTTATCACAAACAGACGAAAGTGTTACAGTAGCGAAATCAAAAAAAAAAATCACATCTACACAAGGCCTTAGAGGAGTACATGTTAATGGCACACAAAAAAAAACCAATATTTTTGGTCAAAATTTGGGCACCGGTTATACGAAAAATGGATGACCCACTTGCCGTGCATCGACCGGCGACGAACTATACGCGATTGCATTGACCTTCCGTTTTCAGTAGTCTGTAGCATTAAAATCCGTGCGCATCAGGGCCGAAAGCAGCTAGACTACACATCGCCACCAGTAGTAGTCCTGTTCCGTCTAACTGACTGTTCATACCCTTCTTTATTTCCATGTGTACACGTCCACGTACACTTACGCAGGTACACACGTTTCGTGCCTAGAATTGATTATTTAATAAAACATGCTTTTCTTACTGATTGTGACATTTTGGGTTGCAACACTAGCACGGAATAGTCTTGAAGCATAATTTAATCTATAAAGGATTTTTTTAAAAAATATTTTCTTATGAATCAGAATCTAATCTACCAGCCAAACGATTTTATGAAGTGATTTTGATTTATCATGAGAAATATTTCCTGAGAGAATCTAGATCCAAAGCATTTCTACAAAAATCTGTCAAACATACCAGAAGTGAATGATTGTTTTTCTATTCGGTTCATCGAATCTTGTATTGAATATTCGAGATCTATGAGAAGGCCTTTAATTTCATTATTAAGTACTCACTACACCAGAAAAAGGCCTCTATTTTTTATCATCGATAGTATATCACCGCCGGTCCGTGGCTTCAATCGGTTGTGATAGACTATTATTGTCAGTTCGTGGCTTGAAGCACTAGTGATAGATTATTTTGCACAAAACAAATTCATAAATTGTGCATATGGTGTAGAATCAAAATATTTTTTTAGTAAATTTCAAAACTCTAAGAAACATAATTTTGGAGTATAAATGATTTCAAAGGAAAAGTTATCAACTAAAAGTCGTAGATCTGTAAAAGCTCTAAAACTTTGATACATAGCTTTGTTTTCATCCGATATTAGTATGTTACATGCTAGCTTCAAAAAACTAACTAGACATAGATCTCAAGTTAAATGTCAATACACATTCTCTAATCCATAACAACTCAAAACTAAACCTTTTGTTGGAAAAATGCACCAAGGAGCCGTATTTCATTCCACAAACATGTCAAGTGAGCATTTGTGTCTACAAAATACAATTTCTAACTAACAAGCCCACAATGGGTCCATTGCATAATTTGCATTGGTGGTCACCAGAAAATGAACATTTTTACTTGGCGGTTAAAAGAAGAGGTGAACTAGCACTGCGTCCATTTCTCAATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATGAAGTATGTATTAATTAGTTTGTCTCAAATATAATATATTTAATGCCTTTTTTTGTGGTCACAGTATTAATTGCAGAGGTCGCCCACAAGATATTGTCTGCTCTAATTTCTAGCTAGCTAGCTAGGCCGGCAAGGTCGATCAGAAAGCAATCCGCGACCAATCGACTATATATACTATATTCTTTTACATAATATAAGGGGTCACGTGCACTGTTCACGTGTCTTATTCTCGTCCGCTTCATTGCCGAATATGCTACAACACTAAATAAAATGCTCCCTCCATCTTAGTACATAATATAAGGCGTAACTATACTTTTTGCTTATGTCTTGTAATATAAGACATGCTCTTTCACAAACATAAAAACAGTACTGCTAGTGTTATGAGAGAAATTAAATACGTTTTTGATCTATATACCGGGACAGAGAGAGTATATATATATGAGGATTCCTCTGCATGCGGGTAAAAGAAGAACGAAGATACTGATCGAGATATCGAGCTGCTGGTATACTAGAACAAGCCTGCTACTAGTAGCTTGCATACACATGAGATGCCTTCATGGGGTGATGGGGATTAATCCGGCGGCCGGCGGCGAGCGAGAAACCTAGGCACAGCATAGGTAACTCCGACAGAGTCCGATGTCGGACCGGCCGGCGGTGCATGGTGAATATAAAGCTAGCTTATCCTATCTGTCCGATCCGGCCGGGCAAGTGGCCTCGCAACTTGTTGGTACTACACAGAGGATCAGAGAGGATATATATGCATGGCAGCTCCGGCGGCATCTCTCACTTAATTGCTGGTACCAGCTGCCTGCCTTCGAGATCCAGACCCAGCAGAGCTCTATTTATTTGCTACTTCCTCTCTCATATAGGCCTTTTGTTTAGTTCCAAAAAATTTTGAAAAATCGACACTGTAGCACTTCGTTTGTATTTGACAAATATTGTTCAATCATGGACTAACTAGGTTTAAAAGATTCGTCTCGTCAATTCCGACCAAACTGTGCAATTAGTTTTTATTTTTGTCTATATTTAATACTTCATGCATGCGTCTAAAGATTCGATGTGACGGGGAATCTAAAAAATTTTATAGAATTTTTTGGGAAGTAAACAAGGCCAAGTAGAAAGTAAATATAAGTCGGATTAGCCAGACCTAATGGAACTTTGAAATGAGTTTTATTCTAATTTAATGATGTAATACAGTACGTATTAGTAATCACTGGTAGAGAAAAGAATATCATATCAACATATTATTCCTGGTTAGAATTGAGTCCGGGACTAGTGGGACTGGAGGATGTGTTATCTACTTAGAGAGTGTTTGGTTGGGGGTATTAAAGTTTAACAGCTGTAGTATTTTTGTTTTATTCAGCAATTAGTGTCCAGTTATGAACTAACTAGGCTTAAAAGATTCATCTCGCAAGTTGCTCTTTTGCTATATTTTTAGTTTCGTAAATAGTTTATATTTAGTACTTCATACATGTGTCCAAACATTCAATGTTAACATGAGTTAAACTTTACCCGAAAAAACTAAACGGGGCCTTAGAATAGTTATCTACTCTAAAATGGATGTGTTATTATAAATTTGACAACAATGATGAACTAGAGATTTTATATAATTTGTCGTTAAATCAGAATAGTTTCATTTAGAACAGATCTAATACGCAGTGGAGCCGTATTTTTCTAACCAAATCGGATGAACTAAATAAAATACTAGGTAGATAATGCATGTTGTAGGTGCAGCTAATGGAGCATATGTATGCATGTATGTATATCTCCGTTCTTAATTAAATTATAAGTTGTTATTTTTTCTAGGTTCATAATTTTTTACTAACTCTTAAAATAATACTCCTATATATATTTATCTGCTTAACAAAAAACTATATATGTACACTAGAATAACAGATACTAAGGGATTATGTATATTTATTTATAAAATACTTAGATGAGGAAGTACACGCTGTATATATATACGTGCGTACAGATCGAAAACGAGAACAAGGAAGGCATCATGTGATCCGGTGCGGCGAGCGAGGTGCGTGCTCCTCGTGGTGCATTACAAATTTAGGACGACTGCACGGCAACGCATCGCATCGCATCGCATGGGGGAACAATTAAAAAATTGCATAAACACGCAGCTGATATTGAGCCACGTCGTCGATCAGGTCGTCGTATATCGGATCGTCGTCAACACGTCGCATGATCCATATATATTCTTTTTTGACAACCTTGCATTGCTTGCAATATACTCACATCACTAAACACATGTACTTTGATTACCATTAACACTTTTGTTCCGATAAGCATATATAGTAAAAGTATAAACTCGGACTGCGGGGAATAATAGCATACATCCCTTGAAAATTTTACACTAAGGCCTTGTTTAGTGCTAAATTTTTTTTTACAAAATAGATACCGTAGTATTTTCGTTTGTATTTGATAAATATTGTTCAATCATAGACTAACTAGAACCTTGTTTAGTTCACTCCAAAAAACCAAAAACTTTTCAAGATTCCCCGTCACATCGAATCTTGCGGCACGTGCATAGATAATGAAATATAGTCAAAAATAAAAACTAATTAAATAGTGCATAGATAATGAATCTTGCGGCACCGCCGTGCGCGTGCACGGCCAGCATGGAGGGCGGCGTGGCGGGCAGCCGGCTGACCTGGGCACGGGATGCGGATGCGGCGAGCGTGCGCAGCTGATGCGGCTAGAGTTTGTTGCTTCCTTCTGCCCCTATCCGGTACTTGTACCATACATTCTTTTGAATGATGCAATTCTCATATGGTTTCACTAAACCTTAACCAACCAAAAGTCCTGGTGCACCTTTATTTGGAAAAATGTTCTATATCACTTAAATTTTTGCTATTAAAACTCAGACCTAAAATATCTAAGGGTATGTTTGGGAGGGATTCTTCAAACTGTTTCCTTAGAAATAGAAGAAATTATCAAACAAAAACTAACAACAAGGAGAGAAATTGGATAATCCCTACCAACGATATGTGATCTAGATTCTCGTGAAACACTCTGCAGAACAATTCTATGATAAATCATAATCAGCTTATGAGACCATTTGGCTAATAGGCTACTAGATTCTTAGAAAAAAATGTTTTTTTTTATAAAAAAACATATACCTAAAAATCTGCTTAAAAATCAGAACTTTTAGGCAAAGTACAACATAAACATAACCCAACTGAACTCGTTGTACTTAAAAATGATTAACTGGAAAAACACGTTGGAGAGAAGCCCTTGAGATGGCCTAGAGAATTCTGGCGAAACTAGCCGGTGGCTATTTCGATTCGATGAAGGAAAGGTGAAACGTTTCAACTCAGAATCCCCGCGCCGCAATCGCGCTTGTTACCGTTCATCGTGTTTCTCGTGAGAAAGAGGTCTGATTGTAGCATCAAATTAAACTCTCTATCTCACAACAGAGAATCAGTTCCTCACACAAAAGTGTCTCTTCAAGAATCAACCAAGACACAATCGGCCGGGGAGCTAGGCCTTGTTTAGTTCCCAAAAAATTTTGCAAATTTTTTCAGATTTCCCGTCACATCGAATCTTTAGACGCATGCATGGAGTATTAAATATAGATGAAAATAAAAACTAATTGCATAGTTTGGTCGGAATTGATGAGACGGATCTTTTGAGCCTAGTTAATCCATGATTGGACAATATTTTTCACAAACAAACGAAAATGCTACAGTCTCTGTTTTGCAAAATATTTTGGAACTAAACAAGGCCCTAATAAAGTCAACACAAGAGAAACCCCATAGCTCAAGTCAAGTGTCTCCAAGGGTATACCGTAGGCAGGCAGGAGAAGCTAGCTGAGCAAGAAGAAGCATCAACCAAAACAGTAAGAGAAGATATTTTCTCTCAGCACACCATAGCTTAGCATGGCAAAAGTATCATCATCTGCCACCGCAGCACATCCATGGCATTTTTTTTTCTAGACCTCGTGTTTTTATTCCTACTACTGCTAGGAGGTGGCTGGGTGGGGTGGAGAAATGACAAATCCAGCAGGCAGCCACACCCACCAAATCGGACGATGCAGGGTGCCCAAATCAGGAAGGATTTTAAGGTTAACCGGCTGCCACCGGCCACCGCCGGTGACCCCTCCGATCCTTGTCAACACCCCCCACCCCACTCTCTCTTCTCTTCTCTTCCATCCATCCATCCATCTATATATTACCCGCCTCTTTTTCTCCCTCCCTCCGCCCCACCTCCTTCTTCTTCCTCAGCTCCGTCGCCCCCCGCCACCGCCGCGCCCGCCCGCCGCGCCGGAGCTCGATCTCCGAAACGCCCTCGGCCGCCTTACCTGTAAAAGCCCCACCTTTAGCTAGCTACCTCTTCTTCCCCTAGCTTTGCTAGCTAGCTAGAGAAGGAAGGAAACGGTAGCTGACAGACAGTTGATTGGTGGTCCGGGGTACGGTGTGTTCTTCAGTCGTGAAGCGACCGTACAGTGGGCTAGGGCGTTCTCTCCGGGTTGCGTGCAGGATGGTCGTCAAGGATCGGGAGTGAGGAGGAGGGCAGCTCTGCTCGTGGTCGTGGAGGCTAAATGTACCGCAAGAACGACTCGGCACTCTCCTGCTTCTACCTCTTCCTCCTCTGGTTCTTCTTCTTGAAATAGACCACCAGCTCGCCAGGTAGCTACTAGCTCTAGCAGTAGCAGCCCAGTTGCTTGCGATTGGGGATCGGTCGTCGCCGCCCGGGCTCTCCTCTCAGCTCGTTTGTTTGCAAGTCGTTGGAAGCTAGCTAGCTAGACGAGGGAGATGGAGCTGGATCTGAACGTGGCCGAGGTGGCGCCGGAGAAGCCGGCGGCGGCGATGGAGGCGAGCGACTCGGGGTCGTCGGAGTCGTCGGTGCTGAACGCGGAGGCGGCGTCCGCGAGCGCCGGGGCCGCCGCGCCGGCGGAGGAAGGCTCCAGCTCGACGCCGGCCGTGCTCGAGTTCAGCATCCTCAGG"; //sorghum
    std::string seq2 = "GGTTTCAAGATCTTTTGAGGTTGCTTGGGAGTCTCCAAAGTTGTGGACGACCCCAAGAAGTTTGTATCACCCGCTCTTTGAGCTAAGATAAGAAGAGATTGCCTTGACCTTTGTGGTCGGCTTTTGGAGGATTAGGGTTGAAAAAGACCCGGCCCTTTGTGGGCTCCTCAACGGGGAGTAGGACACCTTTGTGGTGTGGCCGAACCTCGGATTAAATCTTGTGTTTTGTGTTCTTGCTTGTTTTAACTTGAGATATTTTTACTTGCAAGATTGACATAAAGATTTTTCCGTAGAGTTTATCTAGAAGTGAAATAGCCTTTACATATTCGTCGTTTCATCCTAACTTAGCTATATCTTACTTCTCACAAGAACTCGCGAAATACTTTTTCGAGTAGATTTTTAGTCTAAAGTTTATAAAACAGTTGGATTGGTTCAGAGTGGTTCAACTTGCACTCGTAGCTGTTGAACCGGCCTGCACCAGTCGAACTGTGTATTTAACAATCTTTCGAAGTTTGTTGTGAAAATTTTCATATTAGCCTATTCACCCCCCCCCCTCTAGGCTACTTTCAATTGGTATCAGAGCTTCGTACCTCGTTTTATGCTTAACCGCGTGAGGAAACGATCATGTCTACTCAACGGGACCATGTGGATCCTCTCCTCGAGGAAATCCCCATCACATCTTCCGGTGAGGAAGTGGACACCAAGGTCCTCGACCTCGCCATGAAGATTGTCGAGAAGATGTTCCTCAAAATGAAAGAGGAAGATGCTAAGAAAAAGGTCGAAGAAGAAGAATCAAGAAGAAAGGCCGAAGAAGATAAAGGTAAAGGAAAATTTGATTACAATGATGATCTAGTGGATCTTTTGGTGTCGAAGGTGTTGAGCAAGGTAAGTCTCAACACCGAGGGATCATCTACCAAAAGCAAAGGTAACGAATTTAGTAAAGTTCAATTCGATTACTCTAGAAATTTTATCCACAACTTCTCTTCCGCCCCACTTGGAAAGTTACCAACTCTTAGTGAGTTGAACTATGACGAGTGGGCCGACAAGATGAAGTCGCATTTAATCGGTGTGCAGCCTAGTCTTTGGGAGATTGTTAATGTAGGTATGTATAAGCCCGACCAAGGAGAAGAGATGACTCCGGAAATGATGCAACATGTTCATCGCAATGCTCAAGCAGTGAGCATTATTAAAGGAAGTCTTTGTCCGGAAGAATACCAGAAAGTTCAAGGAAGAGAAGATGCCCGTGACATTTGGAATATTCTTAAGATGTCACATGAAGGAGATCCCAAAGCTAAGAGACATAGAGTTGAAGCTTTGGAGAGTGAGCTTGCAAGATATGATTGGACAAAGGGTGAGTCGCTTCAATCACTCTTTGACCAGTTGATGGTGTTAGTAAACAAAATAAGAGTGCTAGGGAGTGAAGATTGGAGTGACTCCAAGGTCACAAGATTATGCATGAGAGCATATAAAGAAAAGGATAAGAGTCTTGCAAGGATGATTAGGGATCGTGATGACTATGAGGATATGATGCCTCATCAATTATTTGCAAAGATTCAACAACACGAGTCCGAAGAAGCCCCTATCAAGTCAAGGGACTCTCATGCCTTGATCACTAATGAACAAGACAACCCCAAGAAGAACAAAGACCACAAAGCAAAGAAAGTGGTCGAGACCTCAAGTGATGAAGATAGCTCAAGTGATGAAGACACAGCTATGTTCATCAAAACATTCAAGAAATTTGTAAGGAAAAATGACAAGTTCCAAAGGAAAGGAAAGAAGAGGGCATGCTATGAATGTGGCCAAACCGGTCATTTCATAGCGGATTGTCCTAACAAGAAGGAACAAGAAGCTAAGAAGGAATACAAGAAGGACAAGTTCAAAAAGGGAGGCAAGACCAAGGGATACTTCAAGAAGAAGAAATATGGTCAAGCCCATATTGGTGAAGAATGGAACTCCGATGAAGAGAGTTCTAGCTCCGAGGAAGAGGAAGTGGTGGCAAACGTGGCCATCCAATCTACATCAAGCTCGCAACTCTTCACCAACCTACAAGACGACTCCTACACTCCAACTTGCCTCATGGCAAAAGGAGATAAGGTAACCTTATTTAGTAATGATTTTCCAAATGATGATGAACAAATTGCCATGAAAAATAAAATGATTAAAGAATTTGGCTTGAATGGATACAATGTTATCACCAAATTCATGGAGAAGCTAGATAAAAGAAAAGCAACTCTTGATGCTCAAGAAGACTTGCTTATCCTTGAAAAGGAAAGAAACCTAGAGCTTCAAGAATTGCTTCCCAATAAAGATGAAATGCTAGATGTCTTGACTAAGGAAGTATCTTTAGTCAAGATAACTATAGAGAATAAAGATAAAGAAGTAATTAATATGAAAACCTCTATAGCTAATCTTGCAAATGAAAAGAATGCACTTGAAACAAGCATGTTAAGCTTGAATGCTCAAAATCAAGAACTTCAAGTGCAACTTGAAAATTGCAAGAACATCAATGTCTCATCTTTAGTTATTGAATCTAAGTCTAGCTCCTTAAATGATAATTTTTGCAAACATTGTGCCAAATATCATGCTTCTTGTTGTCTAACTAACCATGCAAGGAAGAATAGCCCACAGGTGAAGGTCAAAGGAATTTTGAAAAGATGCTCTAGCAATGATGGGTTAAAGAAAGTTGAACCCAAGTACAAGTCCCTAAAGCCCAACAATGGAAGAAGGGGGCTTGGGTTCAACTCATCCAAGGAAAACCCTAGCACAGTGCATAAGGGGTGGAGATCCCCCAAGTTCATAGAGGGAACCACCCTATATGATGCCTTGGGGAGGATTCACTCCTCAAATGACAAGTCACCTCAAGTAAAGGTAAACTTGAGTTCCACAAAGAGTAAGATGAAGGAAGTGGGATCCTCAAGTGGACAACAATTTAATGCTCCCATTTCTCACTCTTATCTTTGTGATTATATGTTGACTTGGGATTTAGGGAAATTGGTTGTCAAATATGTGGGTGCCTACACTAAAATAAAAGTAATGAAAAGAAGTGTGTGGGTACCCAAGGCTATAACTAACACTGTAGGACCCAATTCAATTTGGGTACCTAAAAGCATAGCCTAAAGTTGTTTTGCAGGTCTACTTCTCTGGTGGGTCAAGTTGGGTGCTTGACAGTGGATATACAAATCACATGACCGGGGAGAAAGACATGTTTCATTCATTGCAACTAACTCAAGAAGCACAAGAAATTGTGTTTGGAGATAGTGGCAAGAGTAGGGTGATTGGTATTGGTAAAATTCCTATCTCTGACCAACAATCACTTTCAAATGTTTTATTGGTAGATTCTTTAAGCTACAATTTGTTGTCCATTTCACAACTTTGTGGAATGGGTTATAATTGTTTATTTTCTGATGTGGATGTGAAGATCCTTAGAAGGGAGGACTCCTCAGTTGCCTTTACGGGTCTCTTGAAGGGCAAGCTTTATCTTGTTGATTTCACAACAAGTAAAGTGACGCCTGAGACTTGTTTAGTGGCAAAGTCCGACAAGGGTTGGCTATGGCATCGCCGGCTAGCCCATGTCGGTATGAGGAATTTGGCCAAACTTCAAAAGGATAATCACATCATTGGACTAACAAATGTTGTATTTGAGAAAGATAGGATTTGTGGCGCATGCCAAGCAGGAAAGCAACATGGAGTCCCGCATCAATCAAAGAATGTGGTCACAACAAAGAGGCCATTGGAGCTTCTTCACATGGACCTCTTCGGACCTGTGGCTCATTGGTGGTAGTAAGTATGGTTTAGTCATTGTTGATGATTTTTCTCGATTCACCTGGGTTTTCTTTTTAAGTGATAAAGGTGAAACTCAAGAGATATTGAAGAAATTCATGAGGAGAGCTCAAAATGAGTTTGAGCTCAAAATCAAGAAAGTGAGAAGTGATAATGGGACGGAATTCAAGAACACAGATGTCGAAGAATTTTTAGGCGAAGAAGGAATCAAACATGAGTTCTCAGTTCCTTACACTCCACAACAAAATGGTGTTGTGGAGAGAAAGAACCGAACTTTAATTGAAGCAGCTAGAACCATGTTGGATGAGTACAAGACACCTGTCAATTTTTGGGCAGAGGCGGTCAACACCGCCTGTCATGCAATCAACCGTCTCTATCTTCATAAGATCTACAAAAAGACTGCTTATGAGCTTCTCACTGGTAACAAACCTAAAGTTGATTATTTTAGAGTATTTGGTTGTAAGTGTTTTATTCTTAACAAGAAAGTCAAGAGCTCAAAGTTTGCTCCTAGAGTGGACGAGGGTTTCTTGCTTGGTTATGCATCAAATGCGCATGGATATCGCGTTTTCAACAATACCACCGGTCTTGTTGAAATAGCGATAGACGTGACATTTGATGAGTCTAATGGCTCGCAAGGGCATGTTTCTAATGGCACTGCAGGAAATGAAGAACTACCTTGTGAGGCCATAAAGAAACTTGCAATAGGTGAAGTGAGACCTCAAGAAAAGAATGATGAGGAAGGAACTTTGTGGATGACCAATGAGGTTGTTGATGTGAGTGCAAAGGTGGTGGGTGACAAATCCTCCACCCAAGCAAACCCATCAACCTCAAGTCATCCAAGTCTTGAAGAAAATCATCAATTCCAAAGGATGCCAACTATGGTAGAAGATGAACACGAAAGTGTGGATGGTGAAGTGCCTCTTGATCAAGTGAATGATAAGGAGGAGCAAATATAAAGACAACCATCAGTGCCTCATCCTAGAGTCCATCATACCATTCAAAGGGACCATCTGGTGGACATCATCCTGGGTAGCATCAGGAGAGGGGTAACAACTCGATCTCGTTTAGCTAATTTTTGTGAATTTTACTCGTTTGTTTCCTCTCTTGAGCCACTTAAGGTTGAAGAAGCATTGGGTGATCCGGATTGGATAATCGCCATGCAAGAGGAGTTGAACAACTTCACCCGGAATGAAGTCTGGTCCTTAGTCCAAAGACCCAAAAAAATGTGATTGGGACTAAATGGGTCTTTAGGAACAAACAAGATGAACATGGCGTAGTTACAAGAAACAAAGCACGGTTGGTTGCCCAAGGCTACACTCAAGTGGAAGGACTTGATTTTGGTGAAACATATGCGCCGATAGCAAGGTTAGAGTCAATTAGAATATTAATTGCCTATGCTACTAACCATGATTTCAAGCTATATCAAATGGATGTCAAGAGCGCATTTCTAAATGGACCACTACAAGAGAGGGTCTATGTGGAGCAACCACCGGGCTTTGAAGATCCAAAGAAGCCAAACCATGTTTATTTACTTCACAAAGCACTCTACGGGCTTAAACAAGCCCCTAGAGCCTGGTATGACTGCCTTAAAGATTTTTTTAATTAAGAATGGGTTTACAATAGGAAAAGCTGACTCTACATTATTTACTCGCAAAGTTGACAATGAATTATTTGTGTGCCAAATATATGTTGATGACATTATATTTGGTAGTACTAATGAAAAATTTTGTGAGGAGTTTAGCAAAGTAATGACAAACAGGTTTGAGATGTCTATGATGGGCGAGCTTAAATACTTCCTGGGATTTCAAGTCAAACAGCTCAAGGAGATATCGAATTGTTTCATGTGAGCACCGAAAATCAACTAGCCGATATCTTCACAAAACCCCTCGATGAGACTAGGTTTTGCTTTCTTAGGAGTGAATTAAATATCTTGGATTCTCGTAACTTAACTTGAAAAAGGACACACAGATTGTTTGATTCACTACTCTTGATTGATGATTTTAAAGCTTTGTGATGATCTTTCAAAACTTGGTGGAAATGTTAAAATTCGAATGAATTTTAGTGCTAAGATTTTTTTTTAAAAAAACTTTTCACTTGAAGCTCAACTGGCTTGGCTTTCTGAACTTTCTAGTTCAACTGGAGGAGGCCAGCTGAACTTCGGAGCTGAACTTCAGCTCGGCTGCTTTTAACCAGATTTTACCAGGACCATCCTGGTCAAAACCGGTTGAACCGGCATTGAAGACCGGTTGAACCGGTCTCTGACCGGTCTGTCCGCGCAGGTCTGTCCGCGCAGGCCTGTCCGCGCAGGCCTGTCAGCTGCAGCTTTTTCAGTTGCACAGTATTTGTCTTTGTTTTTTCCGCGCAGCAAGCCATGCAGATTTATTATTTTCTTCTCTGGCCGGTCAGAGCCGAGCAGCAGCGTCAGGCGCGCGGAATATTTTTTTGAACCGGTCGCATCCTGCAGGGCGCGCAAATTTTTTATTTAGGGCTGTCCGGTCACAGCCAGGCATGCAGTGTCACGGCGCCATGCAGAATTTTTTTTGGTCACAACCATGCTCTAGTGTCAGGGCGCGTAGTTTTTGGTTTGCCCGGTCGCATCCGGTCTGCTGCAGGGCGCGCAGAAATTGTCTTTTGGCCTGTCCGGTCGCTGTTCACTCACCTGTGTCAGCCCACGCAACAACTTTTTTTGTGACCGATCGCATCTGCACAAAGCTTTGTCAGCCCGAGCAGCAAGAATTTTTTTTTGTTTTTCTGTTCAATCACTTCCGCGCAGAAGTGGCAGGGCGGCAGGGCGCTCAGAAATTTGTTTTCTTTTGGCCTGTCCGGTCGCCTCCACGCAGCAGTCTGCTCGCACAAAATTATTTTTATTGGGGCCACGTAGCTTTTTCTTTTGCCGTCCGGTTCCATCTGCGCAGAACTTTTTGACTTATGTGGCGCATGGAATCTCTGCAGTACACTTTTTCAAAATCCTGCGTCAGAAATTTGAATTTTTTGAATTTTTGGCGTTTTCCCTCTGCCGTTTCTGGCTCTCTGACCGCACTGGTCAGAACCAGTTGAACTGGCCTATGGGGGGCGGTTGAACCGGTCCTGAGGACCAGCTCAACTGCCTATATATATTTGACACTCCCCTTCTCTTCTCTCACTCAAACCACTCTCGCCTTAGTCTCTTCGCAACATTTAGAGCCTTCCCTCTCACTCACCTAAACCGGTTCACATTTCCAAAACCTCTTCAATCATGGTGCGTCGTCGCAACCCTTCCGTGATTGAGTTTAGCAGCGACTCGGACAAGATTCAGGAGGAATCTCTCGTCCCCTCTCCTCTGCCTCGTCGCTCCCTCTCTAAGAGAGGACGTGGCTCTGGTAGGATGGATGGAGAGGGTTCTTCCATGCCCTCGCAGCTGACTAGCGGGAGACGCCAAAAGGTTGGCGCCTCTCGACCTCGTCGTAGCACGTCTTCCTCTGGGTATGCTCCCTGCCAGGAGGAGGATGGGGCTTTCGACGACTTCATCGACCTTCCGGCGCCAGATGCTCTCTACGTCACCGGACTGCACATGCATCCGGCGAACGTCACCCGCGGTCCACATGAACCTCCTGTCGACTTCTCTCTGAAGGGAATAGACACCCTTCAGCGCCTGAGGTTCACCAACCCTACTCTCACTCCTCGTCATCCTGGTGTCCAGGACCCTCGGTTCTGGAATCTTTTCCAGGCTGACTTCTACAATTCTGTCATTCTTTCTAAGAAACACCCAGTCGTCAGACATAGGTTTATAGACTGGGAGGGTTGTGAGAATATGGGAGATACTGACATGACGTCAGCTCTCCAAAACTTAGAGAGAAAGGGGTTACAAAACATTATGACCATGGAATACCCCTGGAACGACGAGGTGGTAGCTCAGTTCTATGCCACCCTCTGGATAAAGGAAGTAGATGAGGAAGCTGATGGGTATGATTACCCAGTCATGTATTTCTTCATCAAGGGTAACTGGCACAAAGTCAGTTATCGTCGTTTTGCCCATATTCTAGGCTTCTCTGATGCAGATATTCTGGGCAGCAACATGTGAGTGCATAACATTAGGCTGCCCATGCGGGAAGAGACGAAATTTATTCATGTCTCTACTGAGCGGGAGTTCTGGACGACCGCTAACATGCATAGGTACTACAGGTACCTTAATTCTTTGTCCAGGATGACCGTTCTCCCGATGGGCGGAAATCAAATGAACATCCTTGGGGAGAGTCAAGTCTTCCTCCTCCTCCTTGCTCCAAACAGTCTCGCCCGCATCAATGTGTTTGACATGATCTGGGAAGAGATCGTTCGTGCTTCCTAGTCTCCTCTCAGAGGGTGTCTTCATGCACCCTTTATCATGAAGATGATCGAAGTGGTCACCCAGACCCACTTCGAAAAACCCGTCAAACACTCCCGCTACGTACCCTACTGGGTTGATTCCTCCAACCTAGCTGCTCGTACGAAGCGGGCTCCTTCAGGGTCTGGAGGCCCTGCCTCCTCAACTGAGCCACCTCCCCAGCCTCCCCATCATCCATTTTCCTCTCGCCCCGCCGCTGCCTCGCGTCCCTCTCCTGGCCGCGGATCTCCCATGCCACGTGCTCGTGGTCGTGGTCGAGGTCGTGGCCGAGGGATGGGTGCTCGCTTGGTCCATGGGTTTGCGGCGTTTTTCTCCATGTGCCGCAACATTTCTGCTGATGTCCATGAGGTGGCACGGCGTCAGCGGGAGACTGACGACAATCTTCGTCGTCAAGCTTCCACCATGGGTATGCCATTTGTCCCTCGCTCGCCTGACGTGCCTCTCCATCCTCCTCCCCCGGAGATCAATGAGTGGCACCAGCAGGCCTACGGGGTGCCATTTATGTCTGCAGACGATGAAGAAGAAGAAGCCTATGTTGATGATCGCGAGGAGTTTGCCCCGCCTCCCTACCATGGGGATCCGGGTCAATCATCCTCTCACCCTCCGCCTCCCTACCCGGGTCCGGGACCCATTCCTGGTGACCCTAGGCCATCCGGCTCTGGTTACCATCACCACCCTCCACCTCCTCCTGATCATCGTACATCTTCGGCATCTGTTTTTTCTACGGCTGATGAGCCACCTCGGGACTCCACATTTGAGGAGGACGTGATGAGGACCCTCTTCCCTGACTACACTTCATACCCTCCATCATATCCTCCGTCATATCCTCCTCCTGGAGGATATTGATTCATGTATCCTGATCTCCCTCTTCGTTTTTGGTACTTGTTGCCAAAAAGGGGGAGAAAATCCAACTTCTAGTTTGTGGTTCCTTCGGTTTTCTTTTGGATTATTTTTTGTTTTGGCCACTGTTGTGGCTTGTATCCTTTTACCATTGTATGAAATAAAACTGTTGGGTTGTTATTCAAAGTAATATGATGGTCTTCTGGCCATCATCTTTTGCCCTATCAGTTTGTAATGATTACTTTGCTATGAAGTAAAAATGGTGTGTATTTAGAGATAGTCATGTGCATCACAACTCCTATTATGACACTCTTGCACCCCTCAGATGGTTGTATTCAGTATGGTTTTTGTCACTATCTCTAAATGTGTTGCTCACATGACATATTATGTCAAAACATTGGATTCACAATCGTGCACACATTTAGGGGGAGCCATCTACATACAATCATTCAAAAGAATTTTTCTATTGCAAATGTATATCTTTTTGACCTTAATTGTTTTGGCATCAATCACCAAAAAGGGGGGGATTGTAAGTGCAATCAACCCCAAATGAGGGTTTTGGTGATTTAATGACAAAACAGATAAGGGTACTAACAGTTTTTGTTCTCAGTGTATGATTAGAAGGAAACAAGAACTTAAGGAAGCACCAAGGACTTCATGCAACGGACGTATAAAAATTTATGTAAAATCTTGTGTTGCATGATGTAATTCTTCCTTGTTCTTTTCCGCACAGGTATAGGTGTTTGCTATAGCTCAAGTCCTCTAATTCAAAAGCTATTTTTGAAAACCAAAAATCACTTTAACATTATTTGGTTGACCATGGTTAGGGTTGAGAAACCTAGAGTGTTCTTGCTGAAAAGCAGCTGAACTTCTTCAACTGCAGCTGAGCTTTCCAGCTGAACTTAACTTCAGCTGTGATGAGCTTCTTCAGCTGCAGCTGAGCTTTTCAGCTGAGCTGGACTTCAGCTGAGTTGAACTTTCAACTCCAGCTGAACTTCAACTTCAGCTGTAGACCTCTTTGACCATACCTAGCTGAACTTGCCCCCGGGCAGCTGAGCTGGGGTTTTCTCTCCAAAACTCTCTGGTCTAGACCAGCTGAACTGGTCCTGGGGGCAGCTCAACTGGTCTCTGACCCTCTGACTTCTGATCAGCAGTCTGCCAGTCAGACTGGCAGTTAGGTCGCCAGGGGTGTCAGGGGGCGGTTCAACCGCCCTAGGGGGCGGTTCAACCGGTCCTGGCCTGACCCGCATGCAGGTCAGTCTGCCAGTCAGTCTGGCAGTCAGACTGGCAGACAGTCTGCTGACCTGTCAGGGGGCGGTTCAACCGCCCTAGGGGGCGGTTCAACCGGTTTTCAGCGGGGTTTTTGGTCAAACGGCTAGTTTTTGAGCCGTGACTATAAATAGAACCCTCTCTCTCTCTTCTTCAATACGGCTGAACACTCACTCATTTATTGCTCACCTAACTTGAGACAAAGCACTCATCTCTCCATCTCTCTCACACTTAAATCTTCCTCAAGATTCAACCTTGTTGGAGGAGCTCTTGGGTGTGAGGTTTGTGCTTGTGCTCACGATTTGGTTTTCCTCTCCCATCTCTTTTACAAAGACTTGAGCTTGTGCAAACCCCTCTTGTGCTTTCTTGGTGTTCTTGAAACCCTAGGTTTCAAGATCTTTTGAGGTTGCTTGGGAGTCTCCAAAGTTGTGGACGACCCCAAGAAGTTTGTATCACCCGCTCTTTGAGCTAAGATAAGAAGAGATTGCCTTGACCTTTGTGGTCAGCTTTTGGAGGATTAGGGTTGAAAAAGACCCGGCCCTTTGTGGGCTCCTCAACGGGGAGTAGGACACCTTTGTGGTGTGGCCGAACCTCGGATTAAATCTTGTGTCTTGTGTTCTTGCTTGTTTTAACTTGAGATATTTTTACTTGCAAGATTGACATAAAGATTTTTCCGTAGAGTTTATCTAGAAGTGAAATAGCCTTTACATATTCGTCGTTTCATCCTAACTTAGCTATATCTTACTTCTCACAAGAACTCGCGAAATACTTTTTCGAGTAGATTTTTTGTCTAAAGTTTATAAAACAGTTGGATTGGTTCAGAGTGGTTCAACTTGCACTCGTGGCTGTTGAACCGGCCTGCACCAGTCGAACTGTGTATTTAACAATCTTTCGAAGTTTGTTGTGAAAATTTTCAGATTAGCCTATTCACCCCCCCTCTAGGCTACTTTCATTATTGGTTTATTGTGAACCTTTGGCACCTGTAAAACTTATAGTCTAGAGCAAACTAGTTAGTCCAATTATTTGTGTTGGGCAATTCAACCACCAAAATCAATTAGGAAAAGGTGTAAGCCTATTTACCTTTCATATAGCCTAATAACCTACATTGTTTAATCCTAGTTTAGTTTGCTGAACTTTGTTTTCAGAAAATTTGCTCTAAACATTTTAAGATTTATGCTAGTTTGACAAATAAAAAGTAAAAATACTCTCTCAAAAGTCAAAACAAATGCAGAGTTATATTACATAATTTGTATAATTTTACGACTCATAAAATAATATAACAGTGTCTCCAGTCATATCATTTTGAATAAAGTTGTGTAAATCATCAATATTTTTATTTTTCAGATATTATTTAGAAGTTGAAAAATAGCAGACTCAAACTTTTGTATAGGGCCTTCATGTAAATTTTCGCTTGGGGCCCCCAAAATGTCAGGACCGGCCCTGCATGGCACGGTTCGTAGTGGATCGCAAGCTTAGACTGAGGTTATTTGTTCGTTTTCGTCCCTATCGCCATGCATCGTTTCTCCCTCTCCCTTTCTTGTTCGTCTTCTCTATTATCTCTTTAATGAGGGAAACACGATATGTGGGTAAAGATGTGCAACATCAACAGAGTGTAAACTGGTATATCAATTGTGCTCACGACCATGAGCGGTTCGGACCCTCACATTAGTAAATTATGGAATTAAATTTAACACGTCATGCATTGCATTGTGGGAACTGTTATTTGTGATCTCTTATTTAATTGGGTTAGTATTTACTTATACTTAGAAATTGCTCAATTTTTTTTGACCAACTTAATTAAAACAATGCTCAGCCTCGTTCATTTTCTTGGTAAGCCTTACATTGCACTTATTGAGCCCCAACCATAACTAAGCTCATACATATTTTTTCCCACACTTACTGAGCGATGATCGTATGGGAGCTCACCCTTGCTATAATCACACCAGGTCAAGAACAGGTACCACTGGAGGAGGATCACCATGAAGATGAGTTCAGAGGGGTCTAGGCCGTCGTCTCCTAGTCAACTATGGTTGTGGAGGATCGTCGTTGTCATATGGTATAATTTATTTAGTTATTTTGTACAAAAAATGTTACCATATGTGTGAAACTTGATTCTGTTACACATGTGAGATGCATCTGGATTTTCCTTAAATTCGGGTGTGACAGTGCAGTAGATATAGTGCACTGAACAGGTTACACAACGAGTGTTGCAGATAAGTTATATACTCGCCGGACATATCCGGTGAGACACTGGATACGACATCGGACAAGACTTCTAAAGATAGTTGCACATGAGTGTGAGAGGTAGGGATAGAAAAAAATAGCTCATGGCACCAGGGCCGACCCTTAGGGTGTGCATGGTATGCGATCGCTCATGGCCCCCAAAGTACGTAGAGCCCCCATTTATCCATCTTATATCCAAATACATATACGTTACAGATAAATTTCTATTCTCATTCTTTTGTTGTGTCTCTTTGATTTTCTGGCAATCATCTCTTCGAGATCGCGCAACAGTCTGCTTGATGACAAATACTAAAATTGACGTGGGAAGGATGGAAGGTGTTTTACACCTGTGCCTGTCGCGTTGCTTAATTCTAGCCTACCGGAGGGCCGACCCTAAGGCTAGCTCCAACAACGTCCTGTAAACGGTACTGAACCCTACATTTACCGCAAACTTCCTTTATACCTCCAACAAGGTCCGGTATATTGTGCTCTCAATTTTAGCGTTCAGTGATGAAGTGCTAAATGTGGCGGACGAGCTGCATCCTCTATACGCTAGTCAATCATCCTTGTGCGAGCATAGACGGATCCCGAGGAGATTTTTTTCATGGAGAAGTATTACTTTTTATGTTTTTCTGTGGATGCTCCACGCATGTGACAAAAAAACAGTAAAATAAGAAGTGTAAAAGTTTAATTTCAGTACAAACGAACGGAATTTTGTTTGCGTGTACGGTAAGTTGTCAAAAATTGTCATAGTGACCAATGTTAATTACGGTTTTTTAAATAGGTACGGTAAGTGAAATTCAAATTATATGCAAAAGGTAGAACAAATATTGTAGTTTATGGATTCTGTTAACACTGTAGATATAGAGGACCGAATTTAGGGGATGTTGCTGGAGATGAAGAAGATATAGAGAACAGAATCTTTTAGAGAGTACTGTAAAGGACAGAGAATATTCTTTTAGAGAACGGAATTTAGGGGACCTTGCTGGAGATAGCCTAATAACCTACATTGTTTAATCCTAGTTTAGTTTGCTGAACTTTGTTTTCAGAAAATTTGCTCTAAGCATTTTAAGATTTATGCTAGTTTGACAAATTAAAAGTAAAAATACTCTCTCAAAAGTCAAAACAAATGCAAAGTTATATTACATAATTTGTATAATTTTACGACTCATAAAATAATGTAACAGTGTCTCCAGTCATATCATTTTGAATAAAGTTGTGTAAATCATCAAGTTTTTTTATTTTTCAGATATTATTTAGAAGTTGAAAAAATAGCAGATTCAAACTTTTGTATAGGGCCTTCATGTAAATTTTCACTTGGGCCCCCGAAATGTCAGAACCGGCCCTGCATGCCACGGTTCGTAGTGGATCGCAAGCTTAGACTGAGCCAAACCTGTTTCTCAAGATTGCCACATGGTTGAGCTAAGCGGAGCTAGCTCGATGCCGATTGCAAGCTGAGCTCAGATTGGTCCAACAGAGGCCAAGTACCTAAGGGCGTGATTGGTTGCCGGATGAAGCCAACCTGGCTCGGCCCCCAGATCCAGGCTCATGGGAGCCTGGTTCCAGTGGTGCAAGGGGGGCGTGATTGGTTGCCGGATAAGCAGAACCAGGCTCAACCCAGCAGTTGATTGGTTGCACAGTTGCATTAGCAAATTCAAACACACAGATGCAAGTAGAGTGATTGGTTGCAGCTGCTTACAAGTGGTATGGTCACCACTGCATACAAGTGGTGAGGTTACCCAAACCAAATGATCAACAAGAGCACAACTATCAGTACAAACTGAGAGAGAATCAAACCATCTTTCCAACTGAAACCAACACGATCAGTCTTAATGACTGATTGTTTAAAAGCAGCAACAGCTTCATCTTTTGAGGAGTAGCTTTTGTAGCATGCTCCTTTGAAACCATGAACTTGAGCATGACAAGAATCCCAGGAATCAAACACACCAACTTTTCTTCCATGGTGAACAACATACCACTTCATGTTGGCTATATAGGAAGAAACAAACATCATTGACATTGCACATTTGCTCAGATCATAAGTAGTATAATGCAGCAAGTCATTGGCATATGGGCATTTGCTCAGATTAAAGAAGCAAGAAAAGCATTCAGAATAGTTCATCATGAAGACAGAAGCCAAATGGTTCATACAGGAAATGGTTCATACAACAAGCCAAAGTGCAAGTCAGCAGGATCAGGTGATCAACTGCACAAAAGAATGCTGGTTAGGAAGCTACCAGTCAGCACTACTTCCTACTTACTCAAAAGAACTAACCAGCCTCAACTCAGAACACTCTAAGTGTAGTAATGTTTGCCCAAGAAAGTCCTCATCCAGAGGACTCTATGTGAGTCAGTCATGTGCACAAATGCAGTGCCATGGTTCTTGTGATCCAGCAGATGTGAGTATGCAACAATCAAAGCTTCTTCAGTGAAGCCAGGCATAAACATAACAGCAGCATAGAGGTCAGGATGCACATCTTGTTGCTTGGTCTCCCTAATTGCTTCAGCCACATTGTTAACAGCAGCTGTCATACCACCAAAGACATGGACATCATCATCAGTTAGCACTGACCTCTTCCTTTTTTCCCAGACCCACTTCCAGAACCACCACTAACAGGCACCTCCTTCTTGCTGCCATCCTCTGCCCCCTTAACATTCTCATCAGACTTAAATGAGCTCTCAGCAAAGTCTGAAGGAGAACCCAGTGGCTCACTGGAGCCCATGGCAAACTTGCCAGTGGCCAGGCCATTGCCAAAGATGACCTCCATTTGCTTGTAGTTCTCTATTGGTTTATTAAGTAGTTCAGCATCCTTAGGGTGAGCCTAGATGGTTATGCAAGAAATGAAGTCAGAATTGTAAAGGAGGTTAGAATGAAATGCAAGTCAGGAACTTAAACCTAACAAACCTTAACATGGCCATTATAGTGTTCTTCTTCTAGGGTGATCATATGGTTGTCCTCATCAAACAGAGCCCCACTAAGCTCTCTAAGTTTGATGATCCTAACCCATCTCTGCCTCCACTTGCGCAGGTGATTATACACCTGAGTACCAGTGACATCATGCCCACTGAACTCCTTGAGGTCCTTTGCAACCTGATTAAGGTGCACCTCCTTGAACCCTTTGTCAGTCCTAATCCCAGTTGAAATCAGCTGACAGAAACGCTTAAGAACAAAAGCAGACATCACATTTGTCCATCTCATGGACTGCCTGGAGGCTGCAGGGTTTGGAACAGACTCAGTACCATCAGCCTGAGTCTCAGGGAGCAACTGGGTCATTTTGTCAGTGTAGGACACAACTAAATCATAGACATCTTGCTCAGTATCCATGTCCTAGTGCAAAAATAGCAAAGACAACTCATATATATGAAGAATAAAGCAATCATCATAGCAAATGCCCATGCCTCATGCAGAAAATGTATCATCACAGCAAATGTACAATGTCTCAGCAATGCCTCATACAGTTAAAAAAAGTGTCATCACAGCTTATGCTACCAAATGCCTCATACAACTTAAAAAAAGTATCATCATAGCTTATGCTACCAAATGCCTCATACAACTTAAAAGGGGGTTCTTCTCATGCGCAACATGTTAACACACTCCACATCATTGAAGTTGTAAATTCTATTCAAGTTTCTGTGCCTCTCTTCATCCATTGCACTCATTGGGGCATAAGTGATCCTAGGTCGTGCCATTTGACACTGGACAAGTCTAGAAACAAAATATGCATACACAGCAGCAACAAGTCCAGCTGCCCTCACTATTAACATGCGCCTCTCATCATCATCATCCATCTCTAGATGCAAAGGAACAGAGCCAAGTATCACATAACAGCTCACCAGCAGTAACATGATGCAGTACTAGATCATAGGCATGGCAACAAAATAAGCTCAGATCTAGGGACCAAATCATAGGCATGGCAACAATAAGCTCAGATTAAATGCTATAAAATCACAGGCATGGCAACTTAAGCTCAGACCATTCAAGAAATGAAGGGGGAAAGGGGAGGGGAAGGAATGGCGTCTGTTATACCGAACAGCAACAGGAAGAAGGACGAAGCACCGAGAGGAAGCAGCTACAAGGCGTGAGGAAGCAGCAATCCGAATGCCCTGTCAAGAAATCAAAAGGGGATTAGAAACCCTAGGTCGCGAAATAGAACAAAAAACGAAGCTCCAAACCTACAAAGCACGAATCCACCACGGATCTGAACAAAGAAGAAGAAATAGGGGGAGGATGGGAGGAACTCGGCACCAGTCACCGGCGGCGACAACAGAGCAGGGCGCGCGGAGGAAGGGGAAGGAGGTGGCGGTATATAACAGCGGTCACGGCGGGAAGGAGGGAGGTAACCACAGCTTTCGCGCCATCTCCCGCCATCTCCCGCGCTCACGCGCGCGCTCGCGCGCGCTCGCCCCCGCGTCCCCGCGCTGCACGCCGCCAGCCTGCATCCGGAGAAACGGCCGATTCGGCCGCATGGGGAGAGCCAGGATCGCGGCTAGGGTTTGCATCGGGAGAGCCAGGCCGGGAGGGAGGAAAACGCAACCGAACAAGTGGGCTTGTATCCCGGGAGCCTGGCCCAGGGTGGATTCGGGCAACCAATCACGCCCTAAGACATAGTCTAACTGGCATCGTGGTTGTTGTGTGTTTCTAGCTAGCCTCTAACCAACCAGTACACGTGAAATTATTGACGTACTCAGACGTGCGATCCGCGAACTTCGTCTGATGAACGCTCGGATAGATATATATATTTTATATGTTTTCGAATATTTGTTGTCCACTAGTTTATTTTTTAACTAAAATATAATGAATAATATTTTTGAATATTTGTTGTCCGCTAATTTATTTTTCAACTAAAATGTGATAAATAAAAAAAAAGAATAGACAAGATAATACAAAGTGAATAAACTAAAGCGAGGTAGCATTGGGCGTGTTCGGCAGGCTGCAAGCCGACACTGTTGCAGCTGTTTGGACTGCTGCAGCTGCAATCCATAGAGAGAAAAATACTGTAGAAGCCGCAGCCGCAGCCGGATTGCAGCCGCAGCAAGCCGCAGCGAACAAGCCGATTGTCGTGCTGCTGCGTGCGTGGGATGTCTCCTGTCTGTCCATACCAAATTACTAATATCTGAAATAAAAAATTATATATGTTGATAAGAAAAAAGGAAAAAAACTATAGTAAAATGTACTCTATGTTGAAATGTGATGCCTAGCCTCGACGACGTTCTTTCCGTGCGTGTCGTCGCTGTGGGGGCACTGCCGATGCGCGTAGATGTTCTTCTCGGGCAGTGGGTGGAGCACGGCGCTCGCCGGAGAAAACATGTTTGAGCCGTGCGTCACGCGGCGACCATGGCCCCTGGCCCTGGGAGTATGTAGTAGTATGCTTACTTGCTACGGAAACGAGCCGGTAGCCCGCCCGTGCCAGTCACGGGACTAGGGATAGCAGTGGGTAAGGTATCCAGTGGGCACATGAAGTCATACCCATATTCACTAGAAAAAAATTAATTTATACTCATACCCACTACCCATCATGATTATAAAATTACGCCATATCCATACCCATTATGGGTAGCGGGTACCCATCGGGTACCCATACCCATTAACATTACAATAAACATGATCATTTGAATCGCAAAAATAAGTACCAGTATAATAAATGTATCGAACCCAACATAATTTTATCACAGATTCAACAATATTCAACATTAAATGGCAATATCATGCCTTCAAATCACACAAAACATAGCATTGATTTTACAAAACTGTGGACAAAACATAGCATTCATTTTATATGAGACATATAAATCCGGACCCACGTGTCATATATTAGCGGGTTCCTTATCAGGTAGCGGGTATGGGGCAACAGAGAAAATAACCACACTCACCGTACCAGATAGGCATAATAAATGACCCAATAAAATACACATGGGTATCAAAATCCGCTATACCCGTGCCATAATATGATTTTTACCCATCGAGTTTCGAGTTTCGGATACCCATTGCCATCCTTACACGGACCGGGGGGTTTGCCGGGCAGCCGGCGCCGCCGTCACGCCATCCGCGCAGGCCGCGGCCCCTCTCACGCCTAGCAGCCGCTGTAGCTGCAGCAGGCCCGCCCGCCCAGTGACACGCTATTTATTCTACTATATTAGTATTACCAGGCCATCCGCGCCCCATTAGCAACTGCGCGCACCTCCTGACCTCCACTCCAGCCGTCTCCTCCCTCCGCGTCCGGCTCCCGGTTCCGTTCCGCTTCCGCGATGGGGCAGCAAGAGGGCCACAGCAACCCCCACCCCGACGAGGCGGCGCTCCCGACGACCACCGCCGCACCGTCGCCGTCGCCGTCGCCGTCGCCGTCGCCGTCCTCGTCCAGGCTCTACACGGCGTGGCTGGTGGCGTCGTGGTACGCGTCCAACATCGGCGTGCTGCTGCTCAACAAGTACCTGCTCTCCGTGTACGGCTTCCGCTTCCCGCTGCTGCTGACCGCGTGCCACATGTCCGCCTGCGCCGTCCTCTCCACGCTCGCCCAGCACGCCTCGCCCCGCCCCCGCAGCAGCTCCTCCCCGCGCTCCCACCGCCAGCTCGCTCGCGTCGCGCTCCTCGGCGCCGTCTTCTGCGCGTCCGTCGTCGCCGGGAACGTCTCGCTCCGCCACCTCCCCGTCTCCTTCAACCAGGCCGTCGGCGCCACCACGCCTTTCTTCACCGCGCTCCTCGCCTACGCCGTCGCCGCCCGCCGCGAGGCCTGCGCCACCTACGCCGCGCTCGTCCCCGTCGTCGCCGGCGTCGCCATCGCCACCGGGGTGAGGACGCTAGGCTAGTAGCAACGCGGCGGATGGGGATTAGAGTTCCTGTCGACACTGCCTCCCGATAAATTTACCTTTTTTTTCGCAATGCAGGGCGAGCCGAGCTTCCACCTCTTCGGATTCGTCATGTGCGTGGCCGCCACGGTGGGGCGCGCGCTCAAGACCGTGCTGCAGGGCATCCTGCTATCGTCCGAGGAGTGAGCACGAGCAATTTCCGTTCCCCCTTCGCTCCCCATTGTCGCTCCACTAGCCGCTGGAAACATTCCCATCCCGCGCAGGGAGAAGATGGACTCCATGGACCTCCTCCGCTACATGGCGCCGGTGGCCGTGCTGCTGCTGGTACCGGCGACGCTGGCCATGGAGCGGGACGCGTTCGGCGTGGTGGCGGGCCTCGCGCGGGAGGACCCCAGCTTCCTCTGGCTCCTGCTCTGCAACTCCTGCCTCGCCTACTTCGTCAACCTCACCAACTTCCTCGTCACAAAGCACACCAGCCCGCTCACTCTCCAGGTACGACGCCGGCCTCCACATGCCGGCCGCCGTCCTGTGACGTACGCACCGATGCCCTTTGACATCCCTCATCCCCCCCGCCCTTCGCACGCAGGTCCTCGGGAACGCCAAAGGCGCCGTCGCCGTCGTCGTCTCCATCCTCATCTTCAGGAACCCCGTCACGGTCGTGGGCATGCTGGGCTACGGGGTCACGGTCGCTGGCGTCGTCTTGTACGGCGAGGCCAAGAAGAGGAGCAAGTGAAGCGACGGACAGGCTCCTGTACCTGTACATTTGCTGTCCGGCCTTGGAATGCAACGCCGGCGGTGATAGCTAGGGAGAAGACAAGCAAGCAAAGCGACAGCTGGGGTGGTTTGTTCGTCGCATCGTTTCAGAAGGAGGTGAACCATACAATACAAGAACTCGTGTCGTAGCACAGAGAGAGAGAGTGTGTGTCGTTTAGCATTTGAACGAGCTTGCCGTCCAAAGCCTACGATTTTAATCCGGCCTGGAGAAGAAGAAAAAAAACGAAGCTAGCTAATTATCAAACCAACACCTAGCATCTGAACATATGTGCAGCAATCCAACTGGGTTTCCATATTTTAATGTTGTCGATGAATTGGTCCTTCCTAAACAACATTTCCGTTTGAACATTTACGAAGAGCCAAAGCCAACGGACTGCATTGTTCACTGAAATGTTTTTTATAGTGTACACTTTCATCCACCTTTTTCCTTAATTACTAAATTAACGGCGTCACCACAAGTTAAGAGTGAGATGCTAAAAGAGATGCCATGAATTATCACATTGTTTTATGTATAAATGTAGTGCCAATTGTGATTATACCCTTCTGGAGCATAACAATTGTGTTTCTACTTTGTGGTTTTACCCTTCTGTTTTGAAAACGAAGCTATCGTTTGTTCCTCTCTCTTCCTCTCTCTGTTAGTGAAAATAATAATTTATGATTAAAAAATAAAAAACAAAGATAAAATCACAATTATAAGTTAAGAGAAGGGTATAATTATAATTAATCTTGAAAATAATTTTGAAAATTTTATGATTAATTATGATTATACCCTCATCGCACCTTGTAATTGTGATTTTGCTCTTTAGTTTTTTTTATCACTAACTCTGATTTTCACCAATGGTGCTAAGCAGTAGGAGCAAACAACGACTTTATTTTCAAAACAAAAGAATAAAATCATAAGTTAAAAAAATGAGGGATAAAATCACAATTGTAAGGTTAAATTAATAATGGTATAATCATAATTGCCTATAAAATGTATGCGTTCACGGAAATGAAGTTATGGCTGTCTCCAGCCATACACGCTAAAGGATCTTGATGCCACTGCCAAACGACTTGTCATTCTTTATTATAGTTTGGTTCTTTATCGTACCAATGTGACCAATGATGCATAGAGAGAGTATAGATTAAAAAGAATGAGATATAAAGGATGAATTTTAAAGAATATTATTGGAGACAGAGAAAATATAGATGACAATATTTTTTTAGAATGTATTGTAAGTGACAGAGAATATCCTTTTAGTGAATGAAATTTAGATGATGTTGCTAGAGATAGTAAACCCAGAGTGCGGTTCATATAGAAAACCGGTTCAACTTAACGTGTCAAATTGGGTTGATCCTTTGCCTTGCAAATGCATTAGAGTTTCCAATGAATCAGATAACGCAAAGATCCATTAATCCTAAGATCATTTCTGAGAGCCTATTAATCGGCACTAAGATTTATCTATATGATGAATAAAGCCATAGCAATATAAGAGCAACTAACTAGTTAACAAAATTGAGTAAAATTTATTTCTGGTCCCTTAACTTGTACGAGTGTGTCAATTTGATCTATAAACTTTTAAAATGTGTCTCTTAAATTCTTAATCGTTGTTTGTCTCTCTCTTAAGTTCAAACGGATTAGAATTTCATATGCAGCCTAATCCTAAACCATTTGGATCTAAAAGAGTTCTGAATAGTGATTAAGACTTTAAAAAATACACTTTAAAAGTTTATGGACGGAGTTGACATGCACGCAAGTTTAGAGGGACCAGAAATATATTTTACTCTAATTTTTTTTGTGTGTGCCGCAAAATAAGGTAAAAATGGCTTGGCGCAGTGGCACAAAAAGCGCGCCAACAAAACATTTCGAATCCTAGCCCCCCTCAACTGCAAAAGCCCCTCGTCGAGTCGTCGCGGCCCTCCTCCATCAATCCATCGTCGATCGATGGCGAACATGCTGGTCTCGGAGATGCGGCTTCCTCTCTACCTCGCCCACCTACTCGCCGCGCGCCGCCTTGACACCGCCAAGGTCAGCGTCCCGATTCCACCGCCTCCGCTAGAACCCCACCTCCCCTAATAACCGCGCGCGCTAGTGCCTGCGCTCTGTTCCTATCCAGGACGTGCTGTCGCTGCCAGAGGTGGAGCTCATGGCCATTCTCGACGCCGGCCTTCCCACCGCACGCGCTGCCATCGCCCACGTCAGCGAGGCCGCCTGCCCGCCCTGCCAGACGGTACCTCTCATCCCCCCTTATGCTCTGCTTGCGACTCGGTGGTCGATTTGGGGTTTAGTTTGATTCTTTGGGTCGGGGTTTACGCTGCATTTACCCGCCCTGCGCAGGCGCTCTCACTTCTGGAGGAGCGCGTGAGGTTAGGAGGCGGCGGCCGGCTGGCCACCATGCTCTGCGGGTTGGACGAGGCATTGGGCGGAGGAATCCCCACGGGGAAGCTCACTGAGGTCGTTGGGCCCTCGGGGATCGGCAAAACGCAGGTGAACTGCTGTCGGAGGTCCAAAGCTTAATTTGCGCATCCCGTAGTATTAGTGGTGCTGTGGGTGCCCGTATTTACAGTGTTCAGCGTACACAAGTGCGCACCAGGTGTTTGTGAATATGTAAATGCTATGCATCTCCATGAGTTTGGTATGCAATACGGAACCAGAAGGAAATGTACCTTGTGTATATGCTTATGTGATGCTGCATCATGCACCTGCTTGTGTTGTGATTTCTTTTCCAATGATCTTGGACATTTTGATCGAATCTACCGATCAATCACCATTTCTTCGAGCCAAATGGTTATCGAGCCAGAACATTTTGGGCTAAATTTATTGGGATTACCTCATTCTGTTGCTTGCCTTGTATGTTTTTTTTTAAATTTGTTGTTATGGTTTTGCGCATCGCAGTTCTGCCTGAAGCTTGCGTTGTTAGCAGCATTGCCGGAATATTATGGAGGTTTAGACGGTCGAGTTGTGTATATTGACACGGAATTCAAGTTCTCTCCAAGGAGGTTGCTGAAGCTAATGGCTACTGTTTACTACTTCCAAGTATCTTATTTAATATTTTTCTTTATTTGAACTTCTAATGTTGCCATTGATATGTTGCAGGATGATTGAGATTGGTCAGAAAAGCTTTCCTCAAATATTTCGACAAGAAGGCTTGGCGCAGAAGGTTCTCATCCTCCACATTTTAATACTTGGTGTTGACACCTTTTCGGAGGCGCCAATCACTCAAAAGAACCAGCGGCGGTGCTCTCTGGTCAGGCGCGGACGGTCCGCGGCTTGGGGCCGGACGGTCCGCGACCTGGCGCGACGTGGCGGTGCTCTCTGGTCAGACGCGGACGGTCCGCGGCACAGGGCCGGACGGTCCGCGACCTGGTGCAGGAGCTCGGGTTCCCTGCCTGACGGCCGGACGGTCCGCGCTCTAGGGCCGGACGGTCCGCGCGTGCGCAGGGGCGGCGGAAGATCGCCGGCGGCGCCTGGATCTCGCTCCCGGGAGGGACCCCGTCGGGGAGGAGAGATCCTAGGAGTTGTCTAGGCTCGGGCCGACCGACCTAGACTCCTCTAATCGACGTAGAGTCGAGGAGAGGCAGAGAATTTGGGGATTGGAATACTAAACTAGGGCTAAACTAGAACTAGACTAGAACTACTCCTAATTGTGCTGAAAATAAATGCGAGATAGAAGTTGTATTGGTTCGATTATTGGGGGTTCAATCGGCCGTAGCCCTTCATCTATATAAAGGGGGAGGTCTGGATCCGCTTCCAACTGATTTCCGAGTTAATCCCGCGGTTTTAGGTAACAAATCCCGCGAGAAACTAGGAACCCTAACTGACTCTGCGCACGCGCGGACCGTCCGCGCCACCACCGCGGACGGTCCGGACCGCGGACCGTCCGGCCTCCGGGCTGGACCGTCCGCACGGTCATTTTGGGTTCCAACATATGCCCCCCTGCCTTTTGGTGAAGGTTGACGAACCAAAAGCATATGAACTAAACCTGATGTAAGTCACCGACTTTTCGATATGGAGATTATTCAATAAAGCACCAATATAAAGGCCGTTTCGGATTGTATCTTTCTCGGCCATGACCATTTGATCAATGGATCAAAAGGAATAGAATGGAGGTGCCCCCCAGTCTGGATAGACGAAGGGACTATACATGTACCATGGATTCATCATCGTGCCATTCCATGTTTGAACAGGATAATATACCGACGATGAGTAAATAGGTGGAAAGTACCCTGGTCTCATAGAATGAATAGGCGATGCTTGTTGTGTCGCCTTTCGGGCCGTTTTGTTTAACTTTTGTTTTAGCAGGTGGCCGGGGTTTCTTTGTTGACCGATCACGTGGAACAGTCTTTTTGCTAGCATTTTTGGAGAGCAACTGATCAAAAGTCGGATCGGCTTTGATCAGCCGATTATATGTGCTTTGACATTGCGTCTTTTTCCTTGCTTTGTGTAGAGGTTGACGCCTTTGGTCATAGGGGGACTGTCCGGCTGAGTTAGCCGGACCGTTTGTCTGAGCACCGGATCGTCCGCACGTAGGTGCCGGACCATCTACGATGCTCGGGCTAGGCTGATGTGTTTGGTTCATTAACTGTGCCTGCCCCCCAGTGTCTTTGGACTTTTTAGCCTTCTCGTCCGAAGCCTTTCGAGCAATCTCCTTTTGTGATATATCTGACATGCGGGGATCGCCAATGACGATACCCTTGCCTTTGCCTTTATCGGCCATTTCGGGCCGAACCAAGACCTTTTTGCATGTGAGTTCTATTGTATTAACAGGAAAGGGTTGTGTGTCTACTTGCATCTCCTGAAAAGCCAATCGACCCTCATTTATGGCCGATTGTATCTGTCGACGAAAAACATTACAATCATTGGTGGCATGAGAAAAGGAATTATGCCACTTACAATAAGCACGCCTCTTTAATTCATCAGGAGGAGGAATAGTATGAGTTAATTTAATGTTGCCGTTTTTCAGTAACTCGTCAAATATTTTATCGCATTTGGCAACATTAAACGTAAACTTAACTTCTTCTTGTCGATTCTTTCGAATCGACTGTAAAGCAGAACACGGTGAAGGTTTAGCCTGCCCTGGCCAAACCAGTTCGGCGGTGTATACATCCGTGGATTCATCGTTCGAGTTATCGTATCCCACTAGATGCATCTTATGGCTAGCCGATTTTGATGTTTTCTTACTTCGGCTTTCACATGTCATAGCCCGCTGGTGCAAGTGCACTAGCGAAAAGAATTGTGTGCCATTTAATTTTTCTTTTAAGTAGGGTCGCAACCCATTGAAAGCTAGCCCTGTTAGTTGTTTTTCTGCGACATGAATCTGAAAGCATCGGTTTCTAGTGTCCCGGAATCTCCGGATATAGTCATTAACCGATTCTTCAGGCCCCTGTCGGACTGAGACTAAGTCAGCCAATTCTAATTCATGTTCTCCTGAGAAGAAGTGTTCATGAAATTTTTGCTCTAATTCTTCCCAAGAATTAATAGAGTTTGGTGGCAAAGCTGCGTACCATGCAAATGCAGTATCAGTAAGGGACAACGAAAATAAGCGAACGCGGTAGGCTTCCCCATCAGCCAATTCTCCTAAGTGTGCTATGAATTGGCTAATATGTTCGTGTGTGCTTTTCCCACCTTCACCAGAAAACTTAGAGAAGTCTGGTATTCTAGTTCCCTGTGGATATGGCACGGTGTCGAATCGGTGGCTATAAGGCTTCCGATACGATTGCCCCGTACCTGACAAACTAACACCGAGTTTGTCCCTGAACATCCCGGCTACCTCGTCTCTGATCCTCTCCGCCATATCTGGCGACCATCCATTGAGTTTGTGGGTGGAAACCTCAGGTTGCCTGACATCACTCTGTCGGCTTTCCTCCCAGGGATGTTTGAGGTGTGCGCTAGGTGTCCTATTAATGTCACCAACTCCTGCCCTATACCTCTCAGGCTCTCTTGTGGCCGAACAGTATTCAGCCCTATGTCCATGAGGTGGTATGGCATGATTATAATGTGTTACTGGTGGGGCACCATAATGTTGCTGCGATGAATGCGGAAACTGTGCGTATTGCGCAGCTCTCGGCTCTGCGTATGCATATCCGGACGGTCCGGCGTAAGAGGCCGGATGGTCCGCGACCTGGCCAAATGGTCCGAAGGTATATCCGGATTGTCCAGCCGTATATGGTGCTACATGGGTAGTCTCGTACCCAGACCGTCTGTCGTAGAGAACTGGACGGTCCGCGATCGGGCCGAATGGTCCAGGGCTGTACCCGGACGGATCGGTCATATATGGTGCGACTTGAGTAGTCTGCGCGTGGACTCGTGTAGAGTCGGATGGTCCAGCGTAATACGTCGGCCGGTTCATGCCATATCCGGACGGTCCGGCGTAAGATGGCGGACGGTCCGCAACATATCCGGACGGTCCGGCGTGATGCACCGGACGGTCCGTGATGGGGCCGAAGTGTTCAGGGTTGTACCCTGATGGTCCGGCCATGTACGGTGCGACCTGAGTGTTTTGTGTGCGGGCCTGCATCTGGTCGGACGGCCCGACGAAATACGCCGGACGGTCCGCGGTATGACCGGACTGTCTGAGATGATGCTCGGACGGTCCGGTCGCGTCCAGTACCTGCCGTGTGCCTAGTGGTTGCGGCTGTGATGGTGACACGAAAGCGTGCATCGGCATACCATATGATGGCTGGGTCAAGGGTGACCTGTTTATAGCCGATGTGTTTGACGTGAGTGGCACTGTATTAGACTCTCGTGATGGAAAGTTAGGAGCAGCTGATTTCTCGTATGCACGTAAATGCATTTTAATAGATTCATCTACATATTGCTTTAACTGATCTCCTCGTTGATCCATGAAAGTCGTAAGAGATAGGTCTGGAGTACTTACAACGGGACCCTGAAGTGGCGGTAGGAGAGATTCCATATCGATCTCCCCTTGACGGACGATCTTCTGGTGGCGATCTACCGTGAAGTGTGACAAGTACTTGTCCGCCGCCTCCTTGCGCCTTTCGGAGAGTTTATGCAGCAGTTCCGCCTCTTCCTTGTCGCGTTGTTCGTTCCATTGCCGCATCACCTCCTTCTCATCGCGCATTACGAGGTCTTCAAAGGGCCTTTGGTCATCAGCCGGGAGCGCCTCCATGGCCGGCTTGATGATGTTGGTAGTGGAGATCTTGGTGGGATCCTTAGAACCGGCCATTTATGGGCCGATTTTTAGCAGATCTAGACACCTAGTCCCCAGCGGAGTCGCCAAAAAGTATGTTGACACCTTTTCGGAGGCGCCAATCACTCAAAAGAACCAGCGGCGGTGCTCTCTGGTCAGGCGCGGACGGTCCGCGGCTTGGGGCCGGACGGTCCGCGACCTGGCGCGACGTGGCGGTGCTCTCTGGTCAGACGCGGACGGTCCGCGGCACAGGGCCGGACGGTCCGCGACCTGGTGCAGGAGCTCGGGTTTCCTGCCTGACGGCCGGACGGTCTGCGCTCTAGGGCCGGACGGTCCGCGCGTGCGCAGGGGCGGTGGAAGATCGCCGGCGGCGCCTGGATCTCGCTCCCGGGAGGGACCCCGTCGGGGAGGAGAGATCCTAGGAGTTGTCTAGGCTCGGGCCGGCCGACCTAGACTCCTCTAATCGACGTAGAGTCGAGGAGAGGCAGAGAATTTGGGGATTGGAATACTAAACTAGGGCTAAACTAGAACTAGACTAGAACTACTCCTAATTGTGCTGAAAATAAATGCGAGATAGAAGTTGTATTGGTTCGATTGTTGGGGGTTCAATCGGCCGTAGCCCTTCATCTATATAAAGGGGGAGGTCTGGATCCGCTTCCAACTGATTTCCGAGTTAATCCCGCGGTTTTAGGTAACAAATCCCGCGAGAAACTAGGAACCCTAACTGACTCTGCGCACGCGCGGACCGTCCGCGCCACCACCGCGGACGGTCCGGACCGCGGACCGTCCGGCCTCCGGGCCGGACCGTCCGCACGGTCATTTTGGGTTCCAACACTTGGGAACCTTCTTCCTCACCAAGAAACAATGCTAGGGAACCTTGCTGATATTTCTTTCAGTTCTTGTAGTATTGTTCTGTGATAGGGCCTCTAGTTGAGTTGGTTAGGTGGTCTGAGTAGCACTCCTTTGGTCCTAAGTTCGAATCCTAATGTGAGCGAATTTCAGTCTAAGGTTAAAAAATGTCAATTGCTGGTTCCCCTGGTTGTGTGCACATGAGATGGACTGACCTATGGAGGGCGGATCATTGTGCAGGGTCTGAGAGGGCTCAAAGAACAAGTAAAGATCCGATCTATAGGGGGTAGACCCTCATGCTGCATGGGGGCCAGCTTTTGTGACCTTTTTTCATTCGGGGCTCCGATTGAGCTTTTTCTTAATATAACACTGTGGGGGCGGTCTTTTCCCTATTGTCTGAGTATTTTTTGTAGTACTGTTCTCATAATCCTGAAATGTCTTGCTTGTCGTTTCATATTCTTTGCTTTTTTGCAGATGGCAGGGAGGATCCTAGTAATGCGACCAACATCTTTAGCTGATTTCACCAAGAGGTATAACTGTTCTTGGTATCAGAGAACTTCATTACCAGCTCTTTGGACTGACATTTCGAGTATTTTACTTAGGCGCATGAAGGACCGAAATTCCTTAAGCTATCTTTGATATTCTAGTACCATGTTGTCGAGTGTTGAAAATTTGGCTTCTGCTCAATCTCGTGCCATTATAGTTTAGATCTTGCGCTTTAGGTTTATGCTCCTGACTAAACTTGAACAGCCCCCCACCCTCCTTTACTGTCTTTCTACAGTTTGGAAGAGATGAAGGTGACTCTTCTCCAGCATGATGTGAAGTTACTCATTGTTGATAGTGTGGCTGCTCTTATGTCGATGTAAGGGACATGGAATTCCTTAGTCTAGTTTTCAATAATTTGCTGTATATTTTAAAAGTTTGAAGTGTATCCTTGTAGACACTGTTCTACATATATATAAGTATAGTCTAGTTTTCAATAACTTGCTGTATATCTTTTGCTGGCAATCAGATCGTCAAAGTATCCAATTTCAGATTGAAGTATCAAGTATGTTTGACTATTGTTGTAGAACAAGAGCCATGGCCACTATGTATACGGACTAGTTTCTGTGCCATCTCACTGCCAGCTACACCTCCAATGCTTGCAATGTTTTGGAGCATTTCACATGTGTAGGCTTTAGGTTTTTGCATAAGCTGTTCCTGAATATTTTTGCTACCGGTTGATGGAAAACTAATCTTCTATTACATCCTGTCTGTAACAGGGAAAATGAGAAGGCTACAGCAGGTTTCAGGCAACACCCTTTAAGATGGAGCCTTTCTTTTCTTAAGTATGTCTCCCAGTCTCTTAATGCTCACATAACTTTCTCTAGATTGTTCTAGACTACTACTAGTTCACCAATATTCATCATGCTGAAATACCTTTTGTTTCCTCTCTTTCTCTTTTCATTGTGCAGGTCTATAGCAGAGTTCTCAAGAATTCCAGTTGTGGTTACAAACCAAGTACGCAGCCAGAGTAATGATGATGGTTACCATTTTTCCTTTGAAGGTTTGGACTTTCTTCCTTCAACAGCTATGTTTTAGATAATTTGTCGAACTAAACTAAAGTATCTGTCAAAATCATCTGCAGTGGACAGAAAGGATGGTAGTAACTGTGCTGAAAGGTTTGATTCTCATCTTGTTGCTGCACTAGGGATTCAGTGGGCTCATGCTGTAACTGTCCGTCTAGTCTTTGAATCCCATTCAGGTCTGTCCTGATATTTTTCTAAGCACAGAAGGGGCTTATGTTCTCTTTTCATCTGCCATTGATAGATCGCAAGTTCGCAACACTCATACCTGGATTTTGGAATTTATATCACATCAGTACCTATTATTAAGATCACATGAATGCTTGAATCCGTAGTTATGATTGCATCACGTAAGTTTTACCTTAGAATCTGACTTCTGTTTTGATGTTTGATTGAAGAGTGCAGAAGCACAACTAAATGCGAGCGGCATTTAAAAACACATCATTTTTTTTTATCTTGCAACATAGCACTACAGCAGTGAGTAACATTTAGTGCTGCGTCATTTATATGGCAACTTTGGCATCACTTGGTAAAAGTTGGCAGGATTTGAATAGCAATGAAAGAAAAATACAAGTATGTTTAGGACTAGTTTGGCAACTATATTTTTCCAAGAGATTCCCATTTTCCCAAGAGAAAATAAATTAATTTTCCTTGAGAAAATAGAAAATCCTTTGGAAAAATGGGGTTGCCAAACTAGCCCTTAGGAGGCTTGCTCTGACAGTAAAATATTTATATGAGTTGCTGAAGATGTGCATTTTGGCTTATTAAAATACTCCCTTTTTGCTATCATGTTTAGTCACTTGCTATTTTACTAACCAAATAGCAATTTCACATTTATCCATTGTATTATTTTAGACGAAATAATTATAAATTGCATTCTCCATAGGCCACAGGTTCATTAAGGTGGCAAAATCACCTATGTCTCCGGCAGTAGCATTTCCATTCGCTGTTGAGTCATCAGGTATTACATTACTAAGTGACGAGAGCATTGATGTGACAGGTCCTGATATCACCTCAATTCGTTGTCAAGGTACAATCTATGTTTGTTCTGCAGATAATGGGCTTAAGCTTGCAGTCTTCTCTTGTTATTGTCCTTCTAACACAGATATTTCTCCAGGGCAAAACGTTCTGGCCTGATAACAATTTGAATGCTACTGTGGCTACTGAAAAAGTGTAATCTTGCTCTTTGTTGTCTACTGGACAATATGAATATATTTAATTTAATTTCAGCGCACTCTCCTCAAAGTTCAGTAGTAGTCTTGCTTCTGTATCCTGTTTTCTCCACACAGTGTAGCATGTAAAAGATACAAGATCATTTTAGCATTGCTGTGACTATAAAAGGGATTAAGTGCTTTGCAATTGTTTGGTTGGGGCTATTAAGGGGATGAAAATGCTTTTGACTGGAAACTGTGGGGTCAGTTACCCACATGACATATTTTTTTTCATCAAATAAAAGGAAATATACAAGGGAGTTTCATGGGCTGGAGGCCCATTAAAAGCGAAAAAGGAAAGAAGCATCATGAGAAACTACTAAGCCAATTACTAAGAAAATGCACATCTGAAGGTTTTACTCTAAGTAAACCTTAACAAATTCAACTTTGAATTGTCTGTTTCATCTTGAAGGTAATATTGGATCATTACCACTCTGTACGTGCCTAGCCTTAAGTGATAACATAGAAGTTGCCAGTCCTAAATATCTTGACCCTAGGTGCTCTTTTGTGTTGCAATTGTGACCTGTAGGTATAGTTCCATCCTGTATAGTTCTGGTACCAATATGACTAACTACAGTGATGATCTGGTTCGTTATTTTGCACTATGAGATGTATAGGAATAAGGCTCACCTGCGTGGATCAAACTGGCTGATGTCAAATCTGGTTCGTTATTTTGCACTATGAGCTGTATAGGAATAAGGCTCACCTGCGTGGATCAAACTGGTTGATGTCAAATCTGGTTACTGTGGGTGGTATGGTCTCTTTTCTGTTCAGTGAGTGGGAAGGGTATTCCTGCTGCTGCCTGCTTTATTCTGCATTAAAAAAACACGAAAACTATGTCAATCTTCTCAGTTACGTTTCTGTGAAACTTCCCTCCTGCAACCGTCCAACTAATTTGAAAAGATCAACTCTCTTCCTTTGAAAGTCTGCAGGTCTCTATCCTTTGACAACAAGTTCTGCTCGAGTGTAACGATTGAGGGTTCCTTCGTTCCACGAGGGCCGAGGAAGATTGAAGTAGCTGGTATTCGGTTGTTGTTTCTTCCCAAGTACTATCCTTTTGTGATAACCTAGCAAATTGATGCAGTTTCGTGTAACTCATGTAGAAAGGTCAAAGCAACAGTTTTTGTTTTTACTGAATTGATTTATTTTGCCTTTGTATGAACTCTGGACTTGTAGGTTTCTTCGCCCTATCTGTACTTGATGATGACACCTGCCGTTTCCGTCTGCCTCAACAAAACACGCGGTACGGCAACTTTTACTCGCGAGCTCTGAGCTAAATTGATGCAGAAAATATGATTTAGTCTCAATCTGTCCCTGTTTTTGTTACGCGTGCGTGACTCTGGTGCCTGATCCGTTACTGTCATTAGATAATGCATAAAAATCTCATAGAATCCCTGATTGGAGCACGTCGTTTTTAAGTGGAGGTGATGTCTGACCTGCGGGTTTAATTGGGTATGGAGGACTAGTGGACTTTCTAAGAATGGACACCAATCCTCCACGTCCCGCGATAGCGTTTTTTATTTATGCCTCGATAACCTCGTATGGTGTGTGTGCGGCACATCATCACGTTTCCGCTTTCCGTTTTTTCCACCCACGGAAAAGAACAACAGGAATTGCTGTCACACCTTGTGGCACCAGCACCAGATGGTCCAAATCCTTTGAATACTTTGTAGGCTTTTTTTTTTGTTCACTCCAATCTAGGACTTGCTCTTCTCGTCTCTTCGGCCGCCCTTCCTGCGTAGAATGGGCAGACTTGTTCGGTAACAAGAGATAAAATTATACCATAGATTTAATCCGGATCTAGATTAGATGAATCTATATTTTACACCTAATCTCATAGTAGGATTAAATCCATCATTTAATCAATATATCTTTCAAAGCTAATTATTCTTTTGATCCATGAATTCTAATATAAACTTGTGCAATATCGTCATGGATTTGTATCGTAACTAATGAGTTATAGATAGAATCAATACCTTTCATAGTTAAAAAACCTTTCCAAACCCTTGGGTTATAATTAATGAGGAGGCTTTGTTCGGTTATTCCTATCTTATGTAGATTGAATGAGATTGAAAAAAAATTATAAAAGATTTTGATTTTGTTTGGGAATTAAACTCACTCAATCCCATCCAATCCACATGTATTGAGAACAAAACGAACAATCCCTGATGTGTTTAACCGAATAAAATAATTTGAATAGACATGGATTATTTTGGGATATGGATTGTGGTATGAATTTAATTCCAATTTATACTAATCCACAGCTGGATTGGTATAAAACGAACAAGGCCGGCGAAGCAACTGAGCTGACTTGTTAGTAGACAGGAAGAATCACTTTTACGGGCGAGGGCCTCTGCCGGTCGCTCGGCAATTATATAAATAATAAATAGAGATCTATATAAACAAAAAGCTCGATGCTTCCTTTCGTGCGATTTGGACGGAACGCTGTCGAGCAGCAGGGCCACTTTGGTGGGTTTATCCTGATTGGACGAAATGGGTATGTCTCTAGAAGCAACATATTTTGATGCGTGATGGGACTCATCATCGACAGGATTATATATGCACGGGGAACAAGCACACCACACCAGTTGCGACTGCAAAGACCTAATTAACATGCCCTGGAAAAAGATAGCAGGGAAATATCCACCACATATATATTTGTATACTACTATGTTAGTAGTGCGTTTCCTCATTTGGCGATAAGCTATTTTTTTTTAAAATTAAAAAAATAGCAACCTGGCCAGTGCTAGCTATTCTGAGTATATGCATGAGATAAAGCATTTAACTGGGGGGCATCTGAATAAACATGAAGACCGGATTCTGTCTAGAAGCAGCAACTTTAATTGTTCCTACTAAAGGGCAGCTGATTCTCATTGTAATCGCACCGCCGTCTCGTCTAGCTAGCTAGGCCTGGTTATTTGTGTAGATGGAATTAGGTATTTGACCTGGAAGTTGGCCGGGGTTGCTCGAGTCAGCTAACCAGAGTTGCCTATTAGCGAAACCTGATTGATTGGGAGGAAGGAAGCGAAACCTGAACCGATCTGATTGATTGGGAGGAAGGAAGCGACCGGGCTGGGGAGGAACCTAGAATAGAATAGGACGAGTAGGACAATCCCCATCTTCCGTAGTAGTAGCAGTAGGCAATTGAACTCCTTTGGGGCATGTTCGGTTACACCAATGTATAATGGGATTGGAGATTAAATCCCTTCATAGTCACTATGAGGACTAGTTTGGAAACTCGAATTCCCTTCGGGATTGAAGAGAATTGAGAGGAAATTAGTTTATTTCCACCTTAGGGCTTGTTCGTTTTGTTCTCAATCCATGTGGATTGGGTGGGATTGGATAGGTTTCAATCCCAAATAAGTCAAAATCTTTCATAATTTTTTTCAATCCCATCCAATCCACACGGGATAGGAATAACCGAACAAGCCCTTAATCCCCTCCAATCCTGAAGGGGATTTGAGGTTCCCAAAGTAGCCCTGAGGGGATTTAATCTCCTACAATCCCCTTCTGAATTGGTGTAACCGAACAAGCCTTTGGGGAGGGGCTTAGTAATAGCGCACACAGGCAACAGCTAGGATCATATCGTGTAGGGACACAAAGTAACCAGTACCGTATCGCTACGATCAATCCAATTCAGTTCTGGACGATAGGAGCCATCCGTTTTCCTTTTCATTCGTTTCTCCCTCTCCCTTTCTTGTTCGTCTTCTCTATTATCTCTTTAAAGAGGGAGTAATACATTTAAATTCAAAAGAGGTAGTACCTTTAAAAAGAGAGAGTAAGAGTATGTTTGGTAGGACTTTGGCTCTTCTAAAAATGGTTCTGACTCCCGCTCCTCTGCAGAGTAGCAGCTCCACAGGAGTCCGTGCTTTTTAACCAATCATTTGGTGAAACGGTGTGGGCATCTAACCTCTCAGTAGTGCCTTGGTAATTAATGACAAATATGTGTGGACTGACGATTTCTCTGATAAAGTAAGGTTGCAGGTTGGTCCATGTACTACATCTAGTCTTCAAATGCTTGATGTTTTGGACATGATCTTGTATATATACCAAGGAACTATTAAATAGGGTTGTTGTTTCTATAATATTTAAATACAAGATGTTAGCTTAGATGGTTGCAACGCTGGCTGTAAAAATTGACACACTGGGTTTTGAGGTCGAGCTCTAGTGTCTTTTTTTTTTTGCAATCGTTCGGGACTGACTGAGCTCTAGCCTCTAGGTGTGCCTTTTTTTCCCGACCTGATCGATTGCACTCGTAGGTTTTTTCACGGTTTCAATGCATGCACACATGTAGTGGGAAGAGGAAGCGTGTCTCTGAGTAATTTGTCATGCCTAGGAGAATATGTGGACCACCTAGATCAAGGGTATCCTAGTCTAAAATGAGCTTCAAATGCTTAAACATAAAGTATTTACATCCAAAGTGAAGACCCCTTAAAGTATTTGAAGACTGGATGGAATACATGTGAGCTACTTAGAGTCTTGAACGACTTTGGTGGAAACCGTTAATACAACTCAAGGCAAGGTAAACAATGGATTTTTCATTTTGCCGGTGATTAGAGAAGTAAATTGATCATGTCAGTAGGCTACATAGTCCTACTATATCAAGAGGCAAACAAACAACACAAATTATTAGGTGGATTATTATATAATCTTGATAACTAGATATCATTATATGATTATATAATCTATAAGGTGGATTATATAATCCTGGAAGTCAACAAATTGTGTCTTAAGCTAAGGCCGTGCGCTCCACATAGTGGACGACCTATGTTATTAGGTAGATGACGTCCAAAATGATTAACCAGAGAACCAATGAAAAATAACAAAATCATATATTATTATTATTATTTGCTCGAAAATAGATTCCGGGGAAGATAATTCCTTTTCTGGATTTGGAAAGAAATAAAGAATCCGGTCTTTATGTAGCACCTTTATTTATTTATTTCCTTACCATCTCTGTAGACAACAACAACAACAAAGAAAAAAAAAGATGAGATCGATATATTCCCGCGCGCGTTCCTTTCTTTATTTGCTGTTTTTTTTTTGTTTTTGGCAAAGTGTCACTGTACCTTTGGCCGCCGCGGAAGGCACTGTCGTTTGACTGGTAGTACATCTGATACTACTACCTACCGGTCCGGACAGAATAATGACACACGCACTGCCGCACTGGCACTGCACTTGCGCGGGAAAGGGAGGGGGCCCCACGTACCTGCGTGCTTTCCATTTCTAGCAGTAGCAGAGGCCGTCCGATCCATTGGCGGCGGCGTTGGTCCTGTGGGCTGCTGCCTGCACCGGCACGGCGGCTGGCCAGTGGCCACCATACCAGTATACCACCGCCCCCCGCCCCGTTTGTGTGGCGGGAAAAACGTACAGTGCGACCACGTCGCCGTCTTCTCTTCCTGTTGGACACGGCAGCTAGCAACGTCGTCATGGCCGCAAGCCCACAACACTGTTGCCTCTTGGCCGTCTCAAGGGACAAGTTGTACCAGCTCCCTCTGCTCACAGCCCGGATAATAATAATAATACTTGTTGTAAACAAAGCTAATAATAATGAAAAATAACATATATATATAAATCTTGATGCCCTTTCCAGAAAAAGCGGGTGCAGGTATGAAAGGTCGTCAAAAGGAGCAGGAGAGGAGATGCGTTGGCGATGAGCTTCCTGCCTTTCCCGTTACGTTTGCCTACCTTCACACATTCTACATTTGTAGAAGATGCCTGAAGTGGAAGGTGCCTAGGCCTAGGTTCACATGTTTATTTATGGAAAGCGACATTTGACAAAGAACGCGATTCCGATAGTGCTCTCACAGTCTCACTTCACATCCATTCCATCCATGTATCATGCGTTGTAGTGACGATAAGCCGCCGGTCGTAGGAGTGATCTACATAGATTTATATGACATACATGTATGTCACGGTACAATGCATGTCTTACATGCACGTACGTACGTAGTTATTGCCGGCGGTAACAGTGTACCTTTGCTACTTTGCCTTTGTCTTGCTTCGTCTCGACGCAAGCTTCACAATGATGCCATGTGTCATTTCTGCTTCTTCTTGCTCCAACATCAACCAGATTCAAAGGTCCAAACTCGAAACCCGTCCGCTTATGGTTTTGATGCCCAAACCACTTGACGTCTCTTTCTGTGCGAGAGTAGTGTATCTGATATATTCTTTCCATGATGTTTGCCAAATAAAACTAGTTCTATTGAGGGATGGGGAATTCCCACTAGTTTAGGTGGTCACGAGGTTAAATACAGCTTCCTCAATTTTGCTTTAACAATATAGTCTTTGTAGATTAATCTTTTTACACCAAGTAAGGACAAATTTAAGAGCGATGTTACGTGTACGTAGATGCACGAGAACTATATATAATATAAGTAATTTTTGGACACAAAGCGATATAATATAGTTTCTGGTTCGTCTAACAAACTAATTAAACGCTAGCTCGGTATCGTCGTTTTATGACGATGATGGATTTGACGTTAATTGCTTCTTTGTGGGGAATTCATCGGATCGTACATTCGTATTCATGTCGAGGAGAGGTGCAAGTGGACTCGATGGATGGGAATGCGCACAGTCTGCCTCAGCTAGAAAAATCGTAGATTTCTACTCCCTCCGCTCCAATTTATAAATTCGTTTAATTTTTACACCAAATTTGACATACTTGTCTTATTCAAAATTTTGCAAAAATATATAAAAAAATCAAAGCTACACTTAAAGTATATTATATGCTAAACGGTATCACACTAAAAAAATTAAATAAATTGTAAATTGAAATGGAGGGACTAGCATTTAAAATCTCGAACTGGAACTGATGATCTAAATAGACTAGCTTTCGTTATCCTGTCTACTACAGGAATGATATTTAGTACCGGTTGAAGGTAAACACTATTATTTTTTAGTATCTTACTTTTTAAATGTCTAGATATGTCTAATTCGGGTCTATTTGGTTGGTTCTCAGAAATTCAATTAAACAAACTAAATCTATTAATTATAGCAAACACTTTGATAATTGTGTAAAATTGAATATGTATGTCTATATACTGCAGAAAAAACGTCACAAAGAGTTTGGTATAGAAAAAAGAAAAAGAAACAAATATTCTTTATCGAGTGTTGAGGGAAGACACTCCGTAAAGTATACTTTTGCCGAGGGTTTGCCAGGTTACAATCGACAAAAAGAAGATAATTCCCAAGTGTCAACGTTTGACACTCGGGAAAGTTAATGGCCGTCAACAATAGACGACTCTTGACGACCCTTACCGAGTGTCGCATTTTGCCGAGAGCTTGGCATTCGGTAAAACCTCTTTTGCCGAGTGTTTTTCTATGTCGAGTGTTCCGCTCTCGACATACGAGATCTGTCCCGAGAGCCTTATTAGTTCGCCGAGAGTGACACTCGACAAAGAATGCTTTGCCGAGTGTCCGATTAATTACACTCGGTAAAGTGCTGAGCACTCATTCGAATTCCGATAGTGATGTCTGATTCGGCTCTCTAGGAAAGTGGCTTGGTCCAAACAGTGGCATAGGGTTTGCTTGACTCGCAGGCAAGTGCCATGCCAGGTCTCCAACCAAAAGAGCCTCTTGCCTGATTTTTTCATCAGGAAGCTCAATCTTAATTCTACACATTCATAAATCGACTGTAAAGATCGTTATCTAAAAAACATAGTATAACATATTTCAGGCATTTAAAATATAGCAAATGAGCAAATCCATGTTTGATATGTGGCCGGCCAGCCAGCCAGCCAGGATGATCCCAGTGATGGGCAAGCCCGATCTTGGCGATAAATGCGTGCGAGCGGCCAATCAAAGCGTACGAGAGAGGTGGAGTGTGTGTGGAGTACACGTACGCGTTTGACGTCCTGTTCATACCTTGCCAAACCTGCACGTACACAGGACTAGAGAGTAGAGAGAGACTCTATCTGTATAATCCATATAAATCATTTTATTAATAAAAAATACCGTTAAAATATTTATATCTTCAAAATAAATTTATTACTAGTCGGTTGCCCGTGCGTTACGACTTACAATAATATCCATGTAAACTATCCACAAAAAAATTTCAAATTTTTTTATTGATTGTATTCGCTCTCCGTATAATATTTTTTGATTTGGCTAACTGATGTTATTGTTTACTCCATCCAATATGTCTTGGTACTAGTTGGTTGTCCGTGCGTTGCGACGGTTTACAACAATAACCACGTAAATTGTCTACTAAAAATTTTCAAATTTTTTTATTGATTGTCTCTGCTCTCCGTATAATATTTTTTTTGATTTGGCTAACTGATGTTATTGTTTACTCCATGCAATATGTCTTGGTACAACACGGTAAATAAAGTGAGCGATTAGAAGAGAGTTCACAACGAGTGACTGAACGAACAGAGATAATAAAATGACATAATTCCACCATACAAATACCAAATAAGAGAAAGTTTGTGAGCTCAAGTTTCAAAAATAAGTCACATGAACTCAAACTTATAGATCAAAATATGGAGTGATTGATAAAGTCAGTCATCAATAAAACTGGATGGGCTCCATATAAATTATGCTACTTTGTAGTAATTTACTGACGTTTAAAACCAACAAATAACCTTTCATTTTTTACTGTTAGTGTGACAAATCATTGTTGCTCCATCCAATTCAGCAACCTCAAACACCATGGAGTCCATTGCGCCAATGTGGTCTCATAAACGACCAAATACTTGCAAGAGCCAAGAACGTCGTGCCGTGCTTGGACTGTAGCCTCGGTCCGTAGTGTTGGCCTGACACGATTTTTTTTTTCAAAAAAAACATATATACATATATACAATTTATATTATGTGGTCTCAGAAACGACCTAATGCTTGTAAGAGCCAAAAACGTTGTGTCGTACTTGGACTGTAGCCTCGGCCCGACACGATTATATATTTTTTTTATTAAAAAGCGTATATACATATATACAATTTATATTTAATATTAAAAATATTTGAACATAATGTTCTACTGGTTAGACAGCTTCATCCAGTGTCTCTCGCCCTTTTTCTATTAGGGCATGGGTTCAAACTCCACATTCTACATTGTTTTTTACATTTTATGTTGATTTAATCAAATGAGCCGATGGGCTAATAGGCTGGCCCGACACAGTCAGCAGGCCGGCATGACGTATCTGTGCCAAAGTTATGGCCCGCGGGCATCTGGTCCGTGTCGGGCCGCCGTTTGACCATCTATAGGCGTGCAACGAATTAATTTGAAATAGCTGTGAGGGGTTATTTGTAAAAAGATGACACGTGACGACCGTTGAAACTGGTGCTTTAAGTATAGTATAGATAACACGGTCAATGAAGTGATCAATTAAAAGAGAGTTCACAACGACTGACTGAAAGAACAAAGATTATAAAATGACATAATTCCACCATACAGAGACCAAATAAAAGAAAGTTTGTGAGCTCAAGTTTCTAAAATAAGTCACATGAACTCAAACTTATAAAAAAGATAGATCAAAATATGGAGTGATTGCTAAAGTCAGGCATCAATAAAAATTGAACGTGCTCCATATAAATTATGCTTCTTCATAGCAATTACTAACGTTTAAAACCAACAAATAACCTTTCATTTTATTGTTAGTGTGACAAATCATTGTTGCTCCATCCAACTCAGCAACCTCAAACACCATGGAGTCCATTGCGTCTCAGAAACGACATAATGCTTGCAAGAGCCAAGAACGTCGTACCGTGCTTGGGCTGTAGCCTCGGCCTGTAGTGCTGGCCCGTCCCTACACGATTATATTTTTTATTTTACAAAAAAACGTATATACATATATATAATTTATATTCAATATTAAAAATATCTGAGCATGATGTTCTACTGGTTAGACAGCTTCACCTAGTGTCTTCTGTCCTTCTTCCATCAGGGCATGGGTCCAAACCCGACCTTTTGCACTATTTTTAACATTTTACGCTGATTTAATTAAATGGCCCACGGGCTAACGGGCTGGTCCGACACAGTCAGCAGGTCGGCATGACGTGTCTGGCCAGAGTTGTGGCCCACGGGCGTCTGGTCCGTGCCGGGCCGCCGTTTGGCCATCTATAGACGTGTAGGTTAATTTGAAATAGCTGTGAGATGTTATTTGTAAAAAAGAATAACACGGGACGACCGTTGAAACTGGTATTTTAAGTATAATACAGATAAAAGATATTTAATATTTAAATCTAATGATATTTATTACGTACACAGGACACACAACGGAAGGACCTTTTTACTTGTTATATATACATGTATACTAGTTAGGTACCATGTGAAAGTCGCCTAGAGGGGGGTGAATAGGGCGAATCTGAAATTTACAAACTTAATCACAACTACAAGCCGGGTTAGCGTTAGAAATATAATTGAGTCCGAGAGAGAGGGTGCAAAACAAATCGCAAGCGAATAAAGAGTGTCATTCTTGCAAACCTACTCCCCGTTGAGGTGGTCACAAAGACCGGGTTTCTTTCAACCCTTTCCCTCTCTCAAACGGTCCCTCGGACCGAGTGAGCTTCTCTTTTCAAATCAATTGGGAACCAAACTTCCCGCAAGGACCACCACACAATTGGTGTCTCTTGCCTTGTTTACAATTGAGATGATCTTAAGAACGAATGAGAAAAAGAAGCAATCCAAGCACAAGAGCTCAAAAGAACACAACAAATCTCTCACACTTAACACTAAAGCTTTTGTGGAATTGGGAGAGGATTTGATCACTTGGGTGTGTCTTGTATTGAATGCCTAGCTCTTGTAAGTAGTTGGAAGGTGGAAAACTTGGATAACTTGAATGTGGGGTGGTTTGGGGTATTTATAGCCCCAACCACCAAACTAGCCGTTTGGTGGAGGCTGCTGTCGCATGGCGCACCAGACAGTGTCCGGTGCGCCAGCCACGTCACCTGGCCGTTGGGTTCCGACCGTTAGAGCTCTGACTTGTGGGCCCGCCTGGCTATCCGGTGGTGCACCGAACAAGTCCTGTAGACTGTCCGGTGTGCCACCCGCGCGTGCTCTGTCCTCTGCGCGCGCTGGCGCGCATTTAATGCGTTGCAGTCGACCGTTGCGCGCGATGTAGTCGTTGCTCCGCTGGCTCACCGGACAGTCCGGTGTGCACCAGACATTTCCGGTGAATTATAGCGGAGCGGAAATCCGAAGCTGGCGAGTTCAGAGTTGCTCTCCTCTGGGGCACCGGACACTGTCCGGTGGTGCACCGGACAGTCCGGTGAATTATAGCGGAGCGCCTCTGCGTTTTCCCGAAGGTGAAGAGTTCAGCTTGGAGTCCCCTGGTGCACCGGACACTGTCCGGTGGCACACCGGACAGTCCGGTGCGCCAGACCAGGGCAGCCTTCGGTTGTCCCTTGCTCTCTTGGTAGAACCCATTTCTTGGTCTTTTTATTGGCTAAGTGTGAACCTTTGGCACATGTATAACTTATAGACTAGAGCAAACTAGTTAGTCCAATTATTTGTGTTGGGCAATTCAACCACCAAAATCTTTTTAGGAAATATGTGTAAGCCTAATTCCCTTTCAATCTCCCCCTTTTTGGTGATTGATGCCAACACAAACCAAAGCAAGTATAGAAGTGTATAATTGAACTAGTTTGCATAATGTAAGTGCAGAGGTTACTTGGAATTGAGCCAATATGAATACTTATAAGATATGCATGGATTGTTTCTTTAATTTTAACATTTTGGACCACACTTGCACCACATGTTTTGTTTTTGCAAATTCTTTTGTGAATCCTTTTCAAAGTTCTTTTGCAAATAGTCAAAGGTAAATGAATAAGATTTTGCAAAGCATTTTCAAGATTTGAAATTTTCTCCCCCTGTTTCAAATGTTTTTCCTTTGACTAAACAAAACTCCCCCTTAATGAAATCCTCCTCTTAGTGTTCAAGAGGGTTTTAAGATGTCAATTTTGAAAATACTACTTTCTCCCCCTTTTGAATACAATAAGATACCAATTTGAAATTTACCAACTGAAAATCATTAGTTTTTAAAATTAGGTAATGGTGCGGTCCTTTTGCTTTGGGCTAATACTTTCTCCCCCTTTGGCATGAATCGCCAAAAACGGATACTTAGAATGAAATATAAGCCCTTTTAACACTACTTTCTCCCCTTTGGCAAATAAAACATGAGTGAAGATTATACCAAAGACGGAGAGTTGCTCGGAGCGACGCCGAAGGATGAGTTATGGAGTGGAAGCCTTTGTCTTCGCCGAAGACACCAATTCCCTTTCAATACACCTATGACTTGGTTTGAAATTCACTTGAGAGCACATTAGTCATAGCATATAAAAGAGACATGATCAAAGGTATACTTATGAGCTATATGTGCAAGACATCAAAAGAAATTCCTAGAATCAAGAATATTTAGCTCATGCCTAAGTTTGTTAAAAGTTTGTTCATCAAGTGGCTTGGTAAAGATATCGGCTAGTTGATCTTTAGTATTAATGTATGCAATCTCGATATCTCCCTTTTGTTGGTGATCCCTTAGAAAATGATACCGAATGGCTATGTGTTTAGTGCGACTATGTTCAATGGGATTATCCGCCATGCGGATTGCACTCTCATTATCACATAGAAGAGGAACTTTGGTTAATTTGTAACCGTAGTCCCTAAGGGTTTGCCTCATCCAAAGTAGTTGCACGCAATAATGACCTGTGACAATGTACTCGGCTTCGGCGGTAGAAAGAGCTACGGAATTTTGCTTCTTTGAAGCCCAAGACACGAAGGATCTTCCCAAGAACTGGCAAGTCCCCGATGTGCTCTTTCTGTTAATCTTACACCCCGCCCAATCGGCATCCGAATAACCAATTAAATAAAAAGTGGATCCCCTAGGATACCAAAGCCCAAACTTAGGAGTATAAACTAAATATCTTAAGATTCGTTTTACGGCCGTAAGGTGAGCTTCCTTAGGGTCGGCTTGGAATCTTGCACACATACATATGGAAAGCATAATATCCGGTCGAGATGCACATAGATAGAGTAAAGAGCCTATCATCGACCGGTATACCTTTTGATCGACGGATTTACCTTCCGTGTCGAGGTCGAGATGCCCATTGGTTCCCATGGGTGTCTTGATGGGTTTGGCATCCTTCATTCCAAACTTGCTTAGAATATCTTGAGTATACTTTGTTTGGCTAATGAAGGTGCCCTCTTGGAGTTGCTTTACTTGAAATCCCAAGAAATACTTCAACTCCCCCATCATAGACATCTCGAATTTTTGTGTCATGATCCTACTAAACTCTTCACATGTAGATTCGTTAGTAGACCCAAATATAATATCATCAACATAAATTTGGCATACAAACAAATCATTTTCAAGTGTTTTGGTAAAGAGTGTAGGATCGACCTTTCCGACTTTGAATCCATTAGTGATAAGGAAATCTCTAAGGCATTCATACCATGCTCTTGGGGCTTGCTTGAGCCCATAAAGCACCTTAGAGAGTTTATAAACATGGTTAGGGTACTCACTATCTTCAAAGCCGAGAGGTTGCTCAACACAGACCTCTTCCTTGATTGGTCCATTGAGGAAGGCACTCTTCACGTCCATTTGATAAAGCTTGAAGCCATGGTAAGTAGCATAGGCCAATAATATACGAATTGACTCAAGCCTAGCTACGGGTGCATAGGTTTCACTGAAATCCAAACCTTCGACTTGGGAGTATCCCTTGGCCACAAGTCGGGCTCTGTTCCTTGTCACCACACCATGTTCATCTTGCGTGTTGTGGAAAACCCACTTGGTTCCTACAACATTTTGATTAGGACGTAGAACCAAATGCCATACCTCATTCCTAGTGAAGTTGTTGAGCTCCTCTTGCATCGCCACCACCCAATCCGAATCTTGTAGTGCTTCCTCTACCCTGTGTGGCTCAATAGAGGAAACACAAGAGTAATGCTCACAAAAATGTGCAACACGAGATCTAGTGGTTACCCCCTTATGAATGTCGCCGAGGATGGTGTCGACGGGGTGATCTCGTTGGATTGCTTGGTGGACTCTTGGGTGTGGCGGCCTTGGTTCTTCATCCTCCTTGTCTTGATTATTTGCATCTTCCCCTTGATCATTGCCGTCATCTTGAGGTGGCTCATCTCTTTGATCTTCTCCTTCATCAACTTGAGCCTCATCCTCATTTTGAGTTGGTGGAGATGCTTGCGTGGAGGAGGATGGTTGATCTTGTGCATTTGGAGGCTCTTCGGATTCCTTAGGACACACATCCCCAATGGACATGTTCCTTAGCGCGATGCACGGAGCCTCTTCACCTATCTCATCAAGATCAACTTGCTCTACTTGAGAGCCGTTAGTCTCATCAAACACAACGTCACAAGAAACTTCAACTTGTCCAGAGGACTTGTTAAAGACTCTATATGCCCTTGTGTTTGAATCATATCCTAGTAAAAAGCCTTCTACAGTCTTAGGAGCAAATTTAGATTTTCTACCTCTTTTAACAAGAATAAAGCATTTGCTACCAAAGACTCTAAAATATGAAATATTGGGCTTTTTACCGGTTAGGAGTTCGTATGATGTCTTCTTGAGGATTCGGTGTAGATATAACCGGTTGATGGCGTAGCAGGCGGTGTTGACCGCCTCGGCCCAAAACCGATCCGAAGTCTTGTACTCATCAAGCATGGTTCTTGCCATGTCCAATAGAGTTCTATTCTTCCTCTCCACTACTCCATTTTGTTGTGGGGTGTAGGGAGAAGAGAACTCATGCTTGATTCCCTCATCCTCAAGGAAGCCTTCAATTTGAGAGTTCTTGAACTCCGTCCCGTTGTCGCTTCTAATTTTCTTGATCCTTAAGCCGAACTCATTTTGAGCCCGTCTCAAGAATCCCTTTAAGGTTTCTTGGGTATAAGATTTTTCCTGCAAAAAGAATACCCAAGTGAAGCGAGAATAATCATCCACAATAACTAGACAGTACTTACTCCCGCCGATGCTTATGTAAGCTATCGGGCCGAATAGGTCCATGTGTAGGAGTTCCAGTGGCCTCTCAGTCGTCATGATGTTCTTATGTGGATGATGAGCACCAACTTGCTTCCCTGCTTAGCATGCGCTACAAATCCTGTCTTTCTCAAAATGAACATTTGTTAATCCTAAAATGTGCTCTCCCTTTAGAAGCTTATGAAGATTCTTCATCCCAACATGGGCTAGTCGGCGGTGCTAGAGCCAACCCATGTTAGTCTTAGCAATTAAGCAAGTGTTGAGTTCAGCTCTATCAAAATCTACCAAGTATAGCTGACCTTCTAACACTCCCTTAAATGCTATTGAATCATCACTTCTTCTAAAGACAGTGACACCTACATCAGTAAATAGACAGTTGTAGCCCATTTGACATAATTGAGAAACGGAAAGCAAATTGTAATCTAAAGAATCTACAAGAAAAATATTGGAAATAGAATGGTTAGGAGATATAGCAATTTTACCCAATCCTTTGACCAAACCTTGATTTCCATCCCCGAATGTGATAGCTCGTTGGGGATCTTGGTTTTTCTCGTAGGAGGAGAACATCTTCTTCTCCCCTGTCATGTGGTTTGTGCACCCGCTATCGATGATCCAACTTGAGCCCCCGGATACATAAACCTACAAAACAATTTTAGTTCTTGACTTTAGGTACCCAAATGGTTTTGGGTCCTTTGGCATTAGACACAAGAACTTTGGGTACCCAAACACAAGTCTTTGACCCCTTGTGCTTGCCCCCAACATATTTGGCAACTACCTTGCCGGATTTGTTAGTTAAGACATAAGAGGCATCAAAAGTTTTAAATGAAATGCTATGGTCATTTGATGCACTAGGAGTTTTCTTCTTAGGCAACTTAGCATGGGTTGGTTGCCTAGAACTAGATGTCTCACCCTTATACATAAAAGCATGATTAGGGCCATGGTGAGACTTCTTAGAATGAATTCTCCTAATTTTGCTCTCGGGATAACCGGTAGGGTATAAAATATAACCCTCGTTATCCTGAGGCATGGGAGCCTTGCCCTTAACAAAGTTGGACAATTTCTTAGGAGGGGCATTAAGTTTGACATTGCCTCCCTGTTGGAAGCCAATGCCATCCTTGATGCCAGGGCATCTCCCATTATAAAGCATACTACGAGCAAACTTAAATTTTTCATTTTCTAGTTCATGCTCGGCAATTTTAGCATCTAATTTAGCTATATGATCATTTTGTTGTTTAATTAAAGCCATGTGATCATGAATAGCATTAACATCAACATCTCTACATCTAGTGCAAATAGTGACATGCTCAGTGGTAGATGTAGAGGGTTTGCAATAATTAAGTTCAACAATCTTAGCACGCAATATATCATTTTTATCTCTAAGATCGGAAATTGTAATATTGCAAACATCTAGTTCTTTAGCCTTAGCAAGCAATTTTTTATTCTCAATCCTAATGCTAGCAAGAGATTTGTTAAATTCTTCAATCTTAGCAAGCAAATCAACATTATCATCTCTAAGATTGGGAATTGAAACATTACAAGCATTTGATTCAACCTTAGCAATTAAACTAGCATTTTCATTTCTAAGGTTGTCAATAGTCTCATGGCAAGTGCTTAGCTCACTAGAAAGTTTTTCACATTTCTCTACTTCTAGAGCGTAAGCATTTTTAACCTTAACATGCTTCTTATTTTCTTTAATAAGGAAGTCCTCTTGGGTGTCCAAGAGATCATCCTTCTCATGGATAGCACTAATCAATTCATTTAGTTTCTCTTTTTGTTGCATGTTGAGGTTGGCAAAAAGAATGCGTAAGTTATCCTCCTCATTGCTAGCATTATCCTCATCACTAGAGGTTTCATATTTTGTGGAGGATCTTGATTTTACCTTCTTCCTTTTGCCGTCCTTTGCCATGAGGCACTTGTGGCCGACGTTGGGGAAGAGAAGTCCCTTGGTGACGGCGATGTTTGCGGCGTCCTCGTCGGAGGAGGAGTCGGTGGAGCTCTCGTCGGAGTCCCACTCCCGACAAACATGGGCATCACCGCCCTTCTTCTTGTAGTACTTCTTCTTCTCCTTTATTCTCCCCTTCTTGTCGTCGCCCCTGTCACTGTCACTAGAAATAGGACATTTTGCTATAAAGTGACCGGGCTTACCACACTTGTAGCAAACCTTCTTGGAGCGGGGCTTGTAATCTTTCCCCCTCCTTTGCTTGAGGATTTGGCGGAAGCTCTTGATGACGAGCGCCATTTCCTCGTTGTCGAGCTTGGAGGCGTCGATGGGTGTTCTACTCGATGTAGACTCCTCCTTCTTCTCCTCCGTCGCCTTAAATGCGACCGGTTGTGCTTCGAACATGGAGGCACCGTCAAGCTCGTTGATCTTCTTTGAGCCTTTGATCATAAACTCAAAACTCACAAAATTTCTGATTACTTCCTCGGGAGTCATTAGTGTATATCTAGGATTGCCACAAATTAATTGGACTTGAGTAGGGTTAAGGAAAACAAGTTATCTTAGAATAACCTTAACCATTTCGTGGTCATCCCATTTTTTGCTCCCGAGGTTGCGCAATTGATTCACCAAGGTTTTGAGCCGGTTGTACATGTCTTGTGGCTCCTCCCCTTGGCGAAGACGGAAGCGACCGAGCTCCCCCTCGATCGTCTCCCGCTTGGTGATCTTGGTTAGCTCATCTCCCTCGTGCGCGGTCTTGAGCACATCCCAAACTTCCTTGGCGCTCTTTAATCCTTGCACCTTGTTGTACTCCTCTCTACTTAGAGAGGCGAGGAGTATAGTTGTGGCTTGGGAGTTGAAGTGCTCGATTTGGGCCACTTCATCCTCATCATAATCTTCATCCCCTACGGATGGTACCTGTACACCAAACTCAACTACATCCCATATACTTTTGTGGAGTGAGGTTAGATGAAATCGCATTAAATCACTCCACCTAGCATAATCTTCACCATCAAACGTTGGTGGTTTGCCTAATGGGACGGAAAGTAAAGGCGTATGTTTAGGAATGCGAGGATAGCATAAGGGGATCTTACTAAACTTCTTGCGCTCATGGCGCTTAGAAGTAATGGATGGCGTGTCGGAGCCGGAGGTGGAAGGTGATGAAGTGTCGGTCTCGTAGTAGACCACCTTCCTCATCTTCTTTTCTTGTCCCCACTCCGATGCGACTTGTGGGAAGAGGGTTTCTTCTCCTTCCCCTTCCCTTTGTTGCGGGACTCTTCCGATGAAGCCTTCCTGTGGCTTGTAGTGGGCTTGTCACCGGTCTCCATCTCCTTCTTGGCGTGAATTCCCGACATCACTTTGAGCGGTTAGGCTCTAATGAAGCACCGGGCTTTGATACCAATTGAAAGTCGCCTAGAGGGGGGGGGGGTGAATAGGGCGAATCTGAAATTTACAAACTTAATCACAACTACAAGCCGGGTTAGCGTTAGAAATATAATCGAGTCCGAGAGAGAGGGTGCAAAACAAATCGCAAGCGAATAAAGAGTGTGACACGCGGATTTGTTTTACCGAGGTTCGGTTCTTGCAAACCTACTCCCCGTTGAGGTGGTCACAAAGACCAGGTCTCTTTCAACCCTTTCCCTCTCTCAAACGGTCCCTCGGACCGAGTGAGCTTCTCTTCTCAAATCAATTGGGAACCAAACTTCCCGCAAGGACCACCACACAATTGGTGTCTCTTGCCTTGTTTACAATTGAGATGATCTTAAGAACGAATGAGAAAAAGAAGCAATCCAAGCGCAAGAGCTCAAAAGAACACAACAAATCTCTCTCACTTAACACTAAAGCTTTTGTGGAATTGGGAGAGGATTTGATCACTTGGGTGTGTCTTGTATTGAATGCCTAGCTCTTGTAAGTAGTTGGAAGGTGGAAAACTTGGATAACTTGAATGTGGGGTGGTTGGGGGTATTTATAGCCCCAACCACCAAACTAGCCGTTTGGTGGAGGCTGCTGTCGCATGGCGCACCGGACAGTCCGGTGCGCCACCGGACAGTGTCCGGTGCGCCAGCCACGTCACCTGGCCGTTGGGTTCCGACCGTTGGAGCTCTGACTTGTGGGCCCGCCTGGCTGTCCGGTGGTGCACCGGACAAGTCCTGTAGACTGTCCGGTGTGCCACCCGCGCGTGCTCTGTCCTCTGCGCGCGCTGGTGTGCATTTATTGAGTTGCAGTCGATCGTTGCGCGCGATGTAGTCGTTGCTCCGCTGGCTCACCGGACAGTCCGGTGTGCACCGGACATGTCCGGTGAATTATAGCGGAGCGGAAATCCGAAGCTGGCGAGTTCAGAGTTGCTCTCCTCTGGGGCACCGGACACTGTCCGGTGGTGCACCGGACAGTCCGGTGAATTATAGCGGAGCACCTCTGCGTTTTCCCGAAGGTGAAGAGTTCAACTTGGAGTCCCCTGGTGCACCGGACACTGTCCGGTGGCACACCGGACAGTCTGGTGCGCCAGACCAGGGCAACCTTCGGTTGTCCCTTGCTCTCTTGGTAGAACCCATTTCTTGGTCTTTTTATTGGCTAAGTGTGAACCTTTGGCACATGTATAACTTATAGACTAGAGCAAACTAGTTAGTACAATTATTTGTGTTGGACAACTCAACCACCAAAATCTTTTTAGGAAATAGGTGTAAGCCTAATTTCCTTTCACCATGCATGCTACGAAATCTATATACATCAAATAATTAGATCTTGTTTGAATAGAATAGATGCACGTTCTAAATATTTCTTTACAATAATGTGTATTCTATTTGATCTGGTCACTTTGCCCGTTCTACCATGTAATCACAATTCCAATCCATGTACCGTCTCTTGTGAGGTGGCTCGGATGACACCAAACCCTAGCCACTGCCACCGCCTCATCTTCCTCCCCTATCCATCTTGCCGCCACCAGAGTAAGCAGTCAAAAGTCCAAAGCGACTCGGGAGGATGGTGGCTGCGGGGTCATTCTAGCGTGGCTCAGCACGGTTTAGCCCGGTCCAACCATCGTGTCTTGTCTATGCCTTAGGTGCGGCACGACGGACGAGACCCGACATGGTCTGGTTTATGTTTTTCCTATTTTCATAAAAATATCTATATTATATTTAGGGTTTAGAGTTTTATCGGGAACATGATTCTTAGGCCCCCGGGGCCCACTATTCAATAGAAGGACATGGAAAAGATCCACCTAGGAAGGGCCCGCCCCGGATAAACCCAAGGCAACGAGCCAAGGGTCCGGCGAGGTGGAGATGCAAGGGGGTTGGGCGGATCAAGGGGAAGTCACTCGAGTGATTCCACATAACTGACCACTGTGCCGCAATAACTGCGCCTGGCACGAGCCATAGATGGTGAGTCGGTCCACTTTCCACTCGTGACCTACGAGCTCCTAGGCATGCGACCCACTCATGACCTATCGAGCCTAAGCCCGAGCGCGATGGCCGCCTGCAGGCCTCACAACACGATGAGACCTAGATTATGGTGGCTAGGGGCCGGTGCAGCACAGGGACGACAACAATGGTGCTGGGGCTCCTCGCCGCCCGCCATAAAGACCAAGAGTAGGGACACCTGGGTGACCAAGGCATATACGGACAAACGCAGCCTTACCAGGGTAGAACACTCAGTTTTACCACAACAAAGTAGCTCCTTCTCAAGGAAGCTATCGGTGGCCGACCCCCATCTTGGCATATATAAGGAGGAGGTTAGCCAGGAGACTAGGGAGATACACACTCGACAAGACCAGAAGACCTAGATGTGCGGGATGGATATCCCCTGGGTCACCACGGAACAAATACACTTCATAGGCCGCCTATAAGGCCCATGATGGATCGCCTATCGGACGAGTCGTGGGCCTCAATAGGGTGCAACCAAAGCTGAGTAAAGCGATCACGCACTGAGTCAAGGAGATCCTCGAAAGCGTGTCGCATGCATGCAGTCGTGTGCAAGCCTCATGGGGAACCAACCGAACACTAAAAGACGCAACACGTCTGTAGCGGAAGTGGTTAGACGCGTAAAGGACTCTGATTCGTATAGGGTTATCTCTTCAGCCAACCATATAAGGCTGAGGTATGCCTACAAAAACATTAGACTCCATTAGACCTGTCAGAGAAAAAAACTCCAGCCGAGTGGCAAAGTTCACCCACCTAATTCTTTGTACAGGGAGAGCCTAGGGTTTGCCTGCCGAGTCAGATGAACACAAAAACACACAAGGGTTTTGAGAGTGGTTCATGCCACCGTAACGTAATACCCTACTCTACTTGTTAGTTGTATTGCTGAGAGCTTGTGAAAACCTAGAATTGATTGGAGATGCCTCCTAGGGATTCAAAAGCGTGTGTGTGCCTATGAGCTCTCTCGGCCTGAACAACTTCAGGTGGGTCTGCTTTCTACCGAGTGCACCCTCCCTTTTATACCCTAAGGGAGGGGCCTACATGGCACCTGAGCCCCGACAAGTGGGCCTAGGGATACAAAGAGATCACATAGTAATGCTAGTGGATCGAGGTCTTCTCCGTTCGGTCTTGTTCTCAGCCCTGTCTTATTAGCGTGTGGGCAAAGGTTGCGGCATAGGTGAGTGTAAGGAGAATACAAACGGTGTCGTCCCGTCCCCTCATGCTTCGTTAACACATGAGGAATGGTGGATGTCATAATGAAGGGCACCCAAGCCACGCATTGAGTCAATGGTCTTGCTTTGGTAACGCGTGGGTTTCTGACGTGGCACGGGACTTACAGGCCGCGCATGCCACCTATCACCGGCGTGAACAGTCGCGCGCTCTGCCCACCTTAAACGCTGGTGGCCCCCTGCTAGAGGCTTACGTCAGGTTCTTCCCGATGGCTTACGTCACACGCCATAAGACATGTGGCAGCACCGAGTCTTCACACGAGCAGCGAGCGACCAGCCATGAGTGTATATAGACATGGTACATCCGGACCACCACTTGTTGGGGTGTGGAAGGTGCTCGAAGGAGTCTAGGACCCTCTGGCATGGGTCTGGACGACATGTGGCAGTATTGGACTTCTTTAAGTCGAAGCTACATGTAGCATACGTGAGAACACCTAGAGGGGGGGGGGTTGAATAGGTGATCCTGCAAAACTCAACTCTAAACAACACAAACTTGATTTGTAATGAGTGTTAGTGGAAGCTAAAACCAAGTTGGTATCAGTAGAGATGGAGAGAGGAGAACTCTTCACTTGATTGCTCCTTTGAAATCTGTATCAAACTTAGGAGCAATAATACAAGTGAATATGAGAACTCAAGAAGACGATAATCGCAATAAAGTAAATGGACAAGAGACACAACGATTTTATCCCATGGTTGGACACATGCCTACTCCACATTGTGGTGAACTCCTTCGGTCAAGGGTTTCACTCAACCCATCTCAAGTGATCCAAAGATCAAACTTGAGTACCACAGTTGTCTTCCTTATATCAATTCCCGGTTGTGAGAAATCTCCACAAATTAGAGTCTCTCATGCCTTACATAATTGATCTTAATCAAACCCGGACTAAGGATGGGAGTAGAAACACACACAAGAACACAATCGCAACAACACACACACAAGCCAAGAAGAGAGCGTAAGAATCAAAACAGGAGAGTCACAACTCAAGAATGTGCTCAAATATCTATCTAATGAATTAATTGCGTGGTCGCTTGTGTACCGCTCCATGGCACCTAGGGGGTCTTTTATATCCCAAAGGGACCTAGGAGCCATTGGAACTCAATTTGTAAGGCTCTGGTTGCCTTCTGTCTGCGGGTGCATTGGACTGTCCGGTGCACACCAGACATTGAACAGTGGCCAATTTATTTCCATGTTTGGTTGAGCTGACCGTTCAGTCAGCCGTTGGAGCCTTCGATCGCCTGGCACACCAAACTATCCGGTTCACACTGGACTGTCCGGTGGTGCCTCTCGAACCTTGTCCCAGTAGACCTGGCCGCCACATGGTCGCCATTGATCGCGCGCCGAGCGTTGGCTCTGGCGAGGGCGCTGCCCGCCTGGCGCACCGGACATGTCCGATGCACATCAAACTGTTCGGTGAATTATAGCTGGAGAGCCTCGTTGTCTTCCTGAGAGCGACCTGTTTGGCCAGCCGAAGACTGGCCAATCAGCCTAGGCACCGGACTATCCGGTGCACACCAGATTGTGCGGTGAGGCCCAGACTAGTCCAACTTTGGCCTAAATTAGCTAGTATCCTCCACTTCAATTTGTATTAACTTGGAGAGTTTCCTAGCACTTGGATGAACATTGCTAGCCACTCAAATAATTTTCTAAGTGCTAGAACCTTACCTTTGTTCTTTCCATATGATTTGCAATATACACACATATGCACTTTTCCCATGTTTCTCAATTTTTCTTTCAAATACATGTTCTTATGCTCATACATCAAGATGTTAGCTAAGATAATGTGTTGAGCATTTAATCACCAAAACAAAATAGAAATGACCCAAGGGAACATTCCCCTTTCAATACACAGGTGTCCGGGACCATTCATAAGGGGTTCGACCCCCCGAGACAACTGAAGGGATTCCCCATACTTGGAGCACATGGTGGTTTCGGAGCCCCACTACTAGAGATCCGGAGAGCCCTATGTGGAAGCTTCCCCCACTTGAGGCCAATGGCCTCGTCGGATGGGCCCTTATTTTGAGCTAGGGTGGACCATCTGCGGCTATGCGACCAGAGATAGTCACAAGGGTCCCGCCTTGGTCTAGCAGGAGTGGATATGGCCTCGTCGAACAGGCCCTTATTTTGAGCTAGGGTGGACCATCTGCGGCTATGCGACCAGGGATAGTCACAAGGGTCCCGCCTTGGTCCAGCAGGAGTGGATACCCTAGTCTTGAGGTACCAAAAGTGGCCCCTGGTCGACCTCGAGTGAGGCGGTGAACTCGCAGGTGGGGCCATCTCTGTATCTTGTTCCCAATGTGCTGACGACAAGTATACATTTGGAACGTACTAGGAAGCTCAGTAGGTTAGGTTTTTCGGGGACATAACCCTACTAGACTAGACCTTGTTGCCTGCTGTAAAAGCAAGTTTTTACTCTATCGAGACCTAGTTATCCAGACCGTCCGAGAGTAATATTGTACCAGTGGCCGACCTAGGTATCGCTCAAGGTCGGTGTATAGTAGCGCACCCTCGGGACCCTTACTTCTATTGTATTCAATCCATTGAATCATGCATCTAAGTCTTGCTGAACCTGGTTTTGGGACCAACCCTTGAATCAGTACGAGGTTTAGACTTCTAGGTACCTGGTCCAGAAAGTACTAAGTCAACTATTCTAGGGTGGTGGAAGATGTAGGCTTACAGTACGGGACCAGGCTAAGCGACTACCTAGCTCCGGATCACCGTGAGAATGTGGACTTTTCCCTTAGACCAAGTTCCTAGTAATCCAATCCTCATTCATCATCAAGTCGGGTCCAAAACCAAGTGTAAGGCTAGATACATGACTCCGATTTCATGATTGCAACATAAGTCTAGACCCACAGAAGGGGTAGAGAAAAACAACATGTGTAAACAGGAAAATAAATCGATGGATGAGAACAAATTCTTCCTTTATAACTAAATACATGAGATCTATGCAGTAAATACCTATGACCGGGTTTATAAGGGCTGTCCTCGCTGGGCTGGTCCACAAATGGTGCCAGTGGCTGAGCTTATGAGGGTCGGACCTTACCGGTTGGGCCACAAATAAATGCTATCCTATTCTAGAGGAACAATACATTCTCTTGCTTTTTTAGACAGAAGGTATGTGGCTGACCTATGATGCATGAGCATTGCCTACCCTTAGACACGAAGGGGAAGCAACACCCCGTTCTGCCTTGGCGATAGCCGAAGGCAGTGCTTCACGGTGGTCATTAGCAGCTACCTCAAACGCCATGGAGCTCGGTCTACTGGTCGACGCAGCAAACCCCGAATGCCATGGGCTTGTGCCATGGCCGACGAGTTGGGTCTAGCGCCACGAAGTCCTGTAGTAGGTTAGGTCTTGGCGCTGCCCTTGAAAACCTACCAGTGCGCGTGGTGGGCCCCTAGGCCCCCACGGCTGGTAGGCCACCTCTAGCGTGATGAAGTTGACAAAGGCGTAGCCCTTGTTGGCCTGCCTTCTAAAATCCATTAGCAGGTAGAGGAAATCACGGCACGACTTCCCATATCTAGCCAAATCATTCCAGGAAGGAGGGCCAACATAGGTTAGGTGACGAGAAGCACCCGCCTCACCCTCACCAGTTGTTGGAGGGAGAATTTCATCACCTAGTTAGGAGAATCATCTTCGAGATTCTCCCTACCAGAACGAGGAAATGACTCTTAGAGGAAAACATCTCTAGTTCCCCCATCCCTCTAGTGTAAAATGTCAGAGACGTAGACTCCGGACGGGTGATGGAATGCACCAGCTGGAGGCGCTAGAGGGAGGGCGCAAGAACTTACTAGGGCTTAGCATCCCCTACCTGCATGCCAATGTCGGAGAGAAAACTCCAGACGGGTAAGGAACTTCAGCTTAGAGTTTGGCGCTAGAAGGAGGGTGACTCTCATGAACACTAAGAGATTTTTCAAGATGGATGGAGATGCCTTGGGTGTTTCCCCTATGGAGTCGTCTTGACCTAAATAACTTTAGGTCGGTCTACTTTCTACCGAGTGCACCCTCCCTTTTATAGCATAAGGGAGGGGCCTACATGGTGCCTGAGCCCCAATAAGTGGGCCAAGGGATACAAAGAGATCATACAGTAATGCTAGTGGACCGAGGTTTTCTCAATCCGGTCTTGTTCTCTGCCCGGTGTTGTTGGCATGTGGGCAAAGGCTACGACATAGCTTGCGACATAGGTGAGCATAGGGATAATACAGACAGTGTCATCCTGTCCCCTCCGCCTATTCTCTCATACCTCGGTAACGCGCGAGAAAAGTGGTGGTCATAATGAAGGGAACCCAAGCGACGCATTGAGTCAATGGTCTTGCCTTGGTAACACGTGGGTTTTTGATGCAGTGCAGGCCTTGCAGGGTGTGTGCCACCTATCACTGACGTGAGAACAGTCGTGCACTCTGCCTCACCTTAAATGTAGGTGCCCCTGCTTGGGGCTTACGTTAGGTTCTTCCCGACGACTTACGTCACACGCCATAAGACACGTGGCAGCACCAGGTCTTCACCCCAAAAGGGATCGACTAACCATGAGTGTATGTATACATGGCCTAGCTAGACTCGTAGGCCACTTTCATCTTCTACCTTCAACAATCGCCTGTAATCCTCAAACTCAGAGCAATACAACCACTAGAAACTCACAGGACGTAGGGCATCACGTATCTAGGTGGCTCGAACCTGTATAATCCGTTGTGTTCCATCGTGCCTGGCACGAACAATCGATTTGAGATCAATGCCACCATCCTCCCAAAAGCACGATATGGTATACCCCGTAGTGCACCGTCGGCACCAAACACCGACAGCTGGCGCACCATGTAGGGGGTGCAACGTCGATTCACAAAAGCTCCATGCCCATCATCTTCAACATCGTCTTCAACTTCGCCACCGGGGGTCTTCCACTCCGAACCCCTCTCATACGTTGCTGGCCGAAGGGGCGCTCTGCACCGCATAGCGAACATGGAGAAGGACGCCCTGAGTTCCCACATCGCCACCTCCATGGAAAGGAGGCCACCACCGAGATCTGGGACAGCCCGAGTGACGGTCGTGACTAAACCTCGGTTGGTAACGGCTCAACCCACGCTGACACATCATCGGACCCTAGGGCCCCGATTTAACCGTGTTGTATCAGGCCACCTCCGACTGGCAAAAAAGTGCTTGGTCCCCTTAGCACCTCCGCCTCGGCCAGTCGTTACGTCGTTCCTAACGTGGTGGACCCTATTATCCACCTCACCGACCCAGTCCTAGACCAGGATCATTCCGAAGAAGGACGGGAATAAGGTCTATGCCGCCCCTAGGGCCGGCAAAGCGCGACCGGCTGAAGTCACAATCACCCGACCCTAGGACGTAGGGTATTACACTACTAATATATTAAAGCTAAGTTAAAGTAGTAATCATGAAAAAATATCGTGTTAAAGTCTCGTACCAATGAAAACTAAAAATATACTTAAAATTATATTCATATTAGTGACAACAATTTTATGTAAGTTGTTACTTTAGATAATTTTATATTTTTAAATGAAACATTAGGATATATCTATATTGATAGCATAACTATAAATACTAATTTTCATATATCACCCTAGTAGTACTGTATAGTTTTTCTGATAAGAAAATCTAACATCTTATACAATGAAGTATATGAAACGTTCAATTCAGACTACATGTGATCTGTTGGACCGAGCTATATAGACTAGCATGAGCCACGACCATTTTAGCATAACACGACACGATAATAACCGTACCCGGCCGACACGAGGCCCATCATAGATCATGCTTGGTCTTGCGCACAAGGCCAAGCACGACCTGTTTGTAGTTGTGGCCCATTTAACTGGATGCTAGGCACAGATAGCCCATCAAATCGAGCACGAGTGCGACCCATTTGTGGGCCGACCTTACTTGCTCCAATCTAGCCCACCATATATACAACGATGACCGACTCAACCCCATGTTTGAGCTATTTTCGAGCTCTATTACAAAAAAACGGAGAAAATGAACAAAAAAATAGAGAAAATAAATTATTTTTGAACTCTATTACAAGACAGGCTCTAAGCTGTCCGACACGGTAGCACGACACGACCCATTTAAAACTAGTCACTACTGGAGAAGGCTTCTTTGCCGAGTGTCATAAACACTCGGCAAGGCCTCAAAAGCACTCACCATGAATATCATTTCTACATTCATTGTTTTCCATTTTTTGTGTGACATACAGTTCAAATTTGAATTATCCAAAAAAATCAACTAAATGAACTAAATGAAATAGGTATAGCAAACACTTTGATAAATGTGCTAAAATTGAACATGTAGGTGCAGAGATTAGAGACTGTAGGAACACAATAGAACATTTTTTGGGAAAAAAAGAAAATAAACAAATAAACTTTCCCGAGTGTCAAGAAAGACACTCGGCAAAAAGAATTCTTTGCAGAGTGTTTGTTTGTGCCACTCGGCAAAGACTTGCTTTGTCGAGTGTCTTCTGTAGGCACTCGACATAGTTAACGTTGGTTGCTGGTGGCGCTTTGCTAAGTGTCCTATATATGTTACCGAGAGTTGAACACTCGAGTGTATTTTTATACCGAGTGTTCAGTTCTTGACAAACATGTTATTTGCCGAGCGTGTTTGTTTGCTGAGCGTTGTACTGGATAAAGTATATCTTTGTCGAGTGCCCGTGTATAGACACTCGACAAACCACCTAACACTCAAACCCACTATTTATCTGTACCGTGTCTTGGGTCCGATTAAATTGGCATGACTCATGACGGGCTCGAGCTGTGGCGGACCTATTGGACGAAAGTGTTATGCCCTCACGAATGTGCTCCAAAGTCAGAAGTTCTGACACGGTCGGAAGAAGTTTCGACATAGCAAAACTCGTAAGAAACCTAAAAAACCTTCCTTGGTGCCTCGGATTGTCCCCACTCCCCAGTGGGATAAAACTGTCGGAGTTTTGACAGAAGTCGAGACTTCTCACGAAGCTAAAAAAAGACCTTATTCAGACCGGCCTGCCTGATCTGTTGGATTTTGATCCTAATTATTTTATACTCATGAAAATATGTCTTGAAAAGGTTTTCTATCGCGATAAGGCCACCTCTTTTATATACACTAAAGAATCATGACCCATTGAGTAATCTAACGTCCAATCCCTCTTTTTTCCGCCTATTAGGCTCAGCCCAGCAGAGTTCCTATTTATATTCGCATGCGTAAAATAGGTAAAGAATGACTATTGTTTACCGTTCCAACAGAGCTCGTATCCGTCTTTGCAAATTGCGAGGGGCTCGCGACGCTTCCCTTTTCTGCGCTTCCTCTCTCCTCCCCCTAAATCCAAATCAAACAACCGTCCCATAGTACGGCTGTCTCGCCCCTCTTCTGTCTGCAGCGGCAGTTGCAAATGCGAGTTGGAATAGGAAGAGCCGCTGGAGATAAATAAATAGATAGGAAGAGTTTTGTGTGGATCCTACGTACTGTACCTCCTAAAATAGAAATGCGAATCCGATATCTACTAGAGTCGGTCTTACTTCTTTTGCCCTAGTCTCTTCCACCCGCTATATGTTGTTATTCATCGCTTGATCTAACGATCAAGAATATAATCCAATGATGTCCTTGTCCCGATGGGACAACCTCATATATATATAAACTCTTTGTGAAACGAAACTATATTTTTGGATAAATCTAGAGCTCACTAATTAATCATTGGGCCCACACGTCGTTTATGTATATCTTGCATTGCTATATGCTATAGCTACGCACGCACGCAAGCTAAGGTCGAAAGCAATTGATCCGCGACCATCATGCACTCATGGTTGGTTAATAATAACTTGACCATCATGCATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATTAAGGGTCACGTGCACTGTTCACGTGTCCTATTCTCGTCCGCTTCATTCATTGGCAAATTAGACTACCTACCTAGCTGCCAGCTAGTTTATTCATCGTAATTAATATATCTAATACTTACGTAGTATGTATATAATACTCCTATATATGTAGTATAATATATGAGAGGCCGATCGATCCTCTGGCAAAATTAAAAGAAGAAGATAACTGATGGATCGAGCTGGTGAGTGATCGATCTCGACGTCGGACCTAGCACTACTAGCATGGCCGCCCAGGCACGACTGATCGATCTCCGATCCGCCGGCGGCGAGAAACCTAGGCAGGCAGCACCAAATAGACGTTTGGATTGGATGGAAGAAAAGGAGAAGGCGGCATCATATTGCATGTGAACCGGTGCGCCCAGCGAGGCGCGTGCGTCTCGTGCTCCATTACAAAACAAATCTAGGACCCAGCGCAACTGTACCGCCCACATCGGGGAACAATTAAAAAAACTGCATAAAATGTGCAGATGATGGCGAGCCACGTCGTTGGTCGGTCGTCCTTGTGGATCGTTTGGCTTTAATCTATACTACCGGTTTCTAAAAATGTTTGTACATTATTTTTTTATCTTTATGAATCACAGCTGTTGCAGTACGTAGTCTATACTATACTTAAAGCACCAGTTTCAACGGTCGTCGTGCGTCATTTTTTTGCAAAAAAACCCTTATACTTCTTTCAAATGAACCCGCCATCCAATCTTTCTCATCCCATCTCTAAATGACCTAAAACCTAAACCCGTCCTTGCCACAGCCGCCACCCCCTCCTATCCCATTCACCGTCGCCACCTCCCTTCCCACCGCCGACTCCGCCCCTCCACTGGAGGCCGAAGCGCATATTCTCCCTTACTGCATCCACCCCCATCCTCACGACACGAATCTGCAATGAACCCCTCTCTAGGGTTCCTTCCTGCCACCTTGCCCCACGCCATCCAGAGGTTGTAAGGCTTGGGCGTCCACGCGTGGCTACTTCTTCCGTGGCGGCGGGAAGCACGACAAGGTGCCGTCGACTTGCTCCGCTCCCACGCATCCCTCGACCACGTGCTCGCCCGCGACGACAACAGTGTCGCCGCCGCAAGCACCTCAGGAAGGGTCCCATTGTGCACTTCCTCCTTGCCATCAATTTACATATCCCCTCTCCCTCTCGCAGGTCCCTCGCGACTAAGGCGTCACAGATCAAGGAGTAGTTCGACTGGCTGGACGCATTCTTCGATGGTGATGTTGCCTGCAAGTTGGTTAGGAGTGTGTATTAGATGAAGCTGCTCCTGCTAATTCGTTTTCGCCGCGGAGTGGCTGCTCCGCCCCGCCAGCGTCGACTGGCTCTGCTCCCACGTGCTCGCCCGTGGCGACAACCGCGTCGCCGCCACCTTCTGTCGCGCGCGCCTATGGAAGGACCCAACTGTGCACTTCCTCCTTGCTGTTAATCTCCAGGTGCGACCCCACCGCTTGGCCCTAATCGCATCACATCTAGATTCCACCCTGGGATAAGGAAGCACGTGGGGCCCATGAGCAGTGAAAAGAGATTCGTTTTTTTCTTTTGTGTGTCCGGATGACCATGCTTCTGATTCTGCAGGTGTGAGGCCATCTGGGCGCGTACAATGTTGTGTTCTACCGCTTCATCCACGGTGATGACGCGTACCACAACGCGTGCTTCAAGATCGCCAACCTGAGACGACGAGCCGGGATGGGTTTCAGGGTAAGCTCCGTGAAGGTGGACAACCATAGCCAACAACAAGACCATGTCAGCGACTAGGTGGGCAGGTTGATGAGTTGCCCGCAGCGAGACAGTGGGGATGAACACTTCAAACATAGTTTTGCAGTTAGCTTTACCACAAACTATCTGGTGGCCACGTATGCCAATAATCATGTTTCTAGCATATGCCAGCTCGTGATGACAAGTTCTTGTTTATGTAGCACATGGTCACATGCTCATCTTTTAGTTTTTGGCATGCCACTAACTGGTGTCTCAAGGTGGACTGGATTGTTAAGTGGAGCAGTTGAGTCTTACAACAAGCCATGAGAAGCAATGGCCTATACATCTGGTCAGTTTTTCTATACTATCTTGTTCTTTTTACCACTATGATTTGTCAACTAAAGCTAATTCTTTCCTAGGTTGATCTTCCCAGGCACATTGATTTGAATTGCATAGGCAACTAGGTTGTTGAGTCTGTTGTCAAGAAGAAGCTTAAGTAGGAACCTTTTCCTTGCCATTCCTCACTGCATGCTAATTTCCTCAGATGAATAATGTAGGTTTCAATGTTGAAATGCAGGTTTCAATGTTGGACTGTAGATAGGGGTGTTAACGGACCGTGATCCAATTGTTTCTTCACAATTTGTTAGGGCACTAGATTAATTTTAGTTCAAAAATAAATAGAAATAGAGCCGATCCTAATCTGATCTGATCCTTAAATTTTATATTGTAAAATTTAGAGACCATTGTCACCCCTAACTGTAGAGTACAAGAACATTAGCTTCACTGTCTGGGGTGTGAGGGGTCAGGACAAGGCATGATCCCAACATCATTTTGTATTACCATTTACATTATATAATGTTATCTAATTTCTTGTTCTTACATAATCCAAGAGCCATCATTCTTTTTAATTTTGTATTTGGTTTCTTAGACCGTGATCTGATTGTTAACAAATTACATTATGGCATCACATGAATAATGTAGCATGCTTCAGAGACTTGCATAGATCATGGAGCAGACTTGGTTATGAGCTCGAAAAAGTGATGACCTGGATGTGTTGTTCATGCAACAATGAGAAACTCAATGCCTCCACAAATGCCTACCAACTTGGGTAACTGAATATTTCACTTTTTCATTCTAGAACATCAAACCTCAGTTGCTTAAATGAGCTTATGTTTGTCATTTGGCTTCTATCAAAAGGTTGTTTTGTAAGTACTTTGTTTCAGCTATGATCTGTGAGGTAGTTTATTATTCTAAACGATTGCAGTTGGTGCTTTGGTTCAACTGATTAATGATCGGATGCTATCTCATTTAATGTCCTTTTCTGAGATAAACTGATTATATTGGTGACTTGCAGGTAGTTGCCAATTGCGTACATGCATTGCAACAGATCTGGGCTCTAGAAGCTGCTAACTTGGAAGCAGAGAGATTGAGGCATTGCATAGCAAACCATTGGTATTTTAATTTTACCTATTAAACAAATGAGTTTCAATCTGGCCTGGAGCCTTTTCTATTAAATAATGCCAATGCGAATTTTCTTGAGTTCTCTATTTTTATCTGCGGGATCAAAGAGTTTAGTGAATGGCCACAGTGCATGTTCTTGAGTTGGCATCAAAGTTCCTTCCATATAATAATAACAAAATATTTGACATAATGAATCTGTTGGAGGATAGACTTCAACATGCAAATGGAGCCGTTGTTTTGGCAACAATTAAAGTTTTTCTCCATTCCCATGACCGTTGTTCATCAACAGGTGCGCATGTTACTCCTAGAGGGGTGGTCTGGATCATGCTGCATTGTTGCTATTACTTCTCAAAAGCATTTCTTTCTTGTTCTAGAATATGAGATTTAGTATATGAGCACCTTCTGGTGAATGGTGATGCATAAAAACAAAATCTTATTAAGATTAAGATCCCATAACGGACTAGCACATGACATCACAACAACAACACACTTGGTCCATAGATAATGTTGTATTTGGAGGTTATTTTACCTGATTTTTTTAATTCCATGGACTCTTTATCTCTGTGGATTTCTTGGTTAATTTTTGTCTTATTTTTGTGTTTGATTTATCTCATCAACCTTTGGTGGAAGATCCATCAAATCTTTTTTAGCTCCTCAGTTTTTGGTTTATAAACGTATTAAAGCACCCTTGCTTACTCTGGTGGGTACTGGCAGTCCTAAATAGTCATAGTCTGTGCTATGCCATCTGCATCTTCTGGTGGAAGCGCATGCAACAAGGTTGTTTCAAAGATCAAGCAGAAGCTGTTAGCAGCAAACAAGATGTTTCAAACATTGGAGAATGCACTGAATACCACATTTTCCATGGACCATGGCTACATATGAGGAACAAAGTCTAAGAGAAAGCTTTACGTATTTATTTTATTCATGTTGGTAATGCAAAACTGGATTTAAACAAAATTGTACCGCTAAAAATATGCAATCCCTATCTGTTATATGTACATGGAATTACTACAGTTTTCACGTCCATGGTTCCACAGTCATCGTGTGACATCATCCTTTTCTTTCACTCTTTATATAAAAATAAGGTAAAATTATGTATAAATGAGGCTTAAACTTTGATTATGTTTCATATTCACCTAACAAATAGAACACATACTTATCTTTCATTAAAATAAAGTCTATTAATGTGATATATCGAAACCGTAGCAATGCACGGGCATTTGACTAGTTTAAATAGATAGCCGACTTTAATCTAAAGTGAACAGTCCTTCGTCGTCGACACGCCACATGATTTCTTTTTTTATTTTTTTGCAACCTGAAAAATATTTTCAATCAACAAGCACACTATTAAGGTCTCCTTTTTGGTAGGGGGTTTAGTGGTCCCTAAAAAATAGTCTTGGCTCTATCAAAATTAAAGCTGCTTTTTTTTTCTAGCGGGGCTAAGAGCATCTCCAATAGTTATGTAAACAAATCCCTAAATATGAATTTAACAAGTTTTCAAGAAAAACTCTTCTCCAACAGTTATGTAAATCAGTTCCCTAAATTTACAAACTCCCTACAAAACCCAATTTCATCCCTACATTTACCAACCTTGTTGGGACTCCCTAAATCAATGTTGGTAATAAATCCATCGTGAAAATCACACTTTTGCCATTGTTTTGCGTTGTCCCACTGTCCTGTTGACCACAAATCCATCTGGAGCGCCCAGATGTCTTGATGTCGTCACTGCTACACAACCTCTTACGTCGAGGCGTCGTCAGAGCACATCCTGTGATGTCTTCATGCCATCACAGTGACGGTCTGATCTCTCATCGGTAACTTCGACCTCTCAAGACATCAGAAAAGGTTGTTGATATTGATTTCTACCGCCGACTCAGTTAAAAATAGAAAGCATGCTCGCTTATTTTATGTTGATATAGACAACTTTTAAATTAGCAACTAATAAATAGGGAACTGTTGAAGCACATACATTTTTCTACTTCTTAAAGTATTTTAGCAACTTCTTAAATCTATAATTTAGGGTGCTAAATTTATATAACTATTGGAGATACTCTAAATATATAAAATAATGAGCGTACGTACGTACGACGTACGACTCGACTCAAGTGATATACTAATGTACGATCGGAATTATCTAATGGTTCATTACTTATTATTTCATAACAAAAGGCAATTGCCAGTTATAGATTAATGTTGGTCGTACTTCAATCGTACTACTGAGTCGATCGTACGTGTGAAGTACTTCTCATTTTTTTTTTTGCAAAACGACTTTTATGCTGATTGTTTTTACCAATAAAAACCCTGTCCGCTTACATTAGTAATATACTCAGTCCAAACTATAAGTCGTTGCACGTCTATGAACAGGGGAAAAAAATCTAACAACCGTAGGAGGGTCTATTTTTGTAAAATGTAACGATGGTGAGAAGAACACTTTGTAAAATGTACTGAGCTAAAGCTGGACAAGGCAGTTTTTTAGCTTACCTTTTCTAACGCAAACGAGAGAAGTGCTTTATTAAAATTGGTAGCTATTTGGTAGGGCTTCACTAAATATTAATAAAAACCAGAGCAGAATTCACATGTAAGAGCATCTCCATTAGTTATATAAATTTAGCGTCCTAAATTATAGATTTAGGGAGTTGCTAAAACACTTTAAGGAGTAAAACATGCATGTTTGTGCTCCAACAATTTTCTAATTATTCATTTCTAATTATTTAAGTTGTCCATGTTAAATTTTTTTTGCCTACTTTCTATTTTTGGTTGAGCTAACCATAGAAATGTGTATCAACAAGCTTTCTTAATGTTTTGGAGGTCGGATATGTGGATGGAAGGTCGGACGATCTCGGCGACAGGCACACACCATCCACAGAGGTTCTCACAGCGACACCGTAGAGGGAGCACATGAGGTGCTCTAGCGACACCTCAGCATAAGAGGAGCATAAGAGGTTGCCTGGCGGTAGCAGCGTCAAGGCGTGCTTGTCTGGTGCCTCAGCGTCAAATTGTTTGGGGTATTGGAGATGGATTTATGGTCATGGGTACGGTGGGACGGCGCAAACGGTGGTTGAAGTATGCTTTTCACGATAGATTTAGTACCAACATTGATTTAGAGCATTTCCAAGAAGGGCCTAAATAAAGTTGTATTTTTGAAATTTAGGGCTTGGGACCATAAAACGTCCTCCAATAGTGGCCCTATTTTATAAAAAAATCATCAAAAAATTATAAGGCTCGGTCTCTCGTGTCCTAAATATAATAACCTATATATACGAGTTGTATCCTCGTTTTCTTTTAATCTACAACTATTTATTTTCTAATAACATGATTTATTTCCTAAATAGCATTATCTAAGGCCTCATTGTTCAAGTAAAAGTTCTTTTGAGACCCTGAATACTTTAAATGTGGTCTTATTACAATTTTTAGGGCTTTATTTTAGAGTTCCTTGTTGGAGATGCTTATAAAAGTACATTCAGACCATTTGGCCTAACACGACCCAAGCTTGAAAAGGCCCATCACGATTGAATTTTGGGTCGTTTTGGCCTGAAATTATTATGGGTCGTGCTGGGCCGGCCCATGGTAAAACATGTCGGGCTCATTTTGGACAGCCTAAAATAATAAAAGCCTGAAATTCAATTTCTCGCCCGAAATTCACACATAACCCAAAATTTATATCAAAGCCTGAAATTCAAAACACATCTATCATTTAGTAGATATTAACATTTCAAGCTTTATAATTGAAAAAAATCAAAGAACACAAAATCAATAGTTGGTTGACACCTCATCTCCATGACTAGACATACGACTGGTCAAATTAAAGCGAATACTTTTCATGATTAGAAGAAGATGATGTTGCCTATAAATGGGCTATTTCGTGCCTACGTAAACGGGTCGTGCCCATGCCCGTCCGTGGGCTACGACCTCGCCCTAAGCATGGTCCGTACATAATTTCGGGTCGGTTCGACCCGAAATTATTTCGAGTCGTGTCGTGCTTGGGTCATGTTTTTTTCGTGCTTCGGATCAGCCCACCAGGTCTGGCCTAAATGTACATCTCTAGAGATGCTCTTAGGGGGTACCAACAAGGTGGGTAAATGTAGGGATAAAATTTAGTTTTATATGAAGTTTATAAATTTAGGAGCGTGCTTTACATAACTCTTTACAAAACTGTTGGGGAAGATTTTTTTAAATTTCTTAAATTTATATTTAGATTGTTTACATACCTACTGCGTATGGTCTAACCGTCCTAAAATTTCCCAAAGGCCCTGTTTGTTTCCCGTCATTATGAGAAATTGGAATCTAAATAATGGAGTAGACTATTTTTTAGAATATGACATTACATAACTTTCCAAAGTTAATATATAAGTCTATTTTAAATTCATGGGATGAGAGGCGGAAATTGATTCTATAGATTTGCATGCTACTTTTCTAATGTACAACTTATATCATACTTTTGTACTAGCTTCTCTATATCATAAATGTAGTGTATAACTATCTCTCTCGTATGATTTAGAATAATATATAAATTTATTACATATATAAATATATTAACTTAATTAGTTTGTGTTTAAATTATGATTATTAAAATGGAATTCAATTCTAACGAAACAAACAGGGCCTAAATCTTTTCTACCATATATGCAGATGCGGGAGTATAAATAATTTTCCTATACGGAGATACGCACACTAATAAACCGTAGGCAGGCAGGGGAAGCTCAGCTAGAATATAGATATCAATGGGGACCCAGTACTAATAACTCTTGAGGAATTTTTCTACTAAGGTTAGAATATAAGACAAAAAAATATCCTCATAGGTATAGACATGCAGTGGCGGAGGAACCATTATGGCTTAGTATGACATCAACCATACCTAAAATTTCAAAACCTTTATATCTTCTATTATTACCATAAACAGAAATGTGGGGCAATGGATAGGCTACTATGATCTAAGCAGCAGGCAGTGAATTCGATTTGTGTTGCCTGTGTATTTCTTTTTTCCCAGCCATACCTAATTTTTGTTTTGTCCTCCGCTACTGTAGACATGGGACAAAATACACAATCATTGGGTAAGCGGAGAGTTTGGAAACAATCATCCAAACTCGATTACCCATGGATAATATTCACCCGGTAAAATTAGTTGTGCAGCCTAATAGGCAATGAGCCCAACAACTAACCCAACCTCTTCTAACCCTTTTAAGCACTCATCCACCCAACCTAGTTGGTCAACCACCTACCGCTCCGTCCACCAGTCCATGTGTCACTTTGTGTTTTGTTATTTTATATGGTTAAAGAATGTTTTTTGTTACATGCCAGATTTATTGGTAATGTATTGTTGTTGGTTGTCTTGTGGTTGTTAATGTGTAACACCTTAGAATTAGGGGTATAAAATTTCTTTGCTGATATCCACAAATTTTAGGGTGTTACAACTCTCTTAGTTATCTTTCCCCTTTACTTTAGAAGGTGGGTTATTTTCTTTGGAGATGTAAGTATATCCTAGTAGATAAAATCCCAAGGGAGTGTGAAATGTTGCATCATGCTGAACCTCAAATCTTACCTTTGTTTGATGCATAAGCTTGAGGTGCTCTTCTTTTGAACTTGTGGCACACTTTTGAATTTGAATCCAAGGAAAAAAGGAAAAATAAAAACGGAAACAAAATTCAGAACAAAAGAAAAAGATAAAGAAGCATGCGCCCCTCCCAAGACAGCCTTTCGACCCAATTCTGGCCCACCCGGGCAGCCTGCCGCTCACCCACGCCCTGTCGCTGCCCCACTGGCCACAACCGTCAGTCGCACCCCACGCGCGCCTGCACCCGCCAAGCATCACACTGACCTGGACCCGCACGTCAGCCCCCTCCTATTTCAAAACCACCCTGAGGCCCACTCGTCGGCTTCATTTTCTCCGCCTTAACTACTGCCCCACTCTCTCCCCTGCTCGACGCGCCCACGTCGTGCGACGGATCCTGAACCATCCCGTGTGCCGGGCCCCCTTTTGAACCCCCTAACATGCCCGTTCTTGCCCCTCGCCTCATTCCTCACCTCACACCCGCCCCACCTCGCCCCGCTTGGAACCACCCGCGCTCGCACTCGCCAAGACCCCGCGCCGGAGCCCGTTGAGCCGCCGTGTTGCGCCCACCGTCGAGCGCGCCAATCCTCGACACACCCGGAGCACCGCCACGACGTCCATGCTCGAGCTTCACCGTGCTCGCCGTCGAGGGTGAGCGTCTAACCCCCACCCTCTCTTCTCGTCTACTCCCCTCCTACTCCTCTCTCTCGTCATGGCTGCCAGACCAAGCGGCCATGCCACTCGTCGTGCCCCGACATGGGCCAGCGCTCCGTCGGGGGCTTCCCTACCCCACATGCACGCGAGGCCGGGGGTAGACTGTGACCCTGATAGCTAGGGCCTGGCTGGCCATCTCCCCACACAAAACCCCTCCCCCCTCTTTGGGTCGTTGACCGGTGGACCCCACTCGTCATTGTAGTCCATGTGCCACGCTTGCGCCACCGCGCCTGCACCTGCCCCTGCCCCTGTGCCGCTGACCGAAGGGCTTGACCGTGCTGAGGCCCCACCTGCAACCCCCCGCACCAAAACCCCCTCCTTTTGTCCGTGCACTCCTGACCGGTGGGCCCTCCCGGTCCCACACGTCAATCGTCCACTTACCTCCGCGCGCGCAGCAATCGCACGGGATACTGTGGGATTCGTCCGCCCTCAGCCACGTGATCCACCGTGTGGTCGGCTCTAGCCGGTGGCCAGCCGCCACATGTCTGGGCTAGGGGTGGGCACTTTTTTACCATAAACCGAAAACCGAACCGAAATTTAGGTTTTTCGGTAGTTCGGTTCGGTTTCGATTTTTGTATCTAAGAAGTTCGGTTTTCGGTATCATAATCAATTTTCACTGTATACCAGACTGAAATACCAAAAAAACCGAATACCAAACTTTATCAATTCTCAGATTTGACTATTTGATTATGTGAACTAAATGTGTGATGTAATTAAATTGTTATTCACTTATTTGTATGTGATGTATGATGTAACTCTAAATATTTGTACCTATATAATTTTTAATTCTTAAAATTATTTGTAATCTATCATGTAAACTTGTTGTATTTATTGCCTTGAGTATCACTTTAGTATTTAGTTTTTACCGAAAAACCGAAACCGAACTTCTCGGTTTTTCATTTTCTAGAAAACCGATCAATTTCTAATGTCTGGAAAACCGAAGTTTTATAAAACTGAAATAGCGAACTGAAGTTTAAAAAACTGAATGCCTAGCCCTAATCCGGGCCACGCGCACCCGCCCTTATCGATCCGAGTCGCCCATGTCTAATCCAACGGCTGAGAGGCCTTGATACCCTTTTGCGCGACATTTTTGTTAAAGAGGCTCTCGTTTTCTGAGATTCAACCCACCGTCCTTCCACTTTTCAACCGAGCCCCTGGTAAATTGCGAATAGGTCCTCGAATCCTTTATATTAACCATAAATTCATTTTTAAAATATGTTTTAGGCCTAAAACTTGTAAAATTCATAAATTCATCATGTGAACTCCAAATCAAGTGATTCAAATTGCATAATTTTCATTATGTTTTTGTCTACCTGCTAGTAGTATGATCATCTACTGTTTTAATGTTTCAATTAATGGATAACTCCTAGTTTGAAAATAAACTGATTAGAGCCTTAGTTAAAATAATTCAAGGATAGAAAACTCTAGCAAATGTCCTAGTCTTGGGACCTAGCACCATCTAGTGATTAATACATCAATAAATTAACTCTATTTGAACGATATTATAGTGACATAGACTTGGTTAATGAGTAGTAAAGTAAGAAATAACCTTGTCTAGTGATCTACCCAAAACTTAGGGAAACTCCACAAATAATGTTACCATTTCTGTTAAACCTGTGAATATTCTGTTGCATATGTTTGGTGTATTGTTCTTTTGCTATTCCCAAAGTGTGTTGAATGGATGATTGTCTTTGTGTAGACAACGAGCAGTCTGTGTTTCCCGAAGGAGTTGCAGGAGAAGTCCCTCAGAAACAGCTTGGTGAAGGCAAGTGTCCTCTGACCTAATGATGTCCTATTCACATAATAATTCATTTCCCCGCATTACATAATTGAAACATAAGGATTGACTAGCTTTTCATTTACCTTGTCTTTGATAACCTTTTGGGCTATTATGGTTAGCTTTATGCTATCACTCGACTTTAATTAATGAACATGATGAGACTACTTTATGATACAATGTGTTTATTTTGGATATGATGATGAGGTTGTAACACTAAAGGAGGCTCGGGCTATTTTACGAGTAGTTGTCTGTAAGGACCTGTTCGTTGGATGACCACCCGATAAAACAGTGCAACCATGTGGGTGGTATGGGATGGCTTTAGCTAATTGAAAGTCGCCTAGAGAGGGGTGGATAGGCGGAAACTAAAATTTACAACTTTAAACACACTACAAGTCGGGGTTAGCGTTAGAATAAAATCCGAGTCCGGGAGAGAGGGGAAAACAAATCAACCAAGAAAATAAAGCGGATGACACAATGATTTGTTTTACTGAGGTTCGGTTCTATAGAACCTAGTCCCCGTTGAGGTGGTCACAAAGACCGGGTCTCTTTCAACCCTTTCCCTCTCTCAAACGGCCACTTAGACCGAGTGAGGCTTCTTCCTTAATCTCACGGGTCACTTAGACCCCGCAAGGATCACCACACAATTGGTGTCTCTTGCCTCGCTTACAAAGCACTTGAGAGTAAGAAGTGAGAAAGAAAAGAAAGCCAAGCCAAGCAAACAAGAGCAACAAGAAACACAAGTGATCCTCTCACAAGTCCTAATGAACTAGATTTGAATTGGGGACTTTGATCGGATCAGTGGTTTTGATTGGTGTCTTGGAGTGTTGCACTTTGCTCTTGTATTGAATGAGGAGTAGTGAATGCTTGGATGGTTGGAGTGGAGGTGGTTGGGGGTATTTATATCCCTCAACCACCAAAACAACCGTTGGGGGTGGCTGCTGTCGATGGGCGCACCAGACAGTCCGGTGCGCCACCGGACACTATCCGGTGCGCCAGCCACGTCACCCAATTGTTAGGGTCCTAGCGGTTTCGACCGTTGGAGCTTTGTCTTCTTATGGCACCGGACAGTCCGGTGCCGCACCGGACAGGTATTGTTCATTGTCCGGTACGCCTCTGACGGCTGCTCTGACTTCTGTGCAAACTGTCCACACACTATAGCGCTTTGCAGGTGTCCGTTGCAGTCGACCGTTGCGCTGGGAGCTGTTGCTCCGCTGGTGCACCGGACAGTCCGGTGGCACACCGGACAGTCCGGTGAATTATAGCGGAGCGGCGCTGGAGAAACCCGAAGGTGAAGGGTTCATCTTTTTGAACCCTAACTTGATCTTTTTATTGGTTTGTGTTGAACCTTTAGCACTTGTAGAATATATAATCTAGAGCAAACTAGTTAGTCCAAATATTTGTGTTGGGAATCCAACCACCAAAATTATTTATAGGAAAAGGTTAAACCCTATTTCCCTTTCACTAATTAATTGGAAGAACTTGAGATGTAGTCTTCTTCGTCGTCGTGCCGTTAATGGGGTCCTAGCACAGTACTTGCTATACCGAGGTTGGTACCAAGGTTCTTTTGTTTTGCTTTTGTTGGACACCCCATGTGGGGAGGGGTACTATGTTTATCAAACTGTAGAAACCTAATAGGCGACTTTGACCTCTGGAGAATCTTTGTAAATGCTACATAGTGAAACCTTGTTGACTCACCATGGGAGTGTTTAAGGGTTTGATCAACTTATGGCAAAAAGGGGTCACGACTCATGAGTAAAGTGTACGACCTTTGCATAGGGTTAGAAACTGATATATCAGTCATGCTCACAATTAAGAACAGCCTTGGGAGCTCCTTTGATTAGAGATACTGTAAATACATTCATGATGATGGTTTGATGATGGTGCCTCTAATTATGATTTCTGGTATTTTCTCTACGAGGAGGTACTCTTTGGGATAATAAGCTAGGTTTTAAGATAAAATTTGGCTTATATTAATGATTAAAACCTGATAAAGTAAAAGCAACATGCTATCAGCTTAACTCCACATAAAGCTAGTCCATTTTAGCCAAACAAGATATTTGCTAAGTACGTTGATGTGTGCAAAATGGAGAACTTTTATCTTAAAACACCAGGTTGTCCACACTACAACCATTGCTCAAGCGAGGATGAAGGCAACATGAAGAACTTTCAGGAGTTTCTAGACTTCAAGGAGTTTTAAACTAGATTAGTGGTAAACCCCAGTCAGCTGTCTCTGAAGGCCTTATCTTTACTTTGTTCCGCGTTAGCACTTTGATTACTTGTTAAGTTGATGGATACATCATGTTGTAATTGAGAGCACCTAGAGGGGGGGGTGAATAGGTGATCCTGTAAAATTAAACACTAAATAGCCACAAAACTTAGTTATAAAAGTGTTAGTGTGGCTAAGTAGTTGAGAGGCGAGTTCTTGTGGACAAAACAATCACAGAGAAAGCAATCACAAGAGACACACGATTTTATCCCGTGGTTTGGCCAAGTAACACTTGCCTACTTCCATATTGTGGCGTCTTAATGGACGAGGGTTGCACTCAAACCCTTTCAAGTGATCCGATGATCAACTTGAATACCACGTCTTTTCCTTTAGAGTATCTTCCCGATTGTGAGGAATCTCCACAGTTTGGAGTCTCTCGCCCTTATAATTAGGATCACAAAGAAAGCACAGAGTAAGGTCGGGAGAAGCAATGCACACACAACACAATTTTGTAGCACACACACGCACACAAGCCAAGACTTGAGCTCGAAAAGTAGCACAAAGAGTTCACAACTCGAACGTAGCTCAAATCACTAACACAATCGATCAAATGCGCGGAGGCGGAGTGTGGGAGTCTTAGAATGCTTAGTGAATGCTTAGATGTTTCCTCCATGCGTCTAGGGGTCCCTTTTATAGCCCCAAGGCAGCTAGGAGCCGTTGGAGATCAACATGGAAGGCTATCCTTGCCTTCTATCGAGTGGCGCACCAGACAGTCCGGTGTGCCACCGGACAGGTCCTGTAGACTGTCCGGTGCGTGATCTCCTTCCAAATTTGGCATATCCGACTGTTGCTCCTCTGGGCTAATTGGCGCACCAGACAGTTCGATGCACCAGCTGACCGTTGGAGCAGTCCACGTGTCGCGCAAAGATTGTGCGGCCGACCGTTGCTCAGGCAACCGTTGGCTCACCGGATAGTCCGGTGAATTATAGCCGTACGCCGCCGTCGAAACCCGAGAGCGGCGAGTTCACCGCGGACCAGCCTGGCGCACCACCGGACAGTCCGGTGTGCCAGACTGAGCACAAGATTGGCTGCACAGAGCCAAGCTTTTCCCTTTTTCTCCCTTCTTCTTTAGTCACTGTTTCTAGCACTTGGATAACCATGATAGTACATAAAACAATTCACCAAGTCTAGAAACATACCTTTTGCCTTGATTTTCACTTCTCACTTTATTTGGCACATAAGAAATTAATAAAACGTGTTGGGCACTTAATCACCAAAACATTATCGAAATGGCCCAAAGGCACATTTCCCTTTCAATCTCCCCCTTTTTGGTGATTTATGCCAACACAACCAAAAGCAACCAAAAGAAGTGCAACATCAATGCAATTAAGGACCAAATTGTTTTTGAATATAATTTGGCATATTTGGATTGCTCTTTGCCACCACTTGGTTTGTTTTTGCAAATCAAATTCATTTTCCTATCTCTAAGTCAAACACGCTTGTTTGGGCACAAAGAGAGATATTCCAACAGAAAAATTGATCAAGTGTCAAAAACTCCCCCTTTTCCCATAATCAAAATTCTCTCCCACAAGAGACCAAATTTTGCAATAAGAATATTTTGGACACATCAAAAGTTCTAAATCTATTGTTTTCAAAATTCACAAGTGGTTGCTGATCCATTTGCTTTGGCCTTAATTTCTCCCCCTTTGGCATTAAGCAACAAAACGGGATCATTTTTTGCCCTTTAAACCCCATTGCCTCACCAAAAATGTCAATTTAAGAGCAAATGGCAATGAAATTGCAAATATGAACTTGGAAGTAAGTACACTTATACCGGAGTGCAGTGGAAGTCTTTTCAACGTCCAAGTTCATTGTTCACTTTCAATGCACCTTTGAGACTACATCAAATGTACTCAAACAAACATGTTAGTCTCAAAGGGTCAAGTTGTAGCACATCTCCCCCATAAATGTGTGCACCATTTGCATATGGACTTGTGAGGTCCGGGGAGATTTTGTACAACTTGAGCACCACAAATAGATAACACAAAATGCATAAAGTAACATGATCAAAGGCATAGAACACATGTATGCTATAGATCAATCCAAGTTACGTGAATCTAAGACATTTAGCTCACTACGCAACCTGCAAAAAGTTTTCTCATCCAACGGCTTGGTAAAGATATCGGCTAGCTGGTTCTCGGTGCTAATATGGTACACCTCGATATCTCCCTTTTGCTGGTGGTCTCTCAGGAAGTGATGCCGGATGTATATGTGCTTAGTGTGGCTGTGTTCAACAGAATTATCCACCATGCGGATAGCACTCTCATTGTCACATAGAAGTGGGACTTTGCTCAGATTGTAGCCAAAGTCCCGGAGGGTTTGCGTCATCCAAAGTAGTTGTGCGCAACACTGTCCTGCGGCAACATACTCGGCCTCAGCGGTGGATAGGGCAATGGATGTTTGTTTCTTAAAACTCCAGGACACCAGGGACCTTCCTAGAAATTGGCACGTCCCCGATGTGCTCTTCCTATCAACCTTGCACCCGGCATAATCGGAGTCTGAGTATCCAATCAAGTCAAAGGTAGACCCCTTTGGATACCAGATCTCGAAACAAGGCGTAGAAACTAAATATCTAAGAATTCGCTTAACAGCCACAAGGTGACATTCCTTGGGGTCGGATCGATATCTAGCACACATGCATACACTTAGCATAATGTCTGGTCTACTTGCACATAAATAAAGTAATGACCCTATCATAGACCGGTATGCCTTTTGATCAACGGACTTACCTCCTTTGTTGAGGTCAACATGTCCGTCGGTTCCCATCGGAGTCTTCGCAGGCTTGGCGTCCTTCATCCCAAGCCGCTTGAGAAGATCTTGTGTGTACTTCGTTTGGGAGATGAAGGTGCTGTCCTTGAGTTGCTTTACTTGGAACCCAAGGAAGTAGGTCAACTCTCCCATCATCGACATCTCGAACTTTTGTGTCATCACCCTGCTAAACTCCTCACAAGACTTTTGATTAGTAGAACCAAATTTTATGTCATCGACATAAATTTGGCACACAAAAAAATCACCAACACAAGTCTTTGTAAAAAGAGTGGGATCGGCTTTCCCAACCTTGAAGGCATTAGCAATTAAGAAATCTCTAAGGCATTCATACCATGCTCTTGGGGCTTGCTTAAGTCCATAGAGCGCCTTAGAGAGCTTGAACACATGGTTGGGGTACCTGTCATCCTCAAATCCAGGGGGTTGTTCCACGTATACCTCCTCCTTTATTGGCCCATTGAGGAAAGCGCTCTTCACATCCATTTGAAATAGCCTGAAAGAGTGGTGAGCGGCATAGGCTAATAGAATTCGAATAGACTCTAGCCTAGCCACATGAGCAAAAGTCTCCTTGAAATCCAAACCTGCGACTTGGGCATAACCTTTTGCCACAAGTCGAGCCTTGTTTCTTGTCACCACCCCGTGCTCGTCTTGCTTGTTGCGGAACACTCACTTGGTTCCCACAAGATTTTGCTTTGGACGTGGCACCAGACTCCAAACTTCATTTCTCTTGAAGTTATTGAGCTCTTCCTGCATGGCCAACACCCAATCTGGATCCTGCAAGGCCTCTTCTGCCCTGAAAGGCTCAATAGAAGAGACAAACGAGTAATGCTTACAAAAGTTAGCTAATCTTGAGCGAGTAGTTACTCCCTTGCTTATGTCACCCAGAATCTGGTCAACAGGGTGATGTCATTGGATCATTGCTCGGACTTGAGTTGGAGGGGCATGTGGTGCTTCTTCCTCCATCACTTGTTCTTCCTGTGCTCCCCCTTGATCCTGTCCCTCTTCTTGAGGTACTTGTTCACTGTCTTGAGTTGGAGGATGCACATTGTGGAGGAAGATGGCTGATCTTACTCCTGTTGTTCCTATGGTCGCACATCACCAATCGCCATCGTGCGCATTGCGGCCGTTGGAACAACATCTTCATCTATGTCATCAAGATCAACTTGCTCTCTTGGAAAGCCATTAGTATCATCAAATACAACGTCGCTAGAGACTTCAACCAAACCCGATGATTTGTTGAAGACTCTATATGCCTTTGTATTTGAGTCATACCCTAGTAAAAACCCTTCTATTGCCTTGGGAGCAAATTTAGAATGTCTACATTTCTTCACCAGAATGTAACATTTGCTCCCAAATACACGAAAGTAGGAGACATTTGGTTTGTTACCGGTAAGGAGTTCGTAGGAGGTATTCTTGAGAAGGCGATGCAGATAGAGCCTGTTTATGGCGTGGCAAGCTGTGTTCACAGCTTCCGACCAAAACCGATCGGGCATCTTGAACTCTCCAAGCATCGTTCTCACCATGTCGATAAATATCCTGTTCTTCCTCTCTACCACACCATTTTGCTGTGGTGTGTAGGGAGCGGAGAACTCGTGCTTGACGGCTTCCTCCTTAAGATACTCCTCAACTTGCAGATTCTTGAATTCGGACCCGTTGTCGCTTCTTATCTTTTTCACCTTTAGCTCAAATTCGTTTTGAGCCCTCCTTAGAAAGCGCTTTAGGGTCCCTTGGGTTTCTGATTTATCCTGCAAAAAGAACACCCAAGTGAAGCGGGAAAAATCATCAACAATTACAAGACCGTACTTACTTCCCCCGATGCTAAGGTAGGCAATGGGTCCGAAGAGGTCCATGTGAAGAAGCTCCAGTGGTCTTGATGTTGTCATCACATTCTTGCTATGATGAGTGCTTCCCACCTGTTTCCCTGCCTGACATGCTGCACAAGGCCTATCTTTCTCGAAACAAACATTGGTTAGTCTTAACACATGTTCTCCCTTTAGAAGTTTATGAAGGTTCTTCATCCCAACATGTGCTAGACGGCGATGCCACAGCCAGCCCATACCAGTCTTAGATATTAAGCATGCATCTAGATCGGCCTCCTCTTTTGAAAAATCAACTAAGTAGAGTTTGTTGTCTAATACACCCTTAAAAGCTAGTGAACCATCACTCCTTCTAAAGACAGATACATCAACGTTTGTAAAAAGACAATTGTAACCCATGTGACACAATTGACTTACAGACAGAAGGTTGTAACCAAGCGACTCTACTAGAAATACATTCGATATGGAGTGCTCGTTGGTGATGGCTATCTGCCCAAGCCTTTGACCATTCCTTGGTTCCCATCTCCAAAGATGATCGTATCCTGGAAATCCTTGTTTTTGACATAGGAGGTGAACATCTTCTTCTCCCCCGTCATGTGGTTTATGCATCCGATGTCGACAATCCAGCTTGAGCCCCCGGATGCATAAACCTACAAGGAATTTAAGCTTGGGTTTTAGGTACCCAACTCTTGTTGGGTCCTACAAGGTTAGTCAAAATAGTTTTTGGACCCCAGATGCAAGTTTTGTCTCCCTTGCATTTGGATCCCAACTTCCTAGCCACTACTTTTGCATTCTTACATGAAAGAATAAAATAAGTGTTGCAAGCATGAAAAACAGTAGTGGGTTCATTGCACATTTTCCTAGGCACATGAGACACAACATTATTCCTCCTAGGCGCATGAGACACAACATGATTGCGCCTAGGCCTATTTCTACTACAATCATAAGTAGAGCTGGAAGCAATCATAGCATGATCATAGTCTCTATTATTATAAGCATTCCTAGAAAATTTTCTATCATAAATGTAGGCATGACACTTTTGAGAACTACTAACCATAGGGGCCTTCCCTTTCTCCTTGTTGAGAATAGGAGCCTTTTGGCTTGTTAAGTTCTTGGTTTAAACCAAGTCCATCCTTAATTGAAGGGTCACTACACCATGACACAATTTTAGCGTCAAACATCTGACGATAAAAATAGACGTTTTGCGACAGACAGTGTGTCACAAAATATTACCATCACACTATCTGTCGGTATTTATTGTATTGCAGTCAACAAATCTGTCGGCAATAATGTCGAGTTATTGTGGCAGTCCATCTGTCGGTAAAACTGATTATGCGTGGCGCACAGTCTGTCGGAAAAAAGTAGCATCAGATCGTTTGTCGGTAAATGGATGTCATGATCCTATAGGAAAAAAATAAACGAGCGTGAGTACGTATGTAGGAATCAGGAGTCGAACCCGTGACCCTTTGCAAAAAGAGAAGACACAAAACCATTAAACTAGTGATGATTATTCAATCTGCTTTTGGTATTGTTTTTTATAAACATGTACACGTTAGTTATTGAAGGGGATACATCACATCGGCCCAGTAAAAAACCTGGCGCCGTACGTGTCCCGCGCGACTCTAGCAACAGGCGGCGCCGATTAAGGTCGCGCGCGAGGCCTGCACGGACGACCAGGCGGCGGCGGTTGAGAGGGCCGAGGCCGGCGCAGGCGAGCAGACCGGCGCGGGCGAGCAGGCCGGCGGGCCAGCAGACCGGCGCGCCGTCACCCGTGCGTGCGTCCCTAGACCTAACAGCAGCGACTGCGCCCGGCGTGCAGGCTTCTTCATCACCTGCGACGGCGGCCGCGAGCACCAGGCCACCCGCGAGGCGCTCTCCCTCCTCGACTCCGTAAGCTGCTCCCCAGTTTTCCGTATCTCTGTGTTTGTCATTGGCCGGCATGTAGTGTGGCCGTTCCTTGCGATGCACGCCCGCTAGGTGTTCGTCAAAGTGTTCTTGTCCTTTCTTAGAGAAGCAACCTCTCCAGACTGCAGCGAAAGTTTTGCATGGGGAGATGTGATGGCTTGTGATGGCAGCAAAGTTTCGCATGAGAGACGCATGGTGTATGTATAGCGAGCAGGAAAACACGTCCATGAAGCAGACCAACAAGAACTTTGCAGTGGCGGGAGGGAGACGAACCTGTTGTGGCGGTGGCGAGGAAGGCGATCTCCGATTGGCGGCAGTGAACAGGAAATCTCAACTCTGCGATGGAAAGACGTATCGATGGATAGATTGCAGGGGCAATGTGGTCAGCTGAGAAGACGAAGGGAAAGGAGTAGGAGATAATTGAAATAATTTTAAAAAATATTAGACTGCTATTAAGAATGACAATATCACTTCTAAATACTATTAGTGAGACGATGGAGTTCAGAGTACTCTAATAGAGAAGCCTGTTAAAACAAATATTATAATAGCAAAAAAACTCGAAGAAGACAACTGTTAGGCCCCGTTTGTTTAGTTGGAATGCCGGCGCGGGCGAGCAGGCGGCGGCGATTGAGAGGTAGCCTACTAGTGTTAGCTCCCTGAAATTTTGCTTGTGTGTTGGTTGCGATTGGATACCTGGTTTTAGCATTATTAGTAGGGTTGTATAATGTGATGGAACTCCATCTCCGACGCGGTGGCAACCAGCTGCCGAGGGAGGAGCCCTAGCGCGCGGAACAGTAGACGGCGGCGCGCTATCGGTGACGGGGATACAGCGACATCGTCCCCCTCGGCCCCTCCTCCATCTCCGGCGTGGCGGCGCCAGTAGCCGGCGTGCAGGTACGGGGTCGAGTCTCCGGTCTCCTTGTTTTCTGTAGTGTGCTCTGTTGTAGTCTGTAGGCAGGCGCGAATTAGGATTTTGGATTGGGGGGGGGGGACAGGCCGAACAGCGGGACGGGGACGAGACCAATCCTTCGAATATCAACTATCAACATAAAAATGTCTAAAACAACATAGTGCAAAATGATGCAGCGCCAAGATATATAGAACTATTTAATGCCTTGACAAAAGTAAAAATATGCCAATATATGTTTTGAGATAAAAAGACACCTGTTCAATATGTTCCAGCTTCTATCCAATAAATTACATATTTATTAACTTGAATTTAGCTTTCCGAGAACTAGCAAGATCATAGGAATCAATGATGTCATCAGAACTGATAGTCGCAGCCAACTCTTTTTCAATGTAAACAACCAAGTAATCGCGGAGAAAGGCATCTCCAATTTTATTCATAGCTGAAAAAGCTCTTTTGGTGGTCGTTGTTGATACTGTAAGAGTTATGACTAATCTTAGCAACCTTTCCACCATTGGATAAGAAGCACTTTTACCTGTCTTAACTAGCCCACATGTCAAGTCAGCATGGACATCAGCAGCTACCAAAGTGCTAGAATCAACTGGCTTAAAGAAACGATGAAGTGTTGTTCTTGGTCTCTTTTGACTCTTCATATCTTTATTGAACTGAAATTGAAAGCTAAAATATTACCATCACACTATTTGTCACAAATTCGGCTTTTATCTTTGATACTTGTAATGATTCCCCCTTTTCCAGCCTTTTTGCACAAATATTGTGAAGCTGGGAAGGGTATATTTCATCTTATTATTTGTTAAGGTACAACATGGATAAAAGGTGGATGAATGAACCTAGGTATATTTCTGACAATGTTCTATAAACCATTTTTTGTAGTAATCACCTTTGACCAACATATCATTTTTGTACCTAATATTTATAGGCACACGAAGGTCTACATAGATGGGATAGCTGGATTTATTGAGTTTGCGTTTAGCAATTCTATAGGTGGAAATAGAATTCTATGTCCTTGTAGAAAATGTGTTAATTCTCTATGGAAAGATGAAAATGAAGTTCGTGCGCACTTGATATGTGATGGCTTTCAAGAGGGATACAAAAGGTGGACTTATCATGGAGAAACAGACTCTATGCCTATTAATAGTCATGACTATAGTGATATAGAGGCTGCAAACAATTCAAATGAAGATGACATAGCCAATTTGGTTCAAGACATGGCTGGTGGTCTAGATGATAAAGGGGATTTTGAAGTATCTGATGGTTTAGATGAAACTGATGTAGATTTAGAAGCCATCAATCAGTTAGTGAAAGATAATACCCAAGAGCTATACCCAGGGTGTAAGAAGTTCTCAAAGTTGCATTTCTTAATAAGAATTCTTCACTTGAAATTGCTTGGAGGGTGGACAGATAAAAGTTTTGATATGCTACTAGATCTACTAATGGAGGCTTTCCCTGATGGCTTGGCATTGCCAAAAAATTTCAATGAAGTAAAGAAAGTAATTAAGTGTCTAGGACTCGGGTATATCAAAATTGATGCATGTGAAAATTATTGCATTTTGTTTAGGAAGGAGTATGCTAATTGTGATGCATGTCCCACATGTGAAAGGTCACGCTGGAAATAAAAAGGAGAAGTTTAAATGGAAAACGTGTACACAAGGTTCCATGCAAAGTTCTTCGCTATTTTCCAATAAAAAGAAGGCTTCAGAGGTTATTTCTATCATCTAAGACAGCAAGTGACATGAGGTGGCACGATGAGGAGCGTATACAAGATGGATTACTAAGGCATCCTGCAGATTCACCTCTCTGGAAGGATTTTGACCAAAAGCACCCAGTGTTTGCTTCTGACAGTCGCAACATTAGACTTGCTTTAGCTACAGATGGTTTTAATCCGTTTAGGTCCATGAATGTTAGCTATAGCATTTGGCCTGTTATTCTCATCCCATATAACTTTGCACCATGGTTGTGTATGAAGCAAACAAATTTCATTATCTCATTGCTGATTCCGGGCCGTAGATCTCCTGGTAGTGATATAGATGTTTATTTTGAGCCACTAATTGATGATATGCTAGATATGTTTGTTGAAGGGGTGAGAACATACGATTCTGCAAAGCATGAGTATTTTCAATTACGTGCTGCAATAATATGGACCATATCTGATTATCTAGGCCTAGGGTACATTGCAGCATGTACTACATCAGGTGAAGTAGCATGTATTGAATGCCACTCATATACATGTTCTCTACGATTGAAACAAGGTACAAAAAATTGTTATATGGGACACCGTAGATTCCTGGGCCCAAACCATCCATATAGATTTGATCTAGATTCGTTTGATGGTAAAGTTGAGGACAGGCCAACGCCTAGACCACTTTCTGGAGAGGAAGTTTTATTGCAAACAGAGAATATGCACACAGTGCATGGGAAAGGTTGCCAACATAAAAAGAAAGCATAAAAGAAGGAAGGGAGCCTACTATCATATGGAAAAGGAGGTCTATTTTTTTAGGCTTCCATATTGGAAGGATTTAATGTTACGACATAATTTAGATGTCATGCATATCGAGAAGAATGTGTGTGATAACATTGTCAACACCCTCCTAGGCATTCAACGAAAGTCCAAGGATAACTTGAAATATCGTTTGGATCTTCAGTCCCTAGGCATTCGAACTGAGCTTCATCCTATTCTAGTGGGTGACAAGTTATATTTACCACCAGCATCATATACAATGAGTGCGTCCGAGAAGACGTTATTTTGTGAGGTACTGAAAGGAGTAAAATTTCCAGATGGATATGCATCAGATATAAGAAACAATGTGGATGTGAAGAAGAAGAAGCTTGTTGGGCTAAAGAGCCATGATAACCATGTTTTGATACAGTATTTACTCCCACTTGCTGTGAGAAGGATTTTACCAGAAATAGTTAGTGCTGCACTTATTCGTGTGAGCAATTTTTTCAAGCAAATCTATTCGCCAGTTATTCGTATAAGCGATATGCATAAGCTAGAGTCAGAAATAGCTGAGACTCTAAGCATTCTTGAGACTATATTCTTGCCATCCTTTTTTGATATTATGGTACACTTAATGGTTCATCTCCCTACTCAAGTCAGAATTGCTGGCCCTGTTCAATTTCGTAATATGTACCCAGTGGAGAGGTAAATTTTAAAATATTTTAGCAATAAATTAAACTTAGTTTTACTCAAAGGAAGGTTTAACATATATGCCATTTCTCCTTTATAGGTTTTTAATGAGATGTAAGGGCCATGTTCGCACTAGGAGTCACCCAGAGGGGTCAATTGCAGAAAGTTATCTATTTGATGAGAGTCTCACATTCTGTTCTCGCTATTTACATGGTGAAACTAGATTTAACCGGAAAAAAAGAAATGATGATGGCTTCAATGTTGATTTAATCAATACCACCCCTTTCTTCCATAACATAGGTCGTGGATTGGTTGGTAAGCGTAGTGTCACATTAGACCACAAAACATGGCTTCAAGCACATAGATATGTCTTGTTCAACTATGACCATATAGAGCCTTACTTGAAGTAAGTTGAACATTATATTTTATAATCAATATCATCCATTTATTTCATGATTACTAAATTATAATGTTTGTTTAGTAAACATATAGAATACCTTTACTCTGCTGGTCATCAAAATCAAAGAGCAATCAGTCGTTTACACTATGAGACTTTTCACGAGTGGTTTAAGTCACACGTAAGTATCTAGGAAACTGTTTATATAGAGTATAGTTCCACATTATTTTTTTAGTGCAATAACAATCTTCAATATCATCAGGTGGAAACAACGAGTGAGGAAGTACCTGAGGAGATAACAATATTAGCTAAGGAGCCTAACATGGTTGCACATTCTTATGATAGTTATTCTATAAATGGGATAAATTTCCATACACATTCTTATGATGTGGGTAGGTCAGTGCAGTGTAGTGGTGTAGCCCTAGTTACACATGCAACAAGTTTTGATGGGACAAACAATAATAACCCAGTATCAATGAGCAAAACTTATTATGGTGTTATAAAAGACATTCTTGAGCTGAACTACCATCACCAGGGAAAAATAGTATTGTTCAAATGTGATTGGATTGATAACCGAGTACGAGACAAATGGGTCAAAATAGACAAGTTTGGGGTGACAATGGTCAACTTCAAGCATTTATTCAATACCGGTGACAAGGAATTAGATGAGCCATTTATTTTTGCATCTCAGGCAACTCAAGTTTACTATGTGCAAGATCCTATTGATGCTGATTGGTTTGCTGTTTTAAAGTCAAAACAACGTGACATGTATGATATGGAAGAAAGGCTAGACAATAGTACAGAAATAGGTTCTTTCTTGCCAGATTTGGATGCAAACACGCAAGTGAATGTATCTATCGGTGGTAGTTGTGTTAGGACTGATATAGATGGAATAATGGTTGCTGAAAACCAACCTAGAAAGTAAGATATGTTTATATCTGTTTACGCAAAAATGTTTCATTTTACGTTGTAATCTCTTTGTCTTGAATCTTGATTGTAATCTCTTTATGTCTTGATTTAGGAGGAAAGTTTGTAAGAGGAAGAGGAATGGGTAAAGAAACGGAGGTCAATGAATATGAACTACAAAGGGCTGCCAACATAAAAGAAAATTTGAAGCGGATGAGATCTCTTAATCTCCATGTACCTAGTACTATGGTTAGTGAAGAAGCAAATGGATCACATCAAGCCAATAAGAAACATAAGGTTAATGTATTTCCCTATTTAATAAGTTCATGTTTTTTACATTTAATTTAATAGTACCACAATCTATTTAATACTGCTTGTACTTTCTTTCCAATATTTCAGTTTCAGACTTTCAGTTCTGCTGCTACTGCTATATTTCCATTCCAATTTCAGACTTTCATTTCTGCTGCTTCACTGCTTGTAATTTCAGACTCAGTGACCATTGAGTTCACCAGTGACCATTGAAATTCACCAGTCACAAGTAACTAGATTTTTCTGCCATAGCCACATAAATTCACCAGTGACCATTAAAAATGGGCTGTTAGGCCCTTGGGATTTTCAGATGATAAATAGCAGTCAGTGCCAATTTCAGTGTTTGACTCAAACCTTGGAGGGTTCCTGTAGCTTAGTTGTAGGGATTTGGGCAAAGGTATCTTTGGTAAGTGGTAGCTTAATTCAAGAAGAATAAGGGGACTGAGGTAATTAAGCCTGACACTACATAAAAAATTTAAAAGAAAATCAATGCTCATATGGACTGACAGGTAGTGGCAACTGAGGGTTATTCAGTCGCTGAAGAACCATAGGACTGCATGTGGCACCTAGATTGGGGCTCTTTGGTACTCAAACCATGCAGTTTTTGGCTGAAGCTTCTGTAGTAAGTTGGTAGACCATAGACCTTCTTATAGGTTGTATGGCTTGGTGCTACACAAAATGTGAGGAAAGTGATTCCCAATTTGCCCTGTCAGACTGTGCAGATGGAGTTGGTTTCAGATAACAGAAATAAAAACCAGGGCAGGATGCTTGACTTTGGGGCAATGTAACACAAAATTCATGTATTATTAGATATGTAATGTTATAGCAAAGTGGAAGACCAATACTTGGTGAACAAGTTTGGTGTAATGGGCAAGGTCTGTCCTTATGTGGAATTCGAAGGGAAGTGGGGAGGAAATTAATATTTCAGTAGATTGATCAGATTCTAGTTTTTCTTTTTCAAAATGTCTGAAGAATATCCTGATCCCTTATCCTTTCTAGCAGATTCTGGGCAGCTTGTAGATAGAGTTAACTTGGGATTCGAAAGTGGAGTATTGGAAATATTTTTTAACTACCGTGAGTGCAGAGAGGAGATAAGTGTGGGGCATGTTCCTTGATGATGACTCAGTTGCTGATTTGCTTTGTTGAAGATGACATAGCCAATTTACTAACTAGTTTTACCTCAAGGAACTCCACTGTTTTGCATAACTCAAAGTTCATTTTTTTCGACTCATCCAATTCCTGCTCTTTTTGCGTAAGACGACCGTGTAATTCCTGTTCTGTCTCCGCCGTAAGAGGTAATCTGACGTTCTGCCAACTGTTCAGATCAACAGCGCTACTGTTCTTTCACCCTGTGTTCTTCAACCCATTGTCAGCCAAAACATTGCAGCCCCTCCCTCCTGAGATCGCCCCATCCCCGACAGTGATGCTATCCTCAAGTCTGCCTACGCCAACCCTGCTATGCTAAACAGTTTTTTATGGCGCTGATGATTGACAGTATCTCACCAGGAAAAATCGCACCACCACTAACTGCAGTTGACTGAATTAATGTCCAGACGACTGAAATCGAGGAGGGGTGTAGTAAAAGAAGATAAATATTGCTAGCGCACAAACAAAAACACGGAGCACACCTCTAATTGACCAAACAGTGCCTAAGCAGAGGAACATAAGTGGTGAAGGCATCAAGAAGCATAACTCATAAACAAACCTAGTATATACAATTGTAATCAGATGCAATACATGGCAACTGCTGAGTTAATTTTATGAATTCACCCACCAATCACACAATGTAGAAAGCATCACCAAGACATCAAGTAATTTCAAGTAAGCGAAGTAGCGCACTCATCAGCAGAGGATGGACCATTATTAGGTCAAAAAAGGGACACGAGGCCGTGCCGTAGCAGTAACATAAACCAGAAGTCACATCGTTCTTCAGTTGCACGAATATCTGTTTTTATGTCAGCGCCATAAAACGATATTGGTGAAACAAAAACTGCCATGAGATATAACCAGTCTCCGGATCCCTATATTTCTGTAGGCTAATAAGCAACTACTCATAGAAGTCACATCACTCCCAGTTCAGTTATAGAAGCATCACGACGAAAATTTCCAGTACTGCTACGGAACCGATCACCAGAACGTCAGATTACCTCTTACGGCGGAGACAGATCAGCATGACGGCAGATTACCTCTTACTGCGGAGACAGAACATGAATTTTCCAGTACTGCTCTTTTTTCGACTCATCCAATTCCTGTTCTAACATGTAGGTATGCATGTAATTGGAAAATGAAACCATGTAAATAGATTGTATTTACATTTTATTAATTGCCACTATATCACAGTAACTGGGCATTCTAACTGAAACTATAGAAAAATGCAAGGTCACAGCTCATAGGGAAGGGCAGAGAAAAATCAGGTATTGGTTTGTCGGAGAATATCATACCTCAAGCTCTTGATTGACTTTACGCTTCAGATTATTATTAGATTCTTGGGGCCTCAGACTGCATATTTGTTCCAAGTTATTTTCAGAACCGATCTCCAAGTTATTTTCAGCATTTACATCTTTCTGATTAGTTTGATAGGATACCAGCTCATCAAGTTCATCAGGTGTGGCTGGAATAGCTGGGGATAACTTCAGCAAGTTGTTTGCAGCACCGCCTGCACATTCCCTTTTTGAACTGTTAAGAGAGGAATCAGATCCATACAATTGTGAATCGTTGTTTGAATTTGGGAGCAAACTGTTAGGTTGCTTGACTACAAGGCACACGGTATTCTGCACAGCCGCAAACTTACTTTCCTGATAGTGTTCTGGATGTTCACTCTCAATAGATTCTGCATGAAAATTAGCAGACCCATTCAGAACTACTTTATTGGGTAAGCTAGGGGGCAACAGAGAAATACTTCCTGGCAACTGAAAGCCTCCAGCAAAATCTGGTTCTGGTTGTTTCATATAACTGCCTTTACATGGTACATGTTGAGCATCACTACAACTGGCAACGCCCAATAGGTTTCCCACTTGATCATCAGCTGAAACATCACGTGTTGGAATGGCACGGCTGTTCTTGTGCGTTGTTCTAGAAGTACGAGAAGGACGCACTCCACCAGGAAAAACAAACAAAGGGATTTCTTTTCTTTTTACATGGGACACTTCAATGTCCATTCCTTCTTTCTGGTGTTGATATGCACAAATAGTATTCTTAAACTCATTTACAATTCCCCTGATGTCGAATTGCTCAGTTTCTTGAGCTTGGGTAGTTTGTTTTCTCCAAAGACCCATGTTTGCCGTGTGGGTCAGCCATGCACTGCTGGCACTTAATATTACTGTTGTGTGGGTAGTTTATATTACAACTATATTACTGCTTGAAATTTCTGCTGCTTCACTGGCAATTGGTTGACTGCTTTAATTTCTGAACCTTAAGTTTTAATGTTTTCATTTCTTTGAGTGTAGACAAGGACATTAGCAGCTTCAACTAATGGACCTACTTTGTGATCAAGGACTGAAAGGGATATTGTTGATGAAGAATCTGAAAATAGTTACGAAGAATCTGAAAACCCCATCAATGATGAAGGTAAAAATCAAACAAACTTTAAATTTATTGTAGCAATTTTGAAAATAGTTATGAACAAAGGTTTTCATCTGTTTTATCTCCTGTAATAAGCAGAGCAACCTATTGGGCATCGAAATATAAAGAAAGGAAGGGGCCCAACCTTGAAGAAGAACATATATTCAAGTAGTGGTGGGCCTAAAATCTGCATTACCCTTACACTCTGTTTGGATCATTGGAATTGAATTCCATTCTAATAATAGTAATTTAGGCATATATCAATTAAGCTAATTCGGTTTTATGCAAAGTATATTTATATACTATTATTAGCAAGATGTCGGAGATATTTATGTGCTACATTTTTACTATAGGGGAGTGAGACGAAGAGTGTCATGTAAGTTATAGATTAGAAACCAATTCTAGTAATGCATAAAATCATTTCCCACCCTCCACCCCATGAATTTGAGATAAGCTTATATCTGAACTTTAGAAAGTGGTAGAATGTCAAATTCCAAACTAAATAAGTTACTTTATTAAGTGAATTCCAATTCCTCTAAAATGAAGGGATCCAAACGCCCCGTTAATGAATATGGACAACCAGTTGGTAGGAATAGCAAAGAGTTTGGTAATTTCATAGGAACTTTGGTGAGGAAAAATATCCCGATTTCTTGTGAGGATTGGAGACTAGTTGATTCTAAAAAGAAGTTTGAGTTATGGACAAAAGTGAAGGTAATCTGATTTAGTGAGATAATTTGTCTTAGCCTATTGTCAGTAGATTTGGTAAGCAACTACATTTTTGTTGAAACAGGAATACTATGAAGTGGATGAATCTGGTATAGATTATGTTATAACATCTGCAGCAAAAAAGTGGAGAGAGTTTAAGGCTGACCTAAAGAACAAATATTTTGATGAAACATTGAGTTTTAAAGAGCTTATTGCTAAAAGGGATGAAAGGGTGAAAGAATCTGATTGGGAGTGGTTGATAACTTATTGGATGTCTCCAGAAGCCGAGGTACATGCTGATTTGTTTTTAATTCTTTTTTGAGATATAAACCTTCATTAAACTTTTGTATGAATGAATTTGTATTATTGTTTTGATATATTCATATTTATGTTGTAGGTTCGTACAAACAGAGGTAAAGATAACCGTTCAAAGCTGACTATGTCGCATGCAGCTGGCAGCAAGAGTTATGCTCGTGTGGGGCATGATCTGGTAAAAAAAATTATAATCACTAGTTTGTATAATCATTAGTTTAATTTATGAAAAATCTAGTCAATTACTGGACCAAATCATTTTTTCTTGTAGGCCGAGCAACAAGGATGCCCTGCGAGAAGAGATGAAATATATGTTAGAACACACACGCGTAAGAACAAAGAAGGGGATATTGTGCCTCTACCTGGAGCAGAAATATTCATTGTGAGTACACATATTTTTGATTTATGTAATAATTAAAATATAAGTGAACAAAATATAACATATGTGTGGTTCTCTTAAGAATAAATTTGAAGAAGCTGTTGCTGAGAACCCAGAACTGAAGGACAGAAGTATTGAAGATGGTGATTTATATGCACATGTTTTTGGAGAGAAAGAACCAAGAGGTCGTATTCGTGGTTTAGGCTTAGGACCAACTCCACAAGATGTAGGCACCCCTGGAACTCAAATGAAAATATCAACAAAGCTTCAAATGGCATTGCAAGCTCGCAGTCAGTCTGAGCAAGAGGTTAGAGCCTTAAGACAGGATATGAATCAAATGAAAGAAAAAATGGATCAAATTTATCAAATGATGGTAGCAGCTCAAGGGGTGAAACATATAGAAAGCCCATCACAACATGGGTCTAATTCTCGACAGGTGATATTATTGCATATGTTTAGTCAGTGTAGCATATATCGAGCAGTTCGTCGTCTATCCAGTTAATTAGGCTTGTTTCCTTTATCTCACTAACACTTTTTCCCTTTTTATTTTTCTTTAGAACTCAAGGGTGCATCGGTCAGATGAGGTGGCACATGATTACAATGGTAATGACTATTCACAAAATATGGTTGAAGATGACTTGCAGCTTACTAGGCGAGTTGCAAGCACTGTGGGCCGTCTAAGGAATCGTGAAGATGTTAATACTATTGAAGAGCAAAGTCGCAGGCGAGACATTGCAACAATGGTACCTCAAAGAAATCATGGATCACCTCAAAATGCAGACGATAATCATGTAATTTCTCATCTATCACTTCATATTTATGTTTCATTATTTGAACTTTGTAATCTTACTTCGCATATTTTATTAGGTTGGTAAAGAAGTCATATTATATTCTATGATGAGATCAGAAGTTCCTGTTGCAATAGCAAATATTATCTCAACCGACCCAACCACTAAGGTTGCAGATGTTCCTCTTGGAAGGGAGTTCACACAGGTTTTTGTGACTCGCGTGCTAAAAAGGGAAAGTACTCTGCCACGACCATACTTAGGTGTAGAGAGTATGGGGGATGCTCTTTTTATGCCCGTTGCATGGCCAACAAACAAGGTAATACTACAACCTTGTCTTTTATTTATGTAATTAATGGAACAAGCTTCTAATGTTCTGTATTTTTAGATGGGTCGTAACAAAAAATCAACATTGACCCAAGGTTCTACAGCAGCAGGTAAATATTTTGCCTGAAAATATGTTTTCAGACTCTTGAAATAAATGACATCTAGTTTACAGCAGGTAAATAGTTTACAGTTTACAACAGTTTATAGTTTAGTTTTAGACTTGTATTTTTAGTGAAGGCGCCACTTGTTCTGATCATGATGTGTGGACCAAACTTTCTTGAACCTAGCTAAGTCTGAGTCGGTCTTCCAAGTCTAGGGGTCCCTGTCCAAGTTACTTAGTCTAAGTCAAAATTCCTGACTAGACCTAGGTTAGTTTCTTTTCTCTTTTGATGTTTTAGCCATCAATGTACCTGGGATCGATCTTCCTGGGTTCTGTTTTATCAAGATTTGTAACTCAGCCTTTCTGGGATCGATCTTTTATAAAGTTTAGAACCTGAAAGGATTGGTGCTCTGTGTTGGATGGAAGCATAGCCTCAGAACGAGGATAATGGATAGAGAAGCAATTCTGCCTTTTCTTGAAGCAATTCTGCCACATTTTATCGTTTGCAATTTTTCCTGGCTCGTGCAAAGAATCTGTTTGCTTTCTCATACACTATTAGCTTGGTCTGACAGAGTCTTTGTTTGGTTCCCATGATTAGCCAGGCTCCATAGAAACGGATAAAAATCTTGCACAGAGAGAGCCAGGCCCGGCGGAACCTTCTGCACATGCTCAGCCTGGCTAGATTGCAAGCGCACCCAACCAATCATGTGTGCTTGCACCGGAGGAGCCTAGTTGTCCTTGGTAGAGCTAACCAAACGTGACCATAGCTTTTGGATGAGTGAGCTAACTGTTGTCATACTCCTGAGAGGCACTATGGATATATAGGAGTCCCTTGTCATGCTCGAGAGGCTCTATTTTGTGTCCAACAGTGAGCATGCCTCTTTGAGCCAGATCGCCCCCTGTCTCCTAAATCTCCAAAATCTGTGAAGGAGGTTCTAGAGGAACCAAAGTTCAGTTGTTTGATCTGGATGAGTGAGCTAACTGAAGCGACGTCCTTACATAGAAATTCCATAAATTTCATACTATTTTCTTCATATATTTTCTGCATGTTGTCTCTATTTTGTCTTAATTTTTAATGTCCAAAATCATGTCTCTGTTTCGTCCTCTTCAGAACTTATTCGGTGATCTTTTGGCCCATCTTTGCTAATACCGACAAAGTACGCATGCTGCTGTTGTCAAGGTATAGTATGTGCCCACCTAAATGACTTCATAAAACTCCATCTCTCCATTCACCCTTAAACTTTCATCTACAGTACACCGAACTTTGTGATTTACTTTCATCTTTGTGATTTACTTTCATCAGCAGTATCACCGAACTTTTTGCAATTTAATACTTTTATTAAACTTTCAGTTTGAATTTTCTGTTGGGAGTTTATAGACCTGTTGATATATATAGTAATTGACAAATTATCAAAATATGTTGTAGGGGAGGTGCTGGAAGCAAGAATCAGGAGTGCTTGAAGTGTGTGAAAAAATGATTACTCGCATGCATGCTAGCCGTTGGCGTAACATCACCCAGCTGATGGAGAGATTTTAGTTTGCTAGCTAGATGAATCTTTAGACTTTGCATGTAATTGAACTAAAAGTTCACTTCAGAATTTAGAAAACATGTATGCTATTATTTTTTTATATGTTCATCACACAAATTGTGTTCCAAAACAATATACGTTGGAATTCAAAGTCAAAGATTGATCTCTTTCTTTTAGTAATCAAATTTCATATATTCTGAGCTACATATTATTCAAGTTGTCCTGTAATGTGCATCATAATCAGTTTATTATTTGTTTCATTAGCATATATTTGTATACATGTTTTTTTAACTACACATAGGGAATAACATCAGACAATCTGTCGGTATATAGACTTTTAGCGACACACAAGCTGTTGGCACAAGAATCTGTCGGCACAAATTTGTCTGACGGCAAAGATATGTCTGTCGGCAATAATATAAACATTACCGACAGGTAGTCTGTCGGTAGATCTATGTCTATCGGTAAAGGTTATTGCCGACAGGACTTTCGCCGACACACATTATCTGTCGGCAAAGCTCTTTACCGACGGACTGTGCGACGCTATTTAGGTGTTTTACCATCAGACAATCCGTTGGCAATATTGTGTCATGGTGTAGTGGGTGTCTACCAATAGTGTAGGCATCCCTAGCAAATTTCAATTTATTAAAATCATCCTTGCAAGTCTTAAGTTGAGCATTAAGACTTGCCACCTCATTACTTAATTTTGTAATAGAAACAATATGTTCATTACAGGCATCAACATCAAAATCCTTACATCTATTACAAATAAAAACATGCTCTACACATGAACTAGATTTATTAACTTCCTCTAGCTTAGCATTTAAATCATCATTTAAACTCCTTAAACTGGAGACAGATTCATGGCAAGCAGACTACTCAGAGGATAACATTTCATTTCTTTTAATTTCTAGAGCAAAAGATTTTTGAACACTAACAAATTTGTCATGCTCTTCATACAAAATATCTTCTTGCTTCTCTAGCAATCTATCTTTTTCATTTAGAGCATCAATTAATTCATTAATTTTCTGTGCTTTGGCTCTATCTAAACCTTTAAATAAACTAGAGTAATCTACTTCATCATCGGAGGATTCTTCATCACTAGATGAAGTGTACTTGGGTGTGTCCCGAACACGTAACTTTTTCTCCTTGGCCATAAGACATGTATGGCATTCGTTGGGGAAGTGTGAGGATTTGTTGAAAGCCGAGGCAGCAAGTCCTTCGTCGTCAGATTCGGAAGAGGAGCAGTCCGAGTCCCATTCCTTTCCAAGGTGCGCCTCGCCCTTTTCCTTCCTGTAGTTCTTCTTTTTCTCCTTCTTCCCGCTCTTTTCTTGTCCCTGGTCACTATCATTATCGGGACAATTTGCAATAAAGTGACCAGTCTTACCACATTTAAAGCAGGAGCGCTTTCCCCTTGTTTTGTTCTTGTTGGGATACTCCTTGTGTCCTTTAAGCGCGGTCTTGAAGCGCTTGATGGTAAGCGCCATCTCATCCTCGTTGAGGCCAGCAGCCTCCACTTGTGCCACCTTGCTTGGAAGCGCCTCCCTGCTACTGGTTGCCTTAAGAGCAACAGGTTGAGGCTCGTAGACGGGTAGTGGACCGTTTAGAGCATCGTCCACGTATCGAGCCTCCTTCACCATCATGTGGCCGCTCACAAATTTTCCAAGAATCTCCTCGGGCGTCATCTTGGTGTACCTGGGATTCTCGCGAATAAGATTGACAAGATGGGGGTCAATTACAGTAAATGACCTTAGCATGAGTCGGACAACATCATGATCCATCCATCTTGTGCTTCCATAGCTTCGGATCTTGTTGATGAGAGTCTTGAGCCTGTTGTAGGTCTGAGTTGGCTCCTCTCCTCTTATCATCGCGAACCTCCCCAGCTCGCCTTCCACTAGTTCCATCTTTGTGATCATAGTGGAATTGTTTCTCTCATGAGATATCTTGAGGGTGTCCCAGATTTGCTTGGCGTTGTACAAGCCACTCACCTTATTATATTCATCCCTGCAAAGTGAAGCTAGCAAGACAGTAGTAGCTTGTGCATTTTTATGAATTTTCTCATTTATGAAAACGGAATTGTCAGTACTATCAAAATGCATTCCATTTTCAACAATCTCCCAGATACTAGGATGGAGAGAAAATAGATGACGACGCATTTTATGACTCCAGAAAGAGTAGTCTTCTCCATCAAAGTGTGGAGGTTTCCCAAGTGGAATTGATAGCAAATGAGCATTTGAATTAAATGGAATACGAGAATAATCAAATGAATAGTTTTGATTAACCGTTTTCTTCTTTGAAGAGTCGTCGTCGTCGTGTGGTGAAGAAGACGAGGCATCGCTGTCGTAGTAGATGATCTTCTTGATGCGCCTCTTCTTATTCCCGTCTCTCTTCTTTTGACTCGAGCCTGAGTCAGTAGGCTTGTCGTCTCTTGGCTCGTTGAAGACGAACTCCTTTTCCTTGTCGTTGACCACCATCCCCTTTTCTTTAGGATCCATCTCTTCGGGCGGTTAGTCCCTTAGATGAAGAGTACGACTCTGATACCAATTGAGAGCACCTAGAGGGAGGTGAATAGGTGATCCTGTAAAATTAAACACTAAATAGCCACAAAACTTAGTTATAAAAGTGTTAGTGCGACTAAGTAGTTGAGAGGCGAGTTCTTGTGGACAAAACAATCACAGAGAAAGCAATCACAAGAGACATGCGATTTTATCATGTGGTTCAGCCAAGTAACACTTGCCTACTTCCACGTTGTGGCGTCTCAATGGACGAGGGTTGCACTCAATCCCTTTCAAGTGATCTGATGATCAACTTGAATACCACGGTCTTTTCCTTTAGAGTATCTTCCCGATTGCGAGAAATCTCCACAGTTTGGAGTCTCTCGCCCTTACAATTAGGATCACAAAGAAAGCACGGAGTAAGGTTGGGAGAAGCAACGCACACACAACACAATTTCGCAGCACACACACGCACACAAGCCAAGACTTGAGCTCGAAAAGTAGCACAGAGAGTTCACAACTCGAACGGAGCTCAAATCACTAACACAATCGATCAAATGCGAGGAGGCGGAGTGTGGGAGTCTTAGAATGCTTAGTGGATGCTTAGATGTTTCCTCCATGCGCCTAGAGGTCCCTTTTATAGCCCCAAGACACCTAAGAGCCGTTGGAGATCAACATGGAATGCTATCCTTGCCTTCTGTCGAGTGGCGCACCGGACAGGTCCTGTAGATTGTTCGGTGCGCGATCTCCTTCCAAATTTGGCATATCCGACCGTTGCTCCTCTGGGCTGATTGGCGCACCGGACACAGTCCGGTGCACACCGGACAGTCCGGTGCACCAGCTGACCGTTGGAGCAGTCCACGTGTCGCGCGAAGATTGCGTGGCCGACCGTTGCTCAGGCGACCGTTGGCTCACCGGACAGTCCGGTGCACCACCGGACAGTCCGGTGAATTATAGTCGTACGCCGCCGTCGAAACCCGAGAGCGGCGAGTTCACAGTGGACCAGCCTGGCGCACCGGACACTGTCCGGTGCACCATTGGACATTGTCCGGTGCACCACCGGACAGTCCCGTGTGCGAGACCGAGCACAAGATTGGCTGCACAGAGCCAAGCTTTTCCCTTTTTCTCCCTTCTTTTTTAGTCACTGTTTCTAGCACTTGGATAACCATGTTAGTACATAAAACAATTCACCAAGTCTAGAAACATACCTTTTGCCTTGATTTTCACTTCTCACTTTATTTGGCACATAAGAACTTAATTAAACGTGTTGGGCACTTAATCACCAAAACATTATAGAAATGGCCCAAAGGCACATTTCCCTTTCAGTAATAAGTTAATACCTCTTTATATCATTATTTGAACACTGTGCAATGATGTTCATTTATGTAATCGTTGTGTACGTCAGTTCTAATTCTAGCACGTACATGGTTCACATCCAATTTGTCTTCTAAAAACGAATGTGACATAATGTCATATGTATGTGATAATGCTTTTTGTTGGGGTCCTTCGTCTTTCAAAGGTCCTCAAAAACACATTTAACCATTGGTTGTTAGCACATCCTTAAGTGTTGCAGGAGCTTTGGTATTGAATACCTTCGGAGCAGGACATGGAGGAAGACGAAGATGTTAGCTTCGTCATAACAACACAAGGAAACGAAGGCAGAAGTGGAACAAGGCCGGGATATGGTGTTTTCAAGACTCTGTATCCAAAGCAAAAAAGACAGAAAGACGATACTGCCCTTACATAATTTGTAAACTATGTGAACAAGTTTTATGGACATGTTTGTAACTTTACACGAAATTGTACCACCACACTATAAATAGATAAATAGTGCCCTGCATGAGGCGCCTCTTGGGAACAATGAGGAACAACTCTGTATAATCCTTTTTCTTCTAAGTACCTTCGGGTTTTCTCCTCATCAAAAAGCCGAAGGTACTATTGTAAATTCGTTTCATATAAAGAAAGAAATCCCAAGTTGTTTGAGATAAGTAATCTTATCTAGCTTTGTTATAGCCATGTGTGTAATCTTTATCTTTATCCTCTGACAATCCTATATATTATATATAATAACCTTCGTACTTTACTTGGATGTCCCGAAGGACAAACTCTTTAAGTACGAAGGATAACATCTTTTTTAATAATGTGTTGCCTTGTTTTTTATTGTGTACAACAATTAAAAACGAGTGACCAACATTTTCATGTCAGGGTATGGGGACCCATTGGAGACTCGATATCCAAATGAGGATGAGTATATGATGAATCCTATACCTATGATGAGTATAAGTATGAGAATCAGGATGAGTATAACTTCATCAGAATAGGTGCGGGGGCGTCCTTGTGGGCGTGCCTACCGTGCGATCGCACAAGGCCTCCAAAACCATAGGGCCCCAAAATTTATAACAATCTTTATATACAATATAAATAAAAAATATTATTTTATATAAAATATTTACCAACATATAGCATAGAATCGTAAAAGCGTTGAAATCGATCTGTTCTTATTGTTATTCAAACTATTTACCTCCAGCACATTGTAGTCATTAGATAAAAAAGATTGAGATCTTATTGTCACTATCTTAAGAAGACACAGTTAAAAGAGGTAGACAATATGTCAATATGCTTTGATGCAAGTGACCAATTCGTGACGTTGAGTTTCCTCTAAGATTTTGTTTAAAAAATTGCTATGTTGACATTCTAAATTTTATAAAGCAGAGGAGCAAAACTGAGTAAAATCGCATTTAATGATAAAAATGTGGAAAGTGACAAAACTAAGAATACAATTTTAAATAGTCCAATATTTTTTTACTATCTTTTGCACAGGGCCTCTCAACTTGGGAGGACGCTTCTGGGTGTGGGTTTACAAATTCGATGAAAAATTCCCCATTGACAAACGATAGGAGGATATTTTTCTCCCAGCACAAAATAGCATAGCCATAAGGCAACAAGGCATGGCAAAGGATCGTATCATCGTCATCCGAGACCCATTGCTTTCTCTCTCTCTCCTCGTGCTTTCATTACTGGGGTGGGGGTGGAGTGGACCAGTGGAGTGGAGAAATGACAAATCCAGGCCCGCAGGCAGCCCCACCCACCAAATCGGCCGAGCAGGGTGCCCAAATCAGGAAGGATTTTAAGGTTAACCGGCTGCCACCGCCCACCGCCGGTGACCCCAGTCTCTCTTCTATCTATATATTACCCGCCTCCTTTTCTCCTCTCTCTCCGCCCCACCCTCCTTCCTCAGCTCCGTTGCGCACCGCCACCGCCGGCCGGCCAGCCGCCGGAGCACCGAAAGACCCCCGTTCTTTCCTGTAAAAAAAAACCCGCCGCCTTTAGCTAGCTAACCGGTCGTCCTCTTCACCCCCTAGCTTTGCTAGCTCTAGCTAGGAACGAAAGAAATTAAAGGATAACTGAGATTGCTGATTGGTGGTCCGGGTACGGTGTTCTTGAGTCGTGAAGCGACAGTACAGTGGCTAGGGTCGTGCCGCCCCTGCAGTCTCCGGGGTTGCGTGCAGGATGGTCGTCAGGGATCGGGAGTGAGGAGGCATCAGCTCTCGCGGTCGTGGAGCCTAAATGTACCGCAACAACGACTCGGCACTCTCCTGCTTCTACCTCTTCCTCCTCTGGTTCTTCTTCTTGAAGTAGACACCACCAGTTCGCCAGGTAGTTAGCAGCCCAGTTGCGACTGGGGATCGGTGGCGGGCTGCCGCTTGCGAGTTGTAAGCTTGGAGGGGAGGGGAGCAGGAGCAGGAGATGCAGTTGGATCTGAACGTGGCCGAGGCGCCGCCGCCGGTGGAGATGGAGGCGAGCGACTCGGGGTCGTCGGTGCTGAACGCGTCGGAAGCGGCGTCGGCGGGCGGCGCGCCCGCGCCGGCGGAGGAGGGATCTAGCTCAACGCCGGCCGTGCTGGAGTTCAGCATCCTCATCCGG";
    int _open_gap_penalty = -2;
    int _extend_gap_penalty = -1;

    int matchingScore = 3;
    int mismatchingPenalty = -2;

    calculateK( seq1, seq2, _open_gap_penalty, _extend_gap_penalty, matchingScore, mismatchingPenalty);
    ASSERT_EQ(0, 0);
}
