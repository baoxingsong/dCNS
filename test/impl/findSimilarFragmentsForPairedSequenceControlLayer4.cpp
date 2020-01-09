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
#include "../../controlLayer.h"

#include <limits>
#include <iomanip>
#include <cmath>
#include <algorithm>

// most general one
TEST(controlLayer, pairCnsXExtend){ // just to make sure that every line has been analysed

    int32_t seed_window_size = 38;
    int32_t mini_cns_score = 40;
    int32_t step_size = 8;
    int32_t matrix_boundary_distance = 0;

    int32_t matchingScore = 2;
    int32_t mismatchingPenalty = -3;
    int32_t openGapPenalty1 = -4;
    int32_t extendGapPenalty1 = -2;

    double lambda = 0.382291;
    double kValue = 0.006662;

    bool onlySyntenic = false;
    double pvalues = 0.1;

    std::string input = "/home/bs674/105538";
    std::string reference = "reference";
    std::string output = "/home/bs674/105538_many_to_many_seed_xtend.sam";
    std::string output2 = "/home/bs674/105538_many_to_many_seed_xtend";

    int32_t w = 10;
    int32_t xDrop = 20;

    pairCnsXExtend(input, reference, output, matchingScore, mismatchingPenalty, openGapPenalty1,
                   extendGapPenalty1, seed_window_size, mini_cns_score,
                   step_size, matrix_boundary_distance, onlySyntenic, lambda, kValue, w, xDrop, pvalues);
}



//general + two gap penalties

TEST(controlLayer, pairCns2Gap2){

    int32_t seed_window_size = 38;
    int32_t mini_cns_score = 40;
    int32_t step_size = 8;
    int32_t matrix_boundary_distance = 0;

    int32_t matchingScore = 2;
    int32_t mismatchingPenalty = -3;
    int32_t openGapPenalty1 = -4;
    int32_t extendGapPenalty1 = -2;

    int32_t openGapPenalty2 = -45;
    int32_t extendGapPenalty2 = 0;

    int32_t zDrop = 50; // this one should be slightly larger than openGapPenalty2
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

    int32_t w = 10;
    int32_t xDrop = 20;

    pairCns2Gaps(input, reference, output, matchingScore, mismatchingPenalty, openGapPenalty1,
                 extendGapPenalty1, openGapPenalty2, extendGapPenalty2, seed_window_size, mini_cns_score,
                 step_size, matrix_boundary_distance, onlySyntenic, lambda, kValue, w, bandwidth, xDrop, zDrop, pvalues);
}





//weighted + 1 gap penalties
TEST(controlLayer, weighted1Gap){

    int32_t seed_window_size = 38;
    int32_t mini_cns_score = 40;
    int32_t step_size = 8;
    int32_t matrix_boundary_distance = 0;

    int32_t matchingScore = 2;
    int32_t mismatchingPenalty = -3;
    int32_t openGapPenalty1 = -4;
    int32_t extendGapPenalty1 = -2;

    int32_t openGapPenalty2 = -45;
    int32_t extendGapPenalty2 = 0;

    int32_t zDrop = 50; // this one should be slightly larger than openGapPenalty2
    int32_t bandwidth = 20000; // using very large value, 12000 there are some strange result, should be checked very carefully

    double lambda = 0.382291;
    double kValue = 0.006662;

    double pvalues = 0.1;

    std::string scoreFoler ="/Users/bs674/scoreMatrix";

    std::string refFasta = "/Users/bs674/Zea_mays.B73_RefGen_v4.dna.toplevel.fa";

    std::string input = "/Users/bs674/SORBI_3009G024600_5_prime";
    std::string reference = "SORBI_3009G024600_maize_V4_+_8_135910122_136010121";
    std::string query = "SORBI_3009G024600_sorghum";
    std::string output = "/Users/bs674/SORBI_3009G024600_5_prime_weighted_2gap.sam";
    std::string output2 = "/Users/bs674/SORBI_3009G024600_5_prime_weighted_2gap";

    std::string gffFile = "/Users/bs674/Zea_mays.B73_RefGen_v4.42.gff3";

    int32_t w = 10;
    int32_t xDrop = 20;
    bool onlySyntenic = false;

    weighted1Gap(input, reference, output, matchingScore, mismatchingPenalty, openGapPenalty1,
                 extendGapPenalty1, seed_window_size, mini_cns_score,
                 step_size, matrix_boundary_distance, onlySyntenic, lambda, kValue, bandwidth, w, xDrop, zDrop, scoreFoler,
                 refFasta, gffFile, pvalues);
}




//weighted + 2 gap penalty
TEST(controlLayer, weighted2Gaps){

    int32_t seed_window_size = 38;
    int32_t mini_cns_score = 40;
    int32_t step_size = 8;
    int32_t matrix_boundary_distance = 0;

    int32_t matchingScore = 2;
    int32_t mismatchingPenalty = -3;
    int32_t openGapPenalty1 = -4;
    int32_t extendGapPenalty1 = -2;

    int32_t openGapPenalty2 = -45;
    int32_t extendGapPenalty2 = 0;

    int32_t zDrop = 50; // this one should be slightly larger than openGapPenalty2
    int32_t bandwidth = 20000; // using very large value, 12000 there are some strange result, should be checked very carefully

    double lambda = 0.382291;
    double kValue = 0.006662;

    double pvalues = 0.1;




    std::string input = "/home/bs674/SORBI_3009G024600_5_prime";
    std::string reference = "SORBI_3009G024600_maize_V4_+_8_135910122_136010121";
    std::string query = "SORBI_3009G024600_sorghum";
    std::string output = "/home/bs674/SORBI_3009G024600_5_prime_weighted_1gap.sam";
    std::string output2 = "/home/bs674/SORBI_3009G024600_5_prime_weighted_1gap";

    std::string gffFile = "/home/bs674/Zea_mays.B73_RefGen_v4.42.gff3";


    int32_t w = 10;
    int32_t xDrop = 20;
    std::string scoreFoler ="/home/bs674/scoreMatrix";
    std::string refFasta = "/home/bs674/Zea_mays.B73_RefGen_v4.dna.toplevel.fa";
    bool onlySyntenic = false;
    weighted2Gaps(input, reference, output, matchingScore, mismatchingPenalty, openGapPenalty1,
                  extendGapPenalty1,  seed_window_size, mini_cns_score,
                  step_size, matrix_boundary_distance, onlySyntenic, lambda, kValue, w, bandwidth, xDrop, zDrop, scoreFoler,
                  refFasta, gffFile, pvalues);
}



std::string input = "/media/bs674/panAndAssemblyfi/sorghum/and_cns_sorghum_maize_v2_corrected_bug_kmer/16154";
std::string reference = "reference";


// cut and link and 1 gap
TEST(controlLayer, cutAndLinkdApproach){ // just to make sure that every line has been analysed

    int32_t seed_window_size = 38;
    int32_t mini_cns_score = 40;
    int32_t step_size = 8;
    int32_t matrix_boundary_distance = 0;

    int32_t matchingScore = 2;
    int32_t mismatchingPenalty = -3;
    int32_t openGapPenalty1 = -4;
    int32_t extendGapPenalty1 = -2;

    int32_t openGapPenalty2 = -45;
    int32_t extendGapPenalty2 = 0;

    int32_t zDrop = 50; // this one should be slightly larger than openGapPenalty2
    int32_t bandwidth = 20000; // using very large value, 12000 there are some strange result, should be checked very carefully

    double lambda = 0.382291;
    double kValue = 0.006662;

    double pvalues = 0.1;

    std::string output = "/home/bs674/16154_cutAndLinkdApproach_1gaps.sam";
    std::string output2 = "/home/bs674/16154_cutAndLinkdApproach_1gaps.o";


    int32_t w = 10;
    int32_t xDrop = 20;
    std::string input = "/Users/bs674/Zm00001d008450_3_prime";
    std::string reference = "reference@Zm00001d008450";
    std::string refFasta = "/Users/bs674/B73_v4_k20_9.fa";
    std::string  queryFasta = "/Users/bs674/sorghum_v4_k20_5.fa";
    bool onlySyntenic = false;

    cut1Gap(input, reference, output, matchingScore, mismatchingPenalty, openGapPenalty1,
            extendGapPenalty1, seed_window_size, mini_cns_score,
            step_size, matrix_boundary_distance, onlySyntenic, lambda, kValue, w, xDrop, refFasta, queryFasta, pvalues);

}




// cut and link and 2 gaps
TEST(controlLayer, cutAndLinkdApproach2gap){ // just to make sure that every line has been analysed

    int32_t seed_window_size = 38;
    int32_t mini_cns_score = 40;
    int32_t step_size = 8;
    int32_t matrix_boundary_distance = 0;

    int32_t matchingScore = 2;
    int32_t mismatchingPenalty = -3;
    int32_t openGapPenalty1 = -4;
    int32_t extendGapPenalty1 = -2;

    int32_t openGapPenalty2 = -45;
    int32_t extendGapPenalty2 = 0;
    int32_t zDrop = 50; // this one should be slightly larger than openGapPenalty2
    int32_t bandwidth = 20000; // using very large value, 12000 there are some strange result, should be checked very carefully

    double lambda = 0.382291;
    double kValue = 0.006662;
    double pvalues = 0.1;

    std::string output = "/Users/bs674/Zm00001d008450_3_prime_cutAndLinkdApproach_2gaps.sam";
    std::string output2 = "/Users/bs674/Zm00001d008450_3_prime_cutAndLinkdApproach_2gaps.o";

    int32_t w = 10;
    int32_t xDrop = 20;

    std::string input = "/Users/bs674/Zm00001d008450_3_prime";
    std::string reference = "reference@Zm00001d008450";
    std::string refFasta = "/Users/bs674/B73_v4_k20_9.fa";
    std::string  queryFasta = "/Users/bs674/sorghum_v4_k20_5.fa";
    bool onlySyntenic = false;


    cut2Gaps(input, reference, output, matchingScore, mismatchingPenalty, openGapPenalty1,
             extendGapPenalty1, openGapPenalty2, extendGapPenalty2, seed_window_size, mini_cns_score,
             step_size, matrix_boundary_distance, onlySyntenic, lambda, kValue, w, bandwidth,
             xDrop, zDrop, refFasta, queryFasta, pvalues);
}





// pairwise sequence alignment to multiple sequence alignment
TEST(controlLayer, multCns){ // just to make sure that every line has been analysed

    int32_t minimumNumberOfSpecies = 2;
    uint32_t mini_cns_size = 7;
    double outputWithMinimumLengthPercentage = 0.8;

    bool onlySyntenic = false;


    std::string referenceGenomeFile = "/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.dna.toplevel.fa";
    std::string output = "/media/bs674/1_8t/AndCns/msa2";

    std::map<std::string, std::string> samFiles;
    samFiles["sorghum"] = "/media/bs674/1_8t/AndCns/sorghum/and_cns_sorghum_maize_V2_smaller_kmer_frequency/5_uniq_2.sam";
    samFiles["setaria"] = "/media/bs674/1_8t/AndCns/Setaria_italica/result/5_uniq_2.sam";
    samFiles["miscanthus"] = "/media/bs674/1_8t/AndCns/Miscanthus_sinensis/result/5_uniq_2.sam";
    samFiles["sugarcane"] = "/media/bs674/1_8t/AndCns/sugarcane_tareploid/result/5_uniq_2.sam";
    samFiles["1013"] = "/media/bs674/1_8t/AndCns/1013Chrysopogonserrulatus/result/5_uniq_2.sam";
    samFiles["1025"] = "/media/bs674/1_8t/AndCns/A1025_08May2019/result/5_uniq_2.sam";
    std::cout << "test line 312 " << std::endl;
    getCnsForMultipleSpecies (  onlySyntenic, output,
            //std::map<std::string, std::string> & sequences, /*species, fastaFile*/
                                samFiles,/*species, sameFile*/
                                referenceGenomeFile, minimumNumberOfSpecies, mini_cns_size,
                                outputWithMinimumLengthPercentage);

}





// pairwise sequence alignment to multiple sequence alignment
TEST(controlLayer, mapCNSToGenome){ // just to make sure that every line has been analysed

    std::string input = "/media/bs674/1_8t/AndCns/mapCnsToReads/sorghum/mapreadsToCns/cns_seq.fa";
    std::string reference = "/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.dna.toplevel.fa";
    std::string output = "/home/bs674/mapCNStoRef.sam";

    int32_t matchingScore = 2;
    int32_t mismatchingPenalty = -3;
    int32_t openGapPenalty1 = -4;
    int32_t extendGapPenalty1 = -2;

    int32_t openGapPenalty2 = -45;
    int32_t extendGapPenalty2 = 0;


    mapCNSToGenome(input, reference, output, matchingScore, mismatchingPenalty, openGapPenalty1,
                   extendGapPenalty1, openGapPenalty2, extendGapPenalty2);

}
