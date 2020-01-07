//
// Created by bs674 on 9/25/19.
//

#include "mapCNSToGenome.h"


// the second algorithm, continue from the output of first result
// using smith-waterman approach to find seeds and using x-drop to extend the seeds
// extend the x-drop result using a 2-piece gap cost approach
PairedSimilarFragment mapCNSToGenome (int8_t * seq1, int8_t * seq1_rev_com, int8_t * seq2, int8_t * seq2_rev_com, int32_t & length1, int32_t & length2,
                                                    const int32_t & _open_gap_penalty1, const int32_t & _extend_gap_penalty1,
                                                    const int32_t & _open_gap_penalty2, const int32_t & _extend_gap_penalty2,
                                                    const int32_t & matchingScore,
                                                    const int32_t & mismatchingPenalty, const Scorei & m){

    int32_t maxScore;
    int32_t endPosition1;
    std::cout << "mapCNSToGenome line 19, length1:" << length1 << " length2:"<< length2 << std::endl;
    mapCnsToGenome(seq1, seq2, length1, length2, _open_gap_penalty1, _extend_gap_penalty1,
                                         _open_gap_penalty2, _extend_gap_penalty2,
                                         maxScore, endPosition1, m);

    int32_t length11 = 1 + endPosition1;
    int32_t  endPosition1_rc;
    std::cout << "mapCNSToGenome line 23, length11:" << length11 << " length2:"<< length2 << std::endl;
    mapCnsToGenome(seq1_rev_com+(length1 - 1 - endPosition1), seq2_rev_com, length11, length2, _open_gap_penalty1, _extend_gap_penalty1,
                                                 _open_gap_penalty2, _extend_gap_penalty2,
                                                 maxScore, endPosition1_rc, m);
    std::cout << "mapCNSToGenome line 30" << std::endl;
    int32_t start1 = endPosition1 - endPosition1_rc; //0 based
    int32_t length21 = endPosition1_rc + 1;

    std::vector<uint32_t> cigar;

    std::cout << "length21:" << length21 << " length2:" << length2 << std::endl;
//    Matrix T (length21 + 1, length2 + 1);
//    std::vector<uint32_t> cigar = mapCnsToGenome(seq1+start1, seq2, length21, length2, _open_gap_penalty1, _extend_gap_penalty1,
//                                                 _open_gap_penalty2, _extend_gap_penalty2,
//                                                 maxScore, endPosition1, m, T);


    double p=0;
    PairedSimilarFragment pairedSimilarFragment(start1, endPosition1 + 1, 1, length2, maxScore, cigar, p, p);
    return pairedSimilarFragment;
}
