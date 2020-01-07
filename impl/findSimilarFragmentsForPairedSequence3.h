//
// Created by Baoxing song on 2019-01-02.
//

#ifndef SONG_CNS_FINDSIMILARFRAGMENTSFORPAIREDSEQUENCE_H
#define SONG_CNS_FINDSIMILARFRAGMENTSFORPAIREDSEQUENCE_H

#include <algorithm>
#include "sequenceAlignment.h"
#include "../model/model.h"
#include <vector>
#include "./fasta.h"
#include <iostream>
#include <sstream>
#include <cassert>
#include <math.h>
#include <iostream>
#include "calculateLambda.h"
#include <seqan/align.h>
#include <seqan/align_extend.h>
#include <seqan/sequence.h>
#include <limits>




// do not use this function anymore, since it could not deal with copy number very well

std::vector<PairedSimilarFragment> findSimilarFragmentsForPairedSequence ( int8_t * seq1, int8_t * seq1_rev_com,
                                                                           int8_t * seq2, int8_t * seq2_rev_com, int32_t & length1, int32_t & length2, int32_t & windowsSize,
                                                                           const int32_t & mini_cns_size, const int32_t & mini_cns_score, const int & matrix_boundary_distance,
                                                                           const int & _open_gap_penalty, const int & _extend_gap_penalty, const int & matchingScore,
                                                                           const int & mismatchingPenalty, const Scorei & m, const int32_t & step_size,
                                                                           std::string & seq1_string, std::string & seq2_string, const double & pvalues,
                                                                           const double & lambda, const double & kValue, const int & lDiag, const int & uDiag, const int & xDrop);


std::vector<PairedSimilarFragment> findSimilarFragmentsForPairedSequence ( int8_t * seq1, int8_t * seq1_rev_com,
                                                                           int8_t * seq2, int8_t * seq2_rev_com, int32_t & length1, int32_t & length2, int32_t & windowsSize,
                                                                           const int32_t & mini_cns_score, const int & matrix_boundary_distance,
                                                                           const int & _open_gap_penalty1, const int & _extend_gap_penalty1,const int & matchingScore,
                                                                           const int & mismatchingPenalty, const Scorei & m, const int32_t & step_size,
                                                                           std::string & seq1_string, std::string & seq2_string, const double & pvalues,
                                                                           const double & lambda, const double & kValue, const int & lDiag, const int & uDiag, const int & xDrop);

std::vector<PairedSimilarFragment> findSimilarFragmentsForPairedSequence ( int8_t * seq1, int8_t * seq1_rev_com,
                                                                           int8_t * seq2, int8_t * seq2_rev_com, int32_t & length1, int32_t & length2, int32_t & windowsSize,
                                                                           const int32_t & mini_cns_score, const int & matrix_boundary_distance,
                                                                           const int & _open_gap_penalty1, const int & _extend_gap_penalty1,
                                                                           const int & _open_gap_penalty2, const int & _extend_gap_penalty2, const int & matchingScore,
                                                                           const int & mismatchingPenalty, const Scorei & m, const int32_t & step_size,
                                                                           std::string & seq1_string, std::string & seq2_string, const double & pvalues,
                                                                           const double & lambda, const double & kValue, const int32_t & zDrop, const int32_t & bandwidth,
                                                                           const int & lDiag, const int & uDiag, const int & xDrop);



std::vector<PairedSimilarFragment> findSimilarFragmentsForPairedSequence_wighted ( int8_t * seq1, int8_t * seq1_rev_com,
                                                                                   int8_t * seq2, int8_t * seq2_rev_com, int32_t & length1, int32_t & length2, int32_t & windowsSize,
                                                                                   const int32_t & mini_cns_score, const int & matrix_boundary_distance,
                                                                                   const int & _open_gap_penalty1, const int & _extend_gap_penalty1, const int & _open_gap_penalty2,
                                                                                   const int & _extend_gap_penalty2, const int & matchingScore, const int & mismatchingPenalty,
                                                                                   const Scorei & m, const int32_t & step_size, std::string & seq1_string, std::string & seq2_string,
                                                                                   const double & pvalues, const double & lambda, const double & kValue, const int32_t & zDrop,
                                                                                   const int32_t & bandwidth, int & _lDiag, int & _uDiag, const int & xDrop, Score & score,
                                                                                   int16_t * weight, int16_t * weight_rev);

std::vector<PairedSimilarFragment> findSimilarFragmentsForPairedSequence_wighted_1gap ( int8_t * seq1, int8_t * seq1_rev_com,
                                                                                        int8_t * seq2, int8_t * seq2_rev_com, int32_t & length1, int32_t & length2, int32_t & windowsSize,
                                                                                        const int32_t & mini_cns_score, const int & matrix_boundary_distance,
                                                                                        const int & _open_gap_penalty1, const int & _extend_gap_penalty1, const int & _open_gap_penalty2,
                                                                                        const int & _extend_gap_penalty2, const int & matchingScore, const int & mismatchingPenalty,
                                                                                        const Scorei & m, const int32_t & step_size, std::string & seq1_string, std::string & seq2_string,
                                                                                        const double & pvalues, const double & lambda, const double & kValue, const int32_t & zDrop,
                                                                                        const int32_t & bandwidth, int & _lDiag, int & _uDiag, const int & xDrop, Score & score,
                                                                                        int16_t * weight, int16_t * weight_rev);


/*
std::vector<PairedSimilarFragment> findSimilarFragmentsForPairedSequence ( int8_t * seq1, int8_t * seq1_rev_com,
        int8_t * seq2, int8_t * seq2_rev_com, uint32_t & length1, uint32_t & length2, uint32_t & windowsSize,
        const uint32_t & mini_cns_size, const int32_t & mini_cns_score, const int & matrix_boundary_distance,
        const int & _open_gap_penalty, const int & _extend_gap_penalty, const int & matchingScore,
        const int & mismatchingPenalty, const Scorei & m, const uint32_t & step_size,
        std::string & seq1_string, std::string & seq2_string, const double & pvalue, const double & lambda, const double & kValue, const int & lDiag, const int & uDiag, const int & xDrop);

std::vector<PairedSimilarFragment> findSimilarFragmentsForPairedSequence ( int8_t * seq1, int8_t * seq1_rev_com,
           int8_t * seq2, int8_t * seq2_rev_com, uint32_t & length1, uint32_t & length2, uint32_t & windowsSize,
           const int32_t & mini_cns_score, const int & matrix_boundary_distance,
           const int & _open_gap_penalty1, const int & _extend_gap_penalty1,
           const int & _open_gap_penalty2, const int & _extend_gap_penalty2, const int & matchingScore,
           const int & mismatchingPenalty, const Scorei & m, const uint32_t & step_size,
           std::string & seq1_string, std::string & seq2_string, const double & pvalues,
           const double & lambda, const double & kValue, const uint32_t & zDrop, const int32_t & bandwidth,
           const int & lDiag, const int & uDiag, const int & xDrop);
*/
/*

std::vector<PairedSimilarFragment> findSimilarFragmentsForPairedSequence ( int8_t * seq1, int8_t * seq1_rev_com,
           int8_t * seq2, int8_t * seq2_rev_com, uint32_t & length1, uint32_t & length2, uint32_t & windowsSize,
           const uint32_t & mini_cns_size, const int32_t & mini_cns_score, const int & matrix_boundary_distance,
           const int & _open_gap_penalty1, const int & _extend_gap_penalty1,
           const int & _open_gap_penalty2, const int & _extend_gap_penalty2, const int & matchingScore,
           const int & mismatchingPenalty, const Scorei & m, const uint32_t & step_size,
           std::string & seq1_string, std::string & seq2_string, const double & pvalues,
           const double & lambda, const double & kValue, const uint32_t & zDrop, const int32_t & bandwidth);

*/
/*
std::vector<PairedSimilarFragment> findSimilarFragmentsForPairedSequence ( int8_t * seq1, int8_t * seq1_rev_com,
           int8_t * seq2, int8_t * seq2_rev_com, uint32_t & length1, uint32_t & length2, uint32_t & windowsSize,
           const int32_t & mini_cns_score, const int & matrix_boundary_distance,
           const int & _open_gap_penalty1, const int & _extend_gap_penalty1,const int & matchingScore,
           const int & mismatchingPenalty, const Scorei & m, const uint32_t & step_size,
           std::string & seq1_string, std::string & seq2_string, const double & pvalues,
           const double & lambda, const double & kValue, const int & lDiag, const int & uDiag, const int & xDrop);

std::vector<PairedSimilarFragment> findSimilarFragmentsForPairedSequence_wighted ( int8_t * seq1, int8_t * seq1_rev_com,
           int8_t * seq2, int8_t * seq2_rev_com, uint32_t & length1, uint32_t & length2, uint32_t & windowsSize,
           const int32_t & mini_cns_score, const int & matrix_boundary_distance,
           const int & _open_gap_penalty1, const int & _extend_gap_penalty1, const int & _open_gap_penalty2,
           const int & _extend_gap_penalty2, const int & matchingScore, const int & mismatchingPenalty,
           const Scorei & m, const uint32_t & step_size, std::string & seq1_string, std::string & seq2_string,
           const double & pvalues, const double & lambda, const double & kValue, const uint32_t & zDrop,
           const int32_t & bandwidth, int & _lDiag, int & _uDiag, const int & xDrop, Score & score,
           int16_t * weight, int16_t * weight_rev);
std::vector<PairedSimilarFragment> findSimilarFragmentsForPairedSequence_wighted_1gap ( int8_t * seq1, int8_t * seq1_rev_com,
            int8_t * seq2, int8_t * seq2_rev_com, uint32_t & length1, uint32_t & length2, uint32_t & windowsSize,
            const int32_t & mini_cns_score, const int & matrix_boundary_distance,
            const int & _open_gap_penalty1, const int & _extend_gap_penalty1, const int & _open_gap_penalty2,
            const int & _extend_gap_penalty2, const int & matchingScore, const int & mismatchingPenalty,
            const Scorei & m, const uint32_t & step_size, std::string & seq1_string, std::string & seq2_string,
            const double & pvalues, const double & lambda, const double & kValue, const uint32_t & zDrop,
            const int32_t & bandwidth, const int & _lDiag, const int & _uDiag, const int & xDrop, Score & score,
            int16_t * weight, int16_t * weight_rev);
            */
#endif //SONG_CNS_FINDSIMILARFRAGMENTSFORPAIREDSEQUENCE_H
