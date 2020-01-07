//
// Created by Baoxing song on 2019-01-02.
//


#include "findSimilarFragmentsForPairedSequence.h"

#include <iostream>

#include <seqan/align.h>
#include <seqan/align_extend.h>
#include <seqan/sequence.h>


void x_extend(uint32_t & start1, uint32_t & end1, uint32_t & start2, uint32_t & end2,
        uint32_t & maxScore, seqan::CharString seq1, seqan::CharString seq2, const int & _open_gap_penalty, const int & _extend_gap_penalty, const int & matchingScore,
                               const int & mismatchingPenalty, const double & pvalue, seqan::Align<seqan::Infix<seqan::CharString const>::Type> & align){
    //std::cout << "before " << start1 << " " << end1 << " " << start2 << " " << end2 << " " << maxScore << std::endl;
    --start1;
    --start2;
    seqan::Score<int> sc(matchingScore, mismatchingPenalty, _extend_gap_penalty, _open_gap_penalty);

    resize(rows(align), 2);
    assignSource(row(align, 0), infix(seq1, start1, end1));
    assignSource(row(align, 1), infix(seq2, start2, end2));
    int score = globalAlignment(align, sc);
    // the following diagonals.
    //int lDiag = 0, uDiag = (end1-start1+1)*2;
    int lDiag = 0, uDiag = 20;
    // Set the x-Drop value to 5.
//    int xDrop = 5 < score/2 ? 5 : score/2;
//    std::cout << "xDrop " << xDrop << std::endl;

    int xDrop = 20;
//    std::cout << "xDrop " << xDrop << std::endl;

    seqan::Tuple<unsigned, 4> positions = { {start1, start2, end1, end2} };

    /*
     *
     * We keep track of the best alignment score, denoted T , detected for a grid point
lying on an antidiagonal before the current one. If we discover that S(i, j) , T ¡ X , where X ¶ 0 is
user-speci ed, then we set S(i, j) to ¡ 1 so as to guarantee that S(i, j) will not play a role in subsequent
evaluations of S.

     */

    score = extendAlignment(align, score, seq1, seq2, positions, seqan::EXTEND_BOTH, lDiag, uDiag, xDrop, sc);

    start1 = clippedBeginPosition(row(align, 0)) + 1;
    end1 = endPosition(row(align, 0));
    start2 = clippedBeginPosition(row(align, 1)) + 1;
    end2 = endPosition(row(align, 1));
    maxScore = score;
}


std::vector<PairedSimilarFragment> findSimilarFragmentsForPairedSequence ( int8_t * seq1, int8_t * seq1_rev_com,
                   int8_t * seq2, int8_t * seq2_rev_com, uint32_t & length1, uint32_t & length2, uint32_t & windowsSize,
                   const uint32_t & mini_cns_size, const int32_t & mini_cns_score, const int & matrix_boundary_distance,
                   const int & _open_gap_penalty, const int & _extend_gap_penalty, const int & matchingScore,
                   const int & mismatchingPenalty, int8_t m[][5], const uint32_t & step_size,
                   std::string & seq1_string, std::string & seq2_string, const double & pvalues){

    if( windowsSize > 2* mini_cns_size ){
        std::cerr << "the windows size is larger than 2*mini_cns_size, this might missing some similar fragments" << std::endl;
    }
    if( step_size*2 > mini_cns_score ){
        std::cerr << "step_size* is larger than mini_cns_score, the sequence alignment method would not work properly" << std::endl;
    }
    double lambda = 0.382291;
    double kValue = 0.006662;

    seqan::CharString charString_seq1 = seq1_string;
    seqan::CharString charString_seq2 = seq2_string;

    std::vector<PairedSimilarFragment> pairedSimilarFragments;

    uint32_t startPosition2 = 0;
    uint32_t this_windowsSize = windowsSize > length2 ? length2 : windowsSize;
    uint32_t maxScore;
    uint32_t endPosition1;
    uint32_t endPosition2;
    uint32_t endPosition1_rc;
    uint32_t endPosition2_rc;
    seqan::Align<seqan::Infix<seqan::CharString const>::Type> align;
    double length1_d = double(length1);
    double length2_d = double(length2);
    double eValue;
    double pvalue;
    uint32_t start1;
    uint32_t end1;
    uint32_t start2;
    uint32_t end2;
    std::string alignment1;
    std::string alignment2;
    std::vector<uint32_t> cigar;

    int32_t oldMaxScore;

    while( ((startPosition2 + this_windowsSize) < length2) ) { // could not put == here, since the first one is startPosition2 + 0
        SmithWaterman(seq1, seq2 + startPosition2, length1, this_windowsSize, _open_gap_penalty, _extend_gap_penalty,
                      maxScore, endPosition1, endPosition2, m, false, false);
        while ((this_windowsSize - endPosition2) <= matrix_boundary_distance &&
               ((startPosition2 + this_windowsSize) < length2) && maxScore < mini_cns_score) { // once it reach mini_cns_score, it suggest this region could be aligned, and begin to extend it
            this_windowsSize += step_size;
            this_windowsSize = this_windowsSize < (length2-startPosition2) ? this_windowsSize : (length2-startPosition2);
            SmithWaterman(seq1, seq2 + startPosition2, length1, this_windowsSize, _open_gap_penalty, _extend_gap_penalty,
                          maxScore, endPosition1, endPosition2, m, false, false);
        }

        if ( mini_cns_score <= maxScore ) { //
            oldMaxScore = maxScore;
            //std::cout << "line 221 maxScore " << maxScore << std::endl << std::endl;
            SmithWaterman(seq1_rev_com + (length1 - 1 - endPosition1), seq2_rev_com + ( length2 - (startPosition2 + endPosition2 ) - 1 ), 1+endPosition1,
                                                         endPosition2+1, _open_gap_penalty, _extend_gap_penalty, maxScore, endPosition1_rc, endPosition2_rc, m, true, false);// this is the reverse alignment
            assert(oldMaxScore == maxScore);
            // std::cout << "line 225 maxScore " << maxScore << std::endl << std::endl;
            if( (endPosition1_rc+1)>mini_cns_size && (endPosition2_rc+1)>mini_cns_size && endPosition1>=endPosition1_rc
                && endPosition2>=endPosition2_rc){
                start1 = endPosition1+1-endPosition1_rc;
                end1 = endPosition1+1;
                start2 = startPosition2+endPosition2-endPosition2_rc+1;
                end2 = startPosition2+endPosition2+1;

                x_extend(start1, end1, start2, end2, maxScore, charString_seq1, charString_seq2,
                       _open_gap_penalty, _extend_gap_penalty, matchingScore, mismatchingPenalty, pvalue, align);
                // this is the reference document, I am using for those euqations
                // https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html

                eValue = kValue * length1_d * length2_d * exp(-1.0 * lambda * double(maxScore));
                pvalue = 1.0 - exp(-1.0*eValue);
                if( pvalue <= pvalues ){
                    alignment1="";
                    alignment2="";
                    cigar.clear();
                    for ( size_t i=0; i<length(row(align, 0)); ++i) {
                        uint32_t op = 0;
                        uint32_t length = 1;
                        if( row(align, 0)[i] == '-' ){
                            op = 1;
                        }else if( row(align, 1)[i] == '-' ){
                            op = 2;
                        }
                        if( cigar.size() == 0 || op != (cigar[cigar.size() - 1]&0xf) ){
                            cigar.push_back(length << 4 | op);
                        }else{
                            cigar[cigar.size()-1] += length<<4;
                        }
                        alignment1 += row(align, 0)[i];
                        alignment2 += row(align, 1)[i];
                    }
                    PairedSimilarFragment pairedSimilarFragment(start1, end1, start2, end2, maxScore, cigar, alignment1, alignment2, pvalue, eValue);
                    // here the ordinate put into pairedSimilarFragment is 1 based
                    pairedSimilarFragments.push_back(pairedSimilarFragment);
                    // std::cout << "line 230 maxScore " << maxScore << std::endl << std::endl;
                }

                startPosition2 = end2-1; // here is a 2 base pair overlap
            }else{
                startPosition2 += endPosition2;
            }
            this_windowsSize = (startPosition2 + windowsSize) > length2 ? (length2 - startPosition2) : windowsSize;
        }else{
            startPosition2 += windowsSize/2;
            this_windowsSize=windowsSize;
        }
    }
    // std::cout << "line 238" << std::endl;
    this_windowsSize = length2-startPosition2;
    if( this_windowsSize > mini_cns_size ){
        SmithWaterman(seq1, seq2 + startPosition2, length1, this_windowsSize, _open_gap_penalty, _extend_gap_penalty,
                      maxScore, endPosition1, endPosition2, m, false, false);
        // std::cout << "line 243" << std::endl;
        if ( mini_cns_score <= maxScore ) { //
            oldMaxScore = maxScore;
            SmithWaterman(seq1_rev_com + (length1 - 1 - endPosition1), seq2_rev_com + ( length2 - (startPosition2 + endPosition2 ) - 1 ), 1+endPosition1,
                                                         endPosition2+1, _open_gap_penalty, _extend_gap_penalty, maxScore, endPosition1_rc, endPosition2_rc, m, true, false);// this is the reverse alignment
            assert(oldMaxScore == maxScore);
            if( (endPosition1_rc+1)>mini_cns_size && (endPosition2_rc+1)>mini_cns_size ){
                start1 = endPosition1+1-endPosition1_rc;
                end1 = endPosition1+1;
                start2 = startPosition2+endPosition2-endPosition2_rc+1;
                end2 = startPosition2+endPosition2+1;
                x_extend(start1, end1, start2, end2, maxScore, charString_seq1,
                         alignment2, _open_gap_penalty, _extend_gap_penalty, matchingScore, mismatchingPenalty, pvalue, align);
                eValue = kValue * length1_d * length2_d * exp(-1.0 * lambda * double(maxScore));
                pvalue = 1.0 - exp(-1.0*eValue);
                if( pvalue <= pvalues ){
                    alignment1="";
                    alignment2="";
                    cigar.clear();
                    for ( size_t i=0; i<length(row(align, 0)); ++i) {
                        uint32_t op = 0;
                        uint32_t length = 1;
                        if( row(align, 0)[i] == '-' ){
                            op = 1;
                        }else if( row(align, 1)[i] == '-' ){
                            op = 2;
                        }
                        if( cigar.size() == 0 || op != (cigar[cigar.size() - 1]&0xf) ){
                            cigar.push_back(length << 4 | op);
                        }else{
                            cigar[cigar.size()-1] += length<<4;
                        }
                        alignment1 += row(align, 0)[i];
                        alignment2 += row(align, 1)[i];
                    }
                    PairedSimilarFragment pairedSimilarFragment(start1, end1, start2, end2, maxScore, cigar, alignment1, alignment2, pvalue, eValue);
                    pairedSimilarFragments.push_back(pairedSimilarFragment);
                }
            }
        }
    }
    return pairedSimilarFragments;
}
