//
// Created by Baoxing song on 2019-01-02.
//


#include "findSimilarFragmentsForPairedSequence.h"



std::vector<PairedSimilarFragment> findSimilarFragmentsForPairedSequence_v0 ( int8_t * seq1, int8_t * seq1_rev_com, int8_t * seq2, int8_t * seq2_rev_com, uint32_t & length1, uint32_t & length2,
                                                                           uint32_t & windowsSize, const uint32_t & mini_cns_size, int32_t & mini_cns_score,
                                                                           const int & matrix_boundary_distance, const int & _open_gap_penalty, const int & _extend_gap_penalty){

    int match_score = 1; // assume the match score is one, not changeable currently

    if( windowsSize > 2* mini_cns_size ){
        std::cerr << "you windows size is larger than 2*mini_cns_size, this might missing some similar fragments" << std::endl;
    }
    std::vector<PairedSimilarFragment> pairedSimilarFragments;
    uint32_t startPosition2 = 0;
    uint32_t this_windowsSize = windowsSize > length2 ? length2 : windowsSize;
    uint32_t maxScore;
    uint32_t endPosition1;
    uint32_t endPosition2;
    uint32_t endPosition1_rc;
    uint32_t endPosition2_rc;

    while( ((startPosition2 + this_windowsSize) <  length2) ) {
        SmithWaterman(seq1, seq2 + startPosition2, length1, this_windowsSize, _open_gap_penalty, _extend_gap_penalty,
                      maxScore, endPosition1, endPosition2, false);
        while ((this_windowsSize - endPosition2) <= matrix_boundary_distance &&
               ((startPosition2 + this_windowsSize) < length2)) {
            this_windowsSize += mini_cns_size;

            //todo, check the this_windowsSize is too large?? go boyound the sequence length
            SmithWaterman(seq1, seq2 + startPosition2, length1, this_windowsSize, _open_gap_penalty, _extend_gap_penalty,
                          maxScore, endPosition1, endPosition2, false);
            // z drop is designed for global alignment, maybe not try that
        }

        if ( mini_cns_score <= maxScore ) { //
            int32_t oldMaxScore = maxScore;
            std::vector<uint32_t > cigar = SmithWaterman(seq1_rev_com + (length1 - 1 - endPosition1), seq2_rev_com + ( length2 - (startPosition2 + endPosition2 ) - 1 ), 1+endPosition1,
                                                         endPosition2+1, _open_gap_penalty, _extend_gap_penalty, maxScore, endPosition1_rc, endPosition2_rc, true);// this is the reverse alignment
            assert(oldMaxScore == maxScore);
            if( (endPosition1_rc+1)>mini_cns_size && (endPosition2_rc+1)>mini_cns_size ){
                PairedSimilarFragment pairedSimilarFragment(endPosition1+1-endPosition1_rc, endPosition1+1,
                                                            startPosition2+endPosition2-endPosition2_rc+1, startPosition2+endPosition2+1, maxScore, cigar);
                // here the ordinate put into pairedSimilarFragment is 1 based
                pairedSimilarFragments.push_back(pairedSimilarFragment);
                //std::cout << "second maxScore " << maxScore << std::endl << std::endl;
            }
            startPosition2 += endPosition2;
            this_windowsSize = (startPosition2 + windowsSize) > length2 ? (length2 - startPosition2) : windowsSize;
        }else{
            startPosition2 += windowsSize/2;
            this_windowsSize=windowsSize;
        }
    }
    this_windowsSize = length2-startPosition2;
    if( this_windowsSize > mini_cns_size ){
        SmithWaterman(seq1, seq2 + startPosition2, length1, this_windowsSize, _open_gap_penalty, _extend_gap_penalty,
                      maxScore, endPosition1, endPosition2, false);
        if ( mini_cns_score <= maxScore ) { //
            int32_t oldMaxScore = maxScore;
            std::vector<uint32_t > cigar = SmithWaterman(seq1_rev_com + (length1 - 1 - endPosition1), seq2_rev_com + ( length2 - (startPosition2 + endPosition2 ) - 1 ), 1+endPosition1,
                                                         endPosition2+1, _open_gap_penalty, _extend_gap_penalty, maxScore, endPosition1_rc, endPosition2_rc, true);// this is the reverse alignment
            assert(oldMaxScore == maxScore);
            if( (endPosition1_rc+1)>mini_cns_size && (endPosition2_rc+1)>mini_cns_size ){
                PairedSimilarFragment pairedSimilarFragment(endPosition1+1-endPosition1_rc, endPosition1+1,
                                                            startPosition2+endPosition2-endPosition2_rc+1, startPosition2+endPosition2+1, maxScore, cigar);
                // here the ordinate put into pairedSimilarFragment is 1 based
                pairedSimilarFragments.push_back(pairedSimilarFragment);
                //std::cout << "second maxScore " << maxScore << std::endl << std::endl;
            }
            startPosition2 += endPosition2;
            this_windowsSize = (startPosition2 + windowsSize) > length2 ? (length2 - startPosition2) : windowsSize;
        }
    }
    return pairedSimilarFragments;
}


std::vector<PairedSimilarFragment> findSimilarFragmentsForPairedSequencev_1 ( int8_t * seq1, int8_t * seq1_rev_com, int8_t * seq2, int8_t * seq2_rev_com, uint32_t & length1, uint32_t & length2,
        uint32_t & windowsSize, const uint32_t & mini_cns_size, int32_t & mini_cns_score,
        const int & matrix_boundary_distance, const int & _open_gap_penalty, const int & _extend_gap_penalty){

    int match_score = 1; // assume the match score is one, not changeable currently

    if( windowsSize > 2* mini_cns_size ){
        std::cerr << "you windows size is larger than 2*mini_cns_size, this might missing some similar fragments" << std::endl;
    }
    std::vector<PairedSimilarFragment> pairedSimilarFragments;
    uint32_t startPosition2 = 0;
    uint32_t this_windowsSize = windowsSize > length2 ? length2 : windowsSize;
    uint32_t maxScore;
    uint32_t align_length1;
    uint32_t align_length2;
    uint32_t align_length1_rv;
    uint32_t align_length2_rv;

    while( ((startPosition2 + this_windowsSize) <=  length2) ) {
        SmithWaterman(seq1, seq2 + startPosition2, length1, this_windowsSize, _open_gap_penalty, _extend_gap_penalty,
                      maxScore, align_length1, align_length2, false);
        // here endPosition1 and endPosition2 are actually the number of basepairs being aligned
        while ((this_windowsSize - align_length2) <= matrix_boundary_distance &&
                ((startPosition2 + this_windowsSize) <= length2)) {
            // if the windows size is too small and the maxscore position is not within this region
            //  then enlarge the windowsSize
            this_windowsSize += mini_cns_size;
            this_windowsSize = this_windowsSize < (length2-startPosition2) ? this_windowsSize : (length2-startPosition2);
//            std::cout << "line 34, window size " << this_windowsSize << std::endl;
            SmithWaterman(seq1, seq2 + startPosition2, length1, this_windowsSize, _open_gap_penalty, _extend_gap_penalty,
                      maxScore, align_length1, align_length2, false);
                      // z drop is designed for global alignment, maybe not try that
        }
//        std::cout << "line 39" << std::endl;
        if ( mini_cns_score <= maxScore ) { //
            int32_t oldMaxScore = maxScore;
//            std::cout << "line 42 length1 " << length1 << " length2 " << length2 << " position1 " <<
//                (length1 - 1 - align_length1) << " position2 " << ( length2 - (startPosition2 + align_length2 ) - 1 ) << " end1 " << 1+align_length1
//                << " end2 " << align_length2+1 << std::endl;
//            std::cout << "length2 " << length2 << " startPosition2 " << startPosition2 << " endPosition2 " << align_length2 << std::endl;
            int32_t  seq2_position = ( length2 - (startPosition2 + align_length2) );
//            int32_t  seq2_position;
//            if( length2 > (startPosition2 + endPosition2 +1) ){
//                seq2_position = ( length2 - (startPosition2 + endPosition2 ) - 1 );
//            }else{
//                seq2_position=0;
//            }
            std::vector<uint32_t > cigar = SmithWaterman(seq1_rev_com + (length1 - align_length1), seq2_rev_com + seq2_position, align_length1,
                                         align_length2, _open_gap_penalty, _extend_gap_penalty, maxScore, align_length1_rv, align_length2_rv, true);// this is the reverse alignment
//            std::cout << "line 45 oldMaxScore: " << oldMaxScore << " maxScore: " << maxScore << std::endl;
            assert(oldMaxScore == maxScore); //todo here is a bug here
//            std::cout << "line 46" << std::endl;
            if( (align_length1_rv)>=mini_cns_size && (align_length2_rv)>=mini_cns_size ){
//                std::cout << "line 48" << std::endl;
                PairedSimilarFragment pairedSimilarFragment(align_length1+1-align_length1_rv, align_length1,
                        startPosition2+align_length2-align_length2_rv+1, startPosition2+align_length2, maxScore, cigar);
                // here the ordinate put into pairedSimilarFragment is 1 based
//                std::cout << "line 52" << std::endl;
                pairedSimilarFragments.push_back(pairedSimilarFragment);
//                std::cout << "50 second maxScore " << maxScore << std::endl << std::endl;
            }
//            std::cout << "line 56" << std::endl;
            startPosition2 += align_length2;
            this_windowsSize = (startPosition2 + windowsSize) > length2 ? (length2 - startPosition2) : windowsSize;
        }else{
            startPosition2 += windowsSize/2;
            this_windowsSize=windowsSize;
        }
    }
//    std::cout << "line 59" << std::endl;
    this_windowsSize = length2-startPosition2;
    if( this_windowsSize > mini_cns_size ){
        SmithWaterman(seq1, seq2 + startPosition2, length1, this_windowsSize, _open_gap_penalty, _extend_gap_penalty,
                      maxScore, align_length1, align_length2, false);
        if ( mini_cns_score <= maxScore ) { //
            int32_t oldMaxScore = maxScore;
            int32_t  seq2_position = ( length2 - (startPosition2 + align_length2 - 1) );
//            int32_t  seq2_position;
//            if( length2 > (startPosition2 + endPosition2 +1) ){
//                seq2_position = ( length2 - (startPosition2 + endPosition2 ) - 1 );
//            }else{
//                seq2_position=0;
//            }
            std::vector<uint32_t > cigar = SmithWaterman(seq1_rev_com + (length1 - align_length1), seq2_rev_com + seq2_position, align_length1,
                                                         align_length2, _open_gap_penalty, _extend_gap_penalty, maxScore, align_length1_rv, align_length2_rv, true);
//            std::vector<uint32_t > cigar = SmithWaterman(seq1_rev_com + (length1 - 1 - endPosition1), seq2_rev_com + seq2_position, 1+endPosition1,
//                                                         endPosition2+1, _open_gap_penalty, _extend_gap_penalty, maxScore, endPosition1_rc, endPosition2_rc, true);// this is the reverse alignment
            assert(oldMaxScore == maxScore);
            if( (align_length1_rv)>=mini_cns_size && (align_length2_rv)>=mini_cns_size ){
                PairedSimilarFragment pairedSimilarFragment(align_length1+1-align_length1_rv, align_length1+1,
                                                            startPosition2+align_length2-align_length2_rv+1, startPosition2+align_length2+1, maxScore, cigar);
                // here the ordinate put into pairedSimilarFragment is 1 based
                pairedSimilarFragments.push_back(pairedSimilarFragment);
                std::cout << "73 second maxScore " << maxScore << std::endl << std::endl;
            }
        }
    }
    return pairedSimilarFragments;
}




std::vector<PairedSimilarFragment> findSimilarFragmentsForPairedSequence ( int8_t * seq1, int8_t * seq1_rev_com, int8_t * seq2, int8_t * seq2_rev_com, uint32_t & length1, uint32_t & length2,
                                                                              uint32_t & windowsSize, const uint32_t & mini_cns_size, const int32_t & mini_cns_score,
                                                                              const int & matrix_boundary_distance, const int & _open_gap_penalty, const int & _extend_gap_penalty, const uint32_t & step_size){

    if( windowsSize > 2* mini_cns_size ){
        std::cerr << "the windows size is larger than 2*mini_cns_size, this might missing some similar fragments" << std::endl;
    }
    if( step_size*2 > mini_cns_score ){
        std::cerr << "step_size* is larger than mini_cns_score, the sequence alignment method would not work properly" << std::endl;
    }
    std::vector<PairedSimilarFragment> pairedSimilarFragments;
    uint32_t startPosition2 = 0;
    uint32_t this_windowsSize = windowsSize > length2 ? length2 : windowsSize;
    uint32_t maxScore;
    uint32_t endPosition1;
    uint32_t endPosition2;
    uint32_t endPosition1_rc;
    uint32_t endPosition2_rc;

    while( ((startPosition2 + this_windowsSize) < length2) ) { // could not put == here, since the first one is startPosition2 + 0
        SmithWaterman(seq1, seq2 + startPosition2, length1, this_windowsSize, _open_gap_penalty, _extend_gap_penalty,
                      maxScore, endPosition1, endPosition2, false);
        while ((this_windowsSize - endPosition2) <= matrix_boundary_distance &&
               ((startPosition2 + this_windowsSize) < length2)) {
            this_windowsSize += step_size;
            this_windowsSize = this_windowsSize < (length2-startPosition2) ? this_windowsSize : (length2-startPosition2);
            // std::cout << "line 212 startPosition2: " << startPosition2 << " this_windowsSize: " << this_windowsSize << std::endl;
            SmithWaterman(seq1, seq2 + startPosition2, length1, this_windowsSize, _open_gap_penalty, _extend_gap_penalty,
                          maxScore, endPosition1, endPosition2, false);
            // std::cout << "line 215 endPosition1: " << endPosition1 << " endPosition2: " << endPosition2 << " startPosition2: " << startPosition2 << std::endl;
            // z drop is designed for global alignment, maybe not try that
        }

        if ( mini_cns_score <= maxScore ) { //
            int32_t oldMaxScore = maxScore;
            //std::cout << "line 221 maxScore " << maxScore << std::endl << std::endl;
            std::vector<uint32_t > cigar = SmithWaterman(seq1_rev_com + (length1 - 1 - endPosition1), seq2_rev_com + ( length2 - (startPosition2 + endPosition2 ) - 1 ), 1+endPosition1,
                                                         endPosition2+1, _open_gap_penalty, _extend_gap_penalty, maxScore, endPosition1_rc, endPosition2_rc, true);// this is the reverse alignment
            assert(oldMaxScore == maxScore);
            // std::cout << "line 225 maxScore " << maxScore << std::endl << std::endl;
            if( (endPosition1_rc+1)>mini_cns_size && (endPosition2_rc+1)>mini_cns_size && endPosition1>=endPosition1_rc
                && endPosition2>=endPosition2_rc){
                PairedSimilarFragment pairedSimilarFragment(endPosition1+1-endPosition1_rc, endPosition1+1,
                                                            startPosition2+endPosition2-endPosition2_rc+1, startPosition2+endPosition2+1, maxScore, cigar);
                // here the ordinate put into pairedSimilarFragment is 1 based
                pairedSimilarFragments.push_back(pairedSimilarFragment);
                // std::cout << "line 230 maxScore " << maxScore << std::endl << std::endl;
            }
            startPosition2 += endPosition2;
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
                      maxScore, endPosition1, endPosition2, false);
        // std::cout << "line 243" << std::endl;
        if ( mini_cns_score <= maxScore ) { //
            int32_t oldMaxScore = maxScore;
            std::vector<uint32_t > cigar = SmithWaterman(seq1_rev_com + (length1 - 1 - endPosition1), seq2_rev_com + ( length2 - (startPosition2 + endPosition2 ) - 1 ), 1+endPosition1,
                                                         endPosition2+1, _open_gap_penalty, _extend_gap_penalty, maxScore, endPosition1_rc, endPosition2_rc, true);// this is the reverse alignment
            assert(oldMaxScore == maxScore);
            if( (endPosition1_rc+1)>mini_cns_size && (endPosition2_rc+1)>mini_cns_size ){
                PairedSimilarFragment pairedSimilarFragment(endPosition1+1-endPosition1_rc, endPosition1+1,
                                                            startPosition2+endPosition2-endPosition2_rc+1, startPosition2+endPosition2+1, maxScore, cigar);
                // here the ordinate put into pairedSimilarFragment is 1 based
                pairedSimilarFragments.push_back(pairedSimilarFragment);
                //std::cout << "second maxScore " << maxScore << std::endl << std::endl;
            }
        }
    }
    return pairedSimilarFragments;
}
