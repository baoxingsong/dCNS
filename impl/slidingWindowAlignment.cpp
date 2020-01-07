//
// Created by Baoxing Song on 2019-07-09.
//

#include "slidingWindowAlignment.h"

int32_t min( int32_t & a, int32_t & b){
    return a < b ? a : b;
}

void slidingWindowAlignment ( int8_t ** seqs, std::vector<int32_t> & lengths, std::vector<std::string> & seqNames,
                            int32_t & windowsSize,
                            const int & _open_gap_penalty, const int & _extend_gap_penalty, const int & matchingScore, const int & mismatchingPenalty,
                            const std::string & output, const int32_t & step_size, const Scorei & m){

    int8_t * seqRef = seqs[0];
    int32_t lengthRef = lengths[0];
    int32_t maxScore;
    int i, j;
    int32_t length1 = windowsSize;
    int32_t maxDistance;
    std::ofstream ofile;
    int32_t length2;
    double score;
    int32_t endPosition1, endPosition2;
    int8_t * seq2;

    ofile.open(output);
    for( j=0; j<lengthRef; j+=step_size ){
        int8_t * seq1 = seqRef + j;
        if( (length1 + j) > lengthRef ){
            length1 = lengthRef-j;
        }
        for ( i=1; i< lengths.size(); ++i ){
            seq2 = seqs[i];
            length2 = lengths[i];
            bool returnPosition =false;
            SmithWaterman(seq2, seq1, length2, length1, _open_gap_penalty, _extend_gap_penalty, maxScore, endPosition1, endPosition2, m, returnPosition);

            maxDistance = min(length1, length2) * matchingScore;
            score =  double(maxScore) / double(maxDistance);
            ofile << "" << j << " " << seqNames[i] << " " << score << std::endl;
        }
    }
    ofile.close();
}
