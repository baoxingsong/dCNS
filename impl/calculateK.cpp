//
// Created by bs674 on 6/7/19.
//

#include "calculateK.h"



// I changed to use random sequence to calculate k and lambda values, do not use this function any more
double calculateK( std::string & seq1, std::string & seq2, const int & _open_gap_penalty, const int & _extend_gap_penalty, const int & matchingScore,
                   const int & mismatchingPenalty){
    /*
    seqan::Align<seqan::String<char> > ali;
    resize(rows(ali), 2);
    assignSource(row(ali, 0), seq1);
    assignSource(row(ali, 1), seq2);

    std::cout << "Score = " << localAlignment(ali, seqan::Score<int>(matchingScore, mismatchingPenalty, _extend_gap_penalty, _open_gap_penalty), seqan::DynamicGaps()) << std::endl;
    std::cout << ali;
    std::cout << "Aligns Seq1[" << clippedBeginPosition(row(ali, 0)) << ":" << (clippedEndPosition(row(ali, 0)) - 1) << "]";
    std::cout << " and Seq2[" << clippedBeginPosition(row(ali, 1)) << ":" <<  (clippedEndPosition(row(ali, 1)) - 1) << "]" << std::endl << std::endl;
*/
}



// http://www.bioinfo.rpi.edu/bystrc/courses/biol4540/lecture7.pdf
