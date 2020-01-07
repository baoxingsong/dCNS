//
// Created by bs674 on 6/7/19.
//

#ifndef AND_CNS_CALCULATEK_H
#define AND_CNS_CALCULATEK_H

#include <string>
#include <iostream>


double calculateK( std::string & seq1, std::string & seq2, const int & _open_gap_penalty, const int & _extend_gap_penalty, const int & matchingScore,
                   const int & mismatchingPenalty);

#endif //AND_CNS_CALCULATEK_H
