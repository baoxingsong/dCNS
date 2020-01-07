//
// Created by Baoxing Song on 2019-07-09.
//

#ifndef AND_CNS_SLIDINGWINDOWALIGNMENT_H
#define AND_CNS_SLIDINGWINDOWALIGNMENT_H


#include <algorithm>
#include "sequenceAlignment.h"
#include "../model/model.h"
#include <vector>
#include <map>
#include "fasta.h"
#include <iostream>
#include <sstream>
#include <cassert>
#include <iostream>
#include <limits>


void slidingWindowAlignment ( int8_t ** seqs, std::vector<int32_t> & lengths, std::vector<std::string> & seqNames,
                              int32_t & windowsSize,
                              const int & _open_gap_penalty, const int & _extend_gap_penalty, const int & matchingScore, const int & mismatchingPenalty,
                              const std::string & output, const int32_t & step_size, const Scorei & m);


#endif //AND_CNS_SLIDINGWINDOWALIGNMENT_H
