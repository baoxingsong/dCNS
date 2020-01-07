//
// Created by bs674 on 9/25/19.
//

#ifndef AND_CNS_MAPCNSTOGENOME_H
#define AND_CNS_MAPCNSTOGENOME_H


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
#include <limits>
#include <set>


PairedSimilarFragment mapCNSToGenome ( int8_t * seq1, int8_t * seq1_rev_com, int8_t * seq2, int8_t * seq2_rev_com,
        int32_t & length1, int32_t & length2, const int32_t & _open_gap_penalty1, const int32_t & _extend_gap_penalty1,
        const int32_t & _open_gap_penalty2, const int32_t & _extend_gap_penalty2, const int32_t & matchingScore,
        const int32_t & mismatchingPenalty, const Scorei & m);

#endif //AND_CNS_MAPCNSTOGENOME_H
