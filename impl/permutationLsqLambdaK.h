//
// Created by bs674 on 6/8/19.
//

#ifndef AND_CNS_PERMUTATIONLSQLAMBDAK_H
#define AND_CNS_PERMUTATIONLSQLAMBDAK_H

#include <string>
#include <vector>
#include <map>
#include <cstdlib>
#include <ctime>
#include "fasta.h"
#include "sequenceAlignment.h"
#include "SequenceCharToUInt8.h"
#include "../model/model.h"
#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>

void permutationLsqLambdaK( std::string & referenceFasta, std::vector<std::string> & queryFastas,
                            const int & _open_gap_penalty, const int & _extend_gap_penalty, const int & matchingScore,
                            const int & mismatchingPenalty, int32_t & length, int & permutationTimes, int32_t & seed, bool & removen);
void permutationLsqLambdaKslow( std::string & referenceFasta, std::vector<std::string> & queryFastas,
                                const int & _open_gap_penalty, const int & _extend_gap_penalty, const int & matchingScore,
                                const int & mismatchingPenalty,  int32_t & length, int & permutationTimes, int32_t & seed);
#endif //AND_CNS_PERMUTATIONLSQLAMBDAK_H
