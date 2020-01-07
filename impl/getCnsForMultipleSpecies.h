//
// Created by Baoxing song on 2019-01-09.
//

#ifndef SONG_CNS_GETCNSFORMULTIPLESPECIES_H
#define SONG_CNS_GETCNSFORMULTIPLESPECIES_H

#include <algorithm>
#include "SmithWaterman_bas.h"
#include "../model/model.h"
#include <vector>
#include <map>
#include <set>
#include "./findSimilarFragmentsForPairedSequence.h"
#include "syntenic.h"
#include <iostream>
#include <string>

void getCnsForMultipleSpecies ( int8_t ** seqs, int8_t ** seq_revs, std::vector<int32_t> & lengths,
                                int32_t & windowsSize, int32_t & mini_cns_seed_size, int32_t & mini_cns_score,
                                const int32_t & matrix_boundary_distance, const int32_t & _open_gap_penalty,
                                const int32_t & _extend_gap_penalty, const int32_t & matchingScore,
                                const int32_t & mismatchingPenalty, const bool & onlySyntenic, const std::string & output,
                                std::map<std::string, std::string>& sequences, std::vector<std::string> & seqNames,
                                const int32_t & step_size, std::vector<std::string> & seqs_string,
                                const int32_t & minimumNumberOfSpecies, const Scorei & m, const int32_t & mini_cns_size,
                                const double & outputWithMinimumLengthPercentage, const double & lambda,
                                const double & kValue, const int32_t & w, const int32_t & xDrop);

void getCnsForMultipleSpecies ( const bool & onlySyntenic, const std::string & output,
                                std::map<std::string, std::string> & sequences, /*species, fastaFile*/
                                std::map<std::string, std::string> & samFiles,/*species, sameFile*/
                                std::string & referenceGenomeFile,
                                const int32_t & minimumNumberOfSpecies, const int32_t & mini_cns_size,
                                const double & outputWithMinimumLengthPercentage);

void getCnsForMultipleSpecies ( const bool & onlySyntenic, const std::string & output,
                                //std::map<std::string, std::string> & sequences, /*species, fastaFile*/
                                std::map<std::string, std::string> & samFiles,/*species, sameFile*/
                                std::string & referenceGenomeFile,
                                const int32_t & minimumNumberOfSpecies, const int32_t & mini_cns_size,
                                const double & outputWithMinimumLengthPercentage);

#endif //SONG_CNS_GETCNSFORMULTIPLESPECIES_H
