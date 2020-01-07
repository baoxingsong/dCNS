//
// Created by Baoxing song on 2019-01-02.
//

#ifndef SONG_CNS_SMITHWATERMAN_H
#define SONG_CNS_SMITHWATERMAN_H

#include <string>
#include <stdio.h>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <vector>

std::vector<uint32_t> SmithWaterman(int8_t *seq1, int8_t *seq2, const uint32_t &length1,
                   const uint32_t &length2, const int &_open_gap_penalty, const int &_extend_gap_penalty,
                                    uint32_t &maxScore, uint32_t &endPosition1, uint32_t &endPosition2, bool reverseAlignment, bool returnCigar);
#endif //SONG_CNS_SMITHWATERMAN_H
