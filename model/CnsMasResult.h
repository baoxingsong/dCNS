//
// Created by Baoxing song on 2019-01-08.
//

#ifndef SONG_CNS_CNSMASRESULT_H
#define SONG_CNS_CNSMASRESULT_H
#include <tuple>

#include <list>
#include <stack>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <math.h>
#include <limits.h>

class CnsMasResult {
    private:
        int n;
        double ** matrix;
        std::pair<int, int> * start_ends; // the start and end position of the aligned sequence
    public:
        CnsMasResult(int _n);
        ~CnsMasResult();
        void addSequence(int index, int start, int end);
        void addSimilarity( int index1, int index2, double similarity );
};
#endif //SONG_CNS_CNSMASRESULT_H
