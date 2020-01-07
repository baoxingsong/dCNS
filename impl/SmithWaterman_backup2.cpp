//
// Created by Baoxing song on 2019-01-02.
//

/*
 * This one could only be used to get the maximum score and maximum score position
 * It could not be used to get the sequence alignment, since score matrix is not enough to trace back
 * **/


#include "SmithWaterman_bas.h"

static int8_t m[5][5] = {
        //  A,     C,     N,     G,     T
        {    2,    -1,     0,    -1,   -1} ,   //A
        {   -1,     2,     0,    -1,   -1} ,   //C
        {    0,     0,     -1,     0,    0} ,   //N
        {   -1,    -1,     0,     2,   -1},     //G
        {   -1,    -1,     0,    -1,    2}     //T
};

// for the reverse alignment, the alignment must start from the first bp for both query and target
void reverseSmithWatermanMatrixInitialization( const uint32_t &length1, const uint32_t &length2, int32_t ** E, int32_t ** F, int q, int e){
    int i=0, j=0;
    for (j = 0; j < (length2 + 1); ++j) {
        E[i][j] = q;
        F[i][j] = q + e*j;
    }
    j = 0;
    for (i = 0; i < (length1 + 1); ++i) {
        E[i][j] = q + e*j;
        F[i][j] = q;
    }
}

/**
 * reverseAlignment is used to align the recverse sequence for CNS detection
 * if it is true, the algorithm try hardly to start the alignment from the first base pair for both query and target sequence
 * */

std::vector<uint32_t> SmithWaterman(int8_t *seq1, int8_t *seq2, const uint32_t &length1,
                   const uint32_t &length2, const int &_open_gap_penalty, const int &_extend_gap_penalty,
                                    uint32_t &maxScore, uint32_t &endPosition1, uint32_t &endPosition2, bool reverseAlignment){

    int q = _open_gap_penalty+_extend_gap_penalty;
    int e = _extend_gap_penalty;

    int32_t ** M = new int32_t *[length1 + 1];
    int32_t ** E = new int32_t *[length1 + 1];
    int32_t ** F = new int32_t *[length1 + 1];

    int8_t ** dM = new int8_t *[length1 + 1];
    int8_t ** dE = new int8_t *[length1 + 1];
    int8_t ** dF = new int8_t *[length1 + 1];

    int32_t i, j;
    for (i = 0; i < (length1 + 1); ++i) {
        M[i] = new int32_t [length2 + 1];
        E[i] = new int32_t [length2 + 1];
        F[i] = new int32_t [length2 + 1];
        std::fill_n(M[i], length2+1, 0);
        std::fill_n(E[i], length2+1, 0);
        std::fill_n(F[i], length2+1, 0);

        dM[i] = new int8_t [length2 + 1];
        dE[i] = new int8_t [length2 + 1];
        dF[i] = new int8_t [length2 + 1];
        std::fill_n(dM[i], length2+1, 0);
        std::fill_n(dE[i], length2+1, 1);
        std::fill_n(dF[i], length2+1, 2);
    }
    if( reverseAlignment ){
        //reverseSmithWatermanMatrixInitialization( length1, length2, E, F, q, e);
    }
    endPosition1=0;
    endPosition2=0;
    maxScore = 0;
    for ( i=1; i<=length1; ++i ){
        for ( j=1; j<=length2; ++j ){
            E[i][j] = (q + M[i][j-1]) > (e + E[i][j-1]) ? (q + M[i][j-1]) : (e + E[i][j-1]);
            E[i][j] = E[i][j] > (q + F[i][j-1]) ? E[i][j] : (q + F[i][j-1]);

            F[i][j] = (q + M[i-1][j]) > (q + E[i-1][j]) ?(q + M[i-1][j]) : (q + E[i-1][j]);
            F[i][j] = F[i][j] > (e + F[i-1][j]) ? F[i][j] : (e + F[i-1][j]);

            M[i][j] = (m[seq1[i-1]][seq2[j-1]] + M[i-1][j-1]) > E[i][j] ? (m[seq1[i-1]][seq2[j-1]] + M[i-1][j-1]) : E[i][j];
            M[i][j] = M[i][j] > F[i][j] ? M[i][j] : F[i][j];
            M[i][j] = M[i][j] > 0 ? M[i][j] : 0;
            F[i][j] = F[i][j] > 0 ? F[i][j] : 0;
            E[i][j] = E[i][j] > 0 ? E[i][j] : 0;

            if(!reverseAlignment){

//                F[i][j] = F[i][j] > 0 ? F[i][j] : 0;
//                E[i][j] = E[i][j] > 0 ? E[i][j] : 0;
            }

            if( M[i][j] > maxScore ){ // please do not change > to >=, since we are doing local alignment
                // >= will omit the first similar fragments
                maxScore = M[i][j];
                endPosition1=i; // this means the endPosition1 and endPosition2 is 1 based coordinate
                endPosition2=j;
            }
        }
    }
    //std::cout << "maxScore " << maxScore << std::endl;
    std::vector<uint32_t> cigar;
    if( reverseAlignment ){ // cigar is only needed for reverse alignment
        uint32_t op = 0;
        uint32_t length = 1;
        // trace back begin
        // For speed up and RAM saving purpose, this is just am approximation tracing back implementation
        int ii = endPosition1;
        int jj = endPosition2;
        while (ii>0 && jj>0) {
            if (ii > 0 && jj > 0 &&
                M[ii][jj] == M[ii - 1][jj - 1] + m[seq1[ii - 1]][seq2[jj - 1]] ) {
                --ii;
                --jj;
                op = 0;
            } else if (jj > 0 && M[ii][jj] == E[ii][jj]) {
                --jj;
                op = 1;
            } else if (ii > 0 && M[ii][jj] == F[ii][jj]) {
                --ii;
                op = 2;
            } else if(M[ii][jj]==0){
                //std::cerr << "there is something wrong with smith-waterman algorithm in line 110" << std::endl;
                // here should never run, there is some problem with the code
                break;
            } else {
                std::cout << "M[ii][jj] " << M[ii][jj] << " E[ii][jj] " << E[ii][jj] << " F[ii][jj] " << F[ii][jj] << std::endl;
                std::cout << "there is something wrong with smith-waterman algorithm in line 114" << std::endl;
            }

            if( cigar.size() == 0 || op != (cigar[cigar.size() - 1]&0xf) ){
                cigar.push_back(length << 4 | op);
            }else{
                cigar[cigar.size()-1] += length<<4;
            }
        }
        // trace back end
    }

    for (i = 0; i <= length1; ++i) {
        delete[] M[i];
        delete[] E[i];
        delete[] F[i];
        delete[] dM[i];
        delete[] dE[i];
        delete[] dF[i];
    }
    delete[] M;
    delete[] E;
    delete[] F;
    delete[] dM;
    delete[] dE;
    delete[] dF;
    --endPosition1;// change the endPosition1 and endPosition2 to 0 based coordinate
    --endPosition2;
    return cigar;
}
