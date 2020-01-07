//
// Created by Baoxing song on 2019-01-08.
//

#include "CnsMasResult.h"


CnsMasResult::CnsMasResult( int _n ){
    this->n=_n;
    start_ends = new std::pair<int, int>[n];
    matrix = new double*[n];
    int i, j;
    for( i=0; i<n; ++i ){
        matrix[i] = new double[n];
        start_ends[i] = std::make_pair(0, 0);
        for( j=0; j<n; ++j ){
            matrix[i][j] = INT_MIN;
        }
    }
}
void CnsMasResult::addSequence(int index, int start, int end){
    start_ends[index] = std::make_pair(start, end);
}
void CnsMasResult::addSimilarity( int index1, int index2, double similarity ){
    matrix[index1][index2] = similarity;
    matrix[index2][index1] = similarity;
}
CnsMasResult::~CnsMasResult(){
    for( int i=0; i<n; ++i ){
        delete [] this->matrix[i];
    }
    delete [] matrix;
    delete [] start_ends;
}
