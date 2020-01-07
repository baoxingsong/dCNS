//
// Created by Baoxing song on 2019-01-08.
//


// decide do not implemented finally
// give this task to trimAI https://academic.oup.com/bioinformatics/article/25/15/1972/213148

#include "MsaGraph.h"

MsaGraph::MsaGraph(int n_){
    this->n = n_;
    lengths = new int[n];
    matrix = new bool*[n];
    int i, j;
    for( i=0; i<n; ++i ) {
        matrix[i] = new bool[n];
        for( j=0; j<n; ++j ){
            matrix[i][j] = false;
        }
    }
}

void MsaGraph::addLength(int index, int length){
    lengths[index]=length;
}

void MsaGraph::setOverlap(int index1, int index2){
    matrix[index1][index2] = true;
}

std::vector<int> MsaGraph::getWantedList(){
    std::vector<int> indexShouldKeep;
    int * degrees = new int[n];
    int i, j, degree, max_degree=0;
    for( i=0; i<n; ++i ) {
        degree = 0;
        for( j=0; j<n; ++j ){
            if( matrix[i][j] ){
                ++degree;
            }
        }
        degrees[i] = degree;
        if( max_degree < degree ){
            max_degree=degree;
        }
    }
    for( i=0; i<n; ++i ) {
        if( degrees[n] != max_degree ){ // so this is not a complete graph todo here has not been completed yet

        }
    }
    delete [] degrees;
    return indexShouldKeep;
}
MsaGraph::~MsaGraph(){
    for( int i=0; i<n; ++i ){
        delete [] this->matrix[i];
    }
    delete [] matrix;
    delete [] lengths;
}
