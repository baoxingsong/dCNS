//
// Created by Baoxing song on 2019-01-08.
//

#include "longestPath.h"


void longestPath(int **matrix, std::string *new_mems, int count, int no_seq, std::vector<int *> &LP_Matrix, std::vector<std::string> &LP_MEM){

    std::cout << "Inside LP "<<count<<std::endl;
    int *ml = new int[count]; //MEM length this is the weight???
    for(int i = 0; i < count; i++){
        ml[i] = new_mems[i].length();
    }

    Graph g(count + 2);
    int newCount = count + 2;
    int ** adjMat = (int **) malloc(sizeof(int*)* newCount + 1);
    for(int i = 0; i < newCount; i++){
        adjMat[i] = (int *) malloc(sizeof(int) * 2 + 1);
        for(int j = 0; j < 2; j++){
            adjMat[i][j] = 0;
        }
    }

    int v1 = count, v2 = count + 1;
    std::cout << "count" << count << std::endl;
    for(int p = 0; p < count; p++){

        //v1 -------> p Edge
        if(adjMat[p][0] == 0){
            g.addEdge(v1, p, ml[p]);
            adjMat[p][0] = 1;
        }

        // p -----> v2 Edge
        if(adjMat[p][1] == 0){
            g.addEdge(p, v2, 1);
            adjMat[p][1] = 1;
        }

        for(int q = p + 1; q < count; q++){
            int flag1 = 0, flag2 = 0;
            for(int r = 0; r < no_seq; r++){
                if(matrix[p][r] > matrix[q][r] && matrix[p][r] > (matrix[q][r] + ml[q])) flag1++;
                if(matrix[p][r] < matrix[q][r] && (matrix[p][r] + ml[p]) < matrix[q][r]) flag2++;
            }
            if(flag1 == no_seq ){
                // q ----> p Edge
                g.addEdge(q, p, ml[p]);

                // v1 -----> q Edge
                if(adjMat[q][0] == 0){
                    g.addEdge(v1, q, ml[q]);
                    adjMat[q][0] = 1;
                }
            }
            if(flag2 == no_seq){
                // Edge  p ----> q
                g.addEdge(p, q, ml[q]);

                // Edge q -----> v2
                if(adjMat[q][1] == 0){
                    g.addEdge(q, v2, 1);
                    adjMat[q][1] = 1;
                }
            }
        }
    }
    std::cout<< "Calling longest path ..."<<std::endl;
    g.longestPath(count, matrix, new_mems, count, no_seq, LP_Matrix, LP_MEM);
    std::cout<<"after calling longest path..."<<std::endl;
    delete[] adjMat;
    delete[] ml;
}
