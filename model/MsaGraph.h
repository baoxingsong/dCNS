//
// Created by Baoxing song on 2019-01-08.
//

#ifndef SONG_CNS_MSAGRAPH_H
#define SONG_CNS_MSAGRAPH_H
#include <vector>

class MsaGraph {
    private:
        int n;
        bool ** matrix; // position overlap or non position overlap
        int * lengths;
    public:
        MsaGraph(int n_);
        void addLength(int index, int length);
        void setOverlap(int index1, int index2);
        std::vector<int> getWantedList();
        ~MsaGraph();
};


#endif //SONG_CNS_MSAGRAPH_H
