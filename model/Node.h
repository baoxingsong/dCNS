//
// Created by Baoxing song on 2019-01-08.
//

#ifndef SONG_CNS_NODE_H
#define SONG_CNS_NODE_H


class Node {
    private:
        int v;
        int weight;
    public:
        Node(int _v, int _w);
        int getV();
        int getWeight();
};

#endif //SONG_CNS_NODE_H
