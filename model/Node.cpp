//
// Created by Baoxing song on 2019-01-08.
//

#include "Node.h"

Node::Node(int _v, int _w){
    v = _v;
    weight = _w;
}

int Node::getV(){
    return v;
}
int Node::getWeight()  {
    return weight;
}
