//
// Created by Baoxing song on 2019-01-08.
//

#include "Graph.h"


Graph::Graph(int V){
    this->V = V;
    adj = new std::list<Node>[V];
    revAdj = new std::list<Node>[V];
}

void Graph::addEdge(int u, int v, int weight){ // so this is a directed edge
    Node node(v, weight);
    Node node1(u, weight);
    adj[u].push_back(node); // Add v to uís list
    revAdj[v].push_back(node1);
}

// A recursive function used by longestPath.
void Graph::topologicalSortUtil(int v, bool visited[], std::stack<int> &Stack){
    // Mark the current node as visited
    visited[v] = true;
    // Recur for all the vertices adjacent to this vertex
    std::list<Node>::iterator i;
    for (i = adj[v].begin(); i != adj[v].end(); ++i){
        Node node = *i;
        if (!visited[node.getV()]){
            topologicalSortUtil(node.getV(), visited, Stack);
        }
    }
    // Push current vertex to stack which stores topological sort
    Stack.push(v);
}

// The function to find longest distances from a given vertex. It uses
// recursive topologicalSortUtil() to get topological sorting.
void Graph::longestPath(int s, int **matrix, std::string *new_mems, int count, int no_seq, std::vector<int *> &LP_Matrix, std::vector<std::string> &LP_MEM){
    int u = count;
    //cout<< "u = "<<u<<endl;
    std::list<Node>::iterator i;
    for (i = adj[u].begin(); i != adj[u].end(); ++i){
        //cout<< "I am here..."<<endl;
        //cout<<"i->getV()="<<i->getV()<<" i->getWeight()="<<i->getWeight()<<endl;

    }
    //exit(0);
    std::cout<<"Inside longest path .."<<std::endl;
    int ** mat = (int **) malloc(sizeof(int *) * count);
    for(int k=0; k <count; k++){
        mat[k] = (int *) malloc(sizeof(int)*no_seq + 1);
    }
    for(int i=0; i<count; i++){
        int diff = 0;
        for(int j=0; j<no_seq; j++){
            mat[i][j] = matrix[i][j];
            for(int k=j+1; k<no_seq; k++){
                diff += fabs(matrix[i][j] - matrix[i][k]);
            }
            //mat[k][ind] = fabs(mat[k][2] - mat[k][3]) + fabs(mat[k][2] - mat[k][4]) + fabs(mat[k][3] - mat[k][4]);
        }
        mat[i][no_seq] = diff;
        //k++;
    }

    for(int k=0; k<count; k++){
        for(int j=0; j<no_seq+1; j++){
            std::cout << mat[k][j]<<" ";
        }
        std::cout<<std::endl;
    }
    //exit(0);
    std::stack<int> Stack;
    int dist[V];

    // Mark all the vertices as not visited
    bool *visited = new bool[V];
    for (int i = 0; i < V; i++){
        visited[i] = false;
    }

    // Call the recursive helper function to store Topological Sort
    // starting from all vertices one by one
    for (int i = 0; i < V; i++){
        if (visited[i] == false){
            topologicalSortUtil(i, visited, Stack);
        }
    }
    //cout<<"for testing -- 2\n";
    // Initialize distances to all vertices as infinite and distance
    // to source as 0
    for (int i = 0; i < V; i++) {
        dist[i] = INT_MIN;
    }
    dist[s] = 0;

    // Process vertices in topological order
    //cout<< "testing..while loop"<<endl;
    while (Stack.empty() == false){
        // Get the next vertex from topological order
        int u = Stack.top();
        //cout<<"u="<<u<<endl;
        //cout<<"u->getWeight() "<<u->getWeight()<<endl;
        Stack.pop();

        // Update distances of all adjacent vertices
        std::list<Node>::iterator i;
        if (dist[u] != INT_MIN){
            //cout<<"dist[u] != NINF"<<u<<endl;
            for (i = adj[u].begin(); i != adj[u].end(); ++i){
                //cout<< "I am here..."<<endl;
                //cout<<"i->getV()="<<i->getV()<< "dist[u]="<<dist[u]<<" i->getWeight()="<<i->getWeight()<<endl;
                if (dist[i->getV()] < dist[u] + i->getWeight()){
                    dist[i->getV()] = dist[u] + i->getWeight();
                }
            }
        }
    }
    //cout<<"testing ends...while loop"<<endl;
    // Print the calculated longest distances
    int max = s;
    /*vector<int *> LP_Matrix;
    vector<string> LP_MEM;*/
    //FILE *MEMList = NULL;
    //char temp[100];
    for (int i = 0; i < V; i++){
        //cout<<i<<"--"<<dist[i]<<endl;
        //(dist[i] == NINF)? cout << "INF ": cout << dist[i] << " ";
        if(dist[i] >= dist[s]) max = i;
    }
    //cout<<"The destination vertex is "<<max<<" with score "<<dist[max]<<"\n";

    int flag = 0;
    int start = max;
    int st = max;
    //cout<<max<<"("<<dist[max]<<")";
    std::list<Node>::iterator i1;
    int last = 0;

    while(start != s){
        //cout << "Start:"<<start<<endl;
        int d = -999;
        for(i1 = revAdj[start].begin(); i1 != revAdj[start].end(); ++i1){
            if(dist[i1->getV()] > d){
                d = dist[i1->getV()];
                st = i1->getV() ;
            }
        }
        if(st == s) break;
        //cout<<"TESTING HERE "<<endl;
        //cout<<" <-( ";
        int z=0; int dd = 0; int st1 = 0;
        for(i1 = revAdj[start].begin(); i1 != revAdj[start].end(); ++i1){
            // cout<< " t_s: "<<i1->getV()<<"-"<<dist[i1->getV()]<<" t_e ";
            if(dist[i1->getV()] == d) {
                if(last < 1){
                    if(z == 0){
                        dd = mat[i1->getV()][no_seq];
                        st1 = i1->getV();
                    }else {
                        if(mat[i1->getV()][no_seq] < dd){
                            dd = mat[i1->getV()][no_seq];
                            st1 = i1->getV();
                        }
                    }
                }
                //cout<< "TESTING HERE2 " << endl;
                /* added 10_14 */
                if(last >= 1){
                    if(z == 0){
                        //dd = mat[i1->getV()+1][5];
                        for(int i=0 ; i<no_seq; i++){
                            for(int j=i+1; j<no_seq; j++){
                                //cout<<"start = "<<start<<" i = "<<i<<" j = "<<j<<" i1->getV() = "<<i1->getV()<<endl;
                                int d1 = mat[start][i] - mat[i1->getV()][i];
                                int d2 = mat[start][j] - mat[i1->getV()][j];
                                dd += fabs(d1 - d2);
                                //cout<<" dd = "<<dd<<endl;
                            }
                        }
                        st1 = i1->getV();
                        //cout<<"TESTING HERE3 "<<endl;
                    }
                        //cout<<"TESTING HERE3 "<<endl;
                    else {
                        int ddd;
                        for(int i=0 ; i<no_seq; i++){
                            for(int j=i+1; j<no_seq; j++){
                                int d1 = mat[start][i] - mat[i1->getV()][i];
                                int d2 = mat[start][j] - mat[i1->getV()][j];
                                ddd += fabs(d1 - d2);
                            }
                        }
                        if(ddd< dd){
                            dd = ddd;
                            st1 = i1->getV();
                        }
                    }
                }
                z++;
            }
        }
        //cout<<"--st1:"<<st1<<endl;

        LP_Matrix.push_back(mat[st1]);
        LP_MEM.push_back(new_mems[st1]);
        start = st1;
        last++;
    }
}
