#ifndef _MAINGRAPH_H
#define _MAINGRAPH_H


#include <vector>
#include <set>
#include <map>
#include <list>
#include <iostream>
#include "MxFlowGraph.h"
#include "Groups.h"

using namespace std;


class MainGraph
{
public:

    Nodes nodes; // a node is defined by a list of edges
	double epsilon; // tolerance parameter
	
    // adds an edge to the graph (in both nodes); only intended for use at the start of the 
    // algorithm; will set tension and lambda to 0; flow to the sign (corrected for direction b/c of 2 nodes)
    // and sets the capacity to 1 in the direction of positive sign and infinity in the other (may not be necessary) <------ see that ? 
    void addEdge(const int from, const int to, const int sign, double W); 

    // given a list of nodes, the function generates a graph containing a source and a sink
    // the edges in the graphs are pointers to the edges in the PenaltyGraph object, so that
    // operations in the MaxFlowGraph are automatically stored in the PenaltyGraph as well
    // the map can store all pulls, not only the ones of the current subgraph
    MxFlowGraph* subGraph(const set<int>& subNodes, map<int,Group *>& content, const double& slope);
    
    // copies the pointers to the edges in the subgraph, 
    // in nodePull the pull on each of the nodes will be saved (needed for the source and sink right after)
    void subGraphGetEdges(MxFlowGraph& m, list<pair<int,double> >& nodePull, map<int,Group *>& content, const double& slope);
    // generates the source and sink node, nodePull gets deleted
    void subGraphSourceSink(MxFlowGraph& m, list<pair<int, double> >& nodePull);
    // these two functions are used to create subgraphs on which we test the maximum flow for the groups

    void Construct(int K, vector< vector<double> >& Weights,double eps);
    MainGraph(){};
    // destructor that frees all the edges 
    ~MainGraph();
	
	// print helper
	void printGraph(ostream& outStream);
};

#endif
