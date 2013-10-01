#ifndef _MXFLOWGRAPH_H
#define _MXFLOWGRAPH_H

#include <list>
#include <set>
#include <vector>
#include <map>
#include <queue>
#include <limits>
#include <stdlib.h>

using namespace std;

// constant internal index for source and sink
const int source = 0;
const int sink = 1;
// constant code for hitting times : never and now. 
const double neverSplit=-1;
const double splitNow=-2;

const double infinite = std::numeric_limits<double>::max();

// define a struct that contains all the necessary information on the edges
struct Edge {
    double capacity; // what is the capacity of the node
    double flow; // flow at the node
    double tension; // tension on that node; saved as lambda * t
    double lambda; // lambda when the tension was saved 
	double weights; // constant innitial weights needed for updates
};

struct MaxFlowEdge {
	int from; // from where the edge come from 
    int to; // to what it is pointed
    Edge* edgePtr; // pointer to the edge
	Edge* edgePtrBack; // pointer to the backwards edge, enable to modify the reverse edge easily in mxflow
};

typedef map<int,Edge*> Node; // to
typedef map<int,Node> Nodes; // from

typedef vector<MaxFlowEdge> MaxFlowNode; 
typedef vector<MaxFlowNode> MaxFlowNodes; // vect of vect of edges


class MxFlowGraph
{
	int N ; // nb of nodes <------------------------------------------------------ check when called = nb of nodes + source & sink ?
    MaxFlowNodes G; // vector (1 per node) of vect of mxflow egde
	vector<double> excess; // excess of flow for each node
	vector<int> height,  count; // dist is the height, count number of vertex of this height
	vector<bool> active;
	queue<int> Q; // a queue
    
    // the following variables map internal to external nodes; source and sink always excluded
    map<int,int> nodeMapExtToInt; // maps the external node numbers to the internal node numbers
    vector<int> nodeMapIntToExt; // maps the internal nodes to the external ones
    double slope; // mean value of the pull on each node <----------------------------------------- needed ??????????
    double lambda; // cause different groups are at different lambda. != lambda_group cause this one can be updated to update tau_kl
    double epsilon; // tolerance parameter
	
    friend class MainGraph;
	
	// add edge from, to. A pointer to an edge in maingraph. notation in internal numeration
	void addMxFlowEdge(int from, int to, Edge *edgePtr, Edge *edgePtrBack);
	// add edge from source or to sink 
	void addSourceSink(int from, int to, double capacity);
	
	// create artificial source and sink for mxflow update to be faster
	pair<int,int> addSpecialSourceSink(const vector<double>& overFlow); 
	// delete these sources and sinks
	void removeSpecialSourceSink(const vector<double>& overFlow, const int newSource, const int newSink);
	
	//////////////////////////////////////////////////////////////////////////////////
	/////////////////////////  Push Relabel Functions  ///////////////////////////////
	////////////////////////////////////////////////////////////////////////////////// 
	void Enqueue(int v); // put a vertex in the queue
	void Push(MaxFlowEdge &e); // push excess on one edge
	void Gap(int k); // for gap heuristic
	void Relabel(int v); // the relabel function
	void Discharge(int v); // Discharge vertex v (ie push and relabel until inactive)
	void RunMaxFlow(int s, int t); // give mxflow given source index and sink index
	//////////////////////////////////////////////////////////////////////////////////
	
	bool checkMax(const int sourceNode); // check if the maximum flow is at max capacity
	double validUntil(); // give the lambda for which the maxflow solution ceased to be ok
	vector<int> dfs(int start, bool from); 
	// calculate the distance from the start point following edgePtr if from true and edgePtrBack if false
	
	void updateTension(const double newlambda);
	void updateCapacity(const double newLambda, vector<double>& overFlow);
	void setCapacity();
	void setFlowTo0();
	void setMaxCapacity();
	
	// for destructor
	void deleteAllEdges(int nodeNum);
	
public:
	
    // constructor for MxFlowGraph object
	// does not do much. The edges will be added via a Subgraph call in Maingraph.cpp
    MxFlowGraph(const set<int>& graphNodes);	

    // destructor for MxFlowGraph object
    ~MxFlowGraph();

	// calculate maxflow and return the hitting time lambda it has to be calculated again
	double mxFlowNew(double lambda); // for a group just formed
	double mxFlowUpdate(double lambda); // if we have hitted and need to update
	
	// returns a set of nodes that are reachable from the source with current flow and its complement
    pair<set<int>,set<int> > GetSets(const int sourceNode=source);
	
	// <------------------------------- here needed a function that test nb of nodes in graph and test if possible. 
	
    // prints out the whole graph; used for troubleshooting
	void print_Graph();
	// print active nodes
	void print_active();
    
};

double Dif(double a, double b);

#endif