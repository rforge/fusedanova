#include "MainGraph.h"
#include "MxFlowGraph.h"
#include <Rcpp.h>
#define SIGN(x)( ((x)==0)?(0):( (((x)>0)?(1):(-1)) ) )

/***********************************
Construction & Destruction functions
************************************/

// pseudo-constructor. Our graph is always a clique.  
void MainGraph::Construct(int K, vector< vector<double> >& Weights, double eps){
	epsilon = eps;
	for (int i = 0; i<(K-1); i++){
		for (int j = i+1;j<K;j++){
			addEdge(i,j,-1,Weights[i][j]); // <------------------------------ what happens for equal x ? is sign ok ?
		}
	}
}	

void MainGraph::addEdge(int from, int to, int sign, double W)
{
    Edge* e1 = new Edge();
    Edge* e2 = new Edge();
    
    // initialize some values
    e1->tension=0;
    e1->lambda=0;
    e2->tension=0;
    e2->lambda=0;
	
	e1->weights = W;
    e2->weights = W;
	
    // set the flow (also works as derivative)
    e1->flow = sign * W;
    e2->flow = -sign *W;
    
    // set the capacity 
    if(sign==1){
	
        e1->capacity = W;
        e2->capacity = infinite;
		
    }else if(sign==-1){
	
        e1->capacity = infinite;
        e2->capacity = W;
		
    }else{
	
         e1->capacity = W;
         e2->capacity= W;
    }
    
    nodes[from][to]=e1; // insert the edge 
    nodes[to][from]=e2;
}


MainGraph::~MainGraph()
{
    Node::iterator edgeIt; // iterator over the edges
    Nodes::iterator nodeIt; // iterator over the nodes
    
    for(nodeIt = nodes.begin(); nodeIt != nodes.end(); ++nodeIt)
    {
        for(edgeIt = nodeIt->second.begin(); edgeIt != nodeIt->second.end(); ++edgeIt)
        {
            delete edgeIt->second;
        }
    }
}

/*************************************************************
Functions generating a subgraph on which the maxflow is tested
**************************************************************/

MxFlowGraph* MainGraph::subGraph(const set<int>& subNodes, map<int,Group *>& content, const double& slope)
{

    MxFlowGraph* m = new MxFlowGraph(subNodes);

    list<pair<int,double> >  nodePull; // saves the pull on the nodes; stored with nodeNumbers in internal notation

	m->slope = slope;
	m->epsilon = epsilon; // set the tolerance parameter
	

    subGraphGetEdges(*m, nodePull, content, slope); // nodePull will be changed and so will m

	subGraphSourceSink(*m, nodePull);
	
    return(m);
}


void MainGraph::subGraphGetEdges(MxFlowGraph& m, list<pair<int,double> >& nodePull, map<int,Group *>& content, const double& slope)
{
	map<int,int>::iterator MI; // iterate through the edges in the subgraph
	Node::iterator edgeIt; // iterate through the edges in each node; check for each if in subgraph
	Nodes::iterator nodeIt; // save the position of the node the algorithm is working on
	pair<int,double> foo; // used to save intermediate results
	int fromNodeIntNum, toNodeIntNum; // saves the internal number of the current node
	Edge *edgePtr, *edgePtrBack;
	
	// go through all the nodes in the subgraph
	for(MI = m.nodeMapExtToInt.begin(); MI!=m.nodeMapExtToInt.end(); ++MI){

		// get iterator for the node in MainGraph
		nodeIt = nodes.find(MI->first); // correspond to the "from"

		// foo is used to calculate the pull on the current node; stored in internal notation
		foo.first = MI->second;
		foo.second=0;
		
		fromNodeIntNum = MI->second; 
		
		// cycle through all the edges in the current node
		for(edgeIt=nodeIt->second.begin(); edgeIt!=nodeIt->second.end(); ++edgeIt){
			// for each edge check if it is in the subgraph (targetnode of the edge)
			if(m.nodeMapExtToInt.count(edgeIt->first)){
				
			    // only check edges for which the target node has higher number than origin node
                // edge both ways is included in addMxFlowEdge
				if(edgeIt->first > nodeIt->first){
					// copy the pointer to the edge
					toNodeIntNum = m.nodeMapExtToInt[edgeIt->first];
					edgePtr = edgeIt->second;
					edgePtrBack = nodes[edgeIt->first][nodeIt->first];
					m.addMxFlowEdge(fromNodeIntNum, toNodeIntNum, edgePtr, edgePtrBack);
					m.addMxFlowEdge(toNodeIntNum, fromNodeIntNum, edgePtrBack, edgePtr);

				}
	
			}else{ // count the pull on the node
				foo.second -= SIGN(edgeIt->second->flow) * edgeIt->second->weights; 
				// rewrite weights to avoid error propagation if numerous splits
			}
		}
	
		foo.second = foo.second - (content[MI->first]->nbind) * slope;  // update nodepull with with ngroup*slope
		nodePull.push_front(foo); // total pull on the node
	
	}
}

void MainGraph::subGraphSourceSink(MxFlowGraph& m, list<pair<int, double> >& nodePull)
{
	pair<int,double> foo; // used to save intermediate results
	
	// generate the source and sink node in MxFlowGraph with the appropriate edges
	while(!nodePull.empty()){
		// get the information about the first node
		foo = nodePull.front();
		nodePull.pop_front();

		if(foo.second>0){
			m.addSourceSink(source, foo.first, foo.second);
		}else if(foo.second<0){
			m.addSourceSink(foo.first, sink, -foo.second);
		}
	}
}

void MainGraph::printGraph(ostream& outStream)
{
    Node::iterator edgeIt; // iterator over the edges
    Nodes::iterator nodeIt; // iterator over the nodes
    
    for(nodeIt = nodes.begin(); nodeIt != nodes.end(); ++nodeIt)
    {
        outStream << "Node Number: " << nodeIt->first << endl;
        outStream << "Edges:" << endl;
        for(edgeIt = nodeIt->second.begin(); edgeIt != nodeIt->second.end(); ++edgeIt)
        {
            outStream << "To: " << edgeIt->first << " Cap: " << edgeIt->second->capacity <<
                " Flow: " << edgeIt->second->flow << " Tension: " << edgeIt->second->tension << " Lambda: " << 
                 edgeIt->second->lambda << " Weights: " << edgeIt->second->weights << endl;
        }
        outStream << endl;
    }
    outStream << endl;
}


// int main(int argc, char** argv)
// {
    // MainGraph g;
    
    // g.addEdge(1,2,-1,1);
    // g.addEdge(1,3,-1,1);
    // g.addEdge(1,4,-1,1);
    // g.addEdge(1,5,-1,1);
	
    // g.addEdge(2,3,-1,3);
    // g.addEdge(2,4,-1,3);
    // g.addEdge(2,5,-1,3);
	
    // g.addEdge(3,4,-1,2);
    // g.addEdge(3,5,-1,2);
  
    // g.printGraph(cout);
    
    // set<int> subNodes;
    // subNodes.insert(2);
    // subNodes.insert(3);
    // subNodes.insert(4);
	// subNodes.insert(5);
    
	
	// Group* g2 =0;
	// Group* g3=0;
	// Group* g4 =0;
	// Group* g5 =0;
	// g2 =new Group;
	// g3 =new Group;	
	// g4 =new Group;
	// g5 =new Group;
	
	// g2->nbind = 1;
	// g3->nbind = 1;
	// g4->nbind = 1;
	// g5->nbind = 1;
	
	// g2->slope = 1;
	// g3->slope = 2;
	// g4->slope = 3;
	// g5->slope = 4;
	
	// g2->beta = 1;
	// g3->beta = 2;
	// g4->beta = 3;
	// g5->beta = 4;	
	
	// g2->lambda = 0;
	// g3->lambda = 0;
	// g4->lambda = 0;
	// g5->lambda = 0;
	
	// map<int,Group *> content;
	// content.insert(make_pair(2,g2));
	// content.insert(make_pair(3,g3));	
	// content.insert(make_pair(4,g4));	
	// content.insert(make_pair(5,g5));	
	
	// double slope =0;
	
    // MxFlowGraph *m = g.subGraph(subNodes,content,slope);
    // m->printGraph(cout);
  
	// return 0;
// }





