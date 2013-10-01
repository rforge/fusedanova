#include "MxFlowGraph.h"
#include <Rcpp.h>
#include<utility>

#define Abs(x)    ((x) < 0 ? -(x) : (x))
#define Max(a, b) ((a) > (b) ? (a) : (b))
#define Min(a, b) ((a) < (b) ? (a) : (b))

/**************************************
Construction, add edge, destruction ...
***************************************/

// add edge
void MxFlowGraph::addMxFlowEdge(int from, int to, Edge *edgePtr, Edge *edgePtrBack){

    // int index;
    // MaxFlowEdge e;

    // check index where to put new edge and save it in the mxflowedge
	// index = G[from].size();
    // G[from].resize(index+1);
	// e.from =from;
	// e.to = to;
    // e.edgePtr = edgePtr;
	// e.edgePtrBack = edgePtrBack;
	
	// Rprintf("from : %i, to : %i, index : %i \n", from,to,index);
	
    // G[from][index]=e;

	// modified size => push back : addmxflowedge, addspecialsourcesink
	
	MaxFlowEdge e;
	e.from =from;
	e.to = to;
	e.edgePtr = edgePtr;
	e.edgePtrBack = edgePtrBack;
	
	G[from].push_back(e);
	
}

void MxFlowGraph::addSourceSink(int from, int to, double capacity){

    Edge* e1 = new Edge();
    Edge* e2 = new Edge();
    
    // initialize some values to 0 (they do not really matter)
    e1->tension=0;
    e1->lambda=0;
    e2->tension=0;
    e2->lambda=0;
    e1->flow = 0;
    e2->flow = 0;
    
    // set the capacity; capacity in one direction, 0 in the other
    e1->capacity=capacity;
    e2->capacity=0;
	
	e1->weights=0;
	e2->weights=0;
    
    // add the elements to the graph
    addMxFlowEdge(from, to, e1, e2);
    addMxFlowEdge(to, from, e2, e1);

}

pair<int,int> MxFlowGraph::addSpecialSourceSink(const vector<double>& overFlow){
    int newSource = G.size();
    int newSink = G.size()+1;
	G.resize(G.size()+2);
    for(unsigned int i=0; i< overFlow.size(); ++i)
    {
        if(overFlow[i]>0)
        {
            addSourceSink(newSource, i, overFlow[i]);

        }
        else if(overFlow[i]<0)
        {
            addSourceSink(i,newSink, -overFlow[i]);
        }
    }
	
	N = G.size(); // update size of graph
	
    return(pair<int,int>(newSource, newSink));
}

void MxFlowGraph::removeSpecialSourceSink(const vector<double>& overFlow, const int newSource, const int newSink){
    int numEdges;
    // first go through all other nodes and delete the appropriate edge, which by construction has to be the last one
    for(unsigned int i=0; i< overFlow.size(); ++i)
    {
        if(overFlow[i]!=0)
        {
            numEdges = G[i].size();
            G[i].erase(G[i].begin()+(numEdges-1));
			
        }
    }
    // now delete the appropriate edges
    deleteAllEdges(newSource);
    deleteAllEdges(newSink);
    // now delete the nodes itself; do it with max and min because numbering changes after erasing
    G.erase(G.begin()+Max(newSource,newSink));
    G.erase(G.begin()+Min(newSource,newSink));
	
	N = G.size(); // update size of graph
}


// constructor (part is inline ===>)
MxFlowGraph::MxFlowGraph(const set<int>& graphNodes):G(graphNodes.size()+2,MaxFlowNode(0)), excess(graphNodes.size()+2,0), height(graphNodes.size()+2, 0), count(2*(graphNodes.size()+2)+1,0), active(graphNodes.size(), false), nodeMapIntToExt(graphNodes.size()+2,-1){

	N = graphNodes.size()+2;
	int curIntNode=2;
	set<int>::iterator setIt;
	//set the mapping of the internal to the external nodes

	for(setIt = graphNodes.begin(); setIt!=graphNodes.end(); ++setIt){
		nodeMapIntToExt[curIntNode]=*setIt;
		nodeMapExtToInt[*setIt]=curIntNode;
		++curIntNode;
	}
    
	slope = 0;
	lambda = 0;
	epsilon = 0;
    
}

// destructor
MxFlowGraph::~MxFlowGraph(){

	// delete all edges go to and from the source and sink nodes
	deleteAllEdges(source);
	deleteAllEdges(sink);

	nodeMapExtToInt.clear();
	nodeMapIntToExt.clear();

	G.clear();

}

void MxFlowGraph::deleteAllEdges(int nodeNum){

	MaxFlowNode::iterator edgeIt;
	// go through all the nodes in the node nodeNum and delete them
	for(edgeIt = G[nodeNum].begin(); edgeIt!=G[nodeNum].end(); ++edgeIt){
		delete edgeIt->edgePtr;
		delete edgeIt->edgePtrBack;
	}
	G[nodeNum].clear();
}

/************************************
Functions for update the mxflow graph
*************************************/

void MxFlowGraph::updateTension(const double newlambda){

    MaxFlowNodes::iterator nodeIt;
    MaxFlowNode::iterator edgeIt;
    Edge* e;
    
	for(nodeIt = G.begin()+2; nodeIt != G.end(); ++nodeIt){ // not using source and sink
	
		for(edgeIt = nodeIt->begin(); edgeIt != nodeIt->end(); ++edgeIt){
			// check that it is not pointing to the source or sink
			if((edgeIt->to!=source) && (edgeIt->to!=sink)){
				e = edgeIt->edgePtr;
				e->tension += e->flow *(newlambda - e->lambda);
				e->lambda = newlambda;
			}
		}
	}
	
	lambda = newlambda;

}

void MxFlowGraph::setFlowTo0(){

    MaxFlowNodes::iterator nodeIt;
    MaxFlowNode::iterator edgeIt;
    
    for(nodeIt = G.begin(); nodeIt != G.end(); ++nodeIt) // iterate through all nodes
    {
        for(edgeIt = nodeIt->begin(); edgeIt != nodeIt->end(); ++edgeIt) // iterate through all edges
        {
            edgeIt->edgePtr->flow = 0;
        }
    }

}

void  MxFlowGraph::setMaxCapacity(){

    MaxFlowNodes::iterator nodeIt;
    MaxFlowNode::iterator edgeIt;
    
    for(nodeIt = G.begin()+2; nodeIt != G.end(); ++nodeIt) // iterate through all nodes, but jump over source and sink
    {
        for(edgeIt = nodeIt->begin(); edgeIt != nodeIt->end(); ++edgeIt) // iterate through all edges 
        {
            if((edgeIt->to!=source) && (edgeIt->to!=sink))
            {
                edgeIt->edgePtr->capacity = edgeIt->edgePtr->weights; // 
            }
        }
    }
}

void MxFlowGraph::setCapacity(){

	MaxFlowNodes::iterator nodeIt;
	MaxFlowNode::iterator edgeIt;
	
	for(nodeIt = G.begin()+2; nodeIt != G.end(); ++nodeIt){ // iterate through all nodes, but jump over source and sink
	
		for(edgeIt = nodeIt->begin(); edgeIt != nodeIt->end(); ++edgeIt){ // iterate through all edges
        
			// check that it is not the source or sink
			if((edgeIt->to!=source) && (edgeIt->to!=sink)){
				// set capacity only to weights if maximal tension reached; otherwise infinity
				if((edgeIt->to!=source) && (edgeIt->to!=sink)){
				
					// set capacity only to weights if maximal tension reached; otherwise infinity
					if(Dif(edgeIt->edgePtr->tension, edgeIt->edgePtr->lambda*edgeIt->edgePtr->weights)> epsilon){
						edgeIt->edgePtr->capacity = infinite;
					}else{
						edgeIt->edgePtr->capacity = edgeIt->edgePtr->weights;
					}
				}
			}
		}
	}
}

void MxFlowGraph::updateCapacity(const double newLambda, vector<double>& overFlow){

    // initialize the vector with 0s
    overFlow.assign(G.size(),0);
    
    // adjust the edges if necessary
    MaxFlowNodes::iterator nodeIt;
    MaxFlowNode::iterator edgeIt;
    
	int nodeNum;
	for(nodeIt = G.begin()+2, nodeNum=2; nodeIt != G.end(); ++nodeIt, ++nodeNum){ // iterate through all nodes, but jump over source and sink
	
		for(edgeIt = nodeIt->begin(); edgeIt != nodeIt->end(); ++edgeIt){ // iterate through all edges 
		
			if((edgeIt->to!=source) && (edgeIt->to!=sink)){ // not interested in edge going to source/sink
			
				if(edgeIt->edgePtr->capacity < infinite){
				
					if(edgeIt->edgePtr->tension < edgeIt->edgePtr->lambda * edgeIt->edgePtr->weights - epsilon){ // no flow adjustments necessary	
						edgeIt->edgePtr->capacity=infinite;
						
					}
				}
				else if(edgeIt->edgePtr->capacity > edgeIt->edgePtr->weights){ // ie. infinite
				
					if(edgeIt->edgePtr->tension >= edgeIt->edgePtr->lambda * edgeIt->edgePtr->weights - epsilon){ // adjust the flow and capacity downwards

						edgeIt->edgePtr->capacity = edgeIt->edgePtr->weights;
						
						if(edgeIt->edgePtr->flow > edgeIt->edgePtr->weights){ // correct the flow downward if necessary
						
							overFlow[nodeNum] += edgeIt->edgePtr->flow - edgeIt->edgePtr->weights;
							overFlow[edgeIt->to] -= edgeIt->edgePtr->flow - edgeIt->edgePtr->weights;
							edgeIt->edgePtr->flow = edgeIt->edgePtr->weights;
							edgeIt->edgePtrBack->flow -= edgeIt->edgePtr->weights; // update the other way flow
							// we stock the overflow to try to push it in other edges later
						}
					}
				}
			}
		}
	}
}

/******************************************************************
Functions needed for processing the maximum flow in the given graph.
*******************************************************************/

// if not active but excess >0 put in in the queue
void MxFlowGraph::Enqueue(int v) { 
	if ((!active[v]) && (excess[v] > 0)) {
		active[v] = true; 
		Q.push(v); 
	} 
}

// classical push
void MxFlowGraph::Push(MaxFlowEdge &e) {

	double amt = min(excess[e.from], e.edgePtr->capacity - e.edgePtr->flow); // min(excess, rest of capacity)
	// Rprintf("amt : %.14f \n" , amt);
	if (!(height[e.from] <= height[e.to] || amt == 0)){
		e.edgePtr->flow += amt; // update flow of edge e
		e.edgePtrBack->flow -= amt; // update the other way flow
		excess[e.to] += amt; // update the two excess    
		excess[e.from] -= amt;
		Enqueue(e.to); // if e was not in queue (inactive) and excess >0 put in queue
	}
}
  
// for gap heuristic 
// if there exist k for wich no vertex has this height then you can set height(u) =max(height(u), height(s)+1)
// for each vertex that has height(u)>k. This is beacause any such k is a minimum cut and
// there will not be anymore flow going from any vertex of height>k to  any vertex of height<k  
void MxFlowGraph::Gap(int k) { 
	for (int v = 0; v < N; v++) {
		if (height[v] < k) continue; // if height < k, do not consider this vertex. else:
		count[height[v]]--; // update count 
		height[v] = max(height[v], N+1); // update height 
		count[height[v]]++; // update count
		Enqueue(v); // put v in queue if needed
	}
}

// relabel
void MxFlowGraph::Relabel(int v) {

	count[height[v]]--; // one less v at this height
	height[v] = 2*N; // set height to 2*N ie > any height
	for (unsigned int i = 0; i < G[v].size(); i++) 
		if (G[v][i].edgePtr->capacity - G[v][i].edgePtr->flow > 0)
			height[v] = min(height[v], height[G[v][i].to] + 1); 
			// the height is the min of height of any vertex where we could put some more flow + 1
	count[height[v]]++; // update the count 
	Enqueue(v); // put v in queue (at the begining)
}

void MxFlowGraph::Discharge(int v) {

	for (unsigned int i = 0; (excess[v] > 0) && (i < G[v].size()); i++){
		Push(G[v][i]);
	}
	// for vertex neighboor of current v, while excess(v)>0, push it to vertex 
	if (excess[v] > 0) { // if not everything was passed
		if (count[height[v]] == 1){ // if v was the last of this height => Gap !
			Gap(height[v]); 
		}else{
			Relabel(v);
		}
    } // if everything was passed, no pb
}

// run the max flow algorithm 
void MxFlowGraph::RunMaxFlow(int s, int t) { // main that give the final max flow
	
	// preprocessing 
	count[0] = N-1; 
	count[N] = 1;
	height[s] = N;
	height[t] = 0;
	excess.assign(N,0);
	for (int i = 2; i<N; i++){
		active[i]=false;
		height[i]=0;
	}
	active[s] = active[t] = true;
	
	for (unsigned int i = 0; i < G[s].size(); i++) {
		excess[s] += G[s][i].edgePtr->capacity;
		Push(G[s][i]);
	} // that pushes all from source to neighboors.
	
	while (!Q.empty()) { // while still sth to do :
		int v = Q.front(); // which vertex to consider (in L for book)
		Q.pop(); // vire v de Q
		active[v] = false; // v not active. 
		Discharge(v); // discharge
		// print_Graph();
	}

	// print_Graph();
}

/************************************************
After a run of maxflow, functions for interpeting
************************************************/

bool MxFlowGraph::checkMax(const int sourceNode){
	
    MaxFlowNodes::iterator nodeIt;
    MaxFlowNode::iterator edgeIt;
    
    nodeIt = G.begin()+sourceNode;
    for(edgeIt = nodeIt->begin(); edgeIt != nodeIt->end(); ++edgeIt)
    {
        if(edgeIt->edgePtr->flow < edgeIt->edgePtr->capacity - epsilon)
		// if(Dif(edgeIt->edgePtr->flow , edgeIt->edgePtr->capacity) > epsilon)
        {
            return(false); // not maxed out
        }
    }
    
    return(true); // all edges are maxed out
	// nb if there is no edge linked to source : the result is 0
}

double MxFlowGraph::validUntil()
{
    MaxFlowNodes::iterator nodeIt;
    MaxFlowNode::iterator edgeIt;
    Edge *e;

    double validLambda=infinite;
    double foo, offset;
    
    // go through all edges except those connected to the source or the sink
    for(nodeIt = G.begin()+2; nodeIt != G.end(); ++nodeIt){
	
        for(edgeIt = nodeIt->begin(); edgeIt != nodeIt->end(); ++edgeIt){
		
            // check that it is not pointing to the source or sink
            if((edgeIt->to!=source) && (edgeIt->to!=sink)){
			
                if(edgeIt->edgePtr->flow > edgeIt->edgePtr->weights + epsilon){ // calculate when tension will hit lambda (which is also changing)
				
                    e = edgeIt->edgePtr;
					offset = ((e->lambda * e->weights - e->tension)/(e->flow - e->weights));
					// offset+=epsilon;
					
					if(offset < 0){ // should not happen
						e->tension = e->lambda;
						edgeIt->edgePtrBack->tension = - e->lambda;
						
					}else{
					
						foo = e->lambda + offset;
						validLambda = Min(validLambda, foo);
					}
				}
			}
		}
	}
	
    if(validLambda==infinite)
    {
        validLambda=neverSplit;
    }
    return(validLambda);
}

vector<int> MxFlowGraph::dfs(int start, bool from){

	vector<int> dist(G.size(), G.size()); // initialize the parent vector
	queue<int> next;
	int u; // the number of node that is currently being looked at
	Edge* e;
	MaxFlowNode::iterator edgeIt;
	
	dist[start]=0; 
	next.push(start); // start is the starting point
	
	while(!next.empty()){
	
		u=next.front();
		next.pop();	
		
		// checking node u
		// go through all edges in node u;
        for(edgeIt= G[u].begin(); edgeIt!=G[u].end(); ++edgeIt){
		
			if(from){ // distance from the last node; depends which direction to check capacity
		
				e = edgeIt->edgePtr;
				
			}else{ // to the last node
			
				e = edgeIt->edgePtrBack; // this way the distance when flowing towards start is calculated
			}
     
			if(e->flow < e->capacity-epsilon){ // is reachable
			// if (Dif (e->flow , e->capacity) > epsilon){
				if(dist[edgeIt->to]>dist[u]+1){ // if distance has not been updated yet
                
					dist[edgeIt->to]=dist[u]+1;
					next.push(edgeIt->to);
					
				}
			}
		}
	}
	return(dist);
}

/*********************************
Functions interface with FAgeneral
*********************************/

pair<set<int>,set<int> > MxFlowGraph::GetSets(const int sourceNode){

	set<int> reachable, complement;
	vector<int> dist = dfs(sourceNode,true); // get the distance from the source; nodes.size() means that it can't be reached

	// reachable nodes are all nodes in parent except source and sink
	for(unsigned int i=2; i !=dist.size(); ++i){
	
        if((unsigned int) dist[i] < G.size()){
			// node is reachable; insert node in external notation
			reachable.insert(nodeMapIntToExt[i]); 
		}else{
			complement.insert(nodeMapIntToExt[i]); 
		}
	}
	return(make_pair(reachable, complement));
}


double MxFlowGraph::mxFlowNew(double lambda){
	

	// first update the tensions
	updateTension(lambda);
	
	// print_Graph();
	
	// set the flow to 0 as preparation for calculating tension derivatives
	setFlowTo0();
	
	// constraint the capacity, find an optimal flow to see if the group will ever split
	setMaxCapacity();
	
	RunMaxFlow(source, sink); // call to maxflow with capacities fixed => neversplit if exists
	
	bool flowmax = checkMax(source);
	
	if(flowmax){
	
		return(neverSplit);

	}else{ // we rerun a new maxflow 
		setFlowTo0();
		setCapacity(); 
		RunMaxFlow(source, sink);
		flowmax = checkMax(source);
		if(flowmax){
			return(validUntil());
		}else{
			return(splitNow);
		}
	}
	
	// under is an heuristic from Hoefling. May be adapted later.   
    
	
	//double currentFlow, currentFactor, deltaFactor, deltaFlow, maxFlow;

    // currentFlow=currentFlowFromSource(source);
    // maxFlow = maxFlowFromSource(source);
    // currentFactor=(maxFlow-currentFlow)/currentFlow/2;// a factor of 1 raises the capacity at most by 1;
    // deltaFactor = currentFactor;
    // deltaFlow = currentFlow;
    // setCapacityProportional(currentFactor);
   

    // while(!findMaxFlow(source, sink))
    // {
        // numIter++;
        // deltaFlow = currentFlowFromSource(source) - currentFlow;
        // currentFlow += deltaFlow;
        // ow calculate the new necessary increase in the factor
        // deltaFactor = deltaFactor * (maxFlow-currentFlow)/deltaFlow;
        // currentFactor+=deltaFactor;
        // if(giveOutput)
        // {
	  // printGraph();
	  // Rprintf("Old Flow: %f\nOld Factor: %f\n", currentFlow, currentFactor);
        // }
        
        // if(deltaFlow < tolerance){ // no increase anymore; won't finish
            // have to check first if everything is ok
            // set<int> rfs = reachableFromSource();
            // if(rfs.size() ==size()) {
                // if(giveOutput) {
		  // Rprintf("There is a problem\n");
                // }
                // setFlowTo0();
                // setCapacity();
                // if(findMaxFlow(source, sink)) { // now source nodes at maximum
                    // return(validUntil(giveOutput));
                // }
                // else {
                    // return(splitNow);
                // }
            // }
            // else { // looks like needs to split
                // return(splitNow);
            // }
        // }
        
        // setCapacityProportional(currentFactor);
    // }
    // if(giveOutput)
    // {
        // printGraph();
    // }
    // return(validUntil(giveOutput));


}

double MxFlowGraph::mxFlowUpdate(double lambda){

	// update the tension in the graph
	updateTension(lambda);
	vector<double> overFlow;
    updateCapacity(lambda, overFlow);
	setFlowTo0();
	RunMaxFlow(source, sink);
	bool flowmax = checkMax(source);
	if(flowmax){
		return(validUntil());
	}else{
		return(splitNow);
	}


	// under is an heuristic from Hoefling. May be adapted later.
	
    // update the tension in the graph
    // updateTension(lambda);
    
    // update the capacity, tension update called above
    // vector<double> overFlow;
    // updateCapacity(lambda, overFlow);

	// we first test if the deficient flow can be modified
    // build the new source and sink node
    // pair<int,int> newNodes = addSpecialSourceSink(overFlow);
    // call the maximum flow
	// RunMaxFlow(newNodes.first, newNodes.second); 
    // bool flowmax = checkMax(newNodes.first);

    // remove the new source and sink node
    // removeSpecialSourceSink(overFlow, newNodes.first, newNodes.second);
    
    // if found a new valid flow, return how long it lasts
    // if(flowmax)
    // {
        // return(validUntil()); 
    // }
    // else // there is no valid maximum flow anymore; call a regular maxFlow from 0 to find out where to split
    // {
        // setFlowTo0();
        // RunMaxFlow(source, sink);
		// flowmax = checkMax(source);
		// if(flowmax){
			// return(validUntil());
		// }else{
			// return(splitNow);
		// }
    // }

} 

/****
utils
****/

    
void MxFlowGraph::print_Graph(){

	int edgeNum, nodeNum;
	MaxFlowEdge me;
    
	Rprintf("Group movement: %f. Lambda : %f \n", slope, lambda);
	for(nodeNum = 0; nodeNum != (int) G.size(); ++nodeNum){
	
        // right name for the node; internal to external mapping with special case for source and sink
		if(nodeNum==source){
			Rprintf("Node Number: Source %d\n", nodeNum);
		}else if(nodeNum==sink){
			Rprintf("Node Number: Sink %d\n", nodeNum);
		}else{
			Rprintf("Node Number: %d, %d\n", nodeNum, nodeMapIntToExt[nodeNum]);
		}
		
        // write the excess flow and the distance
		Rprintf("Excess Flow: %.14f Distance: %d\n", excess[nodeNum], height[nodeNum]);
        
		Rprintf("Edges:\n");
		for(edgeNum = 0; edgeNum != (int) G[nodeNum].size(); ++edgeNum){
		
			me = G[nodeNum][edgeNum];
			if(me.to==source){
				Rprintf("To: Source");
			}else if(me.to==sink){
				Rprintf("To: Sink");
			}else{
				Rprintf("To: %d", me.to);
            }
			if(me.edgePtr->capacity==infinite){
				Rprintf(" Cap: Infinite ");
			}else{
				Rprintf(" Cap: %.14f ",me.edgePtr->capacity);
			}
            Rprintf("Flow: %.14f Tension: %.14f Lambda: %.14f", me.edgePtr->flow, me.edgePtr->tension, me.edgePtr->lambda);
            // Rprintf("Diff: %.14f Delta: %.14f Valid: %.14f Diff: %.14f Eligible %.14f", me.edgePtr->tension - me.edgePtr->lambda, ((me.edgePtr->lambda - me.edgePtr->tension)/(me.edgePtr->flow-1)), bar, bar - me.edgePtr->lambda, (me.edgePtr->flow > 1));
			Rprintf("\n");
		}
		Rprintf("\n");
	}
	Rprintf("\n");
}


void MxFlowGraph::print_active(){
	for (int i=0; i < N; i++){
		if(!active[i]){
			Rprintf("0 ");
		}else{
			Rprintf("1 ");
		}
	}
	Rprintf("\n");
}

double Dif(double a, double b)
{
	double c = Abs(a);
	double d = Abs(b);

	d = Max(c, d);

	return d == 0.0 ? 0.0 : Abs(a - b) / d;
}
