#ifndef _WITHSPLIT_H
#define _WITHSPLIT_H

#include <list>
#include <vector>
#include <map>
#include <Rcpp.h> 

#include "MainGraph.h"
#include "MxFlowGraph.h"
#include "Groups.h"

using namespace std;
using namespace Rcpp ;

class FAGeneral{

    Groups myGroups; // list of current groups
    MainGraph myGraph; // a directed clique with capacities for calculating maximum flow. Each Group using subgraph will point to this graph
    Events events;
	int K; // nb of original groups 
	
	vector< vector<double> > Weights; // weight matrix
	vector<double> sl; // the slope of each group at the beginning
	// these three attributes are stocked here for the moment, we'll see later if we can moove them
	
	bool verbose; // print information ?
    int mxSplitSize; // if a group has size > mxSplitSize, it is not check for Size
	double epsilon; 

	// calculate the weights we need 
	void calculateSlope(double*, double*,NumericVector,std::string,double,NumericMatrix);
	void calculateWeights(double*, double*,NumericVector,std::string,double,NumericMatrix); // for the moment basic, merge with calculate slope in some way ? 
	
	// set groups and prefuse if needed
    void initializeGroups(double*,double*);
	
	// initialize events
	void initializeEvents();
    
	// run the algorithm until there is only one group left
    void run();

	double getlambda(Group *g1,Group *g2); // calculate the fusing lambda for 2 groups
	
	void FuseGroups(Events::iterator e_it);
	
	void split(double lambda, Groups::iterator g);
	
	void MxFlowCheck(double lambda, Groups::iterator g, bool newgroup);
	
	void insert_new_merge_event(Group*,Group*,Groups::iterator,Group*); // calculate and insert events in the events list
	
	// cross validation functions
	//void error_cv(Group*, double*, int, double*, double*,double*);
	//void updateError(double*,double*,double*,double*,int,Group*,Group*);
	
	// utilities for self checking
	//void print_groups(Groups); 
	//void print_events(Events);
	
public:

    // constructor 
    FAGeneral(double* x, int K, double* ngroup, NumericVector xv, std::string weights, double gamma, NumericMatrix W, int mxSize, bool verb, double eps);
				
	// destructor  
    ~FAGeneral(){};
	
	Group* GetFinalGroup();
	
	int graphsize;
	
	// utilities for self checking
	void print_groups(); 
	void print_events();
	void print_weights();
	void print_content(Group* g);
	
	void delete_tree(Group*); // remove the tree for memory
	
};

// fused anova functions : append the dataframe(s) result. 
void add_results(Group*,double*,double*,double*,int*,int*);
void add_results_predict(Group*,double*, double*, double*, int*, int*, vector<double>*, double*,vector<double>*, vector<double>*, vector<int>*,int);

// cross validation functions
void error_cv_ws(Group*, double*, int, double*, double*,double*);
void updateError_ws(double*,double*,double*,double*,int,Group*,Group*);

#endif        



