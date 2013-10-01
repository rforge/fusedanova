#include <list>
#include <vector>
#include <map>
#include <Rcpp.h>

using namespace std ;
using namespace Rcpp ;

# define NullGroup ((Group*)0)
# define NullEvent ((Events::iterator)0)

struct Group;

// Doubly-linked list of clusters
typedef std::list<Group*> Groups;

// events are represented as a map, thus e->first is lambda and
// e->second is the pointer to the top group of the merge event.
typedef std::multimap<double,Groups::iterator> Events;
typedef std::pair<double,Groups::iterator> Event;

struct Group {

	int len; // number of classes/individuals in the group 

	int idown; // index of first in the group (meaning smaller coef Beta at lambda = 0).  

	double slope; // the slope of the group
	double beta; // value of the coef at the generation of this group by fusing of the children
	double lambda; // value of the penalty parameter for which the fuse happens
	double nbind; // nb of individuals in group

	Group *childup;  
	Group *childdown; // 2 group fusing with Beta(childup) > Beta(childdown)
	Group *father; 

	Events::iterator nextFusionUp; // enable the easy supression of this event 
};

// tree construction function
Group* maketree(double*,int, double*,double*,double); // general function
vector<double> calculateSlope(NumericVector, NumericVector,NumericVector,std::string,double,NumericMatrix,int); // initialize the slopes
Group* fuse_groups(Groups groups); // iteratively construct the tree of fusing groups
double getlambda(Group *g1,Group *g2); // calculate the fusing lambda for 2 groups
void insert_new_event(Events*,Group*,Group*,Groups::iterator,Group*); // calculate and insert events in the events list
void delete_tree(Group*); // remove the tree for memory

// fused anova functions : append the dataframe(s) result. 
void add_results(Group*,double*,double*,double*,int*,int*,int*);
void add_results_predict(Group*,double*, double*, double*, int*, int*, int*, vector<double>*, double*,vector<double>*, vector<double>*, vector<int>*, vector<int>*,int);
               
// cross validation functions
void error_cv(Group*, double*, int, double*, double*,double*);
void updateError(double*,double*,double*,double*,int,Group*,Group*);

// utilities for self checking
void print_groups(Groups); 
void print_events(Events);
