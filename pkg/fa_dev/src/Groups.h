// Groups structure
#ifndef _GROUP_H
#define _GROUP_H

#include "MxFlowGraph.h"

# define NullGroup ((Group*)0)
# define NullEvent ((Events::iterator)0)

struct Group;
// Doubly-linked list of clusters
typedef std::list<Group*> Groups; 

// events are represented as a map, thus e->first is lambda and
// e->second is an action containing a char : "F" for fusion and "M" for Maxflow test.
// e->second->first is a group iterator, which is the down group if there is a fusion else the group to test the maxflow 
typedef pair<Groups::iterator,char> Action; 
typedef multimap<double,Action> Events;
typedef pair<double, Action> Event;

struct Group {

  int len; // number of classes in the group 
  int nbind; // number of total individuals
  
  set<int> grps; // index of singleton groups in the group
  
  map<int, Group *> content; // usefull to reconstruct groups from split 	 
  
  double slope; // the slope of the group
  double beta; // value of the coef at the generation of this group by fusing of the children
  double lambda; // value of the penalty parameter for which the fuse happens
  
  Group *childup;  
  Group *childdown; // 2 group fusing with Beta(childup) > Beta(childdown)
  Group *fatherup; // only one father if no split, but if split occurence => 2 fathers
  Group *fatherdown;	
  
  MxFlowGraph* m; // maxflowgraph that contains the tension information
  
  Events::iterator nextFusionUp; // enable the easy supression of this event. 
  // contains the lambda for which the group is fusing the one on top of him and an iterator to indicate the group's position 
  Events::iterator nextMxFlow; // same but for split : the lambda for wich it need a check up and the iterator
};

#endif

