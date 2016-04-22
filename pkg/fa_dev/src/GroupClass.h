#ifndef _GROUP_CLASS_H
#define _GROUP_CLASS_H

#include <list>
#include <vector>
#include <map>

# define NullGroup  ((Group*)0)
# define NullFusion ((FusionEvents::iterator)0)

class Group; // needed to declare new typedef

// List of Groups/Cluster
typedef std::list<Group*> GroupList; 

// Each Fusion is pair of two elements
// - f->first  is the value of lambda associated with the fusion
// - f->second is the pointer to the top group of the fusion event
typedef std::pair<double,GroupList::iterator> FusionEvent;

// The set of fusion is a map of fusion events
typedef std::multimap<double,GroupList::iterator> FusionEvents;

// Group are the atomic element of the GroupList
class Group {
 public:

  int nbind        ; // nb of individuals in group
  int len          ; // number of intial classes in the group 
  int idown        ; // index of first in the group (smaller beta at lambda = 0).  
  double slope     ; // the slope of the group
  double beta      ; // value of coefficients at startup by fusion its children
  double lambda    ; // value of the penalty parameter for which the fusion happens

  Group *childup   ; 
  Group *childdown ; // 2 groups fusing with beta(childup) > beta(childdown)
  Group *father    ;  
  
  // contains the lambda for which the group is fusing with the one on
  // top of him and an iterator to indicate the group's position
  FusionEvents::iterator nextFusionUp;
  
  // Constructors
  Group();
  Group(double beta, double lambda, double slope, int nbind, int idown);
  // Destructors
  ~Group();
  
};

#endif        
