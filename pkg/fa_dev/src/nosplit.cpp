#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include "nosplit.h"

using namespace std  ;
using namespace Rcpp ;

# define square(x) ((x)*(x))
# define NEW_EVENT_THRESH 1e-16

vector<double> calculateSlope(NumericVector &xm    ,
			      IntegerVector &ngroup,
			      NumericVector &xv    , 
			      string weights       ,
			      double gamma         ,
			      NumericMatrix &W     ) {
  
  int n = xm.length()  ; 
  int N = sum(ngroup)  ; // number of individuals in all groups
  int sign             ; // some auxiliary variables
  double sum1=0,sum2=0 ;
  vector<double> slopes; // the vector of slopes (output of the function)

  if (weights=="default"){ // easiest
    
    if (N==n){ //slopes are known (clusterpath case)
      for (int i =0; i<n; i++){
	slopes.push_back(n-1-2*i);
      }
    } else { // nk.nl weights (O(n) calculation)
      for (int i=1;i<n;i++){
	sum1 += ngroup[i];
      }
      slopes.push_back(sum1);
      for (int i =0; i<(n-1);i++){
	slopes.push_back(slopes[i]-ngroup[i]-ngroup[i+1]);
      }
    } 
  } else if (weights == "laplace"){ // O(n) calculation
    sum1 = 0;
    sum2 = 0;
    vector<double> slopes1(n),slopes2(n);
    
    for (int i=n-1;i>0;i--){
      sum1 += ngroup[i]*exp(-gamma * xm[i]);
      slopes1[i-1] = sum1;
    }
    slopes1[n-1] =0;
    
    for (int i=1;i<n;i++){
      sum2 += ngroup[i-1]*exp(gamma * xm[i-1]);
      slopes2[i] = sum2;
    }
    slopes2[0] =0;	
    
    for (int i=0;i<n;i++){
      slopes.push_back(slopes1[i]*exp(gamma * xm[i])- slopes2[i]*exp(-gamma * xm[i]));
    }		
    
    // O(n2) calculation for all the rest
    // gaussian and adaptive trick : summing starting by the smallest element for each line
    // the smallest are the one away from the "diagonal"
    // create two vectors of sum and then add them
  } else if (weights == "gaussian"){
    // calculate the positive part
    for(int i=0;i<(n-1);i++){
      sum1=0;
      for(int j =n-1;j>i;j--){
	sum1 += ngroup[j] * exp(-gamma* square(xm[i]-xm[j]));
      }
      slopes.push_back(sum1);
    }
    slopes.push_back(0);
    // calculate the negative part
    
    for(int i=1;i<n;i++){
      sum2=0;
      for (int j = 0;j<i ;j++){
	sum2 += ngroup[j] * exp(-gamma* square(xm[i]-xm[j]));
      }
      slopes[i]-=sum2;
    }
    
  } else if (weights == "adaptive"){
    // calculate the positive part
    for(int i=0;i<(n-1);i++){
      sum1=0;
      for(int j =n-1;j>i;j--){
	sum1 += ngroup[j] /pow(fabs(xm[i]-xm[j]),gamma);
      }
      slopes.push_back(sum1);
    }
    slopes.push_back(0);
    // calculate the negative part
    
    for(int i=1;i<n;i++){
      sum2=0;
      for (int j = 0;j<i ;j++){
	sum2 += ngroup[j] /pow(fabs(xm[i]-xm[j]),gamma);
      }
      slopes[i]-=sum2;
    }

  } else if (weights == "naivettest"){ 
    for (int i=0; i<n; i++){
      sum1 = 0 ;
      sign = -1 ;
      for (int j=0; j<n; j++){
	if (i==j){
	  sign = 1 ;
	}else{
	  sum1+= sign / sqrt(1/ngroup[i]+1/ngroup[j]) ;
	}
      } 
      slopes.push_back(sum1/ ngroup[i]);
    } 	
  } else if (weights == "ttest"){
    for (int i=0; i<n; i++){
      sum1 = 0;
      sign = -1;
      for (int j=0; j<n; j++){
	if (i==j){
	  sign= 1;
	}else{
	  sum1+= sign * sqrt(((ngroup[i]-1)*xv[i] + (ngroup[j]-1)*xv[j])/(ngroup[i]+ngroup[j]-2)) /sqrt(1/ngroup[i]+1/ngroup[j]) ;
	}
      } 
      slopes.push_back(sum1/ ngroup[i]);
    } 	
  } else if (weights == "welch"){
    for (int i=0; i<n; i++){
      sum1 = 0;
      sign = -1 ;
      for (int j=0; j<n; j++){
	if (i==j){
	  sign= 1 ;
	}else{
	  sum1+= sign * sqrt(xv[i]/ngroup[i]+xv[j]/ngroup[j])/(1/ngroup[i]+1/ngroup[j]) ;
	}
      } 
      slopes.push_back(sum1/ngroup[i]);
    } 	
  } else if (weights == "personal"){
    for (int i=0; i<n; i++){
      sum1 = 0;
      sign = -1 ;
      for (int j=0; j<n; j++){
	if (i==j){
	  sign= 1 ;
	}else{
	  sum1 += sign * W(i,j);
	}
      } 
      slopes.push_back(sum1/ngroup[i]);
    } 	
  }
  
  return(slopes);
}

// ==================================================================
// 
// GENERAL FUNCTION TO CONSTRUC THE FUSED-ANOVA TREE
//
// ==================================================================
Group* maketree(NumericVector &x ,vector<double> &slopes, IntegerVector &ngroup){

  // Some interators eneded for getting across the tree
  GroupList::iterator current, next, previous ; 
  
  // ==================================================================
  // CREATE THE INITIAL LIST OF GROUPS
  GroupList myGroups;
  for (int i=0; i<x.size(); i++){
    Group *g = new Group;
    g->beta      = x[i];
    g->lambda    = 0.0;
    g->slope     = slopes[i];
    g->nbind     = ngroup[i];
    g->idown     = i+1;
    g->len       = 1;    
    g->childup   = NullGroup;  
    g->childdown = NullGroup; 
    g->father    = NullGroup;
    // append my list of groups
    myGroups.push_back(g); 
  } // no need to sort because groups are already sorted by ascending beta

  // ==================================================================
  // PREFUSION STEP
  // 2 groups having the same alpha at startup must fuse.
  current = myGroups.begin() ;
  next    = myGroups.begin() ; next++ ;
  while(next != myGroups.end()){
    
    // Merge two groups if they have the same initial coefficients
    if((*current)->beta == (*next)->beta){ 
      
      // Create the new group and set all his attributes
      Group *g  = new Group;
      g->len   = (*current)->len   + (*next)->len  ;
      g->nbind = (*current)->nbind + (*next)->nbind;
      g->idown = (*current)->idown;  // (*current) had smaller beta
      // the new slope is just the weighted mean of the old slopes
      g->slope = ((*current)->nbind * (*current)->slope + (*next)->nbind * (*next)->slope)/(g->nbind); 
      g->beta  = (*current)->beta;
      g->lambda= 0.0;
      g->childdown = *current;
      g->childup   = *next;
      g->father    = NullGroup;

      // set father - TODO : Understand why this can have some influence...
      (*current)->father = g; 
      (*next)->father    = g;
      
      // erase ther current group and replace 
      current  = myGroups.erase(current); // erase the group in current position
      *current = g; // replace value of the second group
      // move to the next group
      next++;     
    } else { // move iterators 
      next++ ; current++;
    }
  }
  
  // ==================================================================
  // GROUP FUSION PROCESS

  // ------------------------------------------------
  // recover the first list of events for each groups
  double lambda ;
  FusionEvents events ;
  current = myGroups.begin();
  next    = myGroups.begin(); next++;
  for(; next != myGroups.end(); current++, next++){
    // calculate initial events list from intersection of all adjacent
    lambda = getlambda(*current, *next);
    if(lambda>0){ 
      (*current)->nextFusionUp = events.insert(FusionEvent(lambda,current));
    }
  }

  // set to 0 here to avoid compiler warnings.
  FusionEvents::iterator my_fusion;
  Group *g=0;
  while(events.size()>0){
    do{
      my_fusion=events.begin(); // my_fusion is automatically the one with lambda min
      lambda = my_fusion->first;
      // current and next are the 2 fusing groups
      current = next = my_fusion->second;
      next    = my_fusion->second; next++;
      // previous is the predecessor of the current fusing group
      previous=my_fusion->second; previous--;
      
      // Create the new group
      g = new Group;
      // and set all hist attributes
      g->len = (*current)->len + (*next)->len;
      g->nbind = (*current)->nbind + (*next)->nbind;
      g->idown = (*current)->idown; // because (*current) had smaller beta
      // new slope is just the ponderate mean of old slopes
      g->slope=((*current)->nbind * (*current)->slope + (*next)->nbind * (*next)->slope)/(g->nbind);
      g->lambda=lambda;
      g->beta = (*current)->beta + (g->lambda - (*current)->lambda)*(*current)->slope; 
      // ie. Beta = old Beta + Delta(lambda)* old slope
      g->childdown= *current;
      g->childup  = *next   ;
      g->father = NullGroup;
      // set father
      (*current)->father = g; // is this useful??? will be erased...
      (*next)->father    = g; 
      
      //set up pointers to neighboring groups
      Group *gtmp = *next;
      next++;      
      
      // We have previous-(*current)-gtmp-next
      // we remove current from the list (erase) and set gtmp=current 
      current = myGroups.erase(current); *current = g;
      
      if(current != myGroups.begin()){ //if the new group != group with min Beta
	insert_new_event(&events,g,*previous,previous,*previous);
      }
      if(next != myGroups.end()){ //if the new group != group with max Beta
	insert_new_event(&events,g,*next,current,gtmp);
      } 
      events.erase(my_fusion);
    } while(events.begin()->first == lambda);
  }
  
  return g;
}

//calculate lambda where g1 and g2 fuse
double getlambda(Group *g1,Group *g2){
  return((g1->beta - g2->beta - g1->slope * g1->lambda + g2->slope * g2->lambda)/(g2->slope-g1->slope));
}

// used for adding/deleting previous and next joining events
void insert_new_event(FusionEvents *eventlist,Group *newg,
		      Group *neighbor_group,GroupList::iterator g_it_up,
		      Group *g_down){
  
  if(g_down-> nextFusionUp != NullFusion){
    (*eventlist).erase(g_down->nextFusionUp); // erase the event that will not happen
  }
  Group *g_up = *g_it_up;
  FusionEvents::iterator new_event;
  FusionEvent e=FusionEvent(getlambda(neighbor_group,newg),g_it_up);
  
  if(newg->lambda -e.first <= NEW_EVENT_THRESH){
    // 2 things : new lambda >= old lambda and add a treshold to prevent round up errors
    new_event=(*eventlist).insert(e);
  }else{
    new_event=NullFusion;
  }
  
  g_up->nextFusionUp = new_event; // contains the lambda for which the group is fusing the one on top of him
}

void print_groups(GroupList myGroups){
  GroupList::const_iterator g_it;
  Rprintf("idown iup beta lambda slope \n");
  for(g_it = myGroups.begin(); g_it != myGroups.end();++g_it){
    Rprintf("%i %i %f %f %f \n",(*g_it)->idown,(int)((*g_it)->idown+(*g_it)->len-1),(*g_it)->beta,(*g_it)->lambda,(*g_it)->slope);
  }
  Rprintf("\n");
}

void print_events(FusionEvents events){
  FusionEvents::const_iterator it;
  for(it = events.begin(); it != events.end();++it){
    Rprintf("%f ",it->first);
    if(it->second == (GroupList::iterator)0) Rprintf("NULL \n");
    else Rprintf("%f \n",(*(it->second))->beta);
  }
  Rprintf("size: %f \n",events.size());
}

void delete_tree(Group *g){
  Group *g1=g->childdown,*g2=g->childup;
  delete g;
  if(g1 != NullGroup)delete_tree(g1);
  if(g2 != NullGroup)delete_tree(g2);
}

void add_results(Group *g,double *beta,
		 double *lambda, double *slope, int *idown,int *iup,int *row){
  
  beta[*row] = g->beta;
  lambda[*row] = g->lambda;
  slope[*row] = g->slope;
  idown[*row] = g->idown;
  iup[*row] = g->idown + g->len -1;
  
  *row +=1;
  
  if(g->childup != NullGroup){ // keep writing
    add_results(g->childdown,beta,lambda,slope,idown,iup,row);
    add_results(g->childup,beta,lambda,slope,idown,iup,row);
  }
}


void add_results_predict(Group *g,double *beta, double *lambda, double *slope, int *idown,int *iup,int *row, vector<double> *pred, double *lambdalist, vector<double> *lambdalistpred, vector<double> *slopepred, vector<int> *idownpred, vector<int> *iuppred,int nlambda){
  
  int k =0;
  
  beta[*row] = g->beta;
  lambda[*row] = g->lambda;
  slope[*row] = g->slope;
  idown[*row] = g->idown;
  iup[*row] = g->idown + g->len -1; 
  
  if (g->father!= NullGroup){ // we have the lambda !
    
    double l = g->father->lambda;
    while (lambdalist[k]>l){ // moove to the right position in the list
      k =k+1;      
      if (k==(nlambda))
	break; 
    }
    
    if (k<(nlambda)){
      while (lambdalist[k]>=g->lambda){ // ok position
	pred->push_back(g->beta + (lambdalist[k]-g->lambda)*g->slope);
	lambdalistpred->push_back(lambdalist[k]);
	slopepred->push_back( g->slope);
	idownpred->push_back(g->idown);
	iuppred->push_back( g->idown + g->len -1);
	k =k+1;      
	if (k==(nlambda))
	  break;
      }
    }
  }else{ // final group : predict further for all 
    while (lambdalist[k]>=g->lambda){ // ok position
      pred->push_back(g->beta + (lambdalist[k]-g->lambda)*g->slope);
      lambdalistpred->push_back(lambdalist[k]);
      slopepred->push_back( g->slope);
      idownpred->push_back(g->idown);
      iuppred->push_back( g->idown + g->len -1);
      k =k+1;
    }
  }
  
  *row +=1;
  
  if(g->childup != NullGroup){ // keep writing
    add_results_predict(g->childdown,beta,lambda,slope,idown,iup,row,pred,lambdalist,lambdalistpred,slopepred,idownpred,iuppred, nlambda);
    add_results_predict(g->childup,beta,lambda,slope,idown,iup,row,pred,lambdalist,lambdalistpred,slopepred,idownpred,iuppred, nlambda);
  }
}


void error_cv(Group *g, double* lambdalist,int lambdalen, double* xtest,int* ngroup, double* error){
  
  int index=0;
  
  if (g->father!= NullGroup){ // we have the lambda
    
    updateError(error, xtest, ngroup, lambdalist,lambdalen, g->father, g);
    
  } else { // final group : calc further err 
    
    // it is over so we can update the end of the lambdalist (case where lambdalist[i] >= final_g->lambda)
    while(lambdalist[index]<g->lambda){
      index+=1;
      if (index==lambdalen){break;}
    } 
    
    if (index!=lambdalen){
      for(int j=(g->idown-1);j<(g->idown -1 +g->len);j++){ // for all starting groups :
	error[index] = error[index] + ngroup[j]* square(g->beta - xtest[j]);
      }// run once, g->slope=0. 
      index +=1;
    }
    while(index<lambdalen){ // now the slope is 0, err[index-1] is touch only once => the rest is constant
      error[index] = error[index-1];
      index+=1;
    }
  }
  
  if(g->childup != NullGroup){ // keep writing
    error_cv(g->childdown, lambdalist, lambdalen, xtest, ngroup,error);
    error_cv(g->childup, lambdalist, lambdalen, xtest, ngroup,error);
  }
}

// allow to add error terms to the err vector for cross validation
void updateError(double* err,double* xtest,int* ngroup,double* lambdalist,int lambdalen, Group* gnew,Group*gold){
  
  int i =0;
  double oldl = gold->lambda;
  double newl = gnew->lambda;
  if (oldl!=newl){ // can be the same for multiple fuse or prefuse : in this case don't calculate
    while(lambdalist[i]<oldl){ // moove i to first interesting lambda value
      i+=1;
      if (i==lambdalen){break;}
    } // lambdalist[i]>=oldl
    
    if (i<lambdalen){
      while(lambdalist[i]<newl){ // we predict on this slope
	for(int j=(gold->idown-1);j < (gold->idown -1 +gold->len);j++){ // for all starting groups :
	  err[i] += ngroup[j]* square(gold->beta + (lambdalist[i]-oldl)*gold->slope - xtest[j]);
	}
	i+=1;
	if (i==lambdalen){break;}
      }
    }
  }
}
