#ifndef _NOSPLIT_H
#define _NOSPLIT_H

#include <list>
#include <vector>
#include <map>
#include <Rcpp.h>
#include "GroupClass.h"

// initialize the slopes
std::vector<double> calculateSlope(Rcpp::NumericVector&,  
				   Rcpp::IntegerVector&,  
				   Rcpp::NumericVector&,
				   std::string, double, Rcpp::NumericMatrix&); 

// general to iteratively construct the tree of fusing groups
Group* maketree(Rcpp::NumericVector&,  std::vector<double>&, Rcpp::IntegerVector&);

// calculate the fusing lambda for 2 groups
double getlambda(Group *g1,Group *g2); 

// calculate and insert fusions in the fusions list
void insert_new_event(FusionEvents*, Group*, Group*, GroupList::iterator, Group*); 

void delete_tree(Group*); // remove the tree for memory

// fused anova functions : append the dataframe(s) result. 
void add_results(Group*,double*,double*,double*,int*,int*,int*);

void add_results_predict(Group*,
			 double*, 
			 double*, 
			 double*, 
			 int*, int*, int*, 
			 std::vector<double>*,  double*,
			 std::vector<double>*, 
			 std::vector<double>*, 
			 std::vector<int>*, 
			 std::vector<int>*,int);
               
// cross validation functions
void error_cv(Group*, double*, int, double*, int*,double*);

void updateError(double*,double*,int*,double*,int,Group*,Group*);

// utilities for self checking
void print_groups(GroupList); 

void print_events(FusionEvents);

#endif        

