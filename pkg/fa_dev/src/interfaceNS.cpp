#include <Rcpp.h>
#include "nosplit.h"

using namespace Rcpp ;
using namespace std  ;

SEXP tree2list(Group *g, int graphSize){
  NumericVector beta(graphSize),lambda(graphSize), slope(graphSize);
  IntegerVector idown(graphSize), iup(graphSize);
  int row=0;
  add_results(g,&beta[0],&lambda[0],&slope[0],&idown[0],&iup[0],&row);
  SEXP D = DataFrame::create(Named("beta",beta),
			     Named("lambda",lambda),
			     Named("slope",slope),
			     Named("idown",idown),
			     Named("iup",iup));
  return List::create(Named("res") = D ,Named("pred")= 0);
}

// function overloading
SEXP tree2list(Group *g, int graphSize, NumericVector lambdalist){
  NumericVector beta(graphSize),lambda(graphSize), slope(graphSize);
  IntegerVector idown(graphSize), iup(graphSize);
  vector<double> pred, lambdalistpred, slopepred;
  vector<int> idownpred, iuppred;
  
  int row=0;
  add_results_predict(g,&beta[0],&lambda[0],&slope[0],&idown[0],&iup[0],&row, &pred, &lambdalist[0],&lambdalistpred, &slopepred,&idownpred,&iuppred,lambdalist.length());
  SEXP res = DataFrame::create(Named("beta",beta),
			       Named("lambda",lambda),
			       Named("slope",slope),
			       Named("idown",idown),
			       Named("iup",iup));
  
  SEXP predict = DataFrame::create(Named("beta",pred),
				   Named("lambda",lambdalistpred),
				   Named("slope",slopepred),
				   Named("idown",idownpred),
				   Named("iup",iuppred));
  
  return List::create(Named("res") = res ,Named("pred")= predict);
}

RcppExport SEXP noSplit(SEXP R_x,SEXP R_xv,SEXP R_ngroup,SEXP R_args){
  
  NumericVector x(R_x), xv(R_xv) ;
  IntegerVector ngroup(R_ngroup) ;
  List args(R_args);
  
  std::string weights = Rcpp::as<std::string>(args["weights"]); 
  double gamma = Rcpp::as<double>(args["gamma"]); 
  bool verbose =  Rcpp::as<bool>(args["verbose"]);
  NumericMatrix W = args["W"];
  NumericVector lambdalist = args["lambdalist"];
  double epsilon = Rcpp::as<double>(args["epsilon"]);
	
  if (verbose){
    Rprintf("Starting on %i groups \n" , x.length());
    Rprintf("gamma = %f \n" , gamma);
    Rprintf("Calculating Slopes");
  }  
  
  vector<double> slopes = calculateSlope(x,ngroup,xv,weights,gamma,W);
  
  if(verbose){
    Rprintf("Fusing groups \n");
  } 
  
  Group *G = maketree(x, slopes, ngroup);
  
  SEXP L ;
  if (lambdalist.length()==0){
    if(verbose){
      Rprintf("Creating result Dataframe.");
    }
    L = tree2list(G,2*x.size()-1);  // tree conversion 
  } else {
    if(verbose){
      Rprintf("Creating result and prediction Dataframe.");
    }
    L = tree2list(G,2*x.size()-1,lambdalist);
  }
  
  delete_tree(G); // delete the old tree for memory
  
  return L ;
}

// we process each dimension individually using this function
RcppExport SEXP noSplitcv(SEXP R_x,SEXP R_xv,SEXP R_ngroup, SEXP R_xtest,SEXP R_ngrouptest ,SEXP R_args){
  
  NumericVector x(R_x), xv(R_xv), xtest(R_xtest);
  IntegerVector ngroup(R_ngroup),  ngrouptest(R_ngrouptest);
  List args(R_args);
  
  std::string weights = Rcpp::as<std::string>(args["weights"]); 
  double gamma = Rcpp::as<double>(args["gamma"]); 
  double epsilon = Rcpp::as<double>(args["epsilon"]);
  NumericMatrix W = args["W"];
  NumericVector lambdalist = args["lambdalist"];
  
  NumericVector error(lambdalist.length());
  
  vector<double> slopes = calculateSlope(x,ngroup,xv,weights,gamma,W);
  
  Group *G = maketree(x, slopes, ngroup);
  
  error_cv(G,&lambdalist[0],lambdalist.length(),&xtest[0], &ngrouptest[0],&error[0]);
  
  delete_tree(G); 
  
  return(error);
}

