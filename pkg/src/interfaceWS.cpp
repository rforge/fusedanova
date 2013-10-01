#include <Rcpp.h>
#include "withsplit.h"

using namespace Rcpp ;


SEXP tree2list(Group *g,int graphSize){
	NumericVector beta(graphSize),lambda(graphSize), slope(graphSize);
	IntegerVector index(graphSize);
	int row =0;
	add_results(g,&beta[0],&lambda[0],&slope[0],&index[0],&row);
	SEXP D = DataFrame::create(Named("beta",beta),
				Named("lambda",lambda),
				Named("slope",slope),
				Named("class",index));
	return List::create(Named("res") = D ,Named("pred")= 0);
}

// function overloading
SEXP tree2list(Group *g, int graphSize, NumericVector lambdalist){

	NumericVector beta(graphSize),lambda(graphSize), slope(graphSize);
	IntegerVector index(graphSize);
	vector<double> pred, lambdalistpred, slopepred;
	vector<int> indexpred;

	int row=0;
	add_results_predict(g,&beta[0],&lambda[0],&slope[0],&index[0],&row, &pred, &lambdalist[0],&lambdalistpred, &slopepred,&indexpred,lambdalist.length());
	SEXP res = DataFrame::create(Named("beta",beta),
				Named("lambda",lambda),
				Named("slope",slope),
				Named("class",index));
 
	SEXP predict = DataFrame::create(Named("beta",pred),
				Named("lambda",lambdalistpred),
				Named("slope",slopepred),
				Named("class",indexpred));

	return List::create(Named("res") = res ,Named("pred")= predict);
}


// we process each dimension individually using this function
RcppExport SEXP withSplit(SEXP R_x,SEXP R_xv,SEXP R_ngroup,SEXP R_args){

	NumericVector x(R_x);
	NumericVector xv(R_xv);
	NumericVector ngroup(R_ngroup);
	List args(R_args);
  
	std::string weights = Rcpp::as<std::string>(args["weights"]); 
	double gamma = Rcpp::as<double>(args["gamma"]); 
	NumericMatrix W = args["W"];
	NumericVector lambdalist = args["lambdalist"];
	bool verbose =  Rcpp::as<bool>(args["verbose"]);
	int mxSize = Rcpp::as<int>(args["mxSplitSize"]);
	double epsilon = Rcpp::as<double>(args["epsilon"]);

	if (verbose){
		Rprintf("Starting on %i groups with no mxflowcheck if size >= %i \n" , x.length(), mxSize);
		Rprintf("epsilon = %.14f \n" , epsilon);
		Rprintf("gamma = %.14f \n" , gamma);
	}
	
	SEXP L=0;

	FAGeneral FA(&x[0], x.length(), &ngroup[0], xv, weights, gamma, W, mxSize, verbose, epsilon);
  
	Group *G = FA.GetFinalGroup();
	int graphSize = FA.graphsize;
	
	if(verbose){
		Rprintf("Recovering solution from tree\n");
	}
	
	if (lambdalist.length()==0){
		if(verbose){
			Rprintf("Creating result Dataframe. \n");
		}
		L = tree2list(G,graphSize);  // tree conversion 
	}else{
		if(verbose){
			Rprintf("Creating result and prediction Dataframe. \n");
		}
		L = tree2list(G,graphSize,lambdalist);
	}

	if(verbose){
		Rprintf("Deleting the tree \n");
	}
	
	FA.delete_tree(G); // delete the old tree for memory

	return L;
	
	return L;
}


// we process each dimension individually using this function
RcppExport SEXP Splitcv(SEXP R_x,SEXP R_xv,SEXP R_ngroup, SEXP R_xtest,SEXP R_ngrouptest ,SEXP R_args){

	NumericVector x(R_x);
	NumericVector xv(R_xv);
	NumericVector xtest(R_xtest);
	NumericVector ngroup(R_ngroup);
	NumericVector ngrouptest(R_ngrouptest);
	List args(R_args);
  
	std::string weights = Rcpp::as<std::string>(args["weights"]); 
	double gamma = Rcpp::as<double>(args["gamma"]); 
	NumericMatrix W = args["W"];
	NumericVector lambdalist = args["lambdalist"];
	bool verbose =  Rcpp::as<bool>(args["verbose"]);
	int mxSize = Rcpp::as<int>(args["mxSplitSize"]);
	double epsilon = Rcpp::as<double>(args["epsilon"]);

	NumericVector error(lambdalist.length());
	
	if (verbose){
		Rprintf("Starting on %i groups with no mxflowcheck if size >= %i \n" , x.length(), mxSize);
		Rprintf("epsilon = %f \n" , epsilon);
		Rprintf("gamma = %f \n" , gamma);
	}

	FAGeneral FA(&x[0], x.length(), &ngroup[0], xv, weights, gamma, W, mxSize, verbose, epsilon);
  
	Group *G = FA.GetFinalGroup();
	
	if(verbose){
		Rprintf("Calculating cv error. \n");
	}
	
	error_cv_ws(G,&lambdalist[0],lambdalist.length(),&xtest[0], &ngrouptest[0],&error[0]);

	if(verbose){
		Rprintf("Deleting the tree \n");
	}
	
	FA.delete_tree(G); // delete the old tree for memory

	return(error);
}

