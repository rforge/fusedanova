#include <Rcpp.h>
//#include <RcppArmadillo.h>
#include "nosplit.h"

using namespace Rcpp ;
//using namespace arma;

SEXP tree2list(Group *g){
	int nrow = 2*g->len -1; // nb infos = K init gpes + K-1 fusions 
	NumericVector beta(nrow),lambda(nrow), slope(nrow);
	IntegerVector idown(nrow), iup(nrow);
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
SEXP tree2list(Group *g, NumericVector lambdalist){
	int nrow = 2*g->len -1; // nb infos = K init gpes + K-1 fusions
	NumericVector beta(nrow),lambda(nrow), slope(nrow);
	IntegerVector idown(nrow), iup(nrow);
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


// we process each dimension individually using this function
RcppExport SEXP noSplit(SEXP R_x,SEXP R_xv,SEXP R_ngroup,SEXP R_args){

	NumericVector x(R_x);
	NumericVector xv(R_xv);
	NumericVector ngroup(R_ngroup);
	List args(R_args);
  
	std::string weights = Rcpp::as<std::string>(args["weights"]); 
	double gamma = Rcpp::as<double>(args["gamma"]); 
	bool verbose =  Rcpp::as<bool>(args["verbose"]);
	NumericMatrix W = args["W"];
	NumericVector lambdalist = args["lambdalist"];
	double epsilon = Rcpp::as<double>(args["epsilon"]);
	
	SEXP L;
  
	if (verbose){
		Rprintf("Starting on %i groups \n" , x.length());
		Rprintf("gamma = %f \n" , gamma);
		Rprintf("Calculating Slopes");
	}  
  
	vector<double> sl = calculateSlope(x,ngroup,xv,weights,gamma,W,x.length());

	if(verbose){
		Rprintf("Fusing groups \n");
	} 
	
	Group *G = maketree(&x[0], x.length(), &sl[0],&ngroup[0],epsilon);

	if (lambdalist.length()==0){
		if(verbose){
			Rprintf("Creating result Dataframe.");
		}
		L = tree2list(G);  // tree conversion 
	}else{
		if(verbose){
			Rprintf("Creating result and prediction Dataframe.");
		}
		L = tree2list(G,lambdalist);
	}

	delete_tree(G); // delete the old tree for memory

	return L;
}


// we process each dimension individually using this function
RcppExport SEXP noSplitcv(SEXP R_x,SEXP R_xv,SEXP R_ngroup, SEXP R_xtest,SEXP R_ngrouptest ,SEXP R_args){

	NumericVector x(R_x);
	NumericVector xv(R_xv);
	NumericVector xtest(R_xtest);
	NumericVector ngroup(R_ngroup);
	NumericVector ngrouptest(R_ngrouptest);
	List args(R_args);
  
	std::string weights = Rcpp::as<std::string>(args["weights"]); 
	double gamma = Rcpp::as<double>(args["gamma"]); 
	double epsilon = Rcpp::as<double>(args["epsilon"]);
	NumericMatrix W = args["W"];
	NumericVector lambdalist = args["lambdalist"];

	NumericVector error(lambdalist.length());

	vector<double> sl = calculateSlope(x,ngroup,xv,weights,gamma,W,x.length());
 
	Group *G = maketree(&x[0], x.length(), &sl[0],&ngroup[0],epsilon);

	error_cv(G,&lambdalist[0],lambdalist.length(),&xtest[0], &ngrouptest[0],&error[0]);

	delete_tree(G); 

	return(error);
}

