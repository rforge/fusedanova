#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include "nosplit.h"

#define square(x) ((x)*(x))
# define NEW_EVENT_THRESH 1e-16
double epsilon = NEW_EVENT_THRESH; // needed because of the C++ round up

vector<double> calculateSlope(NumericVector xm,NumericVector ngroup,NumericVector xv, 
	std::string weights, double gamma, NumericMatrix W, int n){

	vector<double> sl;
	double sign;
	double sum1=0,sum2=0;
	double N = std::accumulate(ngroup.begin(), ngroup.end(), double()); // number of in dividuals in all groups

	if (weights=="default"){ // easiest

		if (N==n){ //slopes are known (clusterpath case)
			for (int i =0; i<n; i++){
				sl.push_back(n-1-2*i);
			}
		}else{ // nk.nl weights (O(n) calculation)
			for (int i=1;i<n;i++){
				sum1 += ngroup[i];
			}
			sl.push_back(sum1);
			for (int i =0; i<(n-1);i++){
				sl.push_back(sl[i]-ngroup[i]-ngroup[i+1]);
			}
		} 
    }else if (weights == "laplace"){ // O(n) calculation
	
		sum1 = 0;
		sum2 = 0;
		vector<double> sl1(n),sl2(n);
	
		for (int i=n-1;i>0;i--){
			sum1 += ngroup[i]*exp(-gamma * xm[i]);
			sl1[i-1] = sum1;
		}
		sl1[n-1] =0;
		
		for (int i=1;i<n;i++){
			sum2 += ngroup[i-1]*exp(gamma * xm[i-1]);
			sl2[i] = sum2;
		}
		sl2[0] =0;	
		
		for (int i=0;i<n;i++){
			sl.push_back(sl1[i]*exp(gamma * xm[i])- sl2[i]*exp(-gamma * xm[i]));
		}		
	
	// O(n2) calculation for all the rest
	// gaussian and adaptive trick : summing starting by the smallest element for each line
	// the smallest are the one away from the "diagonal"
	// create two vectors of sum and then add them
	}else if (weights == "gaussian"){
		// calculate the positive part
		for(int i=0;i<(n-1);i++){
			sum1=0;
			for(int j =n-1;j>i;j--){
				sum1 += ngroup[j] * exp(-gamma* square(xm[i]-xm[j]));
			}
			sl.push_back(sum1);
		}
		sl.push_back(0);
		// calculate the negative part
		
		for(int i=1;i<n;i++){
			sum2=0;
			for (int j = 0;j<i ;j++){
				sum2 += ngroup[j] * exp(-gamma* square(xm[i]-xm[j]));
			}
			sl[i]-=sum2;
		}

	}else if (weights == "adaptive"){
		// calculate the positive part
		for(int i=0;i<(n-1);i++){
			sum1=0;
			for(int j =n-1;j>i;j--){
				sum1 += ngroup[j] /pow(fabs(xm[i]-xm[j]),gamma);
			}
			sl.push_back(sum1);
		}
		sl.push_back(0);
		// calculate the negative part
		
		for(int i=1;i<n;i++){
			sum2=0;
			for (int j = 0;j<i ;j++){
				sum2 += ngroup[j] /pow(fabs(xm[i]-xm[j]),gamma);
			}
			sl[i]-=sum2;
		}

	}else if (weights == "naivettest"){ 
		for (int i=0; i<n; i++){
			sum1 = 0;
			sign = -1.0;
			for (int j=0; j<n; j++){
				if (i==j){
					sign= 1.0;
				}else{
					sum1+= sign / sqrt(1/ngroup[i]+1/ngroup[j]) ;
				}
			} 
			sl.push_back(sum1/ ngroup[i]);
		} 	
	}else if (weights == "ttest"){
		for (int i=0; i<n; i++){
			sum1 = 0;
			sign = -1.0;
			for (int j=0; j<n; j++){
				if (i==j){
					sign= 1.0;
				}else{
					sum1+= sign * sqrt(((ngroup[i]-1)*xv[i] + (ngroup[j]-1)*xv[j])/(ngroup[i]+ngroup[j]-2)) /sqrt(1/ngroup[i]+1/ngroup[j]) ;
				}
			} 
			sl.push_back(sum1/ ngroup[i]);
		} 	
	}else if (weights == "welch"){
		for (int i=0; i<n; i++){
			sum1 = 0;
			sign = -1.0;
			for (int j=0; j<n; j++){
				if (i==j){
					sign= 1.0;
				}else{
					sum1+= sign * sqrt(xv[i]/ngroup[i]+xv[j]/ngroup[j])/(1/ngroup[i]+1/ngroup[j]) ;
				}
			} 
			sl.push_back(sum1/ngroup[i]);
		} 	
	}else if (weights == "personal"){
		for (int i=0; i<n; i++){
			sum1 = 0;
			sign = -1.0;
			for (int j=0; j<n; j++){
				if (i==j){
					sign= 1.0;
				}else{
					sum1 += sign * W(i,j);
				}
			} 
			sl.push_back(sum1/ngroup[i]);
		} 	
	//}else{
	//	Rprintf("Unknown weight parameter formulation. Aborting. \n") ;
	//	exit(0);
	}

	return(sl);

}

// General function
Group* maketree(double* x ,int K, double* sl, double* ngroup, double eps){

	// initialisation 
	Group *g=0,*gi,*gj;
	Groups::iterator g_it, g_it_next;
	Groups myGroups;

	// set epsilon tolerance parameter
	epsilon = eps; 
	
	//create the innitial list of groups 
	for (int i=0; i<K; i++){
		g = new Group;
	
		g->beta = x[i];
		g->lambda = 0.0;
		g->slope = sl[i]; 
		g->nbind = ngroup[i];
		// to avoid weight recalculation we calculate the slopes before pre-fusing
		// not to stock weight matrix, we simply calculate the innitial slopes once
		g->idown = i+1;
		g->len = 1;    

		g->childup = NullGroup;  
		g->childdown = NullGroup; 
		g->father = NullGroup;
	
		myGroups.push_back(g); // append my list of groups
	
	} 

	//print_groups(myGroups);

	// myGroups.sort ? not needed because groups are already sorted by ascending beta
	// indeed : x is sorted in R and the list is appended

	// Prefuse step :
	// if two groups have the same alpha at the begining of the algorithm, they must be fused.
	g_it=g_it_next=myGroups.begin();
	g_it_next++;
	while(g_it_next!=myGroups.end()){
	
		gi =*g_it;
		gj =*g_it_next;
		if(gi->beta == gj->beta){ // if two groups have same initial Beta
			//then merge the 2
		
			// Create the new group
			g = new Group;
			// and set all his attributes
			g->len = gi->len + gj->len;
			g->nbind = gi->nbind + gj->nbind;
			g->idown = gi->idown; // because gi had smaller beta //no need to modify idown it stays the same
			g->slope=(gi->nbind * gi->slope + gj->nbind * gj->slope)/(g->nbind); // new slope is just the ponderate mean of old slopes
			
			g->beta = gi->beta;
			g->lambda=0;
			
			g->childdown = gi;
			g->childup = gj;
			g->father = NullGroup;

			// set father
			gi->father = g; 
			gj->father = g;

			g_it_next++;
			
			g_it=myGroups.erase(g_it); // erase the group in g_it position
			*g_it=g; // replace value of the second group
		
		}else{ // moove iterators 
		
			g_it_next++;
			g_it++;
			
		}
	}

	// print_groups(myGroups);

	// Fusing 
	Group* G = fuse_groups(myGroups);

	return G;
}


// fuse the group and make the tree
Group* fuse_groups(Groups myGroups){
  
	Groups::iterator g_it,g_it_next,g_it_prev;
	Events events;
	// set to 0 here to avoid compiler warnings.
	Group *g=0,*gi,*gj,*g_next, *g_prev=0;
	Event e;
	Events::iterator e_it;
	double lambda;

	// generate the first list of events
	for(g_it=myGroups.begin(),g_it_next=myGroups.begin(),g_it_next++;
		g_it_next!=myGroups.end();
		g_it++,g_it_next++){
		
		g = *g_it;
		g_next=*g_it_next;
		// calculate initial events list from intersection of all adjacent
		lambda= getlambda(g, g_next);
		if(lambda>0){ 
			e=Event(lambda,g_it);
			e_it = events.insert(e);
			g->nextFusionUp = e_it;
		}
	}

	//print_events(events);

	while(events.size()>0){
		do{
			e_it=events.begin(); // e_it is automatically the one with lambda min
			lambda=e_it->first;
			g_it = g_it_next = e_it->second;
			gi= *g_it;
			g_it_next++;
			gj=*g_it_next; // gi and gj are the 2 fusing groups

			//print_events(events);
			//print_groups(myGroups);
		
			// Create the new group
			g = new Group;
			// and set all hist attributes
			g->len = gi->len + gj->len;
			g->nbind =  gi->nbind + gj->nbind;
			g->idown = gi->idown; // because gi had smaller beta
			// new slope is just the ponderate mean of old slopes
			g->slope=(gi->nbind * gi->slope + gj->nbind * gj->slope)/(g->nbind);
			g->lambda=lambda;
			g->beta = gi->beta + (g->lambda - gi->lambda)*gi->slope; 
			// ie. Beta = old Beta + Delta(lambda)* old slope
			g->childdown=gi;
			g->childup=gj;
			g->father = NullGroup;

			// set father
			gi->father = g; 
			gj->father = g;

			//set up pointers to neighboring groups
			g_it_prev=e_it->second;
			g_it_prev--;
			g_it_next++;

			//We have prev-gi-gj-next, we remove gi==g_it from the list (erase)
			//and set gj=g_it 
			g_it=myGroups.erase(g_it);
			*g_it=g;

			if(g_it!=myGroups.begin()){ //if the new group != group with min Beta
				g_prev= *g_it_prev;
				insert_new_event(&events,g,g_prev,g_it_prev,g_prev);
			}
			if(g_it_next != myGroups.end()){ //if the new group != group with max Beta
				insert_new_event(&events,g,*g_it_next,g_it,gj);
			} 
			events.erase(e_it);
		}while(events.begin()->first == lambda);
	}
	
	//print_groups(myGroups);
	return g;

}

//calculate lambda where g1 and g2 fuse
double getlambda(Group *g1,Group *g2){
	return((g1->beta - g2->beta - g1->slope * g1->lambda + g2->slope * g2->lambda)/(g2->slope-g1->slope));
}

// used for adding/deleting previous and next joining events
void insert_new_event(Events *eventlist,Group *newg,
						Group *neighbor_group,Groups::iterator g_it_up,
						Group *g_down){
						
	if(g_down-> nextFusionUp != NullEvent){
		(*eventlist).erase(g_down->nextFusionUp); // erase the event that will not happen
	}
	Group *g_up = *g_it_up;
	Events::iterator new_event;
	Event e=Event(getlambda(neighbor_group,newg),g_it_up);

	if(newg->lambda -e.first <= NEW_EVENT_THRESH){
		// 2 things : new lambda >= old lambda and add a treshold to prevent round up errors
		new_event=(*eventlist).insert(e);
	}else{
		new_event=NullEvent;
	}

	g_up->nextFusionUp = new_event; // contains the lambda for which the group is fusing the one on top of him
}

void print_groups(Groups myGroups){
	Groups::const_iterator g_it;
	Rprintf("idown iup beta lambda slope \n");
	for(g_it = myGroups.begin(); g_it != myGroups.end();++g_it){
		Rprintf("%i %i %f %f %f \n",(*g_it)->idown,(int)((*g_it)->idown+(*g_it)->len-1),(*g_it)->beta,(*g_it)->lambda,(*g_it)->slope);
	}
	Rprintf("\n");
}

void print_events(Events events){
	Events::const_iterator it;
	for(it = events.begin(); it != events.end();++it){
		Rprintf("%f ",it->first);
		if(it->second == (Groups::iterator)0) Rprintf("NULL \n");
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


void add_results_predict(Group *g,double *beta, double *lambda, double *slope, int *idown,int *iup,int *row,
			 vector<double> *pred, double *lambdalist, vector<double> *lambdalistpred, vector<double> *slopepred, vector<int> *idownpred, vector<int> *iuppred,
			 int nlambda){
  
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


void error_cv(Group *g, double* lambdalist,int lambdalen, double* xtest,double* ngroup, double* error){

	int index=0;

	if (g->father!= NullGroup){ // we have the lambda
    
		updateError(error, xtest, ngroup, lambdalist,lambdalen, g->father, g);
    
	}else{ // final group : calc further err 
       
		// it is over so we can update the end of the lambdalist (case where lambdalist[i] >= final_g->lambda)
		while(lambdalist[index]<g->lambda){
			index+=1;
			if (index==lambdalen){break;}
		} // lambdalist[index]>= g->lambda
	
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
void updateError(double* err,double* xtest,double* ngroup,double* lambdalist,int lambdalen, Group* gnew,Group*gold){

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

