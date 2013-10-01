#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include "withsplit.h"

#define square(x) ((x)*(x))

// initialize the group class to have sufficient space such that the largest nodeNumber can still be stored
FAGeneral::FAGeneral(double* x,int K, double* ngroup, NumericVector xv, std::string weights, double gamma, NumericMatrix W, int mxSize, bool verb, double eps):K(K),Weights(K,vector<double>(K)),graphsize(K){
	verbose = verb;
    mxSplitSize = mxSize;
	epsilon = eps;

	// calculate the weights we need
	calculateSlope(x,ngroup,xv,weights,gamma,W);
	calculateWeights(x,ngroup,xv,weights,gamma,W);
	// print_weights();
	
	if (verbose){
		Rprintf("Graph Construction \n");
	}
	
	myGraph.Construct(K,Weights,epsilon);
	// myGraph.printGraph(cout);
	
	// set groups, prefuse and set events 
    initializeGroups(&x[0],&ngroup[0]);

	initializeEvents();
	
	// print_groups(); 
	// print_events();
	
    run();
    
}

void FAGeneral::initializeGroups(double* x, double* ngroup){
	
	// initialisation 
	Group *g=0,*gi,*gj;
	Groups::iterator g_it, g_it_next;
	bool prefuse;
	
	if(verbose){
        Rprintf("Started initializing the %i Groups\n", K);
    }

	//create the innitial list of groups 
	for (int i=0; i<K; i++){
		
		g = new Group;
	
		g->beta = x[i];
		g->lambda = 0.0;
		g->slope = sl[i]; 
		g->nbind = ngroup[i];
		
		// to avoid weight recalculation we calculate the slopes before pre-fusing
		// not to stock weight matrix, we simply calculate the innitial slopes once

		g->grps.insert(i);
		g->content.insert(make_pair(i, g));
		// g->sizes.push_back(ngroup[i]);	
			
		g->len = 1;    

		g->childup = NullGroup;  
		g->childdown = NullGroup; 
		g->fatherup = NullGroup;
		g->fatherdown = NullGroup;
		
		g->m = myGraph.subGraph(g->grps,g->content,g->slope); // update mx flowgraph
		
		myGroups.push_back(g); // append my list of groups

	} 

	//print_groups(myGroups);

	// myGroups.sort ? not needed because groups are already sorted by ascending beta
	// indeed : x is sorted in R and the list is appended

	// Prefuse step :
	// if two groups have the same alpha at the begining of the algorithm, they must be fused.

	if(verbose){
        Rprintf("Started eventual prefuse step. \n");
    }

	prefuse = false;
	g_it=g_it_next=myGroups.begin();
	g_it_next++;
	while(g_it_next!=myGroups.end()){
	
		gi =*g_it;
		gj =*g_it_next;
		if(gi->beta == gj->beta){ // if two groups have same initial Beta
			//then merge the 2
			prefuse= true;
			// Create the new group
			g = new Group;
			// and set all his attributes
			g->len = gi->len + gj->len;
			g->nbind = gi->nbind + gj->nbind;
			
			g->slope=(gi->nbind * gi->slope + gj->nbind * gj->slope)/(g->nbind); // new slope is just the ponderate mean of old slopes
			
			(g->grps).insert(gi->grps.begin(), gi->grps.end());
			(g->grps).insert(gj->grps.begin(), gj->grps.end());
			(g->content).insert(gi->content.begin(),gi->content.end());
			(g->content).insert(gj->content.begin(),gj->content.end());
			
			g->beta = gi->beta;
			g->lambda=0;
			
			g->childdown = gi;
			g->childup = gj;
			g->fatherup = NullGroup;
			g->fatherdown = NullGroup;
			
			g->m = myGraph.subGraph(g->grps,g->content,g->slope); // update mx flowgraph
			
			// delete pointers to mxflow
			delete gi->m;
			delete gj->m;
			
			// set father
			gi->fatherup = g; 
			gj->fatherup = g;
			gi->fatherdown = NullGroup;
			gj->fatherdown = NullGroup;
			
			// update the graphsize (nb here there will be doublon in final dataframe)
			graphsize = graphsize + gi ->grps.size() + gj->grps.size();
			
			g_it_next++;
			
			g_it=myGroups.erase(g_it); // erase the group in g_it position
			*g_it=g; // replace value of the second group
		
		}else{ // moove iterators 
		
			g_it_next++;
			g_it++;
			
		}
	}

	if(verbose && prefuse){
        Rprintf("A Prefuse step occured\n");
    }
	
}

void FAGeneral::initializeEvents(){

 	Groups::iterator g_it,g_it_next;
	Group *g=0,*g_next;
	Event e;
	Action a;
	Events::iterator e_it;
	double lambda;
	
	if(verbose){
        Rprintf("Started initializing Events.\n");
    }

	// generate the first list of events
	for(g_it=myGroups.begin(),g_it_next=myGroups.begin(),g_it_next++;
		g_it_next!=myGroups.end();
		g_it++,g_it_next++){
		
		g = *g_it;
		g_next=*g_it_next;
		// calculate initial events list from intersection of all adjacent
		lambda= getlambda(g, g_next);
		if(lambda>0){ 
			a = Action(g_it,'F');
			e = Event(lambda,a);
			e_it = events.insert(e);
			g->nextFusionUp = e_it;
		}
		// check split only if group size > 1, size groupe < maxsizegroup  and not last group.
		if (((int)g->grps.size()>1)&&((int)g->grps.size() <= mxSplitSize)&&((int)myGroups.size() !=1)){
			MxFlowCheck(0, g_it,true);	// mxflow for split checking	 
		}
	}
 
}

void FAGeneral::run(){
	
	char eventType;
	Events::iterator e_it;
	double lambda;
	
	if(verbose){
        Rprintf("Started main algorithm.\n");
    }
	
	//print_events();
	//print_groups();
	
	while(events.size()>0){
		do{
			
			// print_events();
			// print_groups();
			
			e_it=events.begin(); // e_it is automatically the one with lambda min.
			
			eventType = e_it->second.second; // get the event type "F" or "M"
			lambda = e_it->first;
			
			if (eventType == 'F'){
			
				if(verbose){
					Rprintf("Fusing groups at %f, min %.14f. \n", lambda, (*(e_it->second.first))->beta);
				}
			
				FuseGroups(e_it);

				if(verbose){
					Rprintf("Fusing done. \n", lambda);
				}
				
			}else{ // eventType ="M" for maxflow check
			
				if(verbose){
					Rprintf("Max Flow Check Update at %f. \n", lambda);
				}
				
				MxFlowCheck(e_it->first, e_it->second.first,false); // update split events, do split if needed, etc...
				// nb : we know that size(g)<=mxSplitSize cause if not the event would not exist
			}
			// erase the treated event 
			events.erase(e_it);

			//print_events();
			//print_groups();
			
		}while(events.begin()->first == lambda);
	}
	
	// print_events();
	// print_groups();

	// delete the last maxflow graph 
	delete (myGroups.front())->m ;
	
	fflush(stdout);
	if (verbose){
		Rprintf("Run is complete \n");
		fflush(stdout);
	}
}

void FAGeneral::FuseGroups(Events::iterator e_it){

	Groups::iterator g_it,g_it_next,g_it_prev;
	Group *g=0,*gi,*gj, *g_prev=0;

	double lambda;
	
	lambda = e_it->first;
	g_it = g_it_next = e_it->second.first;

	gi= *g_it;
	g_it_next++;
	gj=*g_it_next; // gi and gj are the 2 fusing groups

	//print_content(gi);
	//print_content(gj);
	
	// Create the new group
	g = new Group;
		
	// and set all hist attributes
	g->len = gi->len + gj->len;
	g->nbind =  gi->nbind + gj->nbind;
		
	g->grps = gi->grps; // insert gi grps	
	(g->grps).insert(gj->grps.begin(), gj->grps.end()); // insert gj grps

	g->content = gi->content;
	(g->content).insert(gj->content.begin(),gj->content.end());
					
	// new slope is just the ponderate mean of old slopes
	g->slope=(gi->nbind * gi->slope + gj->nbind * gj->slope)/(g->nbind);
	g->lambda=lambda;
	g->beta = gi->beta + (g->lambda - gi->lambda)*gi->slope; 
	// ie. Beta = old Beta + Delta(lambda)* old slope
	g->childdown=gi;
	g->childup=gj;
	g->fatherup = NullGroup;
	g->fatherdown= NullGroup;	

	// delete mxflow pointer
	delete gi->m;
	delete gj->m;

	// set to null the next split event (temp)
	g->nextMxFlow = NullEvent;
	g->m = myGraph.subGraph(g->grps, g->content, g->slope); // update mx flow graph

	// set father
	gi->fatherup = g; 
	gj->fatherup = g;
	gi->fatherdown = NullGroup;
	gj->fatherdown = NullGroup;
	
	// update the graphsize
	graphsize = graphsize + gi ->grps.size() + gj->grps.size();
				
	//set up pointers to neighboring groups
	g_it_prev=e_it->second.first;
	g_it_prev--;
	g_it_next++;

	// erase some events that will never happen
	if(gi-> nextMxFlow!= NullEvent){
		events.erase(gi->nextMxFlow); // erase the first mx flow event that will not happen
	} 
	if(gj-> nextMxFlow!= NullEvent){
		events.erase(gj->nextMxFlow); // erase the second mx flow event that will not happen
	} 

	//We have prev-gi-gj-next, we remove gi==g_it from the list (erase)
	//and set gj=g
	g_it=myGroups.erase(g_it);
	*g_it=g;
	
	if(g_it!=myGroups.begin()){ //if the new group != group with min Beta => we delete old merge event and add a new one
		g_prev= *g_it_prev;
		insert_new_merge_event(g,g_prev,g_it_prev,g_prev);
	}	
	
	if(g_it_next != myGroups.end()){ //if the new group != group with max Beta
		insert_new_merge_event(g,*g_it_next,g_it,gj);
	} 

	// if two groups fuse we need to update the split events. 
	if (((int)g->grps.size()<= mxSplitSize) && myGroups.size()!=1){
		MxFlowCheck(e_it->first,g_it,true);		
	}

}

// used for adding/deleting previous and next joining events
void FAGeneral::insert_new_merge_event(Group *newg,
						Group *neighbor_group,Groups::iterator g_it_down,
						Group *g_down){

	// print_events();
	if(g_down->nextFusionUp != NullEvent){	
		events.erase(g_down->nextFusionUp); // erase the event that will not happen
	}

	Group *g_up = *g_it_down;
	Events::iterator new_event;
	Action a = Action(g_it_down,'F');
	Event e=Event(getlambda(neighbor_group,newg),a);

	if(newg->lambda - e.first <= epsilon){
		// 2 things : new lambda >= old lambda and add a treshold to prevent round up errors
		new_event = events.insert(e);
	}else{
		new_event=NullEvent;
	}

	g_up->nextFusionUp = new_event; // contains the lambda for which the group is fusing the one on top of him
}

//calculate lambda where g1 and g2 fuse
double FAGeneral::getlambda(Group *g1,Group *g2){
	return((g1->beta - g2->beta - g1->slope * g1->lambda + g2->slope * g2->lambda)/(g2->slope-g1->slope));
}


void FAGeneral::MxFlowCheck(double lambda, Groups::iterator g,bool newgroup){

	double hitTime;
	
	if(verbose){
		Rprintf("MxFlowCheck \n");
	}
	// calculating the next hitTime : neversplit, splitnow or a double value > current lambda
		
	if (newgroup){ // the group is new
	
		hitTime = (*(*g)->m).mxFlowNew(lambda);
	
	}else{ // max flow result no longer valid
	
		hitTime = (*((*g)->m)).mxFlowUpdate(lambda);
	
	}
	
	// three cases : 
	// this split time is infinite => nothing happens
	// this split time is the current one => split
	// this split time is later : include this in the events
	
	if(hitTime==neverSplit){ // Nothing happens
	
		if(verbose){
			Rprintf("Group will never split.\n");
		}
	
		return;
		
	}else if(hitTime==splitNow){ // make the split and schedule the new events
		
		if(verbose){
			Rprintf("Splitting at %f.\n", lambda);
		}
		
		split(lambda, g);
		
		return;
		  
	}else{ // insert the new event
       
		if(verbose){
			Rprintf("New event scheduled at %f.\n", hitTime);
		}
	   
		Action a = Action(g,'T');
		Event e = Event(hitTime,a);
		Events::iterator e_it = events.insert(e);
		(*g)->nextMxFlow = e_it;
		
		return;
		
	}
	
}


void FAGeneral::split(double lambda, Groups::iterator g){
    
	Group *gdown, *gup=0;
	Groups::iterator g_it_up, g_it_down, g_it_next,g_it_prev;
	set<int>::iterator SI; // iterate through the edges in the subgraph
	
    // save the nodes of the first splitted graph and the second splitted graph
	set<int> splitNodes1, splitNodes2;
	pair<set<int>,set<int> > Sets;
	
	// Get the two sets for the split
	Sets = ((*g)->m)->GetSets();

	// if one of the set is empty : throw error (because it would generate an infinite loop)
	if((Sets.first).empty() || (Sets.second).empty()){

		Rcpp::stop("A splitting group is empty. There are approximation errors. \n Try to increase the epsilon tolerance parameter.");
	}
	
	/////////////////////////////////
	// Reconstruct the two new groups
	/////////////////////////////////
	gdown = new Group;
	gup = new Group;
	
	gup->lambda=lambda;
	gdown->lambda=lambda;
	gup->beta = (*g)->beta + (*g)->slope * (lambda-(*g)->lambda);
	gdown->beta = (*g)->beta + (*g)->slope * (lambda-(*g)->lambda);
	gup->len = 0;
	gdown->len = 0;
	gup->nbind =0;
	gdown->nbind =0;
	gup->slope=0;
	gdown->slope=0;
	
	// slopes, content
	// group connected to source => push >0 => slope >0 => gup
	// the slope can be expresses as the difference between Beta and baricenter divided by lambda
	

	
	for (SI = (Sets.first).begin(); SI!=(Sets.first).end(); ++SI){ 
		(gup->grps).insert(*SI);
		(gup->content).insert(make_pair(*SI, (*g)->content[*SI]));
		gup->nbind += ((*g)->content[*SI])->nbind;
		gup->len +=1 ;
		gup->slope -= ((*g)->content[*SI])->beta * ((*g)->content[*SI])->nbind;
	}	
	gup->slope = (gup->slope/(gup->nbind) + gup->beta  ) / lambda;
	
	// same thing for the other group
	for (SI = (Sets.second).begin(); SI!=(Sets.second).end(); ++SI){ 
		(gdown->grps).insert(*SI);
		(gdown->content).insert(make_pair(*SI, (*g)->content[*SI]));
		gdown->nbind += ((*g)->content[*SI])->nbind;
		gdown->len +=1 ;
		gdown->slope -= (((*g)->content[*SI])->beta) * (((*g)->content[*SI])->nbind);
	}	
	gdown->slope = (gdown->slope/(gdown->nbind) + gdown->beta  ) / lambda;
	
	// update nextFusionUp and nextMxFlow
	gup->nextFusionUp = (*g)->nextFusionUp; // trick to delete correctly the event via insert_event
	gdown->nextFusionUp = NullEvent;
	gup->nextMxFlow = NullEvent;
	gdown->nextMxFlow = NullEvent;

	// update g father 		
	(*g)->fatherup = gup;
	(*g)->fatherdown= gdown;
	
	// update child
	gup->childdown = NullGroup;
	gup->childup = *g; 
	gdown->childdown = NullGroup;
	gdown->childup = *g;
	// init father
	gup->fatherup = NullGroup;
	gup->fatherdown= NullGroup;
	gdown->fatherup = NullGroup;
	gdown->fatherdown= NullGroup;	
	
	// update max flow graphs
	gup->m = myGraph.subGraph(gup->grps,gup->content,gup->slope);
	gdown->m = myGraph.subGraph(gdown->grps,gdown->content,gdown->slope);
	
	// update the graph size 
	graphsize = graphsize + gup ->len + gdown->len;
	
	// insert the new groups in mygroups
	g_it_up = g;
	g_it_down = g;
	
	// insert new groups
	*g_it_up = gup;
	myGroups.insert(g_it_up,gdown);
	
	// set new groups iterator (and reset g_down since we do not know where it is)
	g_it_next = g_it_up;
	g_it_next ++;
	g_it_down = g_it_up;
	g_it_down --;
	
	// generate fuse events
	// note that we do not insert a merge event for the two new groups because then cannot fuse back
	if(g_it_down != myGroups.begin()){ // g_it_prev exists
		g_it_prev = g_it_down;
		g_it_prev--;
		insert_new_merge_event(*g_it_down,*g_it_prev,g_it_prev,*g_it_prev);
	}	
	
	if(g_it_next != myGroups.end()){ //if the new group != group with max Beta
		insert_new_merge_event(*g_it_up,*g_it_next,g_it_up,*g_it_up);
	} 
	
	// print_events();
	// print_groups();
	
	// maxflow check on new groups 
	// don't need to test the max size cause came from a split ie size<previous size 
	if ((int)(*g_it_down)->grps.size()>1){
		MxFlowCheck(lambda, g_it_down,true);
	}
	if ((int)(*g_it_up)->grps.size()>1){	
		MxFlowCheck(lambda, g_it_up,true);
	}
	
}


void FAGeneral::calculateSlope(double* xm,double* ngroup,NumericVector xv, 
	std::string weights, double gamma, NumericMatrix W){
	
	double sign;
	double sum1=0,sum2=0;
	double N = 0;
	for (int i=0; i<K; i++){ 
		N += ngroup[i]; // number of individuals in all groups
	}
	
	if (weights=="default"){ // easiest

		if (N==K){ //slopes are known (clusterpath case)
			for (int i =0; i<K; i++){
				sl.push_back(K-1-2*i);
			}
		}else{ // nk.nl weights (O(n) calculation)
			for (int i=1;i<K;i++){
				sum1 += ngroup[i];
			}
			sl.push_back(sum1);
			for (int i =0; i<(K-1);i++){
				sl.push_back(sl[i]-ngroup[i]-ngroup[i+1]);
			}
		} 
    }else if (weights == "laplace"){ // O(n) calculation
	
		sum1 = 0;
		sum2 = 0;
		vector<double> sl1(K),sl2(K);
	
		for (int i=K-1;i>0;i--){
			sum1 += ngroup[i]*exp(-gamma * xm[i]);
			sl1[i-1] = sum1;
		}
		sl1[K-1] =0;
		
		for (int i=1;i<K;i++){
			sum2 += ngroup[i-1]*exp(gamma * xm[i-1]);
			sl2[i] = sum2;
		}
		sl2[0] =0;	
		
		for (int i=0;i<K;i++){
			sl.push_back(sl1[i]*exp(gamma * xm[i])- sl2[i]*exp(-gamma * xm[i]));
		}		
	
	// O(n2) calculation for all the rest
	// gaussian and adaptive trick : summing starting by the smallest element for each line
	// the smallest are the one away from the "diagonal"
	// create two vectors of sum and then add them
	}else if (weights == "gaussian"){
		// calculate the positive part
		for(int i=0;i<(K-1);i++){
			sum1=0;
			for(int j =K-1;j>i;j--){
				sum1 += ngroup[j] * exp(-gamma* square(xm[i]-xm[j]));
			}
			sl.push_back(sum1);
		}
		sl.push_back(0);
		// calculate the negative part
		
		for(int i=1;i<K;i++){
			sum2=0;
			for (int j = 0;j<i ;j++){
				sum2 += ngroup[j] * exp(-gamma* square(xm[i]-xm[j]));
			}
			sl[i]-=sum2;
		}

	}else if (weights == "adaptive"){
		// calculate the positive part
		for(int i=0;i<(K-1);i++){
			sum1=0;
			for(int j =K-1;j>i;j--){
				sum1 += ngroup[j] /pow(fabs(xm[i]-xm[j]),gamma);
			}
			sl.push_back(sum1);
		}
		sl.push_back(0);
		// calculate the negative part
		
		for(int i=1;i<K;i++){
			sum2=0;
			for (int j = 0;j<i ;j++){
				sum2 += ngroup[j] /pow(fabs(xm[i]-xm[j]),gamma);
			}
			sl[i]-=sum2;
		}

	}else if (weights == "naivettest"){ 
		for (int i=0; i<K; i++){
			sum1 = 0;
			sign = -1.0;
			for (int j=0; j<K; j++){
				if (i==j){
					sign= 1.0;
				}else{
					sum1+= sign / sqrt(1/ngroup[i]+1/ngroup[j]) ;
				}
			} 
			sl.push_back(sum1/ ngroup[i]);
		} 	
	}else if (weights == "ttest"){
		for (int i=0; i<K; i++){
			sum1 = 0;
			sign = -1.0;
			for (int j=0; j<K; j++){
				if (i==j){
					sign= 1.0;
				}else{
					sum1+= sign * sqrt(((ngroup[i]-1)*xv[i] + (ngroup[j]-1)*xv[j])/(ngroup[i]+ngroup[j]-2)) /sqrt(1/ngroup[i]+1/ngroup[j]) ;
				}
			} 
			sl.push_back(sum1/ ngroup[i]);
		} 	
	}else if (weights == "welch"){
		for (int i=0; i<K; i++){
			sum1 = 0;
			sign = -1.0;
			for (int j=0; j<K; j++){
				if (i==j){
					sign= 1.0;
				}else{
					sum1+= sign * sqrt(xv[i]/ngroup[i]+xv[j]/ngroup[j])/(1/ngroup[i]+1/ngroup[j]) ;
				}
			} 
			sl.push_back(sum1/ngroup[i]);
		} 	
	}else if (weights == "personal"){
		for (int i=0; i<K; i++){
			sum1 = 0;
			sign = -1.0;
			for (int j=0; j<K; j++){
				if (i==j){
					sign= 1.0;
				}else{
					sum1 += sign * W(i,j);
				}
			} 
			sl.push_back(sum1/ngroup[i]);
		} 	

	}

}

void FAGeneral::calculateWeights(double* xm,double* ngroup,NumericVector xv, 
	std::string weights, double gamma, NumericMatrix W){
	
	double N = 0;
	for (int i=0; i<K; i++){ 
		N += ngroup[i]; // number of individuals in all groups
	}

	if (weights=="default"){ // easiest

		if (N==K){
			for (int i =0; i<K; i++){
				for (int j = 0; j<K;j++){
					Weights[i][j]=1;
				}
			}
		}else{ // nk.nl weights (O(n) calculation)
			for (int i =0; i<K; i++){
				for (int j = 0; j<K;j++){
					Weights[i][j]=ngroup[i]*ngroup[j];
				}
			}
		} 
    }else if (weights == "laplace"){ // O(n) calculation
		for (int i =0; i<K; i++){
			for (int j = 0; j<K;j++){
				Weights[i][j]=ngroup[i]*ngroup[j]*exp(- gamma  * fabs(xm[i]-xm[j]));
			}
		}

	}else if (weights == "gaussian"){
		for (int i =0; i<K; i++){
			for (int j = 0; j<K;j++){
				Weights[i][j]=ngroup[i]*ngroup[j]*exp(- gamma  * square(xm[i]-xm[j]));
			}
		}

	}else if (weights == "adaptive"){

		for (int i =0; i<K; i++){
			for (int j = 0; j<K;j++){
				Weights[i][j]=ngroup[i]*ngroup[j]/pow(fabs(xm[i]-xm[j]),gamma);
			}
		}
		
	}else if (weights == "naivettest"){ 
	
		for (int i =0; i<K; i++){
			for (int j = 0; j<K;j++){
				Weights[i][j]=1/sqrt(1/ngroup[i]+1/ngroup[j]);
			}
		}
		
	}else if (weights == "ttest"){
	
		for (int i =0; i<K; i++){
			for (int j = 0; j<K;j++){
				Weights[i][j]=sqrt(((ngroup[i]-1)*xv[i] + (ngroup[j]-1)*xv[j])/(ngroup[i]+ngroup[j]-2)) /sqrt(1/ngroup[i]+1/ngroup[j]);
			}
		}
		
	}else if (weights == "welch"){
	
		for (int i =0; i<K; i++){
			for (int j = 0; j<K;j++){
				Weights[i][j]=sqrt(xv[i]/ngroup[i]+xv[j]/ngroup[j])/(1/ngroup[i]+1/ngroup[j]);
			}
		}
		
	}else if (weights == "personal"){
		for (int i =0; i<K; i++){
			for (int j = 0; j<K;j++){
				Weights[i][j] = W(i,j);
			}
		}	

	}

}

Group* FAGeneral::GetFinalGroup(){
return(myGroups.front());
}

void FAGeneral::print_groups(){
	Groups::const_iterator g_it;
	set<int>::const_iterator s_it;
	Rprintf("list beta lambda slope \n");
	for(g_it = myGroups.begin(); g_it != myGroups.end();++g_it){
		for(s_it= ((*g_it)->grps).begin(); s_it!=((*g_it)->grps).end(); ++s_it){
			Rprintf("%i, ",*s_it);
		}
		Rprintf("%f %f %f \n",(*g_it)->beta,(*g_it)->lambda,(*g_it)->slope);
	}
	Rprintf("\n");
}

void FAGeneral::print_events(){
	Events::const_iterator it;
	for(it = events.begin(); it != events.end();++it){
		Rprintf("lambda : %f , ",it->first);
		if(it->second.first == (Groups::iterator)0) Rprintf("NULL \n");
		else Rprintf("beta : %.14f , type : %c \n",(*(it->second.first))->beta,it->second.second );
	}
	Rprintf("size: %i \n",events.size());
}


void FAGeneral::print_content(Group* g){
	Group* curg;
	map<int,Group * >::iterator MI;
	Rprintf("index beta nbind slope lambda \n");
	for (MI = g->content.begin(); MI != g->content.end(); MI++){
		curg = MI->second;
		Rprintf("%i %.14f %i %f %f \n", MI->first, curg->beta, curg->nbind, curg->slope, curg->lambda);
	}

}

void FAGeneral::print_weights(){
	for (int i=0; i<K;i++){
		for (int j = 0;j<K; j ++){
			Rprintf("%f  ", Weights[i][j]);
		}
		Rprintf("\n");
	}
}


void FAGeneral::delete_tree(Group *g){
	Group *gfdown = g->fatherdown;
	if (gfdown!= NullGroup){
		gfdown->childup = 0;
	}
	Group *g1=g->childup,*g2=g->childdown;
	if(g1 != NullGroup)delete_tree(g1);
	if(g2 != NullGroup)delete_tree(g2);
	delete g;
}

void add_results(Group *g,double *beta,
				double *lambda, double *slope, int *index, int *row){

	set<int>::iterator s_it;
				
	for(s_it = (g->grps).begin(); s_it!=(g->grps).end(); ++s_it){
		beta[*row]=g->beta;
		lambda[*row]=g->lambda;
		slope[*row]=g->slope;
		index[*row]= *s_it+1; // because groups are noted from 0 in method 
		*row+=1;
	}

	// if a split has happened we add normally for childup but we don't want to add it again 
	
	if(g->childup != NullGroup){
		if(g->childup->fatherdown != g){ // so we do not pass 2 times at the same point in case of split
			add_results(g->childup,beta,lambda,slope,index,row);
		}
	}
	if (g->childdown != NullGroup){
		add_results(g->childdown,beta,lambda,slope,index,row);
	} 
	
}			


void add_results_predict(Group *g,double *beta, double *lambda, double *slope, int *index,int *row,
			 vector<double> *pred, double *lambdalist, vector<double> *lambdalistpred, vector<double> *slopepred, vector<int> *indexpred,
			 int nlambda){
  
	int k =0;
 	set<int>::iterator s_it;
	
	// add result part
	for(s_it = (g->grps).begin(); s_it!=(g->grps).end(); ++s_it){
		beta[*row]=g->beta;
		lambda[*row]=g->lambda;
		slope[*row]=g->slope;
		index[*row]= *s_it+1; // because groups are noted from 0 in method 
		*row+=1;
	}
 
 
	if (g->fatherup!= NullGroup){ // we have the lambda !
    
		double l = g->fatherup->lambda;
		while (lambdalist[k]>l){ // moove to the right position in the list
			k =k+1;      
			if (k==(nlambda))
				break; 
		}

		if (k<(nlambda)){
			while (lambdalist[k]>=g->lambda){ // ok position
				for(s_it = (g->grps).begin(); s_it!=(g->grps).end(); ++s_it){
					pred->push_back(g->beta + (lambdalist[k]-g->lambda)*g->slope);
					lambdalistpred->push_back(lambdalist[k]);
					slopepred->push_back(g->slope);
					indexpred->push_back(*s_it+1); // because groups are noted from 0 in method
				}
				k =k+1;      
				if (k==(nlambda))
					break;
			}
		}
		
		// commented part following is not necessary since it give the same result !
		
		// same thing for father down if exists, this edge will not be gone through again
		
		// if (g->fatherdown!= NullGroup){ // we have the lambda !
		
			// double l = g->fatherdown->lambda;
			// while (lambdalist[k]>l){ // moove to the right position in the list
				// k =k+1;      
				// if (k==(nlambda))
					// break; 
			// }

			// if (k<(nlambda)){
				// while (lambdalist[k]>=g->lambda){ // ok position
					// for(s_it = (g->grps).begin(); s_it!=(g->grps).end(); ++s_it){
						// pred->push_back(g->beta + (lambdalist[k]-g->lambda)*g->slope);
						// lambdalistpred->push_back(lambdalist[k]);
						// slopepred->push_back(g->slope);
						// indexpred->push_back(*s_it+1); // because groups are noted from 0 in method );
					// }
					// k =k+1;      
					// if (k==(nlambda))
						// break;
				// }
			// }	
		// }
		
		
	}else{ // final group : predict further for all 
		while (lambdalist[k]>=g->lambda){ // ok position
			for(s_it = (g->grps).begin(); s_it!=(g->grps).end(); ++s_it){
				pred->push_back(g->beta + (lambdalist[k]-g->lambda)*g->slope);
				lambdalistpred->push_back(lambdalist[k]);
				slopepred->push_back( g->slope);
				indexpred->push_back(*s_it+1); // because groups are noted from 0 in method 
			}
			k =k+1;
		}
	}
	
	if(g->childup != NullGroup){
		if(g->childup->fatherdown != g){
			add_results_predict(g->childup,beta,lambda,slope,index,row,pred,lambdalist,lambdalistpred,slopepred,indexpred, nlambda);
		}
	}
	if(g->childdown != NullGroup){ // keep writing
		add_results_predict(g->childdown,beta,lambda,slope, index, row,pred,lambdalist,lambdalistpred,slopepred,indexpred, nlambda);
	}
}


void error_cv_ws(Group *g, double* lambdalist,int lambdalen, double* xtest,double* ngroup, double* error){

	int index=0;
	set<int>::iterator s_it;

	if (g->fatherup!= NullGroup){ // we have the lambda
    
		updateError_ws(error, xtest, ngroup, lambdalist,lambdalen, g->fatherup, g);
		
		// comment part below is false cause it counts twice the part.
		
		// if(g->fatherdown!= NullGroup){ // nb, we will never go through this edge again (tested later)
			// updateError_ws(error, xtest, ngroup, lambdalist,lambdalen, g->fatherdown, g);
		// }
    
	}else{ // final group : calc further err 
       
		// it is over so we can update the end of the lambdalist (case where lambdalist[i] >= final_g->lambda)
		while(lambdalist[index]<g->lambda){
			index+=1;
			if (index==lambdalen){break;}
		} // lambdalist[index]>= g->lambda
	
		if (index!=lambdalen){
			for(s_it = (g->grps).begin(); s_it!=(g->grps).end(); ++s_it){ // for all starting groups :
				error[index] = error[index] + ngroup[*s_it]* square(g->beta - xtest[*s_it]);
			}// run once, g->slope=0. 
			index +=1;
		}
		while(index<lambdalen){ // now the slope is 0, err[index-1] is touch only once => the rest is constant
			error[index] = error[index-1];
			index+=1;
		}
	}

	if(g->childup != NullGroup){
		if(g->childup->fatherdown != g){
			error_cv_ws(g->childup, lambdalist, lambdalen, xtest, ngroup,error);
		}
	}
	if(g->childdown != NullGroup){ // keep writing
		error_cv_ws(g->childdown, lambdalist, lambdalen, xtest, ngroup,error);
	}
}


// allow to add error terms to the err vector for cross validation
void updateError_ws(double* err,double* xtest,double* ngroup,double* lambdalist,int lambdalen, Group* gnew,Group*gold){

	int i =0;
	set<int>::iterator s_it;
	double oldl = gold->lambda;
	double newl = gnew->lambda;
	if (oldl!=newl){ // can be the same for multiple fuse or prefuse : in this case don't calculate
		while(lambdalist[i]<oldl){ // moove i to first interesting lambda value
			i+=1;
			if (i==lambdalen){break;}
		} // lambdalist[i]>=oldl

		if (i<lambdalen){
			while(lambdalist[i]<newl){ // we predict on this slope
				for(s_it = (gold->grps).begin(); s_it!=(gold->grps).end(); ++s_it){
					err[i] += ngroup[*s_it] * square(gold->beta + (lambdalist[i]-oldl)*gold->slope - xtest[*s_it]);
				}
				i+=1;
				if (i==lambdalen){break;}
			}
		}
	}
}


