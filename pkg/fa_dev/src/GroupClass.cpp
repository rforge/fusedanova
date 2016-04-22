#include "GroupClass.h"

Group::Group() : beta(0.0), lambda(0.0), slope(0), nbind(1), idown(1), len(1),
		 childup(NullGroup), childdown(NullGroup), father(NullGroup)  {} ;

Group::Group(double beta, double lambda, double slope, int nbind, int idown)  : 
  beta(beta), lambda(lambda), slope(slope), nbind(nbind), idown(idown), 
  len(1), childup(NullGroup), childdown(NullGroup), father(NullGroup) {} ;

Group::~Group() {}

