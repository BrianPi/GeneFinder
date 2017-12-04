#ifndef MUTATION_INFORMATION_H
#define MUTATION_INFORMATION_H

#include<map>
#include<numeric>

class mutationInformation {
	//Holds the relative frequency of base substitution and insertion mutation, with rudimentary support for multi-base insertions
	public:
	double getTsTv() {return (SBPTable[0]+SBPTable[1])/(SBPTable[2]+SBPTable[3]+SBPTable[4]+SBPTable[5]);}
	double getATGC() {return (SBPTable[0]+SBPTable[4])/(SBPTable[1]+SBPTable[3]);}
	vector<double> getSBP();
	double getSBPMeasure() {return accumulate(SBPTable.begin(),SBPTable.end(),0);} //Gets the total amount of mutation double in SBP table
	double getIDSBPRatio() {return (InsertionRate+DeletionRate)/getSBPMeasure();}
	double getIDRatio() {return InsertionRate/DeletionRate;}
	protected:
	//Basepair and insertion mutation rates are kept on an equal scaling relative to themselves
	vector<double> SBPTable={0,0,0,0,0,0};
	double InsertionRate=1;
	double DeletionRate=1;
	double MultipleInsertionDeletionRate=0;
	map<int,double> InsertionSpectrum={};
	map<int,double> DeletionSpectrum={};
};

vector<double> mutationInformation::getSBP(){
	vector<double> returnVector;
	for(int i=0; i<6; i++){
		returnVector.push_back(SBPTable[i]);
	}
	return returnVector;
}

#endif
