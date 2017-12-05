#ifndef MUTATION_PROFILE_H
#define MUTATION_PROFILE_H

#include <fstream>
#include "mutationInformation.h"
#include "rescaleFactors.h"

class mutationProfile : public mutationInformation {
//Reads in the custom mutation profile format and allows adjustments to it
	public:
	mutationProfile(string filename){Load(filename);} //Construct from named profile
	mutationProfile() {}
	void setTsTv(double ratio);
	void setATGC(double ratio);
	void setSBP(double SBPArray[]);
	void setIDSBPRatio(double ratio);
	void setIDRatio(double ratio);
	void Load(string filename); //Load from named file
	void Save(string filename); //Save to named file
	private:
	bool LoadFile(ifstream* stream); //Load data from mutation profile. Returns false on badly formatted input.
	bool LoadSBPID(ifstream* stream, string type); //Load single base pair mutation, insertion, or deletion as specified from ifstream
	void LoadMNID(ifstream* stream); //Load multiple base pair insertion/deletion spectrum from ifstream
};

void mutationProfile::setTsTv(double ratio){ //Set Transition-Transversion to a given ratio
	double sum=getSBPMeasure();
	double desiredDenominator=sum/(1+ratio);
	vector<double> denominator(SBPTable.begin()+2,SBPTable.end());
	rescaleFactors::rescale(&denominator,desiredDenominator);
	SBPTable.erase(SBPTable.begin()+2,SBPTable.end());
	SBPTable.insert(SBPTable.end(),denominator.begin(),denominator.end());
	double desiredNumerator=sum-desiredDenominator;
	vector<double> numerator(SBPTable.begin(),SBPTable.begin()+2);
	rescaleFactors::rescale(&numerator,desiredNumerator);
	SBPTable.erase(SBPTable.begin(),SBPTable.begin()+2);
	SBPTable.insert(SBPTable.begin(),numerator.begin(),numerator.end());
}

void mutationProfile::setATGC(double ratio){ //Set AT to GC mutation to a given ratio
	double sum=getSBPMeasure()-SBPTable[2]-SBPTable[6];
	double desiredDenominator=sum/(1+ratio);
	double desiredNumerator=sum-desiredDenominator;
	vector<double> numerator={SBPTable[0],SBPTable[4]};
	vector<double> denominator={SBPTable[1],SBPTable[3]};
	rescaleFactors::rescale(&denominator,desiredDenominator);
	rescaleFactors::rescale(&numerator,desiredNumerator);
	SBPTable[0]=numerator[0];
	SBPTable[4]=numerator[1];
	SBPTable[1]=denominator[0];
	SBPTable[3]=denominator[1];
}

void mutationProfile::setSBP(double SBPArray[]){
	for(int i=0; i<6; i++){
		SBPTable[i]=SBPArray[i];
	}
}

void mutationProfile::setIDRatio(double ratio){ //Sets Insertion/Deletion ratio for single basepairs
//succeeds trivially when there are no insertions or deletions
	double sum=InsertionRate+DeletionRate;
	DeletionRate=sum/(1+ratio);
	InsertionRate=sum-DeletionRate;
}

void mutationProfile::Load(string filename){
	ifstream inputProfile = ifstream(filename);
	for(int i=0; i<6; i++){
		SBPTable[i]=0;
	}
	InsertionRate=0;
	DeletionRate=0;
	MultipleInsertionDeletionRate=0;
	InsertionSpectrum={};
	DeletionSpectrum={};
	LoadFile(&inputProfile);
}

bool mutationProfile::LoadFile(ifstream* stream){
	while(!stream->eof()){
		string scrapedchars={};
		char newCharacter;
		stream->get(newCharacter);
		while(newCharacter!='\n'&&!stream->eof()){
			if(newCharacter==' '){
				stream->get(newCharacter);
				continue;
			}
			else if(isspace(newCharacter)){ //Any whitespace other than a plain space
				for(int i=0; i<scrapedchars.length(); i++){ //Capitalize tag letters
					scrapedchars[i]=toupper(scrapedchars[i]);
				}
				if(scrapedchars=="MNID"){
					LoadMNID(stream);
				}
				else{ 
					if(!LoadSBPID(stream,scrapedchars)){ //Load stream information into given tag
						return false; //Return false if badly formatted
					}
				}
				break;
			}
			else if(isdigit(newCharacter)){
				return false;
			}
			else if(isalpha(newCharacter)){
				scrapedchars.push_back(newCharacter);
			}
			stream->get(newCharacter);
		}
		if(scrapedchars=="MNID"){
			LoadMNID(stream);
		}
	}
	return true;
}

bool mutationProfile::LoadSBPID(ifstream* stream, string type){
	//The remainder of the line is the number to be loaded
	string number;
	getline(*stream,number);
	switch(type.size()){
		//Four characters used for transition, remove second and fourth for the two character case
		case 4:
			type.erase(1,1);
			type.erase(2,1);
		case 2:
			if(type=="AG"){
				SBPTable[0]+=atof(number.c_str());
			}
			else if(type=="GA"){
				SBPTable[1]+=atof(number.c_str());
			}
			else if(type=="AT"){
				SBPTable[2]+=atof(number.c_str());
			}
			else if(type=="AC"){
				SBPTable[3]+=atof(number.c_str());
			}
			else if(type=="GT"){
				SBPTable[4]+=atof(number.c_str());
			}
			else if(type=="GC"){
				SBPTable[5]+=atof(number.c_str());
			}
			else{
				return false;
			}	
			break;
		case 3:
			if(type=="SNI"){
				InsertionRate+=atof(number.c_str());
			}
			else if(type=="SND"){
				DeletionRate+=atof(number.c_str());
			}
			else{
				return false;
			}
			break;
		default:
			return false;
	}
	return true;
}

void mutationProfile::LoadMNID(ifstream* stream){
	while(!stream->eof()){
		string nextLine;	
		getline(*stream,nextLine);
		size_t tabLoc=nextLine.rfind('\t');
		if(tabLoc==string::npos){
			return;
		}
		double deletionValue=atof(nextLine.substr(tabLoc+1).c_str());
		nextLine.erase(tabLoc);
		tabLoc=nextLine.rfind('\t');
		if(tabLoc==string::npos){
			return;
		}
		double insertionValue=atof(nextLine.substr(tabLoc+1).c_str());
		nextLine.erase(tabLoc);
		string basepairnumstr={};
		for(int i=0; i<nextLine.length(); i++){
			if(!isspace(nextLine[i])){
				basepairnumstr.push_back(nextLine[i]);
			}
		}
		int basepairnumber=stoi(basepairnumstr.c_str());
		MultipleInsertionDeletionRate+=insertionValue+deletionValue;
		pair<int,double> insertion(basepairnumber,insertionValue);
		pair<int,double> deletion(basepairnumber,deletionValue);
		InsertionSpectrum.insert(insertion);
		DeletionSpectrum.insert(deletion);
	}
}

void mutationProfile::Save(string filename){
	ofstream outputProfile=ofstream(filename);
	outputProfile << "A->G\t" << to_string(SBPTable[0]) << endl;
	outputProfile << "G->A\t" << to_string(SBPTable[1]) << endl;
	outputProfile << "A->T\t" << to_string(SBPTable[2]) << endl;
	outputProfile << "A->C\t" << to_string(SBPTable[3]) << endl;
	outputProfile << "G->T\t" << to_string(SBPTable[4]) << endl;
	outputProfile << "G->C\t" << to_string(SBPTable[5]) << endl;
	outputProfile << "SNI\t" << to_string(InsertionRate) << endl;
	outputProfile << "SND\t" << to_string(DeletionRate) << endl;
	if(InsertionSpectrum.size()+DeletionSpectrum.size()>0){
		outputProfile << "MNID" << endl;
		map<int,double>::iterator dit=DeletionSpectrum.begin();
		map<int,double>::iterator iit=InsertionSpectrum.begin();
		int inum;
		int dnum;
		while(iit!=InsertionSpectrum.end()&&dit!=DeletionSpectrum.end()){
			inum=iit->first;
			dnum=dit->first;			
			if(inum==dnum){
				outputProfile << to_string(inum) << '\t' << to_string(iit->second) << '\t' << to_string(dit->second) << endl;
				iit++;
				dit++;
			}
			else if(inum>dnum){
				outputProfile << to_string(dnum) << '\t' << "0" << '\t' << to_string(dit->second) << endl;
				dit++;
			}
			else{
				outputProfile << to_string(inum) << '\t' << to_string(iit->second) << '\t' << "0" << endl;
				iit++;
			}
		}
		if(iit==InsertionSpectrum.end()){
			while(dit!=DeletionSpectrum.end()){
				outputProfile << to_string(dnum) << '\t' << "0" << '\t' << to_string(dit->second) << endl;
				dit++;
			}
		}
		else{
			while(iit!=InsertionSpectrum.end()){
				outputProfile << to_string(inum) << '\t' << to_string(iit->second) << '\t' << "0" << endl;
				iit++;
			}
		}
	}
}

#endif
