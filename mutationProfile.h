#ifndef MUTATION_PROFILE_H
#define MUTATION_PROFILE_H

#include <fstream>

class mutationProfile{ //Reads in the custom mutation profile format. Holds the relative frequency of base substitution and insertion mutations, with rudimentary options for multi-base insertions.
	public:
	void Load(string filename); //Load from named file
	void Save(string filename); //UNIMPLEMENTED!
	mutationProfile(string filename){Load(filename);} //Construct from named profile
	mutationProfile(){}
	double getTsTv() {return (SBPTable[0]+SBPTable[1])/(SBPTable[2]+SBPTable[3]+SBPTable[4]+SBPTable[5]);}
	//void setTsTv(double ratio);
	double getATGC() {return (SBPTable[0]+SBPTable[4])/(SBPTable[1]+SBPTable[3]);}
	//void setATGC(double ratio);
	void setSBP(double SBPArray[]);
	double getIDRate() {return InsertionDeletionRate;}
	//void setIDRate(double rate);
	double getIDRatio() {return InsertionDeletionRatio;}
	//void setIDRatio(double ratio);
	private:
	double SBPTable[6]={0,0,0,0,0,0};
	double InsertionRate=0;
	double DeletionRate=0;
	double InsertionDeletionRate=0;
	double InsertionDeletionRatio=1;
	vector<pair<int,double>> InsertionSpectrum={};
	vector<pair<int,double>> DeletionSpectrum={};
	bool LoadSBPID(ifstream* stream); //Load single base pair mutation, insertion, deletion from ifstream. Returns false on badly formatted input.
	void LoadMNID(ifstream* stream); //Load multiple base pair insertion/deletion spectrum from ifstream
};

void mutationProfile::Load(string filename){
	ifstream inputProfile = ifstream(filename);
	for(int i=0; i<6; i++){
		SBPTable[i]=0;
	}
	InsertionRate=0;
	DeletionRate=0;
	InsertionSpectrum={};
	DeletionSpectrum={};
	LoadSBPID(&inputProfile);
}

bool mutationProfile::LoadSBPID(ifstream* stream){
	while(stream->peek()!=char_traits<char>::eof()){
		string scrapedchars={};
		string nextLine;
		getline(*stream,nextLine);
		for(int i=0; i<nextLine.length(); i++){
			if(nextLine[i]==' '){
				continue;
			}
			else if(isspace(nextLine[i])){
				for(int i=0; i<scrapedchars.length(); i++){
					scrapedchars[i]=toupper(scrapedchars[i]);
				}
				switch(scrapedchars.length()){
					case 4:
						if(scrapedchars=="MNID"){
							LoadMNID(stream);
							break;
						}
						else{ //Four characters used for transition, take first and third and move to case 2.
							scrapedchars.erase(1,1);
							scrapedchars.erase(2,1);
						}
					case 2:
						if(scrapedchars=="AG"){
							SBPTable[0]+=atof(nextLine.substr(i+1).c_str());
						}
						else if(scrapedchars=="GA"){
							SBPTable[1]+=atof(nextLine.substr(i+1).c_str());
						}
						else if(scrapedchars=="AT"){
							SBPTable[2]+=atof(nextLine.substr(i+1).c_str());
						}
						else if(scrapedchars=="AC"){
							SBPTable[3]+=atof(nextLine.substr(i+1).c_str());
						}
						else if(scrapedchars=="GT"){
							SBPTable[4]+=atof(nextLine.substr(i+1).c_str());
						}
						else if(scrapedchars=="GC"){
							SBPTable[5]+=atof(nextLine.substr(i+1).c_str());
						}
						else{
							return false;
						}	
						break;
					case 3:
						if(scrapedchars=="SNI"){
							InsertionRate+=atof(nextLine.substr(i+1).c_str());
						}
						else if(scrapedchars=="SND"){
							DeletionRate+=atof(nextLine.substr(i+1).c_str());
						}
						else{
							return false;
						}
						break;
					default:
						return false;
				}
				break;
			}
			else if(isdigit(nextLine[i])){
				return false;
			}
			else if(isalpha(nextLine[i])){
				scrapedchars.push_back(nextLine[i]);
			}
		}
	}
	return true;
}

void mutationProfile::LoadMNID(ifstream* stream){
	while(stream->peek()!=char_traits<char>::eof()){
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
		pair<int,double> insertion(basepairnumber,insertionValue);
		pair<int,double> deletion(basepairnumber,deletionValue);
		InsertionSpectrum.push_back(insertion);
		DeletionSpectrum.push_back(deletion);
	}
}

void mutationProfile::Save(string filename){
	
}

#endif
