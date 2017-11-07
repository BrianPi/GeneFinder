#ifndef TRANSLATION_TABLE_H
#define TRANSLATION_TABLE_H

#include <fstream>
#include <vector>
#include <string>
#include "codon.h"
#define DEFAULT_TABLE	"transl_tables/transl_table1.txt"
#define DEFAULT_TABLES "transl_tables/transl_table"
using namespace std;

struct codonDef{ //Defines the association between an amino acid or start/stop code and a set of codons
	char aminoCode;
	vector<codon> homologs;
};
static const codonDef nullCodon = {.aminoCode=codon::nullChar}; //Codon for returns to invalid calls.

class translationTable{ //Reads in the NCBI formatted translation table by trans_table integer, or a custom file in that format. Holds the relevant amino acid and syntactic associations.
	private:
	vector<codonDef> AAtable;
	codonDef StartCodon;
	codonDef StopCodon;
	public:
	void Load(string filename); //Load from named file
	//void Save(string filename);
	translationTable(string filename){Load(filename);} //Construct from named file
	translationTable(int NCBInumber){Load(DEFAULT_TABLES+to_string(NCBInumber)+".txt");} //Construct from NCBI number
	translationTable(){Load(DEFAULT_TABLE);} //Construct default table
	char Decode(codon unknownCodon); //Return the amino acid code for a given codon
	codonDef Define(char aminoCode); //Return the definition for a given amino acid code
	codonDef StartCodons(){ return this->StartCodon; } //Return the start codon definition
	codonDef StopCodons(){ return this->StopCodon; } //Return the stop codon definition
	vector<codonDef> Definitions(){ return this->AAtable; } //Return the amino acid definition table
	//codonDef returns may be removed
	bool isStart(codon unknownCodon);
	bool isStop(codon unknownCodon);
};

void translationTable::Load(string filename){
	ifstream inputTable = ifstream(filename);
	string filelines[5];
	for(int i=0;i<5;i++){
		getline(inputTable,filelines[i]);
		if(filelines[i].find("= ")!=string::npos){
			filelines[i]=filelines[i].erase(0,filelines[i].find("= ")+2);
		}
	}
	this->AAtable={};
	for(int i=0;i<filelines[0].length();i++){
		codon basepairs = codon({filelines[2][i],filelines[3][i],filelines[4][i]});
		bool included=false;
		for(int j=0;j<AAtable.size();j++){
			if(filelines[0][i]==AAtable[j].aminoCode){
				AAtable[j].homologs.push_back(basepairs);
				included=true;
				break;
			}
		}
		if(!included){
			codonDef newCodon=codonDef();
			newCodon.aminoCode=filelines[0][i];
			newCodon.homologs.push_back(basepairs);
			AAtable.push_back(newCodon);
		}
	}
	this->StopCodon.aminoCode='*';
	this->StopCodon.homologs={};
	this->StartCodon.homologs={};
	for(int i=0;i<filelines[1].length();i++){
		if(filelines[1][i]!='-'){
			codon basepairs = codon({filelines[2][i],filelines[3][i],filelines[4][i]});
			if(filelines[1][i]==StopCodon.aminoCode){
				StopCodon.homologs.push_back(basepairs);
			}
			else{
				StartCodon.aminoCode=filelines[1][i];
				StartCodon.homologs.push_back(basepairs);
			}
		}
	}
}

char translationTable::Decode(codon unknownCodon){
	for(int i=0;i<AAtable.size();i++){
		for(int j=0;j<AAtable[i].homologs.size();j++){
			if(unknownCodon==AAtable[i].homologs[j]){
				return AAtable[i].aminoCode;
			}
		}
	}
	return nullCodon.aminoCode;
}

codonDef translationTable::Define(char aminoCode){
	for(int i=0;i<AAtable.size();i++){
		if(AAtable[i].aminoCode==aminoCode){
			return AAtable[i];
		}
	}
	return nullCodon;
}

bool translationTable::isStart(codon unknownCodon){
	for(int i=0;i<StartCodon.homologs.size();i++){
		if(unknownCodon==StartCodon.homologs[i]){
			return true;
		}
	}
	return false;
}

bool translationTable::isStop(codon unknownCodon){
	for(int i=0;i<StopCodon.homologs.size();i++){
		if(unknownCodon==StopCodon.homologs[i]){
			return true;
		}
	}
	return false;
}
#endif
