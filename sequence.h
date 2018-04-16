#include "translationTable.h"
#include "mutationProfile.h"

class gene; //Abstract access to genes stored in a sequence.

class sequence{ //Stores FASTA data
	private:
	translationTable translationKey=translationTable();
	mutationProfile mutationKey=mutationProfile();
	std::string rawsequence={};
	bool circular;
	public:
	sequence(string filename, unsigned int n=1, bool loop=false){Load(filename,n);circular=loop;}
	void Load(string filename, unsigned int n=1); //Load nth FASTA sequence in file
	void Save(string filename); //UNIMPLEMENTED!
	vector<gene> getGenes();
	codon operator[](unsigned int index) const;
};

class gene{ 
	public:
	gene(unsigned int start, unsigned int length, sequence* source){startIndex=start;codonLength=length;genome=source;} //Constructs gene from start index and length
	//Consider extending to deal with start/end mutation, noncoding sequences...
	codon operator[](unsigned int index) const;
	unsigned int start(){return startIndex;}
	unsigned int length(){return codonLength;}
	private:
	unsigned int startIndex;
	unsigned int codonLength;
	sequence* genome;
};

codon sequence::operator[](unsigned int index) const{
	string basepairs={};
	for(int i=0; i<3; i++){
		if((index+i)<rawsequence.size()){
			basepairs.push_back(rawsequence[index+i]);
		}
		else if(circular){
			basepairs.push_back(rawsequence[(index+i)%rawsequence.size()]); //Loop around if genome is circular
		}
		else{
			basepairs.push_back(codon::nullChar);
		}
	}
	return codon(basepairs);
}

void sequence::Load(string filename, unsigned int n){
	unsigned int sequence_number=0;
	ifstream inputfile=ifstream(filename);
	while(sequence_number<n){
		string testline;
		getline(inputfile,testline);
		if(testline[0]=='>'){
			sequence_number++;
		}
	}
	string buffer;
	while(getline(inputfile,buffer)){
		if(buffer.length()>0){
			rawsequence.append(buffer);
		}
	}
}

vector<gene> sequence::getGenes(){
	vector<unsigned int> geneStarts;
	vector<gene> returngenes={};
	int limit;
	if(circular){
		limit=rawsequence.size();
	}
	else{
		limit=rawsequence.size()-2;
	}
	for(unsigned int offset=0;offset<3;offset++){ //Loop through all three coding frames
		geneStarts={};
		unsigned int index=0;
		for(;index+offset<limit; index=index+3){
			if(translationKey.isStart(operator[](index+offset))){
				geneStarts.push_back(index+offset);
			}
			else if(translationKey.isStop(operator[](index+offset))){
				for(int i=0; i<geneStarts.size(); i++){
					gene newgene(geneStarts[i],(index+offset-geneStarts[i])/3,this);
					returngenes.push_back(newgene);
				}
				geneStarts={};
			}
		}
		if(circular){
			for(;index+offset<2*rawsequence.size();index=index+3){ //Keep going where loop left off
				//Initializing new genes would be redundant
				if(translationKey.isStop(operator[](index+offset))){
					for(int i=0; i<geneStarts.size(); i++){
						gene newgene(geneStarts[i],(index+offset-geneStarts[i])/3,this);
						returngenes.push_back(newgene);
					}
					geneStarts={};
				}
			}
		}
	}
	return returngenes;
}

codon gene::operator[](unsigned int index) const{
	if(index<=codonLength){	
		return (*genome)[startIndex+3*index];
	}
	else{
		return codon("");
	}
}
