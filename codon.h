#ifndef CODON_H
#define CODON_H

class codon{
	private:
	char bases[3]={nullChar,nullChar,nullChar};
	public:
	static const char nullChar='-';
	codon(std::string basepairs);
	char operator[](unsigned int i) const //Indexed from 1, view first, second or third basepair
		{ if(i<4 && i!=0) return bases[i-1]; return nullChar; }
	std::string basepairs()
		{ return {bases[0],bases[1],bases[2]}; }
	bool operator==(const codon &otherCodon) const {
		for(int i=0;i<3;i++){
			if(bases[i]!=otherCodon[i+1]){
				return false;
			}
		}
		return true;
	}
	//Operator != omitted out of deference to use of nullChar for missing information. While 'ATG' isn't equivalent to 'A-G' it's not nonequivalent either.
};

codon::codon(std::string basepairs){
	if(basepairs.length()<4){ //Don't construct a codon if input string doesn't meet expectations
		for(int i=0;i<basepairs.length();i++){
			if(std::toupper(basepairs[i])=='A'||std::toupper(basepairs[i])=='C'||std::toupper(basepairs[i])=='G'||std::toupper(basepairs[i])=='T'){
				bases[i]=std::toupper(basepairs[i]);
			}
			else if(std::toupper(basepairs[i]=='U')){
				bases[i]='T';
			}
		}
	}
}

#endif
