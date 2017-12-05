#include "sequence.h"
#include <iostream>

int main() { //Test translationTables
	translationTable defaultTable=translationTable();
	translationTable bacteriaarcheaTransTable=translationTable(11);
	cout << defaultTable.Decode(codon("ATG")) << endl;
	cout << bacteriaarcheaTransTable.Decode(codon("GTG")) << endl;
	cout << defaultTable.Decode(codon("AUG")) << endl;
	cout << defaultTable.StartCodons().homologs[2].basepairs() << endl;
	cout << bacteriaarcheaTransTable.StartCodons().homologs[2].basepairs() << endl;
	//Test mutationProfile
	mutationProfile blankProfile();
	mutationProfile coliProfile("EColiMutationProfile.mprof");
	cout << coliProfile.getATGC() << endl;
	cout << coliProfile.getTsTv() << endl;
	coliProfile.setATGC(4.5);
	cout << coliProfile.getATGC() << endl; //Ensure the value is as set
	cout << coliProfile.getTsTv() << endl; //TsTv doesn't change much
	coliProfile.setTsTv(0.01); 
	cout << coliProfile.getATGC() << endl; //An idea of precision limits
	cout << coliProfile.getTsTv() << endl; //...but note how the previously set value changes relatively little
	coliProfile.Save("TestSave.mprof");
	//Test sequence::getGenes()
	sequence tuberculosis("MycobacteriumTuberculosisFASTA.txt");
	cout << tuberculosis[5].basepairs() << endl;
	cout << tuberculosis[10000000].basepairs() << endl; //Shouldn't work, sequence not circular
	tuberculosis = sequence("MycobacteriumTuberculosisFASTA.txt",1,true); //Let's fix that!
	cout << tuberculosis[10000000].basepairs() << endl;
	vector<gene> MTuberculosisGenes = tuberculosis.getGenes();
	cout << MTuberculosisGenes.size() << endl;
}
