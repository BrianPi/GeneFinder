class rescaleFactors{		
	public:
		static vector<double> rescale(vector<double>* vectorPointer, double desiredTotal) {
			double currentTotal=accumulate(vectorPointer->begin(),vectorPointer->end(),0);
			if(currentTotal!=0){
				double multiplier=desiredTotal/currentTotal;
				for(int i=0; i<vectorPointer->size(); i++){
					(*vectorPointer)[i]*=multiplier;
				}
			}
			else{
				for(int i=0; i<vectorPointer->size(); i++){
					(*vectorPointer)[i]=desiredTotal/vectorPointer->size();
				}
			}
			return *vectorPointer;
		}
};
	
