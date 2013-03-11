#ifndef FEATUREENCODER_H_
#define FEATUREENCODER_H_

#include <vector>
#include <string>
#include <assert.h>

#define BASE2UINTLENGTH 256

using namespace std;

// singleton class TODO ensure in C++
class FeatureEncoder
{
	public:
    	FeatureEncoder(int k);
    	//~PelusaOverlapper();
		void encode(string const& input, vector<unsigned int>* features);
		int getNumFeatures();
		
	private:
		int k;
		int numFeatures;
//		string baseString; need for featureToKmer, not yet implemented
		unsigned int shifted; // TODO unsigned below as well?
		unsigned int base2uint[BASE2UINTLENGTH];
	    void initialize();
		unsigned int feature(string const string, int stringIdx);
		unsigned int nextFeature(unsigned int lastFeature, char ch);
		
};
#endif /*FEATUREENCODER_H_*/
