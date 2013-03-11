#include "FeatureEncoder.h"

FeatureEncoder::FeatureEncoder(int k)
	:	k(k),
//		baseString("ACGT"), need for featureToKmer not yet implemented
		shifted(0)
{
	initialize();
}

void FeatureEncoder::initialize()
{
	for (int idx=0; idx < BASE2UINTLENGTH; idx++)
	{
		base2uint[idx] = -1;
	}
	base2uint[(int)'A'] = base2uint[(int)'a'] = 0;
	base2uint[(int)'C'] = base2uint[(int)'c'] = 1;
	base2uint[(int)'G'] = base2uint[(int)'g'] = 2;
	base2uint[(int)'T'] = base2uint[(int)'t'] = 3;

	unsigned int finalSetBit = (unsigned int)(k-1)*2;
	for (unsigned int idx=0; idx < finalSetBit; idx++)
	{
		shifted |= 1 << idx;
	}
	
	numFeatures = 1 << (k*2);
}

void FeatureEncoder::encode(string const& input, vector<unsigned int>* features)
{
	assert(input.length() >= (unsigned int) k);
	features->push_back( feature(input, 0) );
	for(unsigned int featureStart = 1; 
			featureStart < input.length() - k + 1;
			featureStart++)
	{
		unsigned int lastFeature = (*features)[features->size()-1];
		unsigned int feature = 
			nextFeature( lastFeature, input[ featureStart + k - 1 ]);
		features->push_back( feature );
	}
}

unsigned int FeatureEncoder::feature(string const string, int stringIdx)
{
	unsigned int feature = 0;
	for (int seqIdx = stringIdx; seqIdx < stringIdx+k; seqIdx++)
	{
		feature <<= 2;
		feature += base2uint[ (int)string[seqIdx] ];
	}
	return feature;
}

unsigned int FeatureEncoder::nextFeature(unsigned int lastFeature, char ch)
{
	unsigned int newFeature = lastFeature & shifted;
	newFeature <<= 2;
	newFeature += base2uint[ (int)ch ];
	return newFeature;
}

int FeatureEncoder::getNumFeatures()
{
	return numFeatures;
}
