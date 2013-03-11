#ifndef SOAPMAPPINGHIT_
#define SOAPMAPPINGHIT_

#include <stdlib.h>
#include <iostream>

#include "SoapShortHit.h"
#include "SoapHitChain.h"
#include "GusfieldChainer.h"
#include "Rectangle.h"

using namespace std;

class SoapMappingHit
{
  // PORT static variable MIN_CHAIN_LENGTH 8 and MIN_SUPERCHAIN_LENGTH 1
  public:
    string query_id;
    string target_id;
    vector<SoapShortHit*> hits;
	int aligned_length;
	int number_insertions;
	int number_deletions;
	int number_mismatches;
    
    SoapMappingHit(string query_id, string target_id);
      
    void addHit(SoapShortHit * hit);
    SoapHitChain * maximumChain(int diagonalCutoff, 
								int joinCutoff, 
								int minChainLength);
      
    // for debugging
    const string& to_string();

	private:
      string rep;
};




#endif /*SOAPMAPPINGHIT_*/

