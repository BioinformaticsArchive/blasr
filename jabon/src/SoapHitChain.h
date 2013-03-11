#ifndef SOAPHITCHAIN_H_
#define SOAPHITCHAIN_H_

#include <vector>
#include <sstream>
#include <iostream>
#include"SoapShortHit.h"

class SoapHitChain
{
	// PORT JOIN_CUTOFF 6 and DIAGONAL_CUTOFF 3
	
  public:
  	string query_id;
  	string target_id;
  	int query_integer; 
  	int target_integer; // parsed from the query_id and target_id
    int score;
    SoapShortHit min_hit, max_hit;
  	bool positive_strand; // true = "+" false  = "-"
  	vector<SoapShortHit*> hits;

    SoapHitChain(string query_id, string target_id, int score);
    
    void addHit( SoapShortHit *hit );   
    int getLength();
    int compare( SoapShortHit * hit1, SoapShortHit * hit2);
    const string toAmosString();
};

#endif /*SOAPHITCHAIN_H_*/
