#ifndef SOAPSHORTHIT_H_
#define SOAPSHORTHIT_H_

#include <string>
#include <stdlib.h>
#include <iostream>
#include <boost/tokenizer.hpp>

using namespace std;


class SoapShortHit
{
  public:
    string query_id;
    int query_start;
    int query_end;
    int num_hits;
    int query_length;
    int full_query_length; // needed in case we've split up the sequence.
    char target_strand;
    // will probably want to store this as a UID look-up
    // (i.e. use helper class to dynamically assign UIDs to strings)
    string target_id;
    int target_start;
    int target_end;
    int target_length;
    int soapId;
    
    SoapShortHit() : query_id(""), query_start(0), query_end(0),
      num_hits(0), query_length(0), target_strand('+'), target_id(""),
      target_start(0), target_end(0), target_length(0), rep() {
      }
      
    ~SoapShortHit(){
    	// cerr << "Destroying " << this->toString() << endl;
    } 
    void getQueryName( string& name ) const;
    
    bool onPositiveTargetStrand(void);
    bool appendAble(SoapShortHit * otherHit);
    bool append(SoapShortHit * otherHit);
    // for debugging
    const string& toString();
    
    static SoapShortHit* parseLine( const string& line );
    private:
      string rep;
      void reverseTargetCoords(void); 
      static void getPrefix( const string& s, char delim, string& prefix );
      static void getSuffix( const string& s, char delim, string& suffix );
};

#endif /*SOAPSHORTHIT_H_*/
