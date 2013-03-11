#ifndef ALIGNMENT_SAM_ALIGNMENT_H_
#define ALIGNMENT_SAM_ALIGNMENT_H_

#include <string>
#include <sstream>

#include "SAMKeywordValuePair.h"
#include "datastructures/alignment/CigarString.h"

class SAMAlignment {
 public:
  enum SAMAlignmentRequiredFields {
    S_QNAME,
    S_FLAG,
    S_RNAME,
    S_POS,
    S_MAPQV,
    S_CIGAR,
    S_RNEXT,
    S_PNEXT,
    S_TLEN,
    S_SEQ,
    S_QUAL};

  static const char* SAMAlignmentRequiredFieldNames[];

  string qName;
  unsigned int flag;
  string rName;
  unsigned int pos;
  short mapQV;
  CigarString cigar;
  string rNext;
  int pNext;
  int tLen;
  string seq;
  string qual;

  float score;
  int   xs, xe;
  int xl;
  string rg;
  int   as;
  int   xt;
  int   nm;
  int   fi;
  
  SAMAlignment() {
    //
    // Initialize all optional fields.  Required fields will be
    // assigned a value later.
    //
    score = xs = xe = as = xt = nm = fi = xl = 0;
    rg = "";
  }
  bool StoreValues(string &line,  int lineNumber=0) {
    stringstream strm(line);
    vector<bool> usedFields;
    usedFields.resize(S_QUAL);
    fill(usedFields.begin(), usedFields.end(), false);
    string kvPair;
    int i;
    bool parseError = false;
    SAMAlignmentRequiredFields field;
    //
    // Define a temporary mapqv value that gets over a GMAP bug that prints a mapqv < 0.
    //
    int tmpMapQV;
    if (!(strm >> qName)) {
      parseError = true;
      field = S_QNAME;
    }
    else if (! (strm >> flag) ){ 
      parseError = true;
      field = S_FLAG;
    }
    else if (! (strm >> rName) ) {
      parseError = true;
      field = S_RNAME;
    }
    else if (! (strm >> pos) ) {
      parseError = true;
      field = S_POS;
    }
    else if (! (strm >> tmpMapQV)) {
      parseError = true; field = S_MAPQV;
    }
    else if (! (strm >> cigar)) {
      parseError = true; field = S_CIGAR;
    }

    else if (! (strm >> rNext)) {
      parseError = true; field = S_RNEXT;
    }

    else if (! (strm >> pNext)) {
      parseError = true; field = S_PNEXT;
    }
    else if (! (strm >> tLen)) {
      parseError = true; field = S_TLEN;
    }
    else if (! (strm >> seq)) {
      parseError = true; field = S_SEQ;
    }
    else if (! (strm >> qual)) {
      parseError = true; field = S_QUAL;
    }

    mapQV = (unsigned char) tmpMapQV;
    //
    // If not aligned, stop trying to read in elements from the sam string.
    //
    if (rName == "*") {
      return true;
    }

    if (parseError) {
      cout << "Error parsing alignment line " << lineNumber << ". Missing or error in field " << SAMAlignmentRequiredFieldNames[field] << endl;
      exit(1);
    }
    
    //
    // Now parse optional data.
    //
    while (strm) {
      string kvName, kvType, kvValue;
      string typedKVPair;
      if ((strm >> typedKVPair) == 0) {
        break;
      }
      if (TypedKeywordValuePair::Separate(typedKVPair, kvName, kvType, kvValue)) {
        stringstream strm(kvValue);
        if (kvName == "RG") {
          rg = kvValue;
        }
        else if (kvName == "AS") {
          strm >> as;
        }
        else if (kvName == "XS") {
          strm >> xs;
        }
        else if (kvName == "XE") {
          strm >> xe;
        }
        else if (kvName == "XL") {
          strm >> xl;
        }
        else if (kvName == "XT") {
          strm >> xt;
        }
        else if (kvName == "NM") {
          strm >> nm;
        }
        else if (kvName == "FI") {
          strm >> fi;
        }
      }
      else {
        cout <<"ERROR.  Could not parse typed keyword value " << typedKVPair << endl;
      }
    }
  }
};

const char* SAMAlignment::SAMAlignmentRequiredFieldNames[] = { "QNAME", "FLAG", 
                                                               "RNAME", "POS", 
                                                               "MAPQV", "CIGAR", 
                                                               "RNEXT", "PNEXT", 
                                                               "TLEN", "SEQ", "QUAL"} ;


using namespace std;

class SAMPosAlignment : public SAMAlignment {
 public:
  unsigned int qStart, qEnd;
  unsigned int tStart, tEnd;
  int qStrand, tStrand;
};



/* 
 * Write the full one later
 */




  

#endif
