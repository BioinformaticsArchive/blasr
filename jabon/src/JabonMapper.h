#ifndef JABONMAPPER_H_
#define JABONMAPPER_H_

#include <iostream>
#include <iomanip>
#include <ostream>
#include <fstream>
#include <string>
#include <map>
#include <pthread.h>
#include <boost/tokenizer.hpp>

#include "reads.h"
#include "dbseq.h"
#include "align.h"
#include "utilities.h"

#include "JabonException.h"
#include "SoapMappingHit.h"
#include "SoapShortHit.h"

using namespace std;

class JabonMapper
{
  public:
    JabonMapper( const string& queryFile );
    ~JabonMapper();
    void initialize( const string& targetFile );
    void run();
    
    
    
  private:
    RefSeq ref;
    ReadClass read_a;
    string queryFile;
    ifstream fin_a;
    
    bit32_t n_aligned;   //number of reads aligned
    
    void Do_Formatdb();
    void RunProcess();
    void DoSingleAlign();
};

class JabonMapperWorker
{
  public:
    JabonMapperWorker();
    
    void run();
    void postProcessReads();
    void postProcessReadHits(vector<SoapShortHit*> short_hits);

  private:
    SingleAlign a;

};

#endif /*JABONMAPPER_H_*/
