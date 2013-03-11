#ifndef PELUSA_H_
#define PELUSA_H_

#include<string>
#include "PelusaOverlapper.h"

#include<boost/program_options.hpp>
namespace po = boost::program_options;

using namespace std;

class PelusaApp
{
  public:
    PelusaApp( );
    
    int run(int argc, char **argv);
    
  private:
  	bool debug;
  	PelusaOverlapper * overlapper;
    bool processOptions( int argc, char **argv  );
    void printUsage(const po::options_description& desc);
    void setupDefaults();
};

#endif /*PELUSA_H_*/
