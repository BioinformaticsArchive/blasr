#ifndef JABONAPP_H_
#define JABONAPP_H_

#include<string>

#include<boost/program_options.hpp>
namespace po = boost::program_options;

using namespace std;

class JabonApp
{
  public:
    JabonApp( int argc, char **argv );
    
    int run();
    
  private:
    void processOptions( int argc, char **argv );
    void printUsage(const po::options_description& desc);
    void setupDefaults();
    
    bool debug;
    bool shouldExit;
    string queryFile;
    string targetFile;
};

#endif /*JABONAPP_H_*/
