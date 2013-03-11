#include<iostream>
#include "param.h"
#include "JabonApp.h"
#include "JabonMapper.h"
#include "JabonException.h"

using namespace std;

static char *DESCRIPTION = 
  "Jabon: Fast Long-read Mapper\n"
  "Pacific Biosciences (c) 2008\n\n"
  "An alignment algorithm for quickly aligning long reads to a\n"
  "large genomic sequence.";

// since the original SOAP code heavily uses global variables
// we have to keep up the poor design...
Param param;

JabonApp::JabonApp( int argc, char **argv )
{
  shouldExit = false;
  setupDefaults();
  processOptions( argc, argv );
}

void JabonApp::setupDefaults()
{
  param.short_format = true;
  param.num_procs = 2;
  param.chopReads = true;
  param.report_repeat_hits = 2;
  param.max_gap_size = 2;
  param.max_snp_num = 5;
}

void JabonApp::processOptions( int argc, char **argv )
{
  po::options_description desc( "Allowed options" );
  try
  {
    desc.add_options()
      ( "query,q", po::value<string>()->default_value(""), "Query file in FASTA format" )
      ( "target,t", po::value<string>()->default_value(""), "Target file in FASTA format" )
      ( "windowSize,w", po::value<int>()->default_value(25), "Window size for chopping reads" )
      ( "windowStep,d", po::value<int>()->default_value(1), "Window step for chopping reads" )
      ( "seed,s", po::value<int>()->default_value(10), "Seed size for SOAP algorithm" )
      ( "diagonalCutoff,g", po::value<int>()->default_value(3), "Maximum deviation in alignment diagonal to chain two matches" )
      ( "joinCutoff,j", po::value<int>()->default_value(6), "Maximum distance to chain two matches" )
      ( "minChainLength,n", po::value<int>()->default_value(1), "Minimum number of matches per chain" )
      ( "help", "Show this help message" )
      ( "debug", po::bool_switch()->default_value(false), "Enable debugging output" )
      ( "nproc", po::value<int>()->default_value(2), "Number of processors to use" )
    ;
    
    po::variables_map vm;
    po::store( po::parse_command_line( argc, argv, desc ), vm );
    po::notify( vm );
    
    if ( vm.count( "help" ) )
    {
      printUsage( desc );
      shouldExit = true;
      return;
    }
     
    debug = vm["debug"].as<bool>();
    
    if ( vm.count( "query" )==0 || vm["query"].as<string>().length()==0 )
    {
      cerr << "!! No query file specified" << endl << endl;
      printUsage( desc );
      shouldExit = true;
      return;
    }
    else
    {
      queryFile = vm["query"].as<string>();
    }
    
    if ( vm.count( "target" )==0 || vm["target"].as<string>().length()==0 )
    {
      cerr << "!! No target file specified" << endl << endl;
      printUsage( desc );
      shouldExit = true;
      return;
    }
    else
    {
      targetFile = vm["target"].as<string>();
    }
    
    if ( vm.count("nproc")>0 )
    {
      param.num_procs = vm["nproc"].as<int>();
    }
    
    if ( vm.count("windowSize")>0 )
    {
      param.longReadWindowSize = vm["windowSize"].as<int>();
    }
    
    if ( vm.count("windowStep")>0 )
    {
      param.longReadWindowStep = vm["windowStep"].as<int>();
    }
    
    if ( vm.count("seed")>0 )
    {
      param.SetSeedSize( vm["seed"].as<int>() );
    }
    
    if ( vm.count("joinCutoff")>0 )
    {
    	param.joinCutoff = vm["joinCutoff"].as<int>();
    }
    
    if ( vm.count("diagonalCutoff")>0 )
    {
      param.diagonalCutoff = vm["diagonalCutoff"].as<int>();
    }
    
    if ( vm.count("minChainLength")>0 )
    {
      param.minChainLength = vm["minChainLength"].as<int>();
    }
  }
  catch( po::unknown_option u )
  {
    cerr << "!! " << u.what() << endl << endl;
    printUsage( desc );
    shouldExit = true;
  }
}

void JabonApp::printUsage( const po::options_description& desc )
{
  cout << DESCRIPTION << endl;
  cout << endl << desc << endl;
}

int JabonApp::run()
{
  if ( shouldExit )
  {
    return 1;
  }
  
  JabonMapper mapper( queryFile );
  mapper.initialize( targetFile );
  mapper.run();
  
  return 0;
}

int main( int argc, char **argv )
{
  try
  {
    JabonApp app( argc, argv );
    //param.overlapper = true; TODO expose as a parameter? Useful when aligning the same FASTA as target and query.
    return app.run();
  }
  catch( JabonException e )
  {
    cerr << "Caught exception: " << e.getMessage() << endl;
    return 1;
  }
}
