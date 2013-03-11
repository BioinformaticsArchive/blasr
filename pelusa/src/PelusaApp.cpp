#include <iostream>
#include <cmath>
#include "PelusaApp.h"
#include "PelusaOverlapper.h"
#include "PelusaException.h"


using namespace std;

static string DESCRIPTION = 
  "Pelusa: Read overlapper using Bloom filters\n"
  "Pacific Biosciences (c) 2009";


PelusaApp::PelusaApp( )
	: 	debug(false),
	  	overlapper(new PelusaOverlapper())
{

}


bool PelusaApp::processOptions( 
	int argc, 
	char **argv)
{
  po::options_description desc( "Allowed options" );
  try
  {
    desc.add_options()
      ( "query,q", po::value<string>()->default_value(""), "Query file in FASTA format" )
      ( "target,t", po::value<string>()->default_value(""), "Target file in FASTA format" )
      ( "kmerLength,k", po::value<int>()->default_value(8), "Kmer length" )
      ( "bloomWidth,m", po::value<int>()->default_value(2000), "Bloom filter width (should be divisible by 8)" )
      ( "numSegments,s", po::value<int>()->default_value(2), "Number query segments" )
      ( "topColumns,c", po::value<int>()->default_value(10), "Number of top columns to keep" )
      ( "nproc", po::value<int>()->default_value(1), "Number of processors to use" )
      ( "help", "Show this help message" )
      ( "debug", po::bool_switch()->default_value(false), "Enable debugging output" )
    ;
    
    po::variables_map vm;
    po::store( po::parse_command_line( argc, argv, desc ), vm );
    po::notify( vm );
    
    if ( vm.count( "help" ) )
    {
      printUsage( desc );
      return false;
    }
     
    overlapper->debug = vm["debug"].as<bool>();
    
    if ( vm.count( "query" )==0 || vm["query"].as<string>().length()==0 )
    {
      cerr << "!! No query file specified" << endl << endl;
      printUsage( desc );
      return false;
    }
    else
    {
      overlapper->queryFile = vm["query"].as<string>();
    }
    
    if ( vm.count( "target" )==0 || vm["target"].as<string>().length()==0 )
    {
      cerr << "!! No target file specified" << endl << endl;
      printUsage( desc );
      return false;
    }
    else
    {
      	overlapper->targetFile = vm["target"].as<string>();
    }

    if ( vm.count("kmerLength")>0 )
    {
      	overlapper->kmerLength = vm["kmerLength"].as<int>();
    }
    
    if ( vm.count("bloomWidth")>0 )
    {
      	overlapper->setBloomWidth( vm["bloomWidth"].as<int>() );
    }
    
    if ( vm.count("numSegments")>0 )
    {
    	overlapper->numSegments = vm["numSegments"].as<int>();
    }
    
    if ( vm.count("topColumns")>0 )
    {
    	overlapper->topColumns = vm["topColumns"].as<int>();
    }   
    
    if ( vm.count("nproc")>0 )
    {
      	overlapper->numProcs = vm["nproc"].as<int>();
    }
    

  }
  catch( po::unknown_option u )
  {
    cerr << "!! " << u.what() << endl << endl;
    printUsage( desc );
	return false;
  }
  return true;
}

void PelusaApp::printUsage( const po::options_description& desc )
{
  cout << DESCRIPTION << endl;
  cout << endl << desc << endl;
}

int PelusaApp::run(  int argc, char **argv )
{
  	if (!processOptions( argc, argv ))
  	{
    	return 1;
  	}
  	overlapper->run();
	delete(overlapper); // TODO necessary
  	return 0;
}

int main( int argc, char **argv )
{
  try
  {
    PelusaApp app;
    return app.run( argc, argv );
  }
  catch( PelusaException e )
  {
    cerr << "Caught exception: " << e.getMessage() << endl;
    return 1;
  }
}
