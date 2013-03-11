#ifndef PELUSAOVERLAPPER_H_
#define PELUSAOVERLAPPER_H_

#include <iostream>
#include <iomanip>
#include <ostream>
#include <fstream>
#include <string>
#include <map>
#include <set>
#include <pthread.h>
#include <vector>
#include <list>
#include <boost/tokenizer.hpp>
#include <boost/multi_array.hpp>
#include <cmath>

// TODO remove FILE At some point, and these libraries
#include <fstream>
#include <stdio.h>
#include <stdlib.h>


#include "PelusaException.h"
#include "BitArray.h" 	// LICENSE outside library LGPL. Should probably be using bitsetin STL. Ooops.
#include "GeneralHashFunctions.h" // LICENSE outside library CPL (sounds like LGPL)
#include "Fasta.h"
#include "FeatureEncoder.h"

#include <sstream>

#define PELUSA_BYTE_TO_SHORT 256

// typedef boost::multi_array<unsigned long long, 2> CUMULATIVE_SUM_TYPE;
typedef unsigned int uint;
typedef unsigned short ushort;
typedef boost::multi_array<unsigned short, 2> CUMULATIVE_SUM_TYPE;
typedef boost::multi_array<unsigned short, 2> BYTE_TO_SHORT_TYPE;

using namespace std;

class PelusaOverlapper
{
 	public:
    	PelusaOverlapper();
    	~PelusaOverlapper();
    	void run();

		// TODO private these should really be private 
		int debug;
  		string queryFile;
  		string targetFile;
    	int numProcs;
    	int kmerLength;
    	int bloomWidth;
    	int numSegments;
    	int topColumns;
   		int numFeatures;	
	
		const string toString();
		void setBloomWidth(int bloomWidth);
	
		// at least friendly
		int queryRecord(FastaRecord * record, map<string, int>* id2score);
	
		// private:
		void initializeBlooms();
		void populateBlooms();
		void queryBlooms();		
		void addRecordFeatures(FastaRecord * record);
		void findIdPairs(CUMULATIVE_SUM_TYPE cumulativeSums, 
			set< pair<uint, uint> > * idPairs);
		int pairAndedBitSum( pair<uint, uint>, vector<uint>* features);
		
		vector<bit_array_c*> blooms;
		BYTE_TO_SHORT_TYPE * byte2shorts;
		FeatureEncoder * encoder;
		map<pair<uint, uint>, vector<string>* > hashToIds;
		int collisionCount;
		
		
};

class PelusaWorker
{
  	public:
    	PelusaWorker(int id);
    	FILE* queryFh;
    	PelusaOverlapper * overlapper;
    	void run();
		void setQueryFh(FILE* queryFh);
		void setOverlapper(PelusaOverlapper * overlapper);
		void setQueryMutexPtr(pthread_mutex_t * queryMutexPtr);
		void setOutputMutexPtr(pthread_mutex_t * outputMutexPtr);
	
	private:
		int id;
		pthread_mutex_t * queryMutexPtr;
		pthread_mutex_t * outputMutexPtr;
		void outputScores(FastaRecord * record, map<string, int>* id2score);		 // TODO where should this go?

/*  private:
    SingleAlign a;
*/
};


#endif /*PELUSAOVERLAPPER_H_*/
