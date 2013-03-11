#include "PelusaOverlapper.h"

// we have to use global variables until we can switch
// over to a more elegant thread solution
/*
 * RefSeq          *g_ref;
ReadClass       *g_read_a;
ifstream        *g_fin_a;
pthread_mutex_t g_mutex_fin=PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t g_mutex_fout=PTHREAD_MUTEX_INITIALIZER;
// JMS
ref_id_t        g_longReadUid = 0;
bit32_t         *g_n_aligned;
*/

PelusaOverlapper::PelusaOverlapper() 
	: 	debug(false),
		queryFile(""),
		targetFile(""),
		numProcs(1),
		kmerLength(8),
		bloomWidth(2000),
		numSegments(2),
		topColumns(10),
		encoder(NULL),
		collisionCount(0)
{
	
}

PelusaOverlapper::~PelusaOverlapper()
{
	// create our bloom filters
	for (int featureIdx = 0; featureIdx < numFeatures; featureIdx++)
	{
		delete(blooms[featureIdx]);
	}
	delete(encoder);
	for(map<pair<uint,uint>, vector<string>* >::const_iterator it = hashToIds.begin(); it != hashToIds.end(); ++it)
	{
		delete( (*it).second );
	}
	delete(byte2shorts);
}

void PelusaOverlapper::run()
{
	encoder = new FeatureEncoder(kmerLength);
	numFeatures = encoder->getNumFeatures(); // TODO should be in constructor?
	initializeBlooms();
	populateBlooms();
	queryBlooms();
}

void PelusaOverlapper::initializeBlooms()
{
	// create our bloom filters
	for (int featureIdx = 0; featureIdx < numFeatures; featureIdx++)
	{
		blooms.push_back( new bit_array_c( bloomWidth ) );
	}
	
	// populate the byte2ushortptr array which "spreads out" bits in a byte
	// into an eight byte unsigned long long int 
	// this data structure is useful for rapidly summing blooms in the query stage
	int bitsInByte = 8;
	byte2shorts = new BYTE_TO_SHORT_TYPE(boost::extents[256][8]);
	for (int shortIdx = 0; shortIdx < PELUSA_BYTE_TO_SHORT; shortIdx++)
	{	
		for (int bitIdx = 0; bitIdx < bitsInByte; bitIdx++)
		{	
			// TODO make sure we understand 7 - below
			// 7 -  below because we want the zeros place to be in the least significant bit.
			(*byte2shorts)[7 - bitIdx][bitIdx] = (shortIdx >> bitIdx) & 1;
		}
		//cout << byte2ushortptr[byteIdx] << endl;
	}
	
}

void PelusaOverlapper::populateBlooms()
{     
	FILE* file = fopen(targetFile.c_str(), "r");
	while(true)
	{
		FastaRecord * record = new FastaRecord();
		if (!record->parseRecord(file))
		{	
			delete(record);
			break;
		}
		if (record->sequence.length() < (uint) kmerLength)
		{
			cerr << "Sequence of length " << record->sequence.length() << " is too small for k of length " << kmerLength << endl;
			delete(record);
			continue;
		}
		addRecordFeatures(record);	
		
		delete(record);
	}
	fclose(file);
    
    if (debug) cerr << this->toString() << endl; 
    cerr << "Pelusa collision count " << this->collisionCount << endl;
} 

// overload for different hash functions, number etc.
void PelusaOverlapper::addRecordFeatures(FastaRecord * record)
{	
	// calculate hash values
	uint firstHashIdx  = RSHash(record->name) % bloomWidth;
	// uint secondHashIdx = JSHash(record->name) % bloomWidth;
	uint secondHashIdx = JSHash(record->name) % bloomWidth;
	pair<uint, uint> key = 	firstHashIdx < secondHashIdx 				?
		 				 	pair<uint, uint>(firstHashIdx, secondHashIdx) 	:
			 				pair<uint, uint>(secondHashIdx, firstHashIdx);
	
	// if (debug) cerr << "Hash id = " << key.first << " " << key.second << endl; 
	
	// add to our global map for mapping hash values back to record names

	if(hashToIds.count(key) == 0)
	{
		// PelusaOverlapper owns the vector // TODO delete vector in destructor
		vector<string>* strings = new vector<string>;
		hashToIds.insert( pair< pair<uint, uint>, vector<string>*>( key, strings )); 
	}
	else
	{
		collisionCount += 1;
	}
	hashToIds.find(key)->second->push_back( record->name );
	
	// set the bits at the hashed index positions in each feature's bloom filter
	// TODO replace with a stack array for efficiency?
	vector<uint>* features = new vector<uint>; 
	encoder->encode(record->sequence, features);
	if (debug) cerr << "hash=" << key.first << "," << key.second << endl;
	for (uint featureIdx=0; featureIdx < features->size(); featureIdx++)
	{
		// if (debug) cerr << "Feature " << featureIdx << "=" << (*features)[featureIdx] << endl;
		bit_array_c* featureBloom = blooms[ (uint)(*features)[featureIdx] ];
		featureBloom->SetBit(firstHashIdx);
		featureBloom->SetBit(secondHashIdx);	
	}
	
	delete(features);
}

void* threadStarter(void* arg)
{
	void * dummy;
	((PelusaWorker*)arg)->run();	
	return dummy;
}

void PelusaOverlapper::queryBlooms()
{     
	FILE* file = fopen(queryFile.c_str(), "r");
	pthread_mutex_t queryMutex=PTHREAD_MUTEX_INITIALIZER;
	pthread_mutex_t outputMutex=PTHREAD_MUTEX_INITIALIZER;
	
	vector<pthread_t> threads(numProcs);
	vector<PelusaWorker *> workers(numProcs);
	// create threads
	for (int workerIdx=0; workerIdx < numProcs; workerIdx++)
	{
		cerr << "Creating worker " << workerIdx << endl;
		PelusaWorker * worker = new PelusaWorker(workerIdx);
		workers[workerIdx] = worker;
		worker->setOverlapper(this);
		worker->setQueryFh(file);
		worker->setQueryMutexPtr(&queryMutex);
		worker->setOutputMutexPtr(&outputMutex);
		pthread_create(&threads[workerIdx], NULL, threadStarter, (void*)worker);
	}
	
	for (int workerIdx=0; workerIdx < numProcs; workerIdx++)
	{
		pthread_join(threads[workerIdx], NULL);
	}
	
	for (int workerIdx=0; workerIdx < numProcs; workerIdx++)
	{
		delete(workers[workerIdx]);
	}
	
	fclose(file);
}

int PelusaOverlapper::queryRecord(FastaRecord * record, map<string, int>* id2score)
{
	
	if (record->sequence.length() < (uint) kmerLength)
	{
		cerr << "Sequence of length " << record->sequence.length() << " is too small for k of length " << kmerLength << endl;
		return 0;
	}
	vector<uint>* features = new vector<uint>; 
	encoder->encode(record->sequence, features);
	
	// create a 2D array that holds the cumulative sums of all feature bloom filters
	// at any point in the read. this allows us to get the highest hash sums at 
	// between any two bounds in the read.
  	CUMULATIVE_SUM_TYPE cumulativeSum(boost::extents[features->size()+1][bloomWidth]);
  	for (int i=0; i<bloomWidth; i++) cumulativeSum[0][i] = 0;
    for (uint featureIdx = 0; featureIdx < features->size(); featureIdx++)
    {
    	uint bloomIdx = (uint)(*features)[featureIdx];
        for (int bitIdx = 0; bitIdx < bloomWidth; bitIdx++)
        {  			      
            cumulativeSum[featureIdx+1][bitIdx] = (*blooms[bloomIdx])[bitIdx] + cumulativeSum[featureIdx][bitIdx];
        }
    }
    
    if (debug)
    {
    	cerr << " Cumulative sum " << endl;
    	for (uint featureIdx =0 ; featureIdx < features->size()+1; featureIdx++)
    	{
   			for (uint cumIdx = 0; cumIdx < (uint) bloomWidth; cumIdx++)
        	{   
				cerr << (uint)cumulativeSum[featureIdx][cumIdx] << "\t";
        	}
        	cerr << endl;
    	}
    }
    
	//TODO consider keeping the feature vectors for efficiency
	
	// determine good ids pairs from our cumulativeSum array 
	set< pair<uint, uint> > * idPairs= new set< pair<uint, uint> >;
	findIdPairs(cumulativeSum, idPairs);
	
	int numPairs =  idPairs->size();
	if (debug) cerr << "Num pairs found:" << numPairs << endl;	
	for(set< pair<uint, uint> >::iterator it = idPairs->begin(); it != idPairs->end(); it++)
	{
		int score = pairAndedBitSum(*it, features);
		vector<string>* strings = hashToIds.find(*it)->second;
		for(uint stringIdx = 0; stringIdx < strings->size(); stringIdx++)
		{
			// TODO add bitsum
			// TODO maybe just ignore collisions??
			id2score->insert( pair<string, int>((*strings)[stringIdx], score)  );
		}
	}
	delete(features); // TODO put on stack??
	delete(idPairs);
	
	return numPairs;
}

int PelusaOverlapper::pairAndedBitSum( pair<uint, uint> uintPair, vector<uint>* features)
{
	int andBitSum = 0;
	for (uint featureIdx = 0; featureIdx < features->size(); featureIdx++)
	{
		uint bloomIdx = (uint)(*features)[featureIdx];
        int firstBit  = (*blooms[bloomIdx])[ uintPair.first  ];
        int secondBit = (*blooms[bloomIdx])[ uintPair.second ];
        // cerr << "f = " << firstBit << " s = " << secondBit << endl; 
		andBitSum += firstBit * secondBit;
	} 
	// cerr << "size " << features->size() <<  "abs " << andBitSum <<endl;
	return andBitSum;
}


void PelusaOverlapper::findIdPairs(
	CUMULATIVE_SUM_TYPE cumulativeSums, 
	set< pair<uint, uint> > * idPairs)
{
	// TODO allow for more than just first half prefix and second half suffix.
	int numBounds = 3;
	int numReadFeatures = cumulativeSums.shape()[0];
	int bounds[3][2] = { 	{0, numReadFeatures/2}, 
							{numReadFeatures/2+1, numReadFeatures-1}, 
							{0, numReadFeatures-1} };
	
	for (int boundIdx=0; boundIdx < numBounds; boundIdx++)
	{
		if (debug) cerr << "Analyzing bound " << boundIdx << endl;
		list< pair<unsigned short, uint> > * maximizer = 
			new list< pair<ushort, uint> >;
		for (uint cumIdx = 0; cumIdx < cumulativeSums.shape()[1]; cumIdx++)
		{
			// calculate the sum from the first bound to the second bound at this particular
			// cumulative idx
			ushort value =	cumulativeSums[ bounds[boundIdx][1] ][ cumIdx ] 
							      - cumulativeSums[ bounds[boundIdx][0] ][ cumIdx ];
			
			pair<ushort, uint> newPair = pair<ushort, uint>(value, cumIdx);
			maximizer->push_back(newPair);
		}
		// now sort the maximizer to find the chars with the largest sums
		maximizer->sort();
		maximizer->reverse();
		maximizer->resize(topColumns);	
		
		list< pair<ushort, uint> >::iterator firstPairIterator;
  		for ( firstPairIterator=maximizer->begin() ; firstPairIterator != maximizer->end(); firstPairIterator++ )
  		{
  			list< pair<ushort, uint> >::iterator secondPairIterator;
			if (debug) cerr << (uint)(*firstPairIterator).first << " " <<  
				(int)(*firstPairIterator).second<< endl;
			for ( secondPairIterator=firstPairIterator; secondPairIterator != maximizer->end(); secondPairIterator++ )
			{
				if (firstPairIterator == secondPairIterator) continue;

				int firstInt =  (*firstPairIterator).second;
				int secondInt = (*secondPairIterator).second;
				pair<int, int> newPair = firstInt < secondInt ?	
											pair<uint, uint>(firstInt, secondInt) 	:
											pair<uint, uint>(secondInt, firstInt);
				// if so, add them to our vector of idPairs											
				if (hashToIds.count( newPair ) != 0)
				{
					idPairs->insert( newPair );
				}
			}
		}
			
		delete(maximizer);
	}
}


void PelusaOverlapper::setBloomWidth(int bloomWidth)
{
	if ( floor(bloomWidth / 8) != bloomWidth / 8 )
    {
    	throw PelusaException("Bloom width should be divisible by 8!");
    }
	this->bloomWidth = bloomWidth;
}

const string PelusaOverlapper::toString()
{
	stringstream stream;
	stream << "queryFile=" << queryFile << endl;
	stream << "targetFile=" << targetFile << endl;
	stream << "numProcs=" << numProcs << endl;
	stream << "kmerLength=" << kmerLength << endl;		
	stream << "bloomWidth=" << bloomWidth << endl;	
	stream << "numSegments=" << numSegments << endl;
	stream << "topColumns=" << topColumns << endl;
	for (int bloomIdx = 0; bloomIdx < numFeatures; bloomIdx++)
	{	
		stream << bloomIdx << ": ";
		blooms[bloomIdx]->Dump(stream);
		stream << endl;
	}
	
	return stream.str();
	
}





PelusaWorker::PelusaWorker(int id)
	:	id(id)
{
	
}


void PelusaWorker::run()
{
	cerr << "Starting worker " << id << endl;
	while(true)
	{
		FastaRecord * record = new FastaRecord();
		
		// get the input
		pthread_mutex_lock(queryMutexPtr);
		bool wasParsed = record->parseRecord(queryFh)	;
		pthread_mutex_unlock(queryMutexPtr);
		
		// termination condition for all threads
		if (!wasParsed)
		{	
			delete(record);
			break;
		}
		
		// the heavy lifting
		map<string, int>* id2score = new map<string, int>;
		int numHits = overlapper->queryRecord(record, id2score);

		// output
		if (numHits > 0){
			pthread_mutex_lock(outputMutexPtr);
			outputScores(record, id2score);
			pthread_mutex_unlock(outputMutexPtr);
		}
		
		// cleanup
		delete(record);
		delete(id2score); // TODO put both of these on stack?
	}
	cerr << "Finishing worker " << id << endl;
}

void PelusaWorker::setOverlapper(PelusaOverlapper * overlapper)
{
	this->overlapper = overlapper;
}

void PelusaWorker::setQueryFh(FILE* queryFh)
{
	this->queryFh = queryFh;
}

void PelusaWorker::setQueryMutexPtr(pthread_mutex_t * queryMutexPtr)
{
	this->queryMutexPtr = queryMutexPtr;	
}

void PelusaWorker::setOutputMutexPtr(pthread_mutex_t * outputMutexPtr)
{
	this->outputMutexPtr = outputMutexPtr;	
}

// TODO check query parsing parellization. I think it's okay serially, but?

// TODO put in filehandle argument?
void PelusaWorker::outputScores(FastaRecord * record, map<string, int>* id2score)
{
	for(map<string, int>::const_iterator it = id2score->begin(); it != id2score->end(); ++it)
	{
		cout << record->name << "\t" << (*it).first << "\t" << (*it).second << "\n"; 
	}
}


/*

int PelusaOverlapper::queryRecord(FastaRecord * record, map<string, int>* id2score)
{
	
	if (record->sequence.length() < (uint) kmerLength)
	{
		cerr << "Sequence of length " << record->sequence.length() << " is too small for k of length " << kmerLength << endl;
		return 0;
	}
	vector<uint>* features = new vector<uint>; 
	encoder->encode(record->sequence, features);
	
	// create a 2D array that holds the cumulative sums of all feature bloom filters
	// at any point in the read. this allows us to get the highest hash sums at 
	// between any two bounds in the read.
  	CUMULATIVE_SUM_TYPE cumulativeSum(boost::extents[features->size()][bloomWidth/8]);
    for (uint featureIdx = 0; featureIdx  < features->size(); featureIdx++)
    {
        for (int byteIdx = 0; byteIdx < bloomWidth/8; byteIdx++)
        {
        		int bloomIdx = (int)(*features)[featureIdx];
                char byte = blooms[bloomIdx]->get_full_byte(byteIdx);
                // start with this rows byte ...
				cumulativeSum[featureIdx][byteIdx] = byte2ushortptr[(int)byte]; 
        		if (featureIdx != 0) // ... add in the previous row's byte, if one.
        			cumulativeSum[featureIdx][byteIdx] += cumulativeSum[featureIdx-1][byteIdx];
        }
    }
    
    if (debug)
    {
    	cerr << " Cumulative sum " << endl;
    	for (uint featureIdx =0 ; featureIdx < features->size(); featureIdx++)
    	{
   			for (int byteIdx = 0; byteIdx < bloomWidth/8; byteIdx++)
        	{   
				char* chars = (char*)&cumulativeSum[featureIdx][byteIdx];	
        		for (int i = 0; i < 8; i++) cerr << (int)chars[i] << "\t";
        	}
        	cerr << endl;
    	}
    }
    
	//TODO consider keeping the feature vectors for efficiency
	
	// determine good ids pairs from our cumulativeSum array 
	set< pair<int,int> > * idPairs= new set< pair<int,int> >;
	findIdPairs(cumulativeSum, idPairs);
	
	int numPairs =  idPairs->size();
	if (debug) cerr << "Num pairs found:" << numPairs << endl;
	// store away the score (the bitsum of the anded columns) for each pair START

		
	delete(features); // TODO put on stack??
	delete(idPairs);
	
	return numPairs;
}


void PelusaOverlapper::findIdPairs(
	CUMULATIVE_SUM_TYPE cumulativeSums, 
	set< pair<int,int> > * idPairs)
{
	// TODO allow for more than just first half prefix and second half suffix.
	int numBounds = 3;
	int numReadFeatures = cumulativeSums.shape()[0];
	int bounds[3][2] = { 	{0, numReadFeatures/2}, 
							{numReadFeatures/2+1, numReadFeatures-1}, 
							{0, numReadFeatures-1} };
	
	for (int boundIdx=0; boundIdx < numBounds; boundIdx++)
	{
		if (debug) cerr << "Analyzing bound " << boundIdx << endl;
		list< pair<ushort, int> > * maximizer = new list< pair<ushort, int> >;
		for (uint longIdx = 0; longIdx < cumulativeSums.shape()[1]; longIdx++)
		{
			// calculate the sum from the first bound to the second bound at this particular
			// long value
			// TODO this is a BUG, should not subtract as longs.
			unsigned long long value = 	  cumulativeSums[ bounds[boundIdx][1] ][ longIdx ] 
										- cumulativeSums[ bounds[boundIdx][0] ][ longIdx ];
			
			// then place each of the 8 chars in this long value into the maximize list
			ushort* chars = (ushort*)&value;
			for (int charIdx = 0; charIdx < 8; charIdx++) // TODO assumes ulonglong is 8bytes
			{
				// TODO the 8's are ugly!
				pair<ushort, int> newPair = pair<ushort, int>( chars[charIdx], longIdx*8 + charIdx);
				maximizer->push_back(newPair);
			}
			
		}
		// now sort the maximizer to find the chars with the largest sums
		maximizer->sort();
		maximizer->reverse();
		maximizer->resize(topColumns);	
		
		list< pair<ushort, int> >::iterator firstPairIterator;
  		for ( firstPairIterator=maximizer->begin() ; firstPairIterator != maximizer->end(); firstPairIterator++ )
  		{
  			list< pair<ushort, int> >::iterator secondPairIterator;
			if (debug) cerr << (uint)(*firstPairIterator).first << " " <<  
				(int)(*firstPairIterator).second<< endl;
			for ( secondPairIterator=firstPairIterator; secondPairIterator != maximizer->end(); secondPairIterator++ )
			{
				if (firstPairIterator == secondPairIterator) continue;

				int firstInt =  (*firstPairIterator).second;
				int secondInt = (*secondPairIterator).second;
				pair<int, int> newPair = firstInt < secondInt ?	
											pair<int, int>(firstInt, secondInt) 	:
											pair<int, int>(secondInt, firstInt);
				// if so, add them to our vector of idPairs											
				if (hashToIds.count( newPair ) != 0)
				{
					idPairs->insert( newPair );
				}
			}
		}
			
		delete(maximizer);
	}
}*/
