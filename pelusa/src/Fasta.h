#ifndef FASTA_H_
#define FASTA_H_
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <string>
#define MAX_LINE 100
#define MAX_SEQUENCE 1000000

using namespace std;

class FastaRecord {
	
	public:
		string sequence;
		string name;
		string description;
		FastaRecord();
		bool readTitleLine(FILE* fasta_file,   char* name, char* description);
		bool readRawSequence(FILE* fasta_file, char* raw_sequence);
		bool parseRecord(FILE* fasta_file);
};

#endif /*FASTA_H_*/
