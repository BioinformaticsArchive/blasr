#include "../common/FASTAReader.h"
#include "../common/FASTASequence.h"
#include "../common/utils.h"

#include <string>
#include <stdlib.h>
#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char* argv[]) {
		if (argc < 4) {
			cout << "usage: splitContigs in.fa contiglength out" << endl;
			exit(1);
		}
		string inFileName, outFileName;
		inFileName = argv[1];
		int contigLength = atoi(argv[2]);		
		outFileName = argv[3];

		ofstream seqOut;
		CrucialOpen(outFileName, seqOut, std::ios::out);
		FASTAReader reader;
		reader.Init(inFileName);
		FASTASequence seq;
		DNALength curOffset;
		
		while(reader.GetNext(seq)) {
			FASTASequence subseq;
			int i;
			curOffset = 0;
			for (i =0 ; i < seq.length / contigLength + 1; i++ ) {
				subseq.seq = &seq.seq[curOffset];
				subseq.title = seq.title;
				if (curOffset + contigLength > seq.length) {
					subseq.length = seq.length - curOffset;
				}
				else {
					subseq.length = contigLength;
				}
				subseq.PrintSeq(seqOut);
				curOffset += contigLength;
			}
		}
		return 0;
}
