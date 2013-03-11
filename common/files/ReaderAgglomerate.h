#ifndef ALIGNMENT_READER_AGGLOMERATE
#define ALIGNMENT_READER_AGGLOMERATE

#include <cstdlib>

#include "BaseSequenceIO.h"

#include "../FASTAReader.h"
#include "../FASTQReader.h"
#include "../CCSSequence.h"
#include "../SMRTSequence.h"
#include "../Enumerations.h"
#include "../data/hdf/HDFBasReader.h"
#include "../data/hdf/HDFCCSReader.h"
#include "../utils/StringUtils.h"

class ReaderAgglomerate : public BaseSequenceIO {
	FASTAReader fastaReader;
	FASTQReader fastqReader;
	int readQuality;
	int stride;
	int start;
	float subsample;
	bool useRegionTable;
	bool ignoreCCS;
 public:

	//
	// Create interfaces for reading hdf
	//
	T_HDFBasReader<SMRTSequence>  hdfBasReader;
	HDFCCSReader<CCSSequence>     hdfCcsReader;
	vector<SMRTSequence>          readBuffer;
	vector<CCSSequence>           ccsBuffer;
  string readGroupId;

  void SetToUpper() {
    fastaReader.SetToUpper();
  }
	void InitializeParameters() {
		start  = 0;
		stride = 1;
		subsample = 1.1;
		readQuality = 1;
		useRegionTable = true;
		ignoreCCS = true;
	}

	ReaderAgglomerate() {
		InitializeParameters();
	}
	
	ReaderAgglomerate(float _subsample) {
		this->InitializeParameters();
		subsample = _subsample;
	}
	 
	ReaderAgglomerate(int _stride) {
		this->InitializeParameters();
		stride = _stride;
	}
	
	ReaderAgglomerate(int _start, int _stride) {
		this->InitializeParameters();
		start  = _start;
		stride = _stride;
	}

  void GetMovieName(string &movieName) {
    if (fileType == Fasta || fileType == Fastq) {
      movieName = fileName;
    }
    else if (fileType == HDFPulse || fileType == HDFBase || fileType == HDFCCS) {
      movieName = hdfBasReader.GetMovieName();
    }
  }

	bool FileHasZMWInformation() {
		return (fileType == HDFPulse || fileType == HDFBase || fileType == HDFCCS);
	}

	void SkipReadQuality() {
		readQuality = 0;
	}
	
	void IgnoreCCS() {
		ignoreCCS = true;
	}
	
	void UseCCS() {
		ignoreCCS = false;
    hdfBasReader.SetReadBasesFromCCS();
	}

	int Initialize(string &pFileName) {
		if (DetermineFileTypeByExtension(pFileName, fileType)) {
			fileName = pFileName;
			return Initialize();
		}
    return false;
	}
	
	bool SetReadFileName(string &pFileName) {
		if (DetermineFileTypeByExtension(pFileName, fileType)) {
			fileName = pFileName;
			return true;
		}
		else {
			return false;
		}
	}
		
	int Initialize(FileType &pFileType, string &pFileName) {
		SetFiles(pFileType, pFileName);
		return Initialize();
	}

	bool HasRegionTable() {
		switch(fileType) {
		case Fasta:
			return false;
			break;
		case Fastq:
			return false;
			break;
		case HDFPulse:
		case HDFBase:
			return hdfBasReader.HasRegionTable();
			break;
		case HDFCCS:
			return hdfCcsReader.HasRegionTable();
			break;
		}
    return false;
	}
		
	int Initialize() {
		int init = 1;
		switch(fileType) {
		case Fasta:
			init = fastaReader.Init(fileName);
			break;
		case Fastq:
			init = fastqReader.Init(fileName);
			break;
		case HDFPulse:
		case HDFBase:
			//
			// Here one needs to test and see if the hdf file contains ccs.
			// If this is the case, then the file type is HDFCCS.
			if (hdfCcsReader.BasFileHasCCS(fileName) and !ignoreCCS) {
				
				fileType = HDFCCS;
				hdfCcsReader.InitializeDefaultIncludedFields();
				init = hdfCcsReader.Initialize(fileName);
				if (init == 0) return 0;
			}
			else {
				hdfBasReader.InitializeDefaultIncludedFields();
				init = hdfBasReader.Initialize(fileName);

				//
				// This code is added so that meaningful names are printed 
				// when running on simulated data that contains the coordinate
				// information.

				if (init == 0) return 0;
			}
			break;
		}
    readGroupId = "";
		if (init == 0 || (start > 0 && Advance(start) == 0) ){
			return 0;
		};

    string movieName;
    GetMovieName(movieName);
    MakeMD5(movieName, readGroupId, 10);
    
		return 1;
	}

	ReaderAgglomerate &operator=(ReaderAgglomerate &rhs) {
		fileType     = rhs.fileType;
		fileName = rhs.fileName;
		return *this;
	}
	
	bool Subsample(float rate) {
		bool retVal = true;
		while( (rand() % 100 + 1) > (rate * 100) and (retVal = Advance(1)));
		return retVal;
	}

	int GetNext(FASTASequence &seq) {
		int numRecords = 0;
		if (Subsample(subsample) == 0) {
			return 0;
		}
		switch(fileType) {
		case Fasta:
			numRecords = fastaReader.GetNext(seq);
			break;
		case Fastq:
			numRecords = fastqReader.GetNext(seq);
			break;
		case HDFPulse:
    case HDFBase:
			numRecords = hdfBasReader.GetNext(seq);
			break;
    case HDFCCS:
			cout << "ERROR! Reading CCS into a structure that cannot handle it." << endl;
			assert(0);
			break;
		}
    seq.CleanupOnFree();
		return numRecords;
	}

	int GetNext(FASTQSequence &seq) {
		int numRecords = 0;
		if (Subsample(subsample) == 0) {
			return 0;
		}
		switch(fileType) {
		case Fasta:
			numRecords = fastaReader.GetNext(seq);
			break;
		case Fastq:
			numRecords = fastqReader.GetNext(seq);
			break;
		case HDFPulse:
    case HDFBase:
			numRecords = hdfBasReader.GetNext(seq);
			break;
    case HDFCCS:
			cout << "ERROR! Reading CCS into a structure that cannot handle it." << endl;
			assert(0);
			break;
		}
		if (stride > 1)
			Advance(stride-1);
		return numRecords;
	}

	int GetNext(SMRTSequence &seq) {
		int numRecords = 0;

		if (Subsample(subsample) == 0) {
			return 0;
		}
		switch(fileType) {
		case Fasta:
			numRecords = fastaReader.GetNext(seq);
			break;
		case Fastq:
			numRecords = fastqReader.GetNext(seq);
			break;
		case HDFPulse:
    case HDFBase:
			numRecords = hdfBasReader.GetNext(seq);
			break;
    case HDFCCS:
      assert(ignoreCCS == false);
      assert(hdfBasReader.readBasesFromCCS == true);
      numRecords = hdfBasReader.GetNext(seq);
			break;
		}
		if (stride > 1)
			Advance(stride-1);
		return numRecords;
	}

	int GetNext(CCSSequence &seq) {
		int numRecords = 0;
		if (Subsample(subsample) == 0) {
			return 0;
		}

		switch(fileType) {
		case Fasta:
			// This just reads in the fasta sequence as if it were a ccs sequence
			numRecords = fastaReader.GetNext(seq);
			seq.subreadStart = 0;
			seq.subreadEnd   = 0;
			break;
		case Fastq:
			numRecords = fastqReader.GetNext(seq);
			seq.subreadStart = 0;
			seq.subreadEnd   = 0;
			break;
		case HDFPulse:
    case HDFBase:
			numRecords = hdfBasReader.GetNext(seq);
			break;
    case HDFCCS:
			numRecords = hdfCcsReader.GetNext(seq);
			break;
		}

		if (stride > 1)
			Advance(stride-1);
		return numRecords;
	}

	int Advance(int nSteps) {
    int i;
		switch(fileType) {
		case Fasta:
			return fastaReader.Advance(nSteps);
		case HDFPulse:
		case HDFBase:
			return hdfBasReader.Advance(nSteps);
		case HDFCCS:
			return hdfCcsReader.Advance(nSteps);
		case Fastq:
			return fastqReader.Advance(nSteps);
		}
    return false;
	}
	
	void Close() {
		switch(fileType) {

		case Fasta:
			fastaReader.Close();
			break;
		case HDFPulse:
		case HDFBase:
			hdfBasReader.Close();
			//			zmwReader.Close();
			break;
		case HDFCCS:
			hdfCcsReader.Close();
			break;
		}
	}
};

template<typename T_Sequence>
int ReadChunkByNReads(ReaderAgglomerate &reader, vector<T_Sequence> &reads, int maxNReads) {
	T_Sequence seq;
	int nReads = 0;
	while(nReads < maxNReads) {
		if (reader.GetNext(seq)) {
			reads.push_back(seq);
			++nReads;
		}
		else {
			break;
		}
	}
	return nReads;
}

template<typename T_Sequence>
int ReadChunkBySize (ReaderAgglomerate &reader, vector<T_Sequence> &reads, int maxMemorySize) {
	T_Sequence seq;
	int nReads = 0;
	int totalStorage = 0;
	while (totalStorage < maxMemorySize) {
		if (reader.GetNext(seq)) {
			reads.push_back(seq);
			totalStorage += seq.GetStorageSize();
			nReads++;
		}
		else {
			break;
		}
	}
	return nReads;
}
	 
		
			 
#endif
