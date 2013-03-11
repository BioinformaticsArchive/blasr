#include "data/hdf/HDFBasReader.h"
#include "data/hdf/HDFRegionTableReader.h"
#include "data/hdf/HDFRegionTableWriter.h"
#include "data/hdf/HDFBasWriter.h"
#include "data/hdf/PlatformId.h"
#include "utils/StringUtils.h"
#include "utils/RegionUtils.h"
#include "utils/FileOfFileNames.h"
#include "SMRTSequence.h"
#include <string>
#include <iostream>
#include <vector>

void PrintUsage() {
		cout << "usage: writeHDFSubset in out idx1 [idx2 idx3] [-pat p] [-fromto from to]..." << endl;
		cout << " idx is the index of a read (starting at 0)" << endl
				 << " pat p is a pattern to extract all reads with p in their title." << endl
         << " fromto is a range where from < to." << endl;
}
int main(int argc, char* argv[]) {
	
	string inFileName, outFileName;

	if (argc < 3) {
		PrintUsage();
		exit(1);
	}
	inFileName  = argv[1];
	outFileName = argv[2];

	vector<int> readIndices;
	int argi = 3;
	vector<string> patterns;
	vector<int> holeNumbers;
  string regionTableFileName = "";
  int from = 0, to = 0;
	while (argi < argc) {
		if (strlen(argv[argi]) > 0 and argv[argi][0] == '-'){ 
			if (strcmp(argv[argi], "-pat") == 0) {
				patterns.push_back(argv[++argi]);
			}
			else if (strcmp(argv[argi], "-holenumber") == 0) {
				holeNumbers.push_back(atoi(argv[++argi]));
			}
			else if (strcmp(argv[argi], "-regionTable") == 0) {
				regionTableFileName = argv[++argi];
			}
      else if (strcmp(argv[argi], "-fromto") == 0) {
        from = atoi(argv[++argi]);
        to   = atoi(argv[++argi]);
        if (from >= to) {
          cout <<"ERROR. From must be less than to." << endl;
          exit(1);
        }
      }
      else {
        cout <<"ERROR. Bad option " << argv[argi] << endl;
        PrintUsage();
        exit(1);
      }
		}
		else {
			readIndices.push_back(atoi(argv[argi]));
		}
		++argi;
	}
  int index;
  for (index = from; index < to; index++) {
    readIndices.push_back(index);
  }
	std::sort(readIndices.begin(), readIndices.end());
	T_HDFBasReader<SMRTSequence> reader;
  HDFRegionTableReader regionReader;

	HDFBasWriter writer;
  HDFRegionTableWriter regionWriter;
	reader.InitializeDefaultIncludedFields();
	writer.InitializeDefaultIncludedFields();
	writer.IncludeField("HoleNumber");
	writer.IncludeField("HoleXY");

  vector<string> inFiles;
  FileOfFileNames::StoreFileOrFileList(inFileName, inFiles);
  inFileName = inFiles[0];
	reader.Initialize(inFileName);
  RegionTable regionTable;
  if (regionTableFileName != "") {
    regionReader.Initialize(regionTableFileName);
  }
  else {
    regionReader.Initialize(inFileName);
  }
  regionReader.ReadTable(regionTable);
  
  string changeListID;
  reader.GetChangeListID(changeListID);
  
	if (reader.scanDataReader.GetPlatformId() == AstroPlatform) {
		writer.Initialize(outFileName, reader.GetMovieName(), reader.GetRunCode());
	}
	else {
		writer.Initialize(outFileName, reader.GetMovieName(), changeListID);
	}
  regionWriter.Initialize(writer.pulseDataGroup);
  

	int ri;
	int curReadIndex = 0;
	SMRTSequence seq;
	bool printSeq = false;
	ri = 0;
  if (readIndices.size() > 0) {
    reader.PrepareForRandomAccess();
    for (ri = 0; ri < readIndices.size(); ri++) {
      reader.GetReadAt(readIndices[ri], seq);
      writer.Write(seq);

      //
      // Write out region information for the read.
      //
      int low, high;
      FindRegionIndices(readIndices[ri], &regionTable, low, high);
      int regionIndex;
      for (regionIndex = low; regionIndex < high; regionIndex++) {
        regionWriter.Write(regionTable.table[regionIndex]);
      }
    }
    regionWriter.Finalize(regionTable.columnNames,
                          regionTable.regionTypes, 
                          regionTable.regionDescriptions, 
                          regionTable.regionSources
                          );
  }
  else if (patterns.size() > 0) {
    while (reader.GetNext(seq)) {
      printSeq = false;
      if (curReadIndex < readIndices.size() and ri == readIndices[curReadIndex]) {
        ++curReadIndex;
        printSeq = true;
      }
      int p;
      for (p = 0; p < patterns.size(); p++) {
        if (ExactPatternMatch(seq.title, patterns[p])) {
          printSeq = true;
          break;
        }
      }

      for (p = 0; p < holeNumbers.size(); p++) {
        if (seq.holeNumber == holeNumbers[p]) {
          printSeq = true;
          break;
        }
      }

      if (printSeq) {
        cout << "writing " << seq.title << endl;
        writer.Write(seq);
      }
      ++ri;
    }
  }
  
	writer.Flush();

}
