#define __FAST_MATH__

#include "data/hdf/HDFCmpFile.h"
#include "data/hdf/HDFBasReader.h"
#include "data/hdf/HDFPlsReader.h"
#include "data/hdf/HDFCCSReader.h"
#include "data/hdf/PlatformId.h"
#include "datastructures/alignment/CmpFile.h"
#include "datastructures/alignment/CmpAlignment.h"
#include "datastructures/alignment/ByteAlignment.h"
#include "datastructures/reads/BaseFile.h"
#include "datastructures/reads/PulseFile.h"
#include "datastructures/reads/ReadType.h"
#include "utils/FileOfFileNames.h"
#include "utils/TimeUtils.h"
#include "CommandLineParser.h"
#include <map>
#include <set>
#include <string>
#include <algorithm>
#include <assert.h>
#include <numeric>

typedef map<int,int> HoleToStartMap;
typedef map<string, HoleToStartMap*> MovieToHoleMap;
typedef map<string, int> MovieNameToArrayIndex;
typedef map<string, bool> MetricOptionsMap;
using namespace std;

typedef map<string, vector<string> > RequirementMap;

char VERSION[] = "v1.1.0";
char PERFORCE_VERSION_STRING[] = "$Change: 107666 $";



void CapQualityValue(QualityValueVector<QualityValue> &vect, DNALength length, unsigned char maxQualityValue=100) {
  unsigned int i;
  if (vect.data == NULL) {
    return;
  }
  for (i = 0; i < length; i++) {
    vect.data[i] = min(vect.data[i], maxQualityValue);
  }
}


void CapQualityValues(SMRTSequence &seq, unsigned char maxQualityValue = 100) {
  CapQualityValue(seq.qual, seq.length, maxQualityValue);
  CapQualityValue(seq.deletionQV, seq.length, maxQualityValue);
  CapQualityValue(seq.preBaseDeletionQV, seq.length, maxQualityValue);
  CapQualityValue(seq.insertionQV, seq.length, maxQualityValue);
  CapQualityValue(seq.substitutionQV, seq.length, maxQualityValue);
  CapQualityValue(seq.mergeQV, seq.length, maxQualityValue);
}


int CheckCmpFileFormat(  CmpFile &cmpFile) {
  if (cmpFile.readType != ReadType::Standard) {
    cout << "ERROR! Reading pulse information into a cmp.h5 file generated from circular " << endl
         << "consensus called sequences is not supported." << endl;
    exit(0);
  }
  return 1;
}

void BuildRequirementMap(RequirementMap &fieldRequirements) {
	fieldRequirements["StartTimeOffset"].push_back("StartFrame");
	fieldRequirements["StartTimeOffset"].push_back("NumEvent");
	fieldRequirements["StartFrame"].push_back("PreBaseFrames");
	fieldRequirements["StartFrame"].push_back("WidthInFrames");
	fieldRequirements["PulseWidth"].push_back("WidthInFrames");
	fieldRequirements["pkmid"].push_back("MidSignal");
	fieldRequirements["pkmid"].push_back("NumEvent");
	fieldRequirements["IPD"].push_back("StartFrame");
	fieldRequirements["IPD"].push_back("NumEvent");
	fieldRequirements["IPD"].push_back("PreBaseFrames");
	fieldRequirements["IPD"].push_back("WidthInFrames");
	fieldRequirements["IPD"].push_back("NumEvent");
	fieldRequirements["Light"].push_back("MeanSignal");
	fieldRequirements["Light"].push_back("NumEvent");
	fieldRequirements["Light"].push_back("WidthInFrames");
	fieldRequirements["Light"].push_back("NumEvent");
}

void ExclusivelyAdd(const char *value, vector<string> &vect) {
	if (find(vect.begin(), vect.end(), value) == vect.end()) {
		vect.push_back(value);
	}
}

bool AnyFieldRequiresFrameRate(vector<string> &fields) {
  int i;
  for (i = 0; i < fields.size(); i++ ) {
    if (fields[i] == "PulseWidth" or
        fields[i] == "IPD" or
        fields[i] == "Light" or
        fields[i] == "StartTimeOffset" or
        fields[i] == "StartFrame" or
        fields[i] == "PulseWidth" or 
        fields[i] == "PreBaseFrames" or 
        fields[i] == "WidthInFrames") {
      return true;
    }
  }
  return false;
}

template<typename T>
void Free(T* &buf) {
  if (buf != NULL){ 
    delete[] buf;
  }
  buf = NULL;
}

void StoreDatasetFieldsFromPulseFields(MetricOptionsMap &fieldSet,
																			 RequirementMap &fieldRequirements, 
																			 vector<string> &datasetFields) {
	int f;
	int d;
	MetricOptionsMap::iterator optionsIt;
	for (optionsIt = fieldSet.begin(); optionsIt != fieldSet.end(); ++optionsIt) {
		if (optionsIt->second == true) {
			if (fieldRequirements.find(optionsIt->first) == fieldRequirements.end()) {
				ExclusivelyAdd(optionsIt->first.c_str(), datasetFields);
			}
			else {
				for (d = 0; d < fieldRequirements[optionsIt->first].size(); d++) {
					ExclusivelyAdd(fieldRequirements[optionsIt->first][d].c_str(), datasetFields );
				}
			}
		}
	}
}

void ParseMetricsList(string metricListString, MetricOptionsMap &metricOptions) {
	vector<string> metrics;
	Tokenize(metricListString, ",", metrics);
	int m;
	for  (m = 0; m < metrics.size(); m++) {
		if (metricOptions.find(metrics[m]) != metricOptions.end()) {
			metricOptions[metrics[m]] = true;
		}
		else {
			cout << "ERROR! Metric " << metrics[m] << " is not supported." << endl;
			exit(1);
		}
	}
}

void SetDefaultMetricOptions(MetricOptionsMap &metricOptions) {
	metricOptions["QualityValue"]    = true;
	metricOptions["ClassifierQV"]    = true;
	metricOptions["StartFrame"]       = true;
	metricOptions["PulseWidth"]      = true;
	metricOptions["WidthInFrames"]   = true;
	metricOptions["pkmid"]           = true;
	metricOptions["IPD"]             = true;
}

void CreateMetricOptions(MetricOptionsMap &metricOptions) {
	metricOptions["QualityValue"]    = false;
	metricOptions["ClassifierQV"]    = false;
	metricOptions["StartFrame"]       = false;
	metricOptions["PulseWidth"]      = false;
	metricOptions["WidthInFrames"]   = false;
	metricOptions["pkmid"]           = false;
	metricOptions["IPD"]             = false;
	metricOptions["PreBaseFrames"]   = false;
	metricOptions["InsertionQV"]     = false;
	metricOptions["SubstitutionQV"]  = false;
	metricOptions["SubstitutionTag"] = false;
	metricOptions["DeletionQV"]      = false;
	metricOptions["DeletionTag"]     = false;
	metricOptions["StartTimeOffset"] = false;
	metricOptions["WhenStarted"]     = false;
  metricOptions["PulseIndex"]      = false;
  metricOptions["MergeQV"]      = false;
}

void PrintUsage() {
	cout << "  loadPulses - Load pulse information and quality values into a Compare file" << endl;
	cout << "usage: loadPulses movieFile cmpFile [-metrics m1,m2,...] [-useccs] [-byread]" << endl;
	cout << " movieFile may be a movie file or a fofn of movie file names." << endl;
	cout << " metrics m1,m2,... is a comma-separated list (without spaces) of metrics " << endl
       << " to print to the pulse file." << endl;
	cout << " Valid metrics are: " << endl
			 << "    QualityValue, ClassifierQV, MergeQV," << endl
       << "    StartFrame, PulseWidth, pkmid, IPD," << endl
			 << "    WhenStarted, StartTimeOffset, PreBaseFrames," << endl
       << "    InsertionQV, DeletionQV, DeletionTag, SubstitutionQV" << endl
       << "    SubstitutionTag, MergeQV, PulseIndex" << endl;
	cout << "  By default, QualityValue, ClassifierQV, StartFrame, PulseWidth," << endl
			 << "  WidthInFrames, pkmid, and IPD are added" << endl;
  cout << "  -useccs  This option is for older cmp.h5 files that do not have the read type " << endl
       << "    stored.  Newer cmp.h5 files have a read type that indicates the cmp.h5 file " << endl
       << "    has alignments generated from de novo ccs sequences.  Using this flag assuems"<<endl
       << "    ALL alignments in the cmp.h5 file are from ccs sequences, and loads the "<< endl
       << "    quality values from ccs instead of the raw sequence. "<<endl 
       << "  The only metrics that are allowed for de novo ccs sequences are QualityValue, " << endl
       << "  InsertionQV, DeletionQV, and SubstitutionQV" << endl;
  cout <<"   -byread  Reads pulse/base fields by read, rather than reading an entire " << endl
       << "    movie first.  This uses considerably less memory than the defualt mode" << endl
       << "    but is slow." << endl;
	cout << "  Using hdf version " << H5_VERS_MAJOR << "." << H5_VERS_MINOR << "." << H5_VERS_RELEASE << endl;
}

int main(int argc, char* argv[]) {

	string cmpFileName, movieFileName;

	int argi = 3;
	int numMetrics = 8;
	map<string,bool> metricOptions;
	int maxElements = 0;
	//
	// Default is all options are true
	//
	CreateMetricOptions(metricOptions);
	string metricList = "";
  bool useCcs = false;
  bool byRead = false;
  bool failOnMissingData = false;
  CommandLineParser clp;
  bool printVersion = false;

  clp.RegisterStringOption("basFileName", &movieFileName, "The input {bas,pls}.h5 or input.fofn.", true);
  clp.RegisterStringOption("cmpFileName", &cmpFileName, "The cmp.h5 file to load pulse information into.", true);
  clp.RegisterPreviousFlagsAsHidden();
  clp.RegisterStringOption("metrics", &metricList, "The a string delimited list of metrics (with no spaces).The "
                           "valid options are:  QualityValue, ClassifierQV, MergeQV, StartFrame,"
                           "PulseWidth, pkmid, IPD, and Light.");
  clp.RegisterFlagOption("useccs", &useCcs, "Load pulse information for CCS sequences and not raw bases.");
  clp.RegisterFlagOption("byread", &byRead, "Load pulse information by read rather than buffering an entire pls.h5 file.  "
                         "This option will soon be deprecated and on by default.");
  clp.RegisterIntOption("maxElements", &maxElements, "Set a limit on the size of pls/bas file to buffer in.", CommandLineParser::PositiveInteger);
  clp.RegisterFlagOption("failOnMissingData", &failOnMissingData, "Exit if any data fields are missing from the bas.h5 or pls.h5 input that are required to load a metric. Defualt is a warning.");

  clp.SetProgramSummary("Load pulse information such as inter pulse distance, or quality information into the cmp.h5 file."
                        "This allows one to analyze kinetic and quality information by alignment column.");
  clp.ParseCommandLine(argc, argv);

  if (printVersion) {
    cout << VERSION << endl;
    exit(1);
  }

	if (metricList == "") {
		SetDefaultMetricOptions(metricOptions);
	}
  else {
    ParseMetricsList(metricList, metricOptions);
  }


	// 
	// Always read in basecalls since they are used to check the sanity
	// of the alignment indices.
	//
	metricOptions["Basecall"] = true;
	
	//
	// Translate from the metrics to be loaded to the ones that are
	// required to compute them.
	//
	
	vector<string> datasetFields;
	RequirementMap fieldRequirements;
	BuildRequirementMap(fieldRequirements);
	StoreDatasetFieldsFromPulseFields(metricOptions, fieldRequirements, datasetFields);

	
	vector<string> movieFileNames;
	vector<string> fofnMovieNames;

	FileOfFileNames::StoreFileOrFileList(movieFileName, movieFileNames);

	HDFBasReader hdfBasReader;
	HDFPlsReader hdfPlsReader;
  HDFCCSReader<SMRTSequence> hdfCcsReader;

	vector<string> baseFileFields, pulseFileFields;
	int fieldIndex;

	bool useBaseFile = false, usePulseFile = false;

	for (fieldIndex = 0; fieldIndex < datasetFields.size(); fieldIndex++) {
		if (hdfBasReader.ContainsField(datasetFields[fieldIndex])) {
			useBaseFile = true;
			baseFileFields.push_back(datasetFields[fieldIndex]);
		}
	}

	if (maxElements != 0) {
		hdfBasReader.maxAllocNElements = maxElements;
		hdfPlsReader.maxAllocNElements = maxElements;
	}

	//
	// For now, all runs will attempt to use information from a .bas
	// file, since it's assumed that if one has alignments, one has a
	// .bas file.
	//
	useBaseFile = true;
	//
	// Add some default fields.
	//
	hdfBasReader.IncludeField("Basecall");
	hdfBasReader.IncludeField("PulseIndex");

  hdfBasReader.InitializeFields(baseFileFields);

  
	for (fieldIndex = 0; fieldIndex < datasetFields.size(); fieldIndex++) {
		if (hdfPlsReader.ContainsField(datasetFields[fieldIndex])) {
			usePulseFile = true;
			pulseFileFields.push_back(datasetFields[fieldIndex]);
		}
	}
	if (usePulseFile) {
		hdfPlsReader.InitializeFields(pulseFileFields);
	}
	hdfPlsReader.IncludeField("NumEvent");

		
	
	int nMovies = movieFileNames.size();
	int movieIndex;
	MovieNameToArrayIndex movieNameMap;
	//
	// Initialize movies. This accomplishes two tasks.  First, all movie
	// files are opened and initialized, so that if there are data
	// fields missing the program will exit now rather than in the
	// middle of loading pulses.  
	// Next, a list of movie names is created in fofnMovieNames.  The
	// cmp file does not necessarily index movies in the order of the
	// fofn, and so when loading pulses from a movie indexed by a cmp
	// file, one needs to look up the file name of the movie.  This is
	// done by scanning the fofnMovieNames list in order until the movie
	// is found. 

	for (movieIndex = 0; movieIndex < nMovies; movieIndex++) {

    if (!hdfBasReader.Initialize(movieFileNames[movieIndex])) {
      cout << "ERROR, could not initialize HDF file "
           << movieFileNames[movieIndex] << " for reading bases." << endl;
      exit(1);
    }
    else {
      fofnMovieNames.push_back(hdfBasReader.GetMovieName());
      movieNameMap[hdfBasReader.GetMovieName()] = movieIndex;
      hdfBasReader.Close();
    }

		// 
		// The pulse file is optional.  
		//
		if (usePulseFile) {
			if (hdfPlsReader.Initialize(movieFileNames[movieIndex]) == 0) {
				usePulseFile = false;
			}
		}		
	}

	CmpFile cmpFile;
	
	/*
	 * These readers pull information from the same pls file.
	 */
	HDFCmpFile<CmpAlignment> cmpReader;

	if (cmpReader.Initialize(cmpFileName, H5F_ACC_RDWR) == 0) {
		cout << "ERROR, could not open the cmp file." << endl;
		exit(0);
	}
	
  cmpReader.Read(cmpFile);

  string commandLine;
  clp.CommandLineToString(argc, argv, commandLine);
  string versionStr(VERSION);
  AppendPerforceChangelist(PERFORCE_VERSION_STRING, versionStr);
  cmpReader.fileLogGroup.AddEntry(commandLine, "Loading pulse metrics", "loadPulses", GetTimestamp(), versionStr);


	//
	// Group alignment indices by movie so that they may be processed one movie at a time
	// later on.  The movie indices set keeps track of all indices
	// listed in alignment files.  This keeps a reference to all
	// alignments in memory at once.   At the time of writing this, most
	// projects will have at most a few million alignments, and so the
	// size of this structure is modest.
	//

	UInt alignmentIndex;
	map<int, vector<int> > movieIndexSets;

	for (alignmentIndex = 0; alignmentIndex < cmpFile.alnInfo.alignments.size(); alignmentIndex++) {
		movieIndexSets[cmpFile.alnInfo.alignments[alignmentIndex].GetMovieId()].push_back(alignmentIndex);
	}

	vector<float>  computedPulseField;
	string   alignedSequence;
	string   readSequence;
	vector<unsigned char> byteAlignment;
	int m;
	vector<int> baseToAlignmentMap;

	//
	// Load pulses from movies in order they appear in the input fofn.
	//
  int fofnMovieIndex;
  for (fofnMovieIndex = 0; fofnMovieIndex < fofnMovieNames.size(); fofnMovieIndex++) {
    
    if (cmpFile.readType == ReadType::CCS or useCcs) {
      hdfBasReader.SetReadBasesFromCCS();
      hdfCcsReader.Initialize(movieFileNames[fofnMovieIndex]);
    }
    hdfBasReader.Initialize(movieFileNames[fofnMovieIndex]);
		BaseFile  baseFile;
		PulseFile pulseFile;

    if (byRead == false) {
      //
      // Read the entire bas file at once, and then extract values
      // from memory.  This can be faster depending on the chunk
      // size and size of the movie.
      //
      hdfBasReader.ReadBaseFile(baseFile);
      hdfBasReader.Close();
    }
    else {
      //
      // Reads are scanned one by instead of caching all.  It is
      // still necessary to read in some of the datasets entirely,
      // in particular the start positions and hole numbers.
      //


      // This is repeated below for a pulse file.  Since the pulse
      // and base files are separate objects, the scan data is
      // read into each separately.  Somehow later the information
      // should be merged into just one.
      if (hdfBasReader.scanDataReader.fileHasScanData) {
        hdfBasReader.scanDataReader.Read(baseFile.scanData);
      }
      baseFile.readStartPositions.resize(hdfBasReader.nReads+1);
      baseFile.readStartPositions[0] = 0;
      hdfBasReader.GetAllReadLengths(baseFile.readLengths);
      int i;
      assert(baseFile.readLengths.size() + 1 == baseFile.readStartPositions.size());
      for (i = 1; i < hdfBasReader.nReads + 1; i++ ) {
        baseFile.readStartPositions[i] = baseFile.readLengths[i-1] + baseFile.readStartPositions[i-1];
      }
      
      //
      // Although the whole bas file isn't being read in, it is
      // necessary to read in which hole numbers are contained in this
      // bas file since it is possible that the alignment for a
      // particular hole number may be in a different input bas.h5
      // file even if it is the same movie. 
      //
      hdfBasReader.GetAllHoleNumbers(baseFile.holeNumbers);
    }
    set<uint32_t> moviePartHoleNumbers;
    copy(baseFile.holeNumbers.begin(), baseFile.holeNumbers.end(), inserter(moviePartHoleNumbers, moviePartHoleNumbers.begin()));

		
		if (usePulseFile) {
			hdfPlsReader.Initialize(movieFileNames[fofnMovieIndex]);
			hdfPlsReader.IncludeField("NumEvent");
      hdfPlsReader.IncludeField("StartFrame");
      if (byRead == false) { 
        hdfPlsReader.ReadPulseFile(pulseFile);
        hdfPlsReader.Close();
      }
      else {
        if (usePulseFile) {
          pulseFile.pulseStartPositions.resize(hdfBasReader.nReads+1);
          pulseFile.pulseStartPositions[0] = 0;
          hdfPlsReader.GetAllNumEvent(pulseFile.numEvent);
          int i;
          for (i = 1; i < hdfBasReader.nReads + 1; i++ ) {
            pulseFile.pulseStartPositions[i] = pulseFile.numEvent[i-1] + pulseFile.pulseStartPositions[i-1];
          }
          if (hdfPlsReader.scanDataReader.fileHasScanData) {
            hdfPlsReader.scanDataReader.Read(pulseFile.scanData);
          }
        }
      }
		}

    string cmpFileMovieName;

    for (m = 0; m < cmpFile.movieInfo.name.size(); m++) {
      //
      // First find the file name for the movie 'm'
      //
      cmpFileMovieName = cmpFile.movieInfo.name[m];
      int fofnMovieIndex;
      
      if (baseFile.GetMovieName() == cmpFileMovieName) {
				break;
			}
		}

    //
    // If the movie specified in the input.fofn is not found in the
    // cmp file, that indicates something bad is happeing.  Either the
    // input.fofn was not used to generate the cmp.h5 file, or no
    // alignments were found between the input bas.h5 and the
    // reference.  That shouldn't happen.
    // 
		if (m == cmpFile.movieInfo.name.size()) {
			cout << "WARNING: The movie indexed in the compare file " << cmpFileMovieName << " is not listed in the file " << movieFileName << endl;
			continue;
		}
		
		//
		// Open the movie and load its pulses into memory.
		//
		movieIndex = cmpFile.movieInfo.id[m];
		int movieAlignmentIndex;
		float NaN = 0.0/0.0;
    
    UChar missingQualityValue = 255;
    HalfWord missingFrameRateValue    = USHRT_MAX;
    unsigned int missingPulseIndex = UINT_MAX;
    //
    // Since usePulseFile is set when the input file is a pulseFile,
    // and ReadType::CCS becomes the read type when the alignments are
    // ccs, when pulse files are specified for de novo ccs alignments,
    // they will be opened as pulse files.  Since the de novo ccs
    // sequences do not have pulse file information, the auto-reading
    // of pulse files needs to be disabled.  Do that here.
    //
    if (cmpFile.readType == ReadType::CCS or useCcs) {
      usePulseFile = false;
    }


		//
		// Now check the sanity of metric options.
		//

		map<string,bool>::iterator metricIt;
		for (metricIt = metricOptions.begin(); metricIt != metricOptions.end(); ++metricIt) {
			if (metricIt->second == false) {
				continue;
			}
			bool metricMayBeComputed = true;
      if (cmpFile.readType == ReadType::CCS and
          metricIt->first != "QualityValue"  and
          metricIt->first != "DeletionQV" and
          metricIt->first != "SubstitutionQV" and
          metricIt->first != "InsertionQV" and
          metricIt->first != "DeletionTag" and
          metricIt->first != "SubstitutionTag" and
          metricIt->first != "Basecall") {
        cout << "ERROR! The metric " << metricIt->first << " cannot be loaded into de novo ccs alignemnts." << endl;
        //        exit(0);
        metricMayBeComputed = false;
      }
      
			if (metricIt->first == "IPD") {
				//
				// The field requirements for IPD are special. 
				//
				if ((useBaseFile and !hdfBasReader.FieldIsIncluded("PreBaseFrames")) or
						(usePulseFile and (!hdfPlsReader.FieldIsIncluded("StartFrame") and
															 !hdfPlsReader.FieldIsIncluded("WidthInFrames")))) {
					metricMayBeComputed = false;
				}
			}
			else {
				if (fieldRequirements.find(metricIt->first) != fieldRequirements.end()) {
					//
					// There are requirements for this field. Make sure all are
					// present before trying to compute this field.
					//
					int requirementIndex;
					for (requirementIndex = 0; requirementIndex < fieldRequirements[metricIt->first].size(); ++requirementIndex) {
						string requirement;
						requirement = fieldRequirements[metricIt->first][requirementIndex];
				
						if (((useBaseFile == false or ((hdfBasReader.includedFields.find(requirement) == hdfBasReader.includedFields.end() or
 																						hdfBasReader.includedFields[requirement] == false))) and
								 ((usePulseFile == false or (hdfPlsReader.includedFields.find(requirement) == hdfPlsReader.includedFields.end() or
																						 hdfPlsReader.includedFields[requirement] == false))))) {
							metricMayBeComputed = false;
						}
					}
				}
				else {
					//
					// There are no requirements for this field, so it must exist as
					// a datset in either the bas or pls file.
					//
					if ((useBaseFile  == false or ((hdfBasReader.includedFields.find(metricIt->first) == hdfBasReader.includedFields.end() or
																					hdfBasReader.includedFields[metricIt->first] == false))) and
							(usePulseFile == false or (((hdfPlsReader.includedFields.find(metricIt->first) == hdfPlsReader.includedFields.end() or
																					 hdfPlsReader.includedFields[metricIt->first] == false))))) {
						metricMayBeComputed = false;
					}
				}
			}
			if (metricMayBeComputed == false) {
        if (failOnMissingData) {
          cout << "ERROR";
        }
        else {
          cout << "WARNING";
        }
        cout << ": There is insufficient data to compute metric: " << metricIt->first << " in the file " << movieFileNames[fofnMovieIndex] << " ";
        cout << " It will be ignored." << endl;
        if (failOnMissingData) {
          exit(1);
        }
				metricOptions[metricIt->first] = false;
			}
		}

		
		UInt i;
		//
		// This is currently used as a sentinal for showing that an array
		// element does not have a value stored for it, as in deleted
		// bases. 
		//

		
		vector<int> pulseIndexArray;
		vector<unsigned int> statTime;

		if (metricOptions["WhenStarted"]) {
			string whenStarted;
			if (hdfPlsReader.scanDataReader.useWhenStarted == false) {
				cout << "ERROR! Attempting to read WhenStarted from " 
						 << movieFileNames[fofnMovieIndex]
						 << " but the attriubte does not exist." << endl;
				exit(1);
			}
			hdfPlsReader.scanDataReader.ReadWhenStarted(whenStarted);
			
			if (!cmpReader.movieInfoGroup.whenStartedArray.IsInitialized()) {
				cmpReader.movieInfoGroup.whenStartedArray.Initialize(cmpReader.movieInfoGroup.movieInfoGroup, "WhenStarted");
			}

			cmpReader.movieInfoGroup.whenStartedArray.Write(&whenStarted, 1);
		}

    if (AnyFieldRequiresFrameRate(datasetFields)) {
      if (useBaseFile) {
        cmpReader.movieInfoGroup.StoreFrameRate(m, baseFile.GetFrameRate());
      }
      else if (usePulseFile) {
        cmpReader.movieInfoGroup.StoreFrameRate(m, pulseFile.GetFrameRate());
      }
    }
				
		//
		// An index set is a set of indices into the alignment array that
		// are of reads generated by this movie.  Load pulses for all
		// alignments generated for this movie.
		//

		//
		// Movie index sets should be sorted by alignment index. Build a lookup table for this.
		//
		
		std::vector<std::pair<int,int> > toFrom;
		for (movieAlignmentIndex = 0; movieAlignmentIndex < movieIndexSets[movieIndex].size(); movieAlignmentIndex++) {
			alignmentIndex = movieIndexSets[movieIndex][movieAlignmentIndex];
			toFrom.push_back(std::pair<int,int>(cmpFile.alnInfo.alignments[alignmentIndex].GetAlignmentId(), movieAlignmentIndex));
		}
		// orders by first by default.
		std::sort(toFrom.begin(), toFrom.end());

    //
    // Load metrics for alignments from movie 'movieIndex'.
    //
    cout << "loading " <<  movieIndexSets[movieIndex].size() << " alignments for movie " << movieIndex << endl;
		for (movieAlignmentIndex = 0; movieAlignmentIndex < movieIndexSets[movieIndex].size(); movieAlignmentIndex++) {
			alignmentIndex = movieIndexSets[movieIndex][toFrom[movieAlignmentIndex].second];


			//
			// Alignments are groupsd by ref group id then movie id.
			//
			int refGroupId  = cmpFile.alnInfo.alignments[alignmentIndex].GetRefGroupId();
			int movieId     = cmpFile.alnInfo.alignments[alignmentIndex].GetMovieId();
      UInt holeNumber = cmpFile.alnInfo.alignments[alignmentIndex].GetHoleNumber();

      //
      // Since the movie may be split into multiple parts, look to see
      // if this hole number is one of the ones covered by this
      // set. If it is not, just continue. It will be loaded on
      // another pass through a different movie part.
      //
      if (moviePartHoleNumbers.find(holeNumber) == moviePartHoleNumbers.end()) {
        continue;
      }

			//
			// Now locate where this movie is stored.
			//

			if (cmpReader.refGroupIdToArrayIndex.find(refGroupId) == cmpReader.refGroupIdToArrayIndex.end()) {
				cout << "ERROR!  An alignment " << alignmentIndex << " is specified with reference group " << endl
						 << refGroupId << " that is not found as an alignment group." << endl;
				exit(1);
			}
			int refGroupIndex = cmpReader.refGroupIdToArrayIndex[refGroupId];
			
			//
			// Now find the group containing the alignment for this movie.
			//
			if (cmpReader.refAlignGroups[refGroupIndex]->movieIdToIndex.find(movieId) ==
					cmpReader.refAlignGroups[refGroupIndex]->movieIdToIndex.end()) {
				cout << "ERROR!  An alignment " << alignmentIndex << " is specified with movie index " << endl
						 << movieId << " that is not found in the alignment group " << refGroupIndex << endl;
				exit(1);
			}

			int readGroupIndex = cmpReader.refAlignGroups[refGroupIndex]->movieIdToIndex[movieId];
      
			//
			// First do sanity check on the read to make sure the pules and the bases match.
			//

			//
			// Look to see if the output HDF arrays need to be created.
			//
			UInt offsetBegin, offsetEnd;
		
			offsetBegin = cmpFile.alnInfo.alignments[alignmentIndex].GetOffsetBegin();
			offsetEnd   = cmpFile.alnInfo.alignments[alignmentIndex].GetOffsetEnd();
		
			int alignedSequenceLength = offsetEnd - offsetBegin;
			if (alignedSequenceLength >= 0) {
				alignedSequence.resize(alignedSequenceLength);
				byteAlignment.resize(alignedSequenceLength);
			}
	
			//
			// Read the alignment string.  All alignments 
			//
			cmpReader.refAlignGroups[refGroupIndex]->readGroups[readGroupIndex]->alignmentArray.Read(offsetBegin, 
																																															 offsetEnd, 
																																															 &byteAlignment[0]);
		
			//
			// Convert to something we can compare easily.
			//
			ByteAlignmentToQueryString(&byteAlignment[0], byteAlignment.size(), &alignedSequence[0]);


			//
			// Do a sanity check to make sure the pulses and the alignment
			// make sense.  The main check is to see if the query sequence
			// in the alignment is the same as the query sequence in the
			// read. 
			//
		
			//
			// First pull out the bases corresponding to this read.
			//
			int queryStart = cmpFile.alnInfo.alignments[alignmentIndex].GetQueryStart();
			int queryEnd   = cmpFile.alnInfo.alignments[alignmentIndex].GetQueryEnd();

      // Build a map of where 
      CreateSequenceToAlignmentMap(byteAlignment, 
                                   baseToAlignmentMap);


			//
			// Condense gaps in the alignment for easy comparison.
			//
      //
      RemoveGaps(alignedSequence, alignedSequence);
      
			
			//
			// Query the cmp file for a way to look up a read based on
			// coordinate information.  For Astro reads, the coords are
			// based on x and y.  For Springfield, it is read index.  The
			// base files should be able to look up reads by x,y or by
			// index. 
			//
			int readIndex;

          
      if (cmpFile.platformId == Astro) {
        cout << "ASTRO pulse loading is deprecated." << endl;
        exit(0);
      }

      if (baseFile.LookupReadIndexByHoleNumber(holeNumber, readIndex) == false) {
          cout << "ERROR! Alignment has hole number " << holeNumber << " that is not in the movie. " << endl;
          assert(0);
      }

			int readStart, readLength, alignBaseStart, alignBaseEnd, alignBaseLength;
			readStart       = baseFile.readStartPositions[readIndex];
			readLength      = baseFile.readStartPositions[readIndex+1] - baseFile.readStartPositions[readIndex];
			alignBaseStart  = readStart + queryStart;
			alignBaseEnd    = readStart + queryEnd;
			alignBaseLength = alignBaseEnd - alignBaseStart;

			int pulseStart;
			if (usePulseFile) {
				pulseStart      = pulseFile.pulseStartPositions[readIndex];
			}

	    
      //
      // This maps from pulse to a base, since there are more pulses
      // called than bases, and the is one pulse for every base.
      //
			pulseIndexArray.resize(readLength);

      
      SMRTSequence sourceRead;
      unsigned int numPasses;
      //
      // These are not allocated in the regular allocate function
      // since they are only used in loadPulses. (maybe I should
      // subclass SMRTSequence here).
      //
      
      if (byRead) {
        // Read in the data from the bas file if it exsts.
        if (useBaseFile) {
          hdfBasReader.GetReadAt(readIndex, sourceRead);
          if (cmpFile.readType == ReadType::CCS or useCcs) {
            numPasses = hdfCcsReader.GetNumPasses(readIndex);
          }
        }
        // Read in the data from the pls file if it exists.
        if (usePulseFile) {
          hdfPlsReader.GetReadAt(readIndex, sourceRead.pulseIndex, sourceRead);
        }
      }
      else {
        //
        // The entire base/pulse file was read in, so copy data from that into a read
        // For the data used in the read, it is possible to simply
        // reference the data,  but for the pls file it is necessary
        // to copy since there is a packing of data.
        //
        if (useBaseFile) {
          baseFile.CopyReadAt(readIndex, sourceRead);
          if (cmpFile.readType == ReadType::CCS or useCcs) {
            numPasses = hdfCcsReader.GetNumPasses(readIndex);
          }
        }
        if (usePulseFile) {
          //
          // Copy the subset of pulses that correspond to the ones called as bases.
          //
          int i;
          for (i = 0; i < readLength; i++) {
            pulseIndexArray[i] = pulseStart + baseFile.pulseIndex[readStart + i];
          }
          pulseFile.CopyReadAt(readIndex, &pulseIndexArray[0], sourceRead);
        }
      }

      readSequence.resize(queryEnd - queryStart);
      CapQualityValues(sourceRead);
			copy((char*) (sourceRead.seq + queryStart),
					 (char*) (sourceRead.seq + queryEnd),
					 readSequence.begin());
      
			bool stringsMatch = true;
			if (alignedSequence.size() != readSequence.size() or alignedSequence != readSequence) {
				cout << "ERROR, the query sequence does not match the aligned query sequence." << endl;
				cout << "HoleNumber: "<< holeNumber << ", MovieName: " << cmpFileMovieName;
        cout << " ,ReadIndex: " << (int) readIndex << 
				cout << ", qStart: "<< queryStart << ", qEnd: " << queryEnd << endl;
				cout << "Aligned sequence: "<< endl;
				cout << alignedSequence << endl;
				cout << "Original sequence: " << endl;
				cout << readSequence << endl;
				assert(0);
      }

			/*
			 * Compute any necessary data fields.  These usually involve
			 * using differences of pulse indices, pulse widths, etc..
			 * Missing fields are stored as 0's. 
			 */

			vector<float> readPulseMetric;
            vector<float> floatMetric;
      vector<UChar> qvMetric;
      vector<HalfWord> frameRateMetric;
      vector<uint32_t> timeMetric;
			int ungappedAlignedSequenceLength = alignedSequence.size();
			
      floatMetric.resize(alignedSequenceLength+1);
      readPulseMetric.resize(alignedSequenceLength+1);
      qvMetric.resize(alignedSequenceLength+1);
      frameRateMetric.resize(alignedSequenceLength+1);
      timeMetric.resize(alignedSequenceLength+1);

			UInt i;
			UInt pi;

			HDFCmpExperimentGroup* expGroup = cmpReader.refAlignGroups[refGroupIndex]->readGroups[readGroupIndex];

      if (cmpFile.readType == ReadType::CCS or useCcs) {
        if (!cmpReader.alnInfoGroup.numPasses.IsInitialized()) {
          cmpReader.alnInfoGroup.InitializeNumPasses();
        }
        cmpReader.alnInfoGroup.numPasses.WriteToPos(&numPasses, 1, alignmentIndex);
      }
      
			if (metricOptions["StartTimeOffset"] == true) {
				if (!expGroup->startTimeOffset.IsInitialized()) {
					expGroup->startTimeOffset.Initialize(expGroup->experimentGroup, "StartTimeOffset");
				}
        unsigned int readStartTimeOffset = sourceRead.startFrame[queryStart];
				expGroup->startTimeOffset.WriteToPos(&readStartTimeOffset, 1, alignmentIndex);
			}

			if (metricOptions["QualityValue"] == true) {
				if (!expGroup->qualityValue.IsInitialized()) {
					expGroup->qualityValue.Initialize(expGroup->experimentGroup, "QualityValue");
				}
				
				// Store start time normalized to frame rate.
        fill(qvMetric.begin(), qvMetric.end(), missingQualityValue);

				for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
          qvMetric[baseToAlignmentMap[i]] = sourceRead.qual[queryStart + i];
        }
				qvMetric[qvMetric.size()-1] = 0;
				expGroup->qualityValue.WriteToPos(&qvMetric[0], qvMetric.size(), offsetBegin);
			}

			if (metricOptions["InsertionQV"] == true) {
				if (!expGroup->insertionQV.IsInitialized()) {
					expGroup->insertionQV.Initialize(expGroup->experimentGroup, "InsertionQV");
				}
				
				// Store start time normalized to frame rate.
        fill(qvMetric.begin(), qvMetric.end(), missingQualityValue);
				for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
          qvMetric[baseToAlignmentMap[i]] = sourceRead.insertionQV[queryStart+ i];
				}
				qvMetric[qvMetric.size()-1] = 0;
				expGroup->insertionQV.WriteToPos(&qvMetric[0], qvMetric.size(), offsetBegin);
			}

			if (metricOptions["MergeQV"] == true) {
				if (!expGroup->mergeQV.IsInitialized()) {
					expGroup->mergeQV.Initialize(expGroup->experimentGroup, "MergeQV");
				}
				
				// Store start time normalized to frame rate.
        fill(qvMetric.begin(), qvMetric.end(), missingQualityValue);
				for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
          qvMetric[baseToAlignmentMap[i]] = sourceRead.mergeQV[queryStart+ i];
				}
				qvMetric[qvMetric.size()-1] = 0;
				expGroup->mergeQV.WriteToPos(&qvMetric[0], qvMetric.size(), offsetBegin);
			}

			if (metricOptions["DeletionQV"] == true) {
				if (!expGroup->deletionQV.IsInitialized()) {
					expGroup->deletionQV.Initialize(expGroup->experimentGroup, "DeletionQV");
				}
				
				// Store start time normalized to frame rate.
        fill(qvMetric.begin(), qvMetric.end(), missingQualityValue);
				for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
          qvMetric[baseToAlignmentMap[i]] = sourceRead.deletionQV[queryStart+i];
				}
				qvMetric[qvMetric.size()-1] = 0;
				expGroup->deletionQV.WriteToPos(&qvMetric[0], qvMetric.size(), offsetBegin);
			}

			if (metricOptions["DeletionTag"] == true) {
				if (!expGroup->deletionTag.IsInitialized()) {
					expGroup->deletionTag.Initialize(expGroup->experimentGroup, "DeletionTag");
				}
        vector<char> readDeletionTagMetric;
        readDeletionTagMetric.resize(readPulseMetric.size());
				// Store start time normalized to frame rate.
				for (i = 0; i < readDeletionTagMetric.size()-1; i++ ) {
					readDeletionTagMetric[i] = '-';
				}
        readDeletionTagMetric[i] = '\0';
				for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
          assert(baseToAlignmentMap[i] < readDeletionTagMetric.size());
					readDeletionTagMetric[baseToAlignmentMap[i]] = sourceRead.deletionTag[queryStart+i];
				}
				readDeletionTagMetric[readDeletionTagMetric.size()-1] = 0;
				expGroup->deletionTag.WriteToPos(&readDeletionTagMetric[0], readDeletionTagMetric.size(), offsetBegin);
			}

			if (metricOptions["PulseIndex"] == true) {
        
				if (!expGroup->pulseIndex.IsInitialized()) {
					expGroup->pulseIndex.Initialize(expGroup->experimentGroup, "PulseIndex");
				}
				vector<uint32_t> readPulseIndexMetric;
        fill(readPulseIndexMetric.begin(), readPulseIndexMetric.end(), missingPulseIndex);
        readPulseIndexMetric.resize(readPulseMetric.size());
				// Store start time normalized to frame rate.
        assert(readPulseIndexMetric.size() > 0);
				for (i = 0; i < readPulseIndexMetric.size(); i++ ) {
          readPulseIndexMetric[i] = 0;
				}
				for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
					readPulseIndexMetric[baseToAlignmentMap[i]] = sourceRead.pulseIndex[queryStart+i];
				}
				readPulseIndexMetric[readPulseIndexMetric.size()-1] = 0;
				expGroup->pulseIndex.WriteToPos(&readPulseIndexMetric[0], readPulseIndexMetric.size(), offsetBegin);
			}

			if (metricOptions["SubstitutionTag"] == true) {
				if (!expGroup->substitutionTag.IsInitialized()) {
					expGroup->substitutionTag.Initialize(expGroup->experimentGroup, "SubstitutionTag");
				}
				vector<char> readSubstitutionTagMetric;
        readSubstitutionTagMetric.resize(readPulseMetric.size());
				// Store start time normalized to frame rate.
				for (i = 0; i < readSubstitutionTagMetric.size()-1; i++ ) {
          readSubstitutionTagMetric[i] = '-';
				}
        readSubstitutionTagMetric[i] = '\0';
				for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
					readSubstitutionTagMetric[baseToAlignmentMap[i]] = sourceRead.substitutionTag[queryStart+i];
				}
				readSubstitutionTagMetric[readSubstitutionTagMetric.size()-1] = 0;
				expGroup->substitutionTag.WriteToPos(&readSubstitutionTagMetric[0], readSubstitutionTagMetric.size(), offsetBegin);
      }

			if (metricOptions["SubstitutionQV"] == true) {
				if (!expGroup->substitutionQV.IsInitialized()) {
					expGroup->substitutionQV.Initialize(expGroup->experimentGroup, "SubstitutionQV");
				}
				
				// Store start time normalized to frame rate.
        fill(qvMetric.begin(), qvMetric.end(), missingQualityValue);

				for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
          qvMetric[baseToAlignmentMap[i]] = sourceRead.substitutionQV[queryStart+i];
				}
				qvMetric[qvMetric.size()-1] = 0;
				expGroup->substitutionQV.WriteToPos(&qvMetric[0], qvMetric.size(), offsetBegin);
			}

			if (metricOptions["ClassifierQV"] == true) {
				
				if (!expGroup->classifierQV.IsInitialized()) {
					expGroup->classifierQV.Initialize(expGroup->experimentGroup, "ClassifierQV");			
				}
				// Store start time normalized to frame rate.
        fill(floatMetric.begin(), floatMetric.end(), missingQualityValue);

				for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
					floatMetric[baseToAlignmentMap[i]] = sourceRead.classifierQV[i+queryStart];
				}
				floatMetric[floatMetric.size()-1] = 0;
				expGroup->classifierQV.WriteToPos(&floatMetric[0], floatMetric.size(), offsetBegin);
			}

			if (metricOptions["StartFrame"] == true) {
				if (!expGroup->startTime.IsInitialized()) {
					expGroup->startTime.Initialize(expGroup->experimentGroup, "StartFrame");			
				}

        if (useBaseFile) {
          sourceRead.startFrame = new unsigned int[sourceRead.length];
          copy(sourceRead.preBaseFrames, &sourceRead.preBaseFrames[sourceRead.length], sourceRead.startFrame);
          for (i = 0; i < sourceRead.length-1; i++) {
            sourceRead.startFrame[i+1] += sourceRead.widthInFrames[i];
          }
          partial_sum(sourceRead.startFrame, &sourceRead.startFrame[sourceRead.length],  sourceRead.startFrame);
        }
				
				// Store start time normalized to frame rate.
        fill(timeMetric.begin(), timeMetric.end(), missingPulseIndex);
				for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
					timeMetric[baseToAlignmentMap[i]] = sourceRead.startFrame[i+queryStart];
				}
				timeMetric[timeMetric.size()-1] = 0;
				expGroup->startTime.WriteToPos(&timeMetric[0], timeMetric.size(), offsetBegin);
			}

			if (metricOptions["PulseWidth"] == true) {
				if (!expGroup->pulseWidth.IsInitialized()) {
					expGroup->pulseWidth.Initialize(expGroup->experimentGroup, "PulseWidth");			
				}
				// Store start time normalized to frame rate.
        fill(frameRateMetric.begin(), frameRateMetric.end(), missingFrameRateValue);

        //
        // For legacy reasons, it's possible the width in frames is
        // stored in the bas file. If this is the case, use the width
        // in frames there.  Otherwise, use the width in frames stored
        // in the pls file.
        for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
          frameRateMetric[baseToAlignmentMap[i]] = sourceRead.widthInFrames[queryStart + i];
        }
				frameRateMetric[frameRateMetric.size()-1] = 0;
				expGroup->pulseWidth.WriteToPos(&frameRateMetric[0], frameRateMetric.size(), offsetBegin);
			}

			if (metricOptions["PreBaseFrames"] == true) {
				if (!expGroup->preBaseFrames.IsInitialized()) {
					expGroup->preBaseFrames.Initialize(expGroup->experimentGroup, "PreBaseFrames");
				}
				// Compute width in frames normalized to frame rate.
        fill(frameRateMetric.begin(), frameRateMetric.end(), missingFrameRateValue);
				for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
					frameRateMetric[baseToAlignmentMap[i]] = sourceRead.preBaseFrames[i+queryStart];
				}
				frameRateMetric[frameRateMetric.size()-1] = 0;
				expGroup->preBaseFrames.WriteToPos(&frameRateMetric[0], frameRateMetric.size(), offsetBegin);
			}

			if (metricOptions["WidthInFrames"] == true) {
				if (!expGroup->widthInFrames.IsInitialized()) {
					expGroup->widthInFrames.Initialize(expGroup->experimentGroup, "WidthInFrames");
				}
				// Compute width in frames normalized to frame rate.
        fill(frameRateMetric.begin(), frameRateMetric.end(), missingFrameRateValue);

				for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
					if (usePulseFile) {
						frameRateMetric[baseToAlignmentMap[i]] = sourceRead.widthInFrames[i+queryStart];
					}
					else {
						frameRateMetric[baseToAlignmentMap[i]] = sourceRead.widthInFrames[i+queryStart];
          }
				}
				frameRateMetric[frameRateMetric.size()-1] = 0;
				expGroup->widthInFrames.WriteToPos(&frameRateMetric[0], frameRateMetric.size(), offsetBegin);
			}

			if (metricOptions["pkmid"] == true) {

				if (!expGroup->pkmid.IsInitialized()) {
					expGroup->pkmid.Initialize(expGroup->experimentGroup, "pkmid");
				}

				for (i = 0; i < readPulseMetric.size(); i++ ) {
          readPulseMetric[i] = NaN;
				}
				
				for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
          readPulseMetric[baseToAlignmentMap[i]] = sourceRead.midSignal[i+queryStart];
        }
				readPulseMetric[readPulseMetric.size()-1] = 0;
				expGroup->pkmid.WriteToPos(&readPulseMetric[0], readPulseMetric.size(), offsetBegin);
			}

			if (metricOptions["IPD"] == true) {
				if (!expGroup->ipd.IsInitialized()) {
					expGroup->ipd.Initialize(expGroup->experimentGroup, "IPD");
				}
        fill(frameRateMetric.begin(), frameRateMetric.end(), missingFrameRateValue);				

				for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {

					//
					// The IPD is undefined for the first base in a read.
          //
					if (usePulseFile ) {
						if (queryStart == 0 and i == 0) {
              frameRateMetric[baseToAlignmentMap[i]] = 0;
						}
						else {
              frameRateMetric[baseToAlignmentMap[i]] = (sourceRead.startFrame[i+queryStart]  
                                                        - sourceRead.startFrame[i+queryStart-1]
                                                        - sourceRead.widthInFrames[i+queryStart-1]);
						}
					}
					else if (useBaseFile) {
            frameRateMetric[baseToAlignmentMap[i]] = sourceRead.preBaseFrames[i + queryStart];
					}
        }
				frameRateMetric[frameRateMetric.size()-1] = 0;
				expGroup->ipd.WriteToPos(&frameRateMetric[0], frameRateMetric.size(), offsetBegin);			
			}

			
			if (metricOptions["Light"] == true) {
				if (!expGroup->light.IsInitialized()) {
					expGroup->light.Initialize(expGroup->experimentGroup, "Light");
				}
        fill(frameRateMetric.begin(), frameRateMetric.end(), missingFrameRateValue);
				for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
          frameRateMetric[baseToAlignmentMap[i]] = sourceRead.meanSignal[i+queryStart];
          frameRateMetric[baseToAlignmentMap[i]] = (frameRateMetric[baseToAlignmentMap[i]] * 
                                                    sourceRead.widthInFrames[i+queryStart]);
				}
				frameRateMetric[frameRateMetric.size()-1] = 0;
				expGroup->light.WriteToPos(&frameRateMetric[0], frameRateMetric.size(), offsetBegin);			
			}
    
      sourceRead.Free();
      Free(sourceRead.meanSignal);
      Free(sourceRead.maxSignal);
      Free(sourceRead.midSignal);
      Free(sourceRead.startFrame);
      Free(sourceRead.classifierQV);
      Free(sourceRead.widthInFrames);
		}

    if (byRead == true) {
      if (useBaseFile) {
        hdfBasReader.Close();
      }
      if (cmpFile.readType == ReadType::CCS or useCcs) {
        hdfCcsReader.Close();
      }
      if (usePulseFile) {
        hdfPlsReader.Close();
      }
    }
	} // done loading movies


	cmpReader.Close();
}
