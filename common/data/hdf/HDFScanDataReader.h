#ifndef DATA_HDF_HDF_SCAN_DATA_READER_H_
#define DATA_HDF_HDF_SCAN_DATA_READER_H_

#include "HDFGroup.h"
#include "HDFAtom.h"
#include "PlatformId.h"
#include "../../Enumerations.h"
#include "../../datastructures/reads/ScanData.h"
//
// The SanDataReader cannot live outside 

class HDFScanDataReader {
 public:
	bool fileHasScanData, useRunCode;
	HDFGroup scanDataGroup;
	HDFGroup acqParamsGroup;
	HDFGroup runInfoGroup;
	bool initializedAcqParamsGroup, initializedRunInfoGroup;
	bool useWhenStarted;
	HDFAtom<string> whenStartedAtom;
	HDFAtom<unsigned int> platformIdAtom;
	HDFAtom<float> frameRateAtom;
	HDFAtom<unsigned int> numFramesAtom;
	HDFAtom<string> movieNameAtom;
	HDFAtom<string> runCodeAtom;
	//
	// It is useful to cache the movie name in the reader since this is
	// loaded once upon initialization, and may be fetched when loading
	// reads one at a time.
	//
	bool   useMovieName;
	string movieName, runCode;
	PlatformId platformId;
	HDFScanDataReader() {
		//
		// Assume the file is written without a movie name.  This is
		// flipped when a movie name is found.
		//
		useMovieName    = false;
		useRunCode      = false;
		useWhenStarted  = false;
		fileHasScanData = false;
		movieName = "";
		runCode   = "";
    platformId      = NoPlatform;
	}

	int InitializeAcqParamsAtoms() {
		if (frameRateAtom.Initialize(acqParamsGroup.group, "FrameRate") == 0) { return 0; }
		if (numFramesAtom.Initialize(acqParamsGroup.group, "NumFrames") == 0) { return 0; }
		if (acqParamsGroup.ContainsAttribute("WhenStarted")) {
			if (whenStartedAtom.Initialize(acqParamsGroup.group, "WhenStarted") == 0) { return 0; }
			useWhenStarted = true;
		}
		return 1;
	}

	//
	// This is created on top of a file that is already opened, so
	// instead of initializing by opening a file, it is initialized by
	// passing the root group that should contain the ScanData group.  
	// When shutting down, no file handle needs to be closed.
	//
	int Initialize(HDFGroup *pulseDataGroup) {

		//
		// Initiailze groups for reading data.
		//
		
		initializedAcqParamsGroup = false;
		initializedRunInfoGroup   = false;
		if (pulseDataGroup->ContainsObject("ScanData") == 0 or 
				scanDataGroup.Initialize(pulseDataGroup->group, "ScanData") == 0) {
			return 0;
		}
		fileHasScanData = true;

		if (scanDataGroup.ContainsObject("AcqParams") == 0) {
			return 0;
		}
		if (acqParamsGroup.Initialize(scanDataGroup.group, "AcqParams") == 0) {
			return 0;
		}
		initializedAcqParamsGroup = true;

		if (scanDataGroup.ContainsObject("RunInfo") == 0 or
				(runInfoGroup.Initialize(scanDataGroup.group, "RunInfo") == 0)) {
			return 0;
		}
		initializedRunInfoGroup = true;
		if (InitializeAcqParamsAtoms() == 0) {
			return 0;
		}

		//
		// Read in the data that will be used later on either per read or
		// when the entire bas/pls file is read.
		//
		
		if (ReadPlatformId(platformId) == 0) {
			return 0;
		}

		if (runInfoGroup.ContainsAttribute("RunCode") and
				runCodeAtom.Initialize(runInfoGroup, "RunCode")) {
			useRunCode = true;
		}

		//
		// Attempt to load the movie name.  This is not always present.
		//
		LoadMovieName(movieName);
		
		return 1;
	}

	string GetMovieName() {
		LoadMovieName(movieName);
		return movieName;
	}
	
	string GetRunCode() {
		return runCode;
	}

	int Read(ScanData &scanData) {
		// All parameters below are required.
		if (ReadPlatformId(scanData.platformId) == 0) return 0;
		LoadMovieName(scanData.movieName);
		
		if (useRunCode) {
			runCodeAtom.Read(scanData.runCode);
		}
		frameRateAtom.Read(scanData.frameRate);
		numFramesAtom.Read(scanData.numFrames);
		
		if (useWhenStarted) {
			whenStartedAtom.Read(scanData.whenStarted);
		}

	}

	void ReadWhenStarted(string &whenStarted) {
		whenStartedAtom.Read(whenStarted);
	}

	PlatformId GetPlatformId() {
		return platformId;
	}

	int ReadPlatformId(PlatformId &pid) {
		if (runInfoGroup.ContainsAttribute("PlatformId")) {
			if (platformIdAtom.Initialize(runInfoGroup, "PlatformId") == 0) {
				return 0;
			}
			platformIdAtom.Read((unsigned int&)pid);
		}
		else {
			pid = AstroPlatform;
		}
		return 1;
	}

	int LoadMovieName(string &movieName) {
		// Groups for building read names
		if (runInfoGroup.ContainsAttribute("MovieName") and
				movieNameAtom.Initialize(runInfoGroup, "MovieName")) {
			useMovieName = true;
			movieNameAtom.Read(movieName);
			int e = movieName.size() - 1;
			while (e > 0 and movieName[e] == ' ') e--;
			movieName= movieName.substr(0,e+1);
			return 1;
		}
		else {
			return 0;
		}
	}

	void Close() {
		if (useMovieName) {
			movieNameAtom.dataspace.close();
		}
		if (useRunCode) {
			runCodeAtom.dataspace.close();
		}
		scanDataGroup.Close();
		acqParamsGroup.Close();
		runInfoGroup.Close();
	}


};

#endif
