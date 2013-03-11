#ifndef DATASTRUCTURES_READS_SCAN_DATA_H_
#define DATASTRUCTURES_READS_SCAN_DATA_H_

#include "../../Enumerations.h"
#include "../../data/hdf/PlatformId.h"
class ScanData {
 public:
	PlatformId platformId;
	float frameRate;
	unsigned int numFrames;
	string movieName, runCode;
	string whenStarted;
	string GetMovieName() {
		return movieName;
	}

};

#endif
