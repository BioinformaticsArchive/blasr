#ifndef DATASTRUCTURES_READS_PULSE_BASE_COMMON_H_
#define DATASTRUCTURES_READS_PULSE_BASE_COMMON_H_

// 
// This includes values that both pulse and base files must have.
//

#include "ScanData.h"
class PulseBaseCommon {
 public:
	ScanData scanData;
	float GetFrameRate() {
		return scanData.frameRate;
	}
	unsigned int GetNumFrames() {
		return scanData.numFrames;
	}
	string GetMovieName() {
		return scanData.movieName;
	}

};

#endif
