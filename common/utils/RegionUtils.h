#ifndef UTILS_REGION_UTILS_H_
#define UTILS_REGION_UTILS_H_
#include <algorithm>
#include "SMRTSequence.h"
#include "datastructures/reads/ReadInterval.h"
#include "datastructures/reads/RegionTable.h"

bool LookupHQRegion(int holeNumber, RegionTable &regionTable, int &start, int &end, int &score) {
	int regionLowIndex, regionHighIndex;
	regionLowIndex = regionHighIndex = 0;
	regionTable.LookupRegionsByHoleNumber(holeNumber, regionLowIndex, regionHighIndex);
	bool readHasGoodRegion = true;
	int  regionIndex = regionLowIndex;
	while (regionIndex < regionHighIndex and 
				 regionTable.GetType(regionIndex) != HQRegion) {
		regionIndex++;
	}
	
	if (regionIndex == regionHighIndex) {
    start = end = score = 0;
		return false;
	}
	else {
		start = regionTable.GetStart(regionIndex);
		end   = regionTable.GetEnd(regionIndex);
    score = regionTable.GetScore(regionIndex);
		return true;
	}
}


template<typename T_Sequence>
bool MaskRead(T_Sequence &fastaRead, ZMWGroupEntry &zmwData, RegionTable &regionTable) {
	int regionIndex;						 
	int regionLowIndex, regionHighIndex;
	regionLowIndex = regionHighIndex = 0;
	regionTable.LookupRegionsByHoleNumber(zmwData.holeNumber, regionLowIndex, regionHighIndex);
	bool readHasGoodRegion = true;

	DNALength readPos;

	regionIndex = regionLowIndex;
	int lastHQRegionIndex;
	
	int hqRegionStart=0, hqRegionEnd=0, hqRegionScore = 0;
	readHasGoodRegion = LookupHQRegion(zmwData.holeNumber, regionTable, hqRegionStart, hqRegionEnd, hqRegionScore);

	//
	// Mask off the low quality portion of this read.
	//
	for (readPos = 0; (readPos < hqRegionStart and
											 readPos < fastaRead.length); readPos++) {
		fastaRead.seq[readPos] = 'N';
	}
	for (readPos = hqRegionEnd; readPos < fastaRead.length; readPos++) {
		fastaRead.seq[readPos] = 'N';
	}

	//
	// Look to see if there is region information provided, but the entire read is bad.
	//
	if (hqRegionEnd == hqRegionStart) {
		//
		// This read is entirely bad, flag that.
		//
		readHasGoodRegion = false;
	}

	return readHasGoodRegion;
}


template<typename T_Sequence>
bool GetReadTrimCoordinates(T_Sequence &fastaRead,
														ZMWGroupEntry &zmwData,
														RegionTable &regionTable,
														DNALength &readStart,
														DNALength &readEnd,
                            int &score) {

	int regionIndex;						 
	int regionLowIndex, regionHighIndex;
	regionLowIndex = regionHighIndex = 0;
	regionTable.LookupRegionsByHoleNumber(zmwData.holeNumber, regionLowIndex, regionHighIndex);
	bool readHasGoodRegion = true;

	DNALength readPos;

	regionIndex = regionLowIndex;
	int lastHQRegionIndex;
	
	while (regionIndex < regionHighIndex and 
				 regionTable.GetType(regionIndex) != HQRegion) {
		regionIndex++;
	}
	
	if (regionIndex < regionHighIndex ) {
		readStart = regionTable.GetStart(regionIndex);
		readEnd   = regionTable.GetEnd(regionIndex);
    score     = regionTable.GetScore(regionIndex);
		return true;
	}
	else {
		readStart = 0;
		readEnd   = fastaRead.length;
		return false;
	}
}

template<typename T_Sequence>
bool TrimRead(T_Sequence &fastaRead, ZMWGroupEntry &zmwData, RegionTable &regionTable, T_Sequence &trimmedRead) {

	DNALength readStart, readEnd;
	GetReadTrimCoordinates(fastaRead, zmwData, regionTable, readStart, readEnd);
	if (readEnd - readStart > 0) {
		trimmedRead.CopySubsequence((FASTQSequence&)fastaRead, 
																readStart, readEnd);
		// signal that the read has a good region.
		return true;
	}
	else {

		//
		// There is no information for this read. Make it skipped.
		//
		trimmedRead.seq = NULL;
		trimmedRead.CopyTitle(fastaRead.title);
		// signal this read has no good region.
		return false;
	}
}

class CompareRegionIndicesByStart {
public:
	RegionTable *regionTablePtr;
	int operator()(const int a, const int b) const {
		if (regionTablePtr->GetStart(a) == regionTablePtr->GetStart(b)) {
			return (regionTablePtr->GetEnd(a) < regionTablePtr->GetEnd(b));
		}
		else {
			return (regionTablePtr->GetStart(a) < regionTablePtr->GetStart(b));
		}
	}
};

		
int SortRegionIndicesByStart(RegionTable &regionTable, vector<int> &indices) {
	CompareRegionIndicesByStart cmpFctr;
	cmpFctr.regionTablePtr = &regionTable;
	std::sort(indices.begin(), indices.end(), cmpFctr);
  return indices.size();
}

class OrderRegionsByReadStart {
 public:
	int operator()(const ReadInterval &lhs, const ReadInterval &rhs) const {
		return lhs.start < rhs.start;
	}
};

int FindRegionIndices(unsigned int holeNumber, RegionTable *regionTablePtr, int &regionLowIndex, int &regionHighIndex) {
	int regionIndex;						 
  regionLowIndex = regionHighIndex = 0;
	regionTablePtr->LookupRegionsByHoleNumber(holeNumber, regionLowIndex, regionHighIndex);  
  return regionHighIndex - regionLowIndex;
}

int FindRegionIndices(SMRTSequence &read, RegionTable *regionTablePtr, int &regionLowIndex, int &regionHighIndex) {
  return FindRegionIndices(read.zmwData.holeNumber, regionTablePtr, regionLowIndex, regionHighIndex);
}

//
// Collect region indices for either all region types, or just a few specific region types.
//

int CollectRegionIndices(SMRTSequence &read, RegionTable &regionTable, vector<int> &regionIndices,
                         RegionType *regionTypes=NULL, int numRegionTypes = 0) {
  int regionLow, regionHigh;
  int prevNumRegionIndices = regionIndices.size();
  if (FindRegionIndices(read, &regionTable, regionLow, regionHigh)) {
    int i;
    for (i = regionLow; i < regionHigh; i++) {
      if (regionTypes == NULL) {
        regionIndices.push_back(i);
      }
      else {
        int t;
        for (t = 0; t < numRegionTypes; t++) {
          if (regionTable.GetType(i) == regionTypes[t]) {
            regionIndices.push_back(i);
            break;
          }
        }
      }
    }
  }
  return regionIndices.size() - prevNumRegionIndices;
}




template<typename T_Sequence>
int CollectSubreadIntervals(T_Sequence &read, RegionTable *regionTablePtr, vector<ReadInterval> &subreadIntervals, bool byAdapter=false) {
	int regionIndex;						 
	int regionLowIndex, regionHighIndex;
	regionLowIndex = regionHighIndex = 0;
	regionTablePtr->LookupRegionsByHoleNumber(read.zmwData.holeNumber, regionLowIndex, regionHighIndex);
	if (byAdapter == false) {
		for (regionIndex = regionLowIndex; regionIndex < regionHighIndex; regionIndex++) {
			if (regionTablePtr->GetType(regionIndex) ==  Insert) {
				subreadIntervals.push_back(ReadInterval(regionTablePtr->table[regionIndex].row[RegionAnnotation::RegionStart],
																								regionTablePtr->table[regionIndex].row[RegionAnnotation::RegionEnd],
                                                regionTablePtr->table[regionIndex].row[RegionAnnotation::RegionScore]));
			}
		}
	}
	else {
		vector<int> adapterIntervalIndices;
		for (regionIndex = regionLowIndex; regionIndex < regionHighIndex; regionIndex++) {
			if (regionTablePtr->GetType(regionIndex) == Adapter) {
				adapterIntervalIndices.push_back(regionIndex);
			}
		}
		// Sort indices so that the intervals appear in order on the read.
		SortRegionIndicesByStart(*regionTablePtr, adapterIntervalIndices);
		int curIntervalStart = 0;
		int i;
		if (adapterIntervalIndices.size() == 0) {
			subreadIntervals.push_back(ReadInterval(0, read.length));
		}
		else {
			subreadIntervals.push_back(ReadInterval(0, regionTablePtr->table[adapterIntervalIndices[0]].row[RegionAnnotation::RegionStart]));
			for (i = 0; i + 1 < adapterIntervalIndices.size() ; i++) {
				subreadIntervals.push_back(ReadInterval(regionTablePtr->table[adapterIntervalIndices[i  ]].row[RegionAnnotation::RegionEnd],
																								regionTablePtr->table[adapterIntervalIndices[i+1]].row[RegionAnnotation::RegionStart]));
			}
			subreadIntervals.push_back(ReadInterval(regionTablePtr->table[adapterIntervalIndices[adapterIntervalIndices.size()-1]].row[RegionAnnotation::RegionEnd],
																							read.length));
		}
	}
	sort(subreadIntervals.begin(), subreadIntervals.end(), OrderRegionsByReadStart());
}

template<typename T_Sequence>
int RemoveLowQualitySubreads(T_Sequence &read, RegionTable *regionTablePtr, vector<ReadInterval> &subreadIntervals, int highQualityStart, int highQualityEnd ) {
  
  int i;
  for (i = 0; i < subreadIntervals.size();) {
    if (highQualityStart > subreadIntervals[i].end or highQualityEnd < subreadIntervals[i].start) {
      subreadIntervals.erase(subreadIntervals.begin() + i);
    }
    else {
      if (highQualityStart > subreadIntervals[i].start and highQualityStart < subreadIntervals[i].end) {
        subreadIntervals[i].start = highQualityStart;
      }        
    }
  }
}


#endif
