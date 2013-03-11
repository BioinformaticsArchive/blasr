#ifndef ALIGNMENT_SET_TO_CMP_H5_ADAPTER_H_
#define ALIGNMENT_SET_TO_CMP_H5_ADAPTER_H_


#include "data/hdf/HDFCmpFile.h"
#include "AlignmentSet.h"
#include "datastructures/alignment/AlignmentCandidate.h"
#include "datastructures/alignment/ByteAlignment.h"
#include "algorithms/alignment/ScoreMatrices.h"

class RefIndex {
 public:
  string name;
  int id;
  RefIndex() {
    id = 0;
    name = "";
  }

  RefIndex(string n, int i) {
    name = n; id = i;
  }

  RefIndex(const RefIndex &rhs) {
    name = rhs.name;
    id   = rhs.id;
  }
};

bool ParsePBIReadName(string &readName, string &movieName, int &readIndex) {
  vector<string> tokens;
  ParseSeparatedList(readName, tokens, '/');
  if (tokens.size() < 3) {
    movieName = "";
    readIndex = 0;
    return false;
  }
  else {
    movieName = tokens[0];
    readIndex = atoi(tokens[1].c_str());
    return true;
  }
}

template<typename T_CmpFile>
class AlignmentSetToCmpH5Adapter {
 public:

  map<string, RefIndex> refToId;
  map<string, int> knownMovies;
  map<string, int> knownPaths;
  unsigned int numAlignments;
  map<string, int> refNameToIndex;
 
  void Initialize() {
      numAlignments = 0;
  }

  template<typename T_Reference>
  void StoreReferenceInfo(vector<T_Reference> &references, T_CmpFile &cmpFile) {
    int r;
    for (r = 0; r < references.size(); r++) {
      string sequenceName, md5;
      unsigned int id, length;
      sequenceName = references[r].GetSequenceName();
      md5 = references[r].GetMD5();
      length = references[r].GetLength();
      string name;
      id = cmpFile.AddReference(sequenceName, length, md5, name);
      refToId[sequenceName] = RefIndex(name, id);
      refNameToIndex[sequenceName] = r;
    }
  }

  int StoreMovieInfo(string movieName, T_CmpFile &cmpFile) {
    map<string,int>::iterator mapIt;
    mapIt = knownMovies.find(movieName);
    if (mapIt != knownMovies.end()) {
      return mapIt->second;
    }
    else {
      int id;
      id = cmpFile.movieInfoGroup.AddMovie(movieName);
      knownMovies[movieName] = id;
      return id;
    }
  }


  int StorePath(string path, T_CmpFile &cmpFile) {
    if (knownPaths.find(path) != knownPaths.end()) {
      return knownPaths[path];
    }
    else {
      int id = cmpFile.alnGroupGroup.AddPath(path);
      knownPaths[path] = id;
      return id;
    }
  }


  void RemoveGapsAtEndOfAlignment(AlignmentCandidate<> &alignment) {
    int numEndDel = 0, numEndIns = 0;
    if (alignment.gaps.size() > 0) {
      int lastGap = alignment.gaps.size() - 1;
      int g;
      for (g = 0; g < alignment.gaps[lastGap].size(); g++) {
        if (alignment.gaps[lastGap][g].seq == Gap::Target) {
          numEndIns += alignment.gaps[lastGap][g].length;
        }
        else if (alignment.gaps[lastGap][g].seq == Gap::Query) {
          numEndDel += alignment.gaps[lastGap][g].length;
        }
      }
    }
    alignment.qAlignedSeqLength -= numEndIns;
    alignment.tAlignedSeqLength -= numEndDel;
  }

  void StoreAlignmentCandidate(AlignmentCandidate<> &alignment, 
                               int alnSegment,
                               T_CmpFile &cmpFile, int moleculeNumber = -1) {
    //
    // Find out where the movie is going to get stored.
    //
    string movieName;
    int holeNumber = 0;
    bool nameParsedProperly;
    
    nameParsedProperly = ParsePBIReadName(alignment.qName, movieName, holeNumber);
    if ( !nameParsedProperly) {
      cout <<"ERROR. Attempting to store a read with name " << alignment.qName << " that does not " << endl
           << "appear to be a PacBio read." << endl;
      exit(1);
    }
  
    int movieId;
    movieId = StoreMovieInfo(movieName, cmpFile);
  
    map<string, RefIndex>::iterator refToIdIt;
    refToIdIt = refToId.find(alignment.tName);
    if (refToIdIt == refToId.end()) {
      cout <<"ERROR. The reference name " << alignment.tName << " was not found in the list of references." << endl;
      cout << "Perhaps a different reference file was aligned to than " << endl
           << "what was provided for SAM conversion. " << endl;
      exit(1);
    }
      
    string refGroupName = refToIdIt->second.name;
    int    refGroupId   = refToIdIt->second.id;

    if (cmpFile.refGroupIdToArrayIndex.find(refGroupId) == cmpFile.refGroupIdToArrayIndex.end()) {
      cout << "ERROR. The reference ID is not indexed. This is an internal inconsistency." << endl;
      exit(1);
    }
    int    refGroupIndex= cmpFile.refGroupIdToArrayIndex[refGroupId];
    assert(refGroupIndex + 1 == refGroupId);

    string path = "/" + refGroupName + "/" + movieName;
  
    int pathId = StorePath(path, cmpFile);
    int pathIndex = pathId - 1;
  
    vector<unsigned int> alnIndex;
    alnIndex.resize(22);


    RemoveGapsAtEndOfAlignment(alignment);
  
    /*
     * Store the alignment string
     */
    vector<unsigned char> byteAlignment;
    AlignmentToByteAlignment(alignment, 
                             alignment.qAlignedSeq, alignment.tAlignedSeq,
                             byteAlignment);

    unsigned int offsetBegin, offsetEnd;
    cmpFile.StoreAlnArray(byteAlignment, alignment.tName, movieName, offsetBegin, offsetEnd);

    numAlignments++;
    /*    EditDistanceMatrix scoreMat;
    */
    int tmpMatrix[5][5];
    // the 5,5 are indel penalties that do not matter since the score is not stored.
    ComputeAlignmentStats(alignment, alignment.qAlignedSeq.seq, alignment.tAlignedSeq.seq, tmpMatrix, 5, 5);

    /*
      The current AlnIndex column names:
      (0): "AlnID", "AlnGroupID", "MovieID", "RefGroupID", "tStart",
      (5): "tEnd", "RCRefStrand", "HoleNumber", "SetNumber",
      (9): "StrobeNumber", "MoleculeID", "rStart", "rEnd", "MapQV", "nM",
      (15): "nMM", "nIns", "nDel", "Offset_begin", "Offset_end",
      (20): "nBackRead", "nReadOverlap"
    */
    if (moleculeNumber == -1) {
      moleculeNumber = holeNumber * movieId;
    }
    alnIndex[0]  = numAlignments;  // AlnId
    alnIndex[1]  = pathId;        // AlnGroupID
    alnIndex[2]  = movieId;    // MovieID
    alnIndex[3]  = refGroupId; // RefGroupID
    alnIndex[4]  = alignment.tAlignedSeqPos; // tStart
    alnIndex[5]  = alignment.tAlignedSeqPos +  alignment.tAlignedSeqLength; // tEnd
    alnIndex[6]  = alignment.tStrand; // RCRefStrand
    alnIndex[7]  = holeNumber;
    alnIndex[8]  = 0; // SET NUMBER -- parse later!!!!
    alnIndex[9]  = alnSegment; // strobenumber
    alnIndex[10] = moleculeNumber;
    alnIndex[11] = alignment.qAlignedSeqPos; 
    alnIndex[12] = alignment.qAlignedSeqPos + alignment.qAlignedSeqLength;
    alnIndex[13] = alignment.mapQV;
    alnIndex[14] = alignment.nMatch;
    alnIndex[15] = alignment.nMismatch;
    alnIndex[16] = alignment.nIns;
    alnIndex[17] = alignment.nDel;
    alnIndex[18] = offsetBegin;
    alnIndex[19] = offsetEnd;
    alnIndex[20] = 0;
    alnIndex[21] = 0;
    cmpFile.alnInfoGroup.WriteAlnIndex(alnIndex);
  }

  void StoreAlignmentCandidateList(vector<AlignmentCandidate<> > &alignments, T_CmpFile &cmpFile,  int moleculeNumber=-1) {
    int a;
    for (a = 0; a < alignments.size(); a++) {
      StoreAlignmentCandidate(alignments[a], a, cmpFile, moleculeNumber);
    }
  }

  void StoreAlignmentCandidate(AlignmentCandidate<> alignment, 
                               T_CmpFile &cmpFile) {
    StoreAlignmentCandidate(alignment, 0, cmpFile);
  }

};



#endif
