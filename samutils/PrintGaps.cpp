#include "algorithms/alignment/readers/sam/SAMReader.h"
#include "datastructures/alignmentset/SAMAlignment.h"
#include "datastructures/alignmentset/SAMToAlignmentCandidateAdapter.h"
#include "FASTAReader.h"
#include "datastructures/alignment/AlignmentCandidate.h"
#include "algorithms/alignment/printers/StickAlignmentPrinter.h"
#include "CommandLineParser.h"
#include <map>
#include <set>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

class Interval {
public:
  int index;
  unsigned int start;
  unsigned int end;
  Interval(int i, unsigned int st, unsigned int en): index(i), start(st), end(en) {}
  Interval() {}
  bool operator<(const Interval &rhs) const {
    return start < rhs.start;
  }
};


typedef set<Interval> IntervalSet;
typedef map<string, IntervalSet> RegionMap;


void ReadROIFile(istream &roiIn, RegionMap &regions) {
  int lineNumber = 0;

  while(roiIn){
    string chr;
    unsigned int start, end;
    if ( (roiIn >> chr >> start >> end) ) {
      regions[chr].insert(Interval(lineNumber, start,end));
      ++lineNumber;
    }
  }
}


bool FindRegion(RegionMap &regions, string chr, unsigned int start, unsigned int end,
                Interval &foundRegion) {
  if (regions.find(chr) == regions.end()) {
    return false;
  }
  else {
    //
    // Scan this region for a match.  It is expected that start <
    // intv.start and end > intv.end.
    //
    IntervalSet::iterator it;
    Interval tmpIntv;
    tmpIntv.start = start;
    it = regions[chr].lower_bound(tmpIntv);
    while (it != regions[chr].end() and
           (*it).start < end) {
      if ((*it).start >= start and (*it).start < end) {
        foundRegion = (*it);
        return true;
      }
      ++it;
    }
    return false;
  }
}

int main(int argc, char* argv[]) {
  string samFileName, genomeFileName;
  string outFileName;
  string alignmentsFileName = "";
  int minGapLength = 100;
  string gapFastaFileName = "";
  string regionOfInterestFileName = "";
  outFileName = "";
  CommandLineParser clp;
  int minMapQV = 0;
  clp.RegisterStringOption("sam", &samFileName, "Alignments.", true);
  clp.RegisterStringOption("genome", &genomeFileName, "Genome.", true);
  clp.RegisterPreviousFlagsAsHidden();
  clp.RegisterStringOption("out", &outFileName, "Output file. Default to stdout", false);
  clp.RegisterIntOption("minGapLength", &minGapLength, "Minimum gap length", CommandLineParser::NonNegativeInteger, false);
  clp.RegisterStringOption("printAlignments", &alignmentsFileName, "Print human readable alignments.", false);
  clp.RegisterStringOption("printGapFasta", &gapFastaFileName, "Print the sequences of gaps to a file.", false);
  clp.RegisterStringOption("roi", &regionOfInterestFileName, "Regions of interest, comma separated file with 'chr,start,end'."
                           " Only find alignments within these regions.", false);
  clp.RegisterIntOption("minqv", &minMapQV, "Minimum mapping quality value to allow", CommandLineParser::NonNegativeInteger, false);
  clp.ParseCommandLine(argc, argv);
  SAMReader<SAMFullReferenceSequence, SAMReadGroup, SAMPosAlignment> samReader;
  FASTAReader fastaReader;
  samReader.Initialize(samFileName);
  fastaReader.Initialize(genomeFileName);
  AlignmentSet<SAMFullReferenceSequence, SAMReadGroup, SAMPosAlignment> alignmentSet;
  samReader.ReadHeader(alignmentSet)
;
  ostream *outPtr;
  ofstream outFile;
  ofstream alignmentsOut, gapFastaOut;
  if (outFileName == "") {
    outPtr = &cout;
  }
  else {
    CrucialOpen(outFileName, outFile, std::ios::out);
    outPtr = &outFile;
  }
  if (alignmentsFileName != "") {
    CrucialOpen(alignmentsFileName, alignmentsOut, std::ios::out);
  }

  if (gapFastaFileName != "") {
    CrucialOpen(gapFastaFileName, gapFastaOut, std::ios::out);
  }

  ifstream roiIn;
  RegionMap regions;
  bool findRegion = false;
  if (regionOfInterestFileName != "") {
    CrucialOpen(regionOfInterestFileName, roiIn, std::ios::in);
    ReadROIFile(roiIn, regions);
    findRegion = true;
  }
  vector<FASTASequence> references;

  fastaReader.ReadAllSequences(references);
  int i;
  map<string,int> refToIndex;
  for (i = 0; i < references.size(); i++) {
    refToIndex[references[i].title] = i;
  }
  SAMAlignment samAlignment;  

  if (findRegion) {
    cout << "region index "
         << "qname" << " "
         << "tname" << " " 
         << "rgnStart" << " " << "rgnEnd " << " "
         << "gapLength" << " " 
         << "tStart" << " " 
         << "tEnd " << " " 
         << "tGap " << " " << "strand" << endl;
  }
  int index = 0;

  while (samReader.GetNextAlignment(samAlignment)) {
    if (samAlignment.rName == "*") {
      continue;
    }
    vector<AlignmentCandidate<> > convertedAlignments;
    SAMAlignmentsToCandidates(samAlignment, 
                              references, refToIndex,
                              convertedAlignments, true, false);


    int a;
    int g, gi;
    
    for (a = 0; a < convertedAlignments.size(); a++) {
      string textStr, alignStr, queryStr;
      CreateAlignmentStrings( convertedAlignments[a],
														  convertedAlignments[a].qAlignedSeq,
                              convertedAlignments[a].tAlignedSeq,
                              textStr, alignStr, queryStr);
      int window = 10;
      if (textStr.size() < window) {
        continue;
      }
      int j,k;
      for (j = 0; j < textStr.size() - window + 1; j+=window) {
        int nm = 0, ng = 0;
        for (k = j; k < j + window; k++) {
          if (textStr[k] == '-' or queryStr[k] == '-') {
            ng++;
          }
          else {nm++;}
        }
        //        cout << index << " " << j << " " << (1.0*nm)/window << endl;
      }
      ++index;
      if (convertedAlignments[a].mapQV < minMapQV) {
        continue;
      }
      if (findRegion) {
        Interval foundRegion;
        int refIndex = refToIndex[convertedAlignments[a].tName];
        DNALength refLength = references[refIndex].length;
        DNALength tStart, tEnd;
        tStart = convertedAlignments[a].GenomicTBegin();
        tEnd   = convertedAlignments[a].GenomicTEnd();
        
        if (FindRegion(regions, convertedAlignments[a].tName, 
                       tStart, tEnd, foundRegion)) {

          //
          // Now try and look for the gap
          //
          int b;
          unsigned int alnStart = convertedAlignments[a].GenomicTBegin();
          if (convertedAlignments[a].blocks.size() == 0) {
            continue;
          }
          int gapLength = foundRegion.end - foundRegion.start;
          int minGapDiff = 1000000;
          int minGapDiffIndex = 0;
          int minGapDiffG = 0;
          for (b = 0; b < convertedAlignments[a].blocks.size()  - 1; b++) {
            int g;
            
            for (g = 0; g < convertedAlignments[a].gaps[g+1].size(); g++) {
              //              cout << b << " " << convertedAlignments[a].blocks[b].tPos << " g " << g << " " << (int) convertedAlignments[a].gaps[b+1][g].length << endl;
              if ( abs((int) convertedAlignments[a].gaps[b+1][g].length -  gapLength) < minGapDiff) {
                minGapDiff = abs((int) convertedAlignments[a].gaps[b+1][g].length -  gapLength);
                minGapDiffIndex = b;
                minGapDiffG = g;
              }
            }
          }
          b = minGapDiffIndex;
          cout << foundRegion.index << " "
               << convertedAlignments[a].qName << " "
               << convertedAlignments[a].tName << " " 
               << foundRegion.start << " " << foundRegion.end << " " 
               << gapLength << " " 
               << convertedAlignments[a].tAlignedSeqPos + convertedAlignments[a].blocks[b].tPos << " " 
               << convertedAlignments[a].tAlignedSeqPos + convertedAlignments[a].blocks[b+1].tPos << " " 
               << convertedAlignments[a].gaps[b+1][minGapDiffG].length << " " << convertedAlignments[a].tStrand << endl;

          if (alignmentsFileName != "") {
            StickPrintAlignment(convertedAlignments[a], convertedAlignments[a].qAlignedSeq, convertedAlignments[a].tAlignedSeq, alignmentsOut, convertedAlignments[a].qAlignedSeqPos,convertedAlignments[a].tAlignedSeqPos  );
          }
          /*
            if ((convertedAlignments[a].tStrand == convertedAlignments[a].qStrand and 
                 (convertedAlignments[a].blocks[b].tPos + alnStart <= foundRegion.start and
                  convertedAlignments[a].blocks[b+1].tPos + alnStart >= foundRegion.start)) or
                (convertedAlignments[a].tStrand != convertedAlignments[a].qStrand and 
                 ((alnStart + (convertedAlignments[a].TEnd() - (convertedAlignments[a].blocks[b+1].tPos + convertedAlignments[a].blocks[b+1].length - 1)) <= foundRegion.start) and
                  (alnStart + (convertedAlignments[a].TEnd() - (convertedAlignments[a].blocks[b].tPos +1) )) >= foundRegion.start))) {
              unsigned int len = foundRegion.end - foundRegion.start;
              if (convertedAlignments[a].gaps[b+1].size() > 0) {

       
              }
            }
            }
          */
        }
      }
      else {
        for (g = 0; g < convertedAlignments[a].gaps.size(); g++) {
          for (gi = 0; gi < convertedAlignments[a].gaps[g].size(); gi++) {
            if (convertedAlignments[a].gaps[g][gi].length > minGapLength) {
              FASTASequence gap;
              DNALength gapStartPos;
              DNALength gapSeqPos;
              string type = "";
              if (convertedAlignments[a].gaps[g][gi].seq == Gap::Query) {
                gapStartPos = convertedAlignments[a].tPos;
                if (g > 0) {
                  gapStartPos += convertedAlignments[a].blocks[g-1].tPos + convertedAlignments[a].blocks[g-1].length;
                }
              
                gap.CopyTitle(convertedAlignments[a].tName);
                assert(gapStartPos + convertedAlignments[a].gaps[g][gi].length <= convertedAlignments[a].tAlignedSeq.length);       gap.seq = &convertedAlignments[a].tAlignedSeq.seq[gapStartPos];
                gap.length = convertedAlignments[a].gaps[g][gi].length;
                type = "deletion";
              }
              else {
                gapStartPos = convertedAlignments[a].qPos;
                if (g > 0) {
                  gapStartPos += convertedAlignments[a].blocks[g-1].qPos + convertedAlignments[a].blocks[g-1].length;
                }
                gapSeqPos = gapStartPos + convertedAlignments[a].qAlignedSeqPos;
                gap.CopyTitle(convertedAlignments[a].qName);
                assert(gapStartPos + convertedAlignments[a].gaps[g][gi].length <= convertedAlignments[a].qAlignedSeq.length);
                gap.seq = &convertedAlignments[a].qAlignedSeq.seq[gapStartPos];
                gap.length = convertedAlignments[a].gaps[g][gi].length;
                type = "insertion";
              }
              gapSeqPos = (convertedAlignments[a].blocks[g-1].tPos + 
                           convertedAlignments[a].blocks[g-1].length + 
                           + convertedAlignments[a].tAlignedSeqPos);
              string gapSeq;
              gapSeq.insert(0, (char*) gap.seq, gap.length);
              *outPtr << convertedAlignments[a].tName << "\tblasr\t"
                      << type << "\t" 
                      << gapSeqPos << "," << gapSeqPos + gap.length << "\t"
                      << gap.length << "\t";
              if (convertedAlignments[a].qStrand == convertedAlignments[a].tStrand) {
                *outPtr << "+\t";
              }
              else {
                *outPtr << "-\t";
              }
              *outPtr << "0\t";
              *outPtr << "seq " << gapSeq << endl;

              if (alignmentsFileName != "") {
                // hack to get around SAM problem for now
                StickPrintAlignment(convertedAlignments[a], convertedAlignments[a].qAlignedSeq, convertedAlignments[a].tAlignedSeq, alignmentsOut, convertedAlignments[a].qAlignedSeqPos,convertedAlignments[a].tAlignedSeqPos  );
              }
              if (gapFastaFileName != "") {
                gap.PrintSeq(gapFastaOut);
              }
            }
          }
        }
      }
    }
  }
}
