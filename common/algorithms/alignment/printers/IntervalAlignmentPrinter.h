#ifndef INTERVAL_ALIGNMENT_PRINTER_H_
#define INTERVAL_ALIGNMENT_PRINTER_H_

#include "../../../datastructures/alignment/AlignmentCandidate.h"
#include "../../../FASTQSequence.h"

class IntervalAlignmentPrinter {
 public:
  static void Print(AlignmentCandidate<DNASequence,FASTQSequence> &alignment, ostream &outFile) {
      int mapQV = alignment.mapQV;
      int lastBlock = alignment.blocks.size()-1;
			outFile << alignment.qName << " " << alignment.tName << " " << alignment.score << " " 
							<< alignment.pctSimilarity << " " 
							<< alignment.qStrand << " " 
							<< alignment.qAlignedSeqPos + alignment.blocks[0].qPos << " "
							<< (alignment.qAlignedSeqPos 
									+ alignment.blocks[lastBlock].qPos 
									+ alignment.blocks[lastBlock].length) << " "
							<< alignment.qLength << " "
							<< alignment.tStrand << " "
							<< alignment.tAlignedSeqPos + alignment.blocks[0].tPos  << " "
							<< alignment.tAlignedSeqPos + alignment.blocks[lastBlock].tPos + alignment.blocks[lastBlock].length << " "
							<< alignment.tLength << " " << mapQV << " " << alignment.nCells << " " << alignment.clusterScore 
              << " " << alignment.probScore << " " << alignment.numSignificantClusters << endl;

  }
  static void PrintHeader(ostream &out) {
    out << "qname tname score pctsimilarity qstrand qstart qend qseqlength tstrand tstart tend tseqlength mapqv ncells clusterScore probscore numSigClusters" << endl;
  }
};


#endif
