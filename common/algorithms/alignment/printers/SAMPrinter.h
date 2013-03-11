#ifndef ALGORITHMS_ALIGNMENT_PRINTERS_SAMPRINTER_H_
#define ALGORITHMS_ALIGNMENT_PRINTERS_SAMPRINTER_H_

#include <sstream>
#include <stdint.h>
#include "SMRTSequence.h"
#include "datastructures/alignment/AlignmentCandidate.h"
#include "datastructures/alignment/AlignmentContext.h"


#define MULTI_SEGMENTS 0x1
#define ALL_SEGMENTS_ALIGNED 0x2
#define SEGMENT_UNMAPPED 0x4
#define NEXT_SEGMENT_UNMAPPED 0x8
#define SEQ_REVERSED 0x10
#define SEQ_NEXT_REVERSED 0x20
#define FIRST_SEGMENT 0x40
#define LAST_SEGMENT 0x80
#define SECONDARY_ALIGNMENT 0x100
#define NO_PASS_QUALITY 0x200
#define PCR_OR_OPTICAL_DUPLICATE 0x400

namespace SAMOutput {

  enum Clipping {hard, soft, none};

  void BuildFlag(T_AlignmentCandidate &alignment, AlignmentContext &context, uint16_t &flag) {

    /*
     *  Much of the flags are commented out for now since they do not
     *  generate PICARD compliant SAM.  This needs to be worked out. 
     */


    //
    // Without supporting strobe, assume 1 segment per template.
    flag = 0;
    /*
    if (context.nSubreads > 1) {
      flag |= MULTI_SEGMENTS;
      }*/

    //    if (context.AllSubreadsAligned() and context.nSubreads > 1) {
    //      flag |= ALL_SEGMENTS_ALIGNED;
      //    }

    if (alignment.tStrand == 1) {
      flag |= SEQ_REVERSED;
    }
    /*
    if (context.hasNextSubreadPos == false and context.nSubreads > 1) {
      flag |= NEXT_SEGMENT_UNMAPPED;
    }
    if (context.nextSubreadDir == 1) {
      flag |= SEQ_NEXT_REVERSED;
    }
    if (context.IsFirst() and context.nSubreads > 1) {
      flag |= FIRST_SEGMENT;
    }
    */
    else if (context.nSubreads > 1) {
      /*
       * Remember, if you're not first, you're last.
       */
      //      flag |= LAST_SEGMENT;
    }
    if (context.isPrimary == false) {
      flag |= SECONDARY_ALIGNMENT;
    }
  }

  
  // 
  // The aligned sequence is either the sequence from the first
  // aligned base to the last (hard and no clipping), or first high
  // quality base to the last high quality base (soft clipping).
  //
  template<typename T_Sequence>
  void SetAlignedSequence(T_AlignmentCandidate &alignment, T_Sequence &read,
                          T_Sequence &alignedSeq,
                          Clipping clipping = none) {
    DNALength qStart = alignment.QAlignStart();
    DNALength qEnd   = alignment.QAlignEnd();
    //
    // In both no, and hard clipping, the dna sequence that is output
    // solely corresponds to the aligned sequence.
    //
    DNALength clippedReadLength = 0;
    DNALength clippedStartPos   = 0;


    if (clipping == none or clipping == hard) {
      clippedReadLength = qEnd - qStart;
      clippedStartPos   = qStart;
      if (alignment.tStrand == 0) {
        alignedSeq.ReferenceSubstring(read, clippedStartPos, clippedReadLength);
      }
      else {
        T_Sequence subSeq;
        subSeq.ReferenceSubstring(read, clippedStartPos, clippedReadLength);
        subSeq.MakeRC(alignedSeq);
        alignedSeq.deleteOnExit = true;
      }
    }
    else if (clipping == soft) {
      clippedReadLength = read.length - read.lowQualityPrefix - read.lowQualitySuffix;
      if (alignment.tStrand == 0) {
        clippedStartPos = read.lowQualityPrefix;
        alignedSeq.ReferenceSubstring(read, read.lowQualityPrefix, clippedReadLength);
        alignedSeq.deleteOnExit = false;
      }
      else {
        //
        // If the alignment is a reverse complement, the soft clipped
        // string still excludes the low quality prefix and suffix,
        // but since the reverse string is printed, the prefix and
        // suffix are reversed.n
        //

        T_Sequence subSeq;
        subSeq.ReferenceSubstring(read, read.lowQualityPrefix, clippedReadLength);
        subSeq.MakeRC(alignedSeq);
      }
    }
    else {
      cout <<" ERROR! The clipping must be none, hard, or soft when setting the aligned sequence." << endl;
      assert(0);
    }
  }



 void CreateDNAString(DNASequence &seq, 
                      DNASequence &clippedSeq,
                      //
                      // Trimming is used for both hard non-clipping
                      // so it is called trim instead of clip.
                      //
                      int trimFront=0,
                      int trimEnd=0) {
   assert(seq.length - trimEnd >= trimFront);

   clippedSeq.seq    = &seq.seq[trimFront];
   clippedSeq.length = seq.length - trimEnd - trimFront;
 }

 void AddGaps(T_AlignmentCandidate &alignment, int gapIndex,
              vector<int> &opSize, vector<char> &opChar) {
   int g;
   for (g = 0; g < alignment.gaps[gapIndex].size(); g++) {
     if (alignment.gaps[gapIndex][g].seq == Gap::Query) {
       opSize.push_back(alignment.gaps[gapIndex][g].length);
       opChar.push_back('D');
     }
     else if (alignment.gaps[gapIndex][g].seq == Gap::Target) {
       opSize.push_back(alignment.gaps[gapIndex][g].length);
       opChar.push_back('I');
     }
   }
 }

 void CreateNoClippingCigarOps(T_AlignmentCandidate &alignment, 
                                vector<int> &opSize, 
                                vector<char> &opChar) {
    //
    // Create the cigar string for the aligned region of a read,
    // excluding the clipping.
    //
    int b;
    // Each block creates a match NM (N=block.length)
    int nBlocks = alignment.blocks.size();
    int nGaps   = alignment.gaps.size();
    opSize.clear();
    opChar.clear();
    //
    // Add gaps at the beginning of the alignment.
    //
    if (nGaps > 0) {
      AddGaps(alignment, 0, opSize, opChar);
    }
    for (b = 0; b < nBlocks; b++) {
      //
      // Determine the gap before printing the match, since it is
      // possible that the qurey and target are gapped at the same
      // time, which merges into a mismatch.
      //
      int qGap=0, tGap=0, commonGap=0;
      int matchLength = alignment.blocks[b].length;
      if (nGaps == 0) {
        if (b + 1 < nBlocks) {
          qGap = alignment.blocks[b+1].qPos - alignment.blocks[b].qPos - alignment.blocks[b].length;
          tGap = alignment.blocks[b+1].tPos - alignment.blocks[b].tPos - alignment.blocks[b].length;
          int commonGap;
          commonGap = abs(qGap - tGap);
          qGap -= commonGap;
          tGap -= commonGap;
          matchLength += commonGap;
          opSize.push_back(matchLength);
          opChar.push_back('M');
          assert((qGap > 0 and tGap == 0) or (qGap == 0 and tGap > 0));
          if (qGap > 0) {
            opSize.push_back(qGap);
            opChar.push_back('I');
          }
          if (tGap > 0) {
            opSize.push_back(tGap);
            opChar.push_back('D');
          }
        }
      }
      else {
        opSize.push_back(matchLength);
        opChar.push_back('M');
        int g;
        int gapIndex = b+1;
        AddGaps(alignment, gapIndex, opSize, opChar);
      }
    }
    if (alignment.tStrand == 1) {
      reverse(opSize.begin(), opSize.end());
      reverse(opChar.begin(), opChar.end());
    }
  }

 template<typename T_Sequence>
  void SetSoftClip(T_AlignmentCandidate &alignment,
                   T_Sequence &read,
                   DNALength &softClipPrefix, 
                   DNALength &softClipSuffix) {
    DNALength qStart, qEnd;
    int hardClipPrefix = read.lowQualityPrefix;
    int hardClipSuffix = read.lowQualitySuffix;

    qStart = alignment.QAlignStart();
    qEnd   = alignment.QAlignEnd();

    assert(qStart >= hardClipPrefix);
    softClipPrefix = alignment.QAlignStart() - hardClipPrefix;
    assert(alignment.QAlignEnd() + hardClipSuffix <= read.length);
    softClipSuffix = read.length - hardClipSuffix - alignment.QAlignEnd();
  }
 
 template<typename T_Sequence>
  void SetHardClip(T_AlignmentCandidate &alignment, 
                   T_Sequence &read,
                   int &prefixClip,
                   int &suffixClip) {
    //
    // Set the hard clipping assuming the read is in the forward
    // direction.
    //
    prefixClip = alignment.QAlignStart();
    suffixClip = read.length - alignment.QAlignEnd();
    if (alignment.tStrand == 1) {
      //
      // If the read is instead reverse, swap the clipping lengths.
      //
      swap(prefixClip, suffixClip);
    }
  }

  void CigarOpsToString(vector<int> &opSize,
                        vector<char> &opChar,
                        string &cigarString) {
    stringstream sstrm;
    int i, nElem;
    for (i = 0, nElem = opSize.size(); i < nElem; i++) {
      sstrm << opSize[i] << opChar[i];
    }
    cigarString = sstrm.str();
  }

  //
  // Straight forward: create the cigar string allowing some clipping
  // The read is provided to give length and hq information.
  //
  template<typename T_Sequence>
  void CreateCIGARString(T_AlignmentCandidate &alignment,
                         T_Sequence &read,
                         string &cigarString,
                         Clipping clipping=none) {
    cigarString = "";
    // All cigarString use the no clipping core
    vector<int> opSize;
    vector<char> opChar;
    CreateNoClippingCigarOps(alignment, opSize, opChar);

    // Clipping needs to be added

    if (clipping == hard) {
      int prefixClip, suffixClip;
      SetHardClip(alignment, read, prefixClip, suffixClip);
      if (prefixClip > 0) {
        opSize.insert(opSize.begin(), prefixClip);
        opChar.insert(opChar.begin(), 'H');
      }
      if (suffixClip > 0) {
        opSize.push_back(suffixClip);
        opChar.push_back('H');
      }
    }
    if (clipping == soft) {
      //
      // Even if clipping is soft, the hard clipping removes the 
      // low quality regions
      //
      DNALength prefixSoftClip, suffixSoftClip;
      DNALength prefixHardClip, suffixHardClip;
      prefixHardClip = read.lowQualityPrefix;
      suffixHardClip = read.lowQualitySuffix;

      SetSoftClip(alignment, read, prefixSoftClip, suffixSoftClip);



      if (alignment.tStrand == 1) {
        swap(prefixHardClip, suffixHardClip);
        swap(prefixSoftClip, suffixSoftClip);
      }

      //
      // Insert the hard and soft clipping so that they are in the
      // order H then S if both exist.
      //
      if (prefixSoftClip > 0) {
        opSize.insert(opSize.begin(), prefixSoftClip);
        opChar.insert(opChar.begin(), 'S');
      }
      if (prefixHardClip > 0) {
        opSize.insert(opSize.begin(), prefixHardClip);
        opChar.insert(opChar.begin(), 'H');
      }
      
      //
      // Append the hard and soft clipping so they are in the order S
      // then H. 
      //
      if (suffixSoftClip > 0) {
        opSize.push_back(suffixSoftClip);
        opChar.push_back('S');
      }
      if (suffixHardClip > 0) {
        opSize.push_back(suffixHardClip);
        opChar.push_back('H');
      }
    }

    CigarOpsToString(opSize, opChar, cigarString);

  }
  template<typename T_Sequence>
  void PrintAlignment(T_AlignmentCandidate &alignment,
                      T_Sequence &read,
                      ostream &samFile,
                      AlignmentContext &context,
                      Clipping clipping = none,
                      int subreadIndex = 0,
                      int nSubreads = 0) {
    
    string cigarString;
    uint16_t flag;
    T_Sequence alignedSequence;
    CreateCIGARString(alignment, read, cigarString, clipping);
    SetAlignedSequence(alignment, read, alignedSequence, clipping);
    BuildFlag(alignment, context, flag);
    samFile << alignment.qName << "\t" 
            << flag << "\t" 
            << alignment.tName << "\t";   // RNAME
    if (alignment.tStrand == 0) {
      samFile << alignment.TAlignStart() + 1 << "\t"; // POS, add 1 to
                                             // get 1 based
                                             // coordinate
                                             // system   
    }
    else {
      samFile << alignment.tLength - (alignment.TAlignStart() + alignment.TEnd()) + 1 << "\t"; // includes - 1 for rev-comp,  +1 for one-based
    }
    samFile << (int) alignment.mapQV << "\t"// MAPQ
            << cigarString << "\t"; // CIGAR
      
      //
      // Determine RNEXT

    string rNext;
    rNext = "*";
    /*
    if (context.hasNextSubreadPos == false) {
      rNext = "*";
    }
    else {
      if (context.rNext == alignment.tName) {
        rNext = "=";
      }
      else {
        rNext = context.rNext;
      }
    }
    */
    samFile << rNext << "\t"; // RNEXT
    
    DNALength nextSubreadPos = 0;
    /*
    if (context.hasNextSubreadPos) {
      nextSubreadPos = context.nextSubreadPos + 1;
      }*/
    samFile << nextSubreadPos << "\t"; // RNEXT, add 1 for 1 based
                                           // indexing


    DNALength tLen = alignment.GenomicTEnd() - alignment.GenomicTBegin();
    samFile << tLen << "\t"; // TLEN
    // Print the sequence on one line, and suppress printing the
    // newline (by setting the line length to alignedSequence.length
    ((DNASequence)alignedSequence).PrintSeq(samFile, 0);  // SEQ
    samFile << "\t";
    if (alignedSequence.qual.data != NULL) {
      alignedSequence.PrintAsciiQual(samFile, 0);  // QUAL
    }
    else {
      samFile <<"*";
    }
    samFile << "\t";
    //
    // Add optional fields
    //
    samFile << "RG:Z:" << context.readGroupId << "\t";
    samFile << "AS:i:" << alignment.score << "\t";

    DNALength prefixSoftClip=0, suffixSoftClip=0;
    DNALength prefixHardClip=0, suffixHardClip=0;
    if (clipping == soft) {
      SetSoftClip(alignment, read, prefixSoftClip, suffixSoftClip);
    }
    samFile << "XS:i:" << alignment.QAlignStart() + 1 - prefixSoftClip << "\t";
    samFile << "XE:i:" << alignment.QAlignEnd() + 1 - prefixSoftClip << "\t";
    samFile << "XL:i:" << alignment.qAlignedSeq.length << "\t";
    samFile << "XT:i:1\t"; // reads are allways continuous reads, not
                        // referenced based circular consensus when
                        // output by blasr.
    samFile << "NM:i:" << context.nSubreads << "\t";
    samFile << "FI:i:" << alignment.qAlignedSeqPos + 1;

    samFile << endl;

    alignedSequence.FreeIfControlled();
  }
};




#endif
