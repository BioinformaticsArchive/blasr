#ifndef SDP_ALIGN_H_
#define SDP_ALIGN_H_

#include <math.h>
#include "SWAlign.h"
#include "AlignmentUtils.h"
#include "DistanceMatrixScoreFunction.h"
#include "sdp/SparseDynamicProgramming.h"
#include "sdp/SDPFragment.h"
#include "../../tuples/TupleList.h"
#include "../../tuples/DNATuple.h"
#include "../../tuples/TupleMatching.h"
#include "../../tuples/TupleList.h"
#include "../../datastructures/alignment/Path.h"
#include "../../datastructures/alignment/Alignment.h"
#include "../../datastructures/alignment/AlignmentGapList.h"

#define SDP_DETAILED_WORD_SIZE 5
#define SDP_PREFIX_LENGTH 50
#define SDP_SUFFIX_LENGTH 50

template<typename T_QuerySequence, typename T_TargetSequence, typename T_ScoreFn>
int SDPAlign(T_QuerySequence &query, T_TargetSequence &target,
             T_ScoreFn &scoreFn, int wordSize, 
						 int sdpIns, int sdpDel, float indelRate,
						 Alignment &alignment, 
						 AlignmentType alignType=Global,
						 bool detailedAlignment=true,
						 bool extendFrontByLocalAlignment=true) {

  /*
    Since SDP Align uses a large list of buffers, but none are
    provided with this mechanism of calling SDPAlign, allocate the
    buffers on the stack.
  */
	vector<Fragment> fragmentSet, prefixFragmentSet, suffixFragmentSet;
	TupleList<PositionDNATuple> targetTupleList;
	TupleList<PositionDNATuple> targetPrefixTupleList;
	TupleList<PositionDNATuple> targetSuffixTupleList;
	std::vector<int> maxFragmentChain;
  
  return SDPAlign(query, target,
                  scoreFn, wordSize, 
                  sdpIns, sdpDel, indelRate,
                  alignment, 
                  fragmentSet, prefixFragmentSet, suffixFragmentSet, 
                  targetTupleList, targetPrefixTupleList, targetSuffixTupleList,
                  maxFragmentChain,
                  alignType,
                  detailedAlignment,
                  extendFrontByLocalAlignment);
}

template<typename T_QuerySequence, typename T_TargetSequence, typename T_ScoreFn, typename T_BufferCache>
int SDPAlign(T_QuerySequence &query, T_TargetSequence &target,
             T_ScoreFn &scoreFn, int wordSize, 
						 int sdpIns, int sdpDel, float indelRate,
						 Alignment &alignment, 
             T_BufferCache &buffers,
						 AlignmentType alignType=Global,
						 bool detailedAlignment=true,
						 bool extendFrontByLocalAlignment=true, 
             int  noRecurseUnder = 10000) {

  return SDPAlign(query, target, scoreFn, wordSize, 
                  sdpIns, sdpDel, indelRate,
                  alignment,  
                  buffers.sdpFragmentSet,
                  buffers.sdpPrefixFragmentSet,
                  buffers.sdpSuffixFragmentSet,
                  buffers.sdpCachedTargetTupleList,
                  buffers.sdpCachedTargetPrefixTupleList,
                  buffers.sdpCachedTargetSuffixTupleList,
                  buffers.sdpCachedMaxFragmentChain,
                  alignType, detailedAlignment, extendFrontByLocalAlignment, noRecurseUnder);
}

template<typename T_QuerySequence, typename T_TargetSequence, typename T_ScoreFn, typename T_TupleList>
int SDPAlign(T_QuerySequence &query, T_TargetSequence &target,
             T_ScoreFn &scoreFn,
             int wordSize, 
						 int sdpIns, int sdpDel, float indelRate,
						 Alignment &alignment, 
						 vector<Fragment> &fragmentSet,
						 vector<Fragment> &prefixFragmentSet,
						 vector<Fragment> &suffixFragmentSet,
						 T_TupleList &targetTupleList,
						 T_TupleList &targetPrefixTupleList,
						 T_TupleList &targetSuffixTupleList,
						 std::vector<int> &maxFragmentChain,
             // A few optinal parameters, should delete that last one.
						 AlignmentType alignType=Global,
						 bool detailedAlignment=true,
						 bool extendFrontByLocalAlignment=true, 
             int  noRecurseUnder=10000) {

  fragmentSet.clear();
  prefixFragmentSet.clear();
  suffixFragmentSet.clear();
	targetTupleList.clear();
	targetPrefixTupleList.clear();
	targetSuffixTupleList.clear();
	maxFragmentChain.clear();

	/*
		Collect a set of matching fragments between query and target.
		Since this function is an inner-loop for alignment, anything to
		speed it up will help.  One way to speed it up is to re-use the
		vectors that contain the sdp matches. 
	*/

	TupleMetrics tm, tmSmall;
	tm.Initialize(wordSize);

  int smallWordSize = (wordSize < SDP_DETAILED_WORD_SIZE ? wordSize : SDP_DETAILED_WORD_SIZE);
  tmSmall.Initialize(smallWordSize);




  /*
   * Partition the read into a prefix, middle, and suffix.  The prefix
   * and suffix are matched using a smaller word size allowing for
   * higher sensitivity at the ends of reads, which are more likely to
   * be misaligned.
   */
  int prefixLength, middleLength, suffixLength, middlePos, suffixPos; // prefix pos is 0
  prefixLength = min(target.length, (DNALength) SDP_PREFIX_LENGTH);
  suffixLength = min(target.length - prefixLength, (DNALength) SDP_SUFFIX_LENGTH);
  middleLength = target.length - prefixLength - suffixLength;

  DNASequence prefix, middle, suffix;
  DNALength pos = 0;
  prefix.seq = &target.seq[pos];
  prefix.length = prefixLength;
  pos += prefixLength;
  middlePos = pos;
  middle.seq = &target.seq[middlePos];
  middle.length = middleLength;
  pos += middleLength;
  suffixPos = pos;
  suffix.seq = &target.seq[suffixPos];
  suffix.length = suffixLength;
  
	fragmentSet.clear();
  SequenceToTupleList(prefix, tmSmall, targetPrefixTupleList);
  SequenceToTupleList(suffix, tmSmall, targetSuffixTupleList);
	SequenceToTupleList(middle, tm, targetTupleList);

  targetPrefixTupleList.Sort();
  targetSuffixTupleList.Sort();
	targetTupleList.Sort();


  //
  // Store in fragmentSet the tuples that match between the target
  // and query.
  //
  StoreMatchingPositions(query, tmSmall, targetPrefixTupleList, prefixFragmentSet);
  StoreMatchingPositions(query, tmSmall, targetSuffixTupleList, suffixFragmentSet);
	StoreMatchingPositions(query, tm, targetTupleList, fragmentSet); 

  
  // 
  // The method to store matching positions is not weight aware.
  // Store the weight here.
  //
	VectorIndex f;
  
	for (f = 0; f < suffixFragmentSet.size(); f++) {
    (suffixFragmentSet)[f].weight = tmSmall.tupleSize;
  }
	for (f = 0; f < prefixFragmentSet.size(); f++) {
    (prefixFragmentSet)[f].weight = tmSmall.tupleSize;
  }
	for (f = 0; f < fragmentSet.size(); f++) {
		(fragmentSet)[f].weight = tm.tupleSize;
  }


  //
  // Since different partitions of the read are matched, the locations
  // of the matches do not have the correct position because of the
  // offsets.  Fix that here.

	for (f = 0; f < fragmentSet.size(); f++) {
    (fragmentSet)[f].y += middlePos;
	}
	for (f = 0; f < suffixFragmentSet.size(); f++) {
    (suffixFragmentSet)[f].y += suffixPos;
  }
  
  //
  // Collect all fragments into one.
  //
  fragmentSet.insert(fragmentSet.begin(), prefixFragmentSet.begin(), prefixFragmentSet.end());
  fragmentSet.insert(fragmentSet.end(), suffixFragmentSet.begin(), suffixFragmentSet.end());

	if (fragmentSet.size() == 0) {
		//
		// This requires at least one seeded tuple to begin an alignment.
		//
		return 0;
	}

  //
  // Find the longest chain of anchors.
  //
  
	SDPLongestCommonSubsequence(query.length, fragmentSet, tm.tupleSize, sdpIns, sdpDel, scoreFn.scoreMatrix[0][0], maxFragmentChain, alignType);

	//
	// Now turn the max fragment chain into real a real alignment.
	//
		
	int startF;
	Alignment chainAlignment;
	alignment.qPos = 0;
	alignment.tPos = 0;
	Block block;
	vector<int> fragScoreMat;
	vector<Arrow> fragPathMat;

  //
  // Patch the sdp fragments into an alignment, possibly breaking the
  // alignment if the gap between two fragments is too large.
  //

	for (f = 0; f < maxFragmentChain.size(); f++ ){ 
		startF = f;
		// Condense contiguous stretches.
		while(f < maxFragmentChain.size()  - 1 and 
					fragmentSet[maxFragmentChain[f]].x == fragmentSet[maxFragmentChain[f+1]].x - 1 and
					fragmentSet[maxFragmentChain[f]].y == fragmentSet[maxFragmentChain[f+1]].y - 1) {
			f++;
    }
			
		block.qPos = fragmentSet[maxFragmentChain[startF]].x;
		block.tPos = fragmentSet[maxFragmentChain[startF]].y;
		block.length = fragmentSet[maxFragmentChain[startF]].weight + (f - startF);

		chainAlignment.blocks.push_back(block);
	}


	//
	// It may be possible that in regions of low similarity, spurious matches fit into the LCS.  
	// Assume that indels cause the matches to diverge from the diagonal on a random walk.  If they 
	// walk more than 3 standard deviations away from the diagonal, they are probably spurious. 
	//
	unsigned int b;
	chainAlignment.qPos = 0;
	chainAlignment.tPos = 0;



	for (b = 0; b < chainAlignment.size()-1; b++){ 
		if (chainAlignment.blocks[b].qPos + chainAlignment.blocks[b].length > chainAlignment.blocks[b+1].qPos) {
			chainAlignment.blocks[b].length = (chainAlignment.blocks[b+1].qPos - chainAlignment.blocks[b].qPos);
		}
		if (chainAlignment.blocks[b].tPos + chainAlignment.blocks[b].length > chainAlignment.blocks[b+1].tPos) {
			chainAlignment.blocks[b].length = (chainAlignment.blocks[b+1].tPos - chainAlignment.blocks[b].tPos);
		}
		// the min indel rate between the two chain blocks is the difference in diagonals between the two sequences.
		int curDiag, nextDiag, diffDiag;
		curDiag = chainAlignment.blocks[b].tPos - chainAlignment.blocks[b].qPos;
		nextDiag = chainAlignment.blocks[b+1].tPos - chainAlignment.blocks[b+1].qPos;
		diffDiag = abs(curDiag - nextDiag);

		//
		// It is expected that the deviation is at least 1, so discount for this
		//
		diffDiag--;
		// compare the alignment distances.  
	}

  vector<bool> blockIsGood;
  blockIsGood.resize(chainAlignment.size());
  fill(blockIsGood.begin(), blockIsGood.end(), true);

  /*
   * The hack that allows anchors of different lengths at the front
   * and end of alignments (to increase sensitivity at the ends of
   * sequences) has the side effect that there may be blocks that have
   * zero length.  This shouldn't happen, so to balance this out
   * remove blocks that have zero length.
   */
   
  bool badBlock;
	for (b = 0; b < chainAlignment.size(); b++){ 
    if (chainAlignment.blocks[b].length == 0) {
      blockIsGood[b] = false;
    }
  }
	for (b = 1; b < chainAlignment.size()-1; b++){ 
		// the min indel rate between the two chain blocks is the difference in diagonals between the two sequences.
		int prevDiag = abs(((int)chainAlignment.blocks[b].tPos -   (int)chainAlignment.blocks[b].qPos)  -
                       ((int)chainAlignment.blocks[b-1].tPos - (int)chainAlignment.blocks[b-1].qPos));

    int prevDist = min(chainAlignment.blocks[b].tPos - chainAlignment.blocks[b-1].tPos,
                       chainAlignment.blocks[b].qPos - chainAlignment.blocks[b-1].qPos);

		int nextDiag = abs(((int)chainAlignment.blocks[b+1].tPos - (int)chainAlignment.blocks[b+1].qPos)  -
                       ((int)chainAlignment.blocks[b].tPos -   (int)chainAlignment.blocks[b].qPos));
		
    int nextDist = min(chainAlignment.blocks[b+1].tPos - chainAlignment.blocks[b].tPos,
                       chainAlignment.blocks[b+1].qPos - chainAlignment.blocks[b].qPos);

    if (prevDist * indelRate < prevDiag and nextDist * indelRate < nextDiag) {
      blockIsGood[b] = false;
    }
  }
  for (b = chainAlignment.size(); b > 0; b--) {
    if (blockIsGood[b-1] == false) {
      chainAlignment.blocks.erase(chainAlignment.blocks.begin() + b-1);
    }
  }

	if (chainAlignment.blocks.size() > 0) {
		T_QuerySequence  qFragment;
    T_TargetSequence tFragment;
		Alignment fragAlignment;			
		unsigned int fb;
		if (alignType == Global) {
			//
			// For Global alignment, refine the alignment from the beginnings of the
			// sequences to the start of the first block.
			//
			if (chainAlignment.blocks[0].qPos > 0 and
					chainAlignment.blocks[0].tPos > 0) {
				
				qFragment.seq = &query.seq[0];
				qFragment.length = chainAlignment.blocks[0].qPos;
				tFragment.seq = &target.seq[0];
				tFragment.length = chainAlignment.blocks[0].tPos;
				for (fb = 0; fb < alignment.blocks.size(); fb++) {
					alignment.blocks.push_back(fragAlignment.blocks[b]);
				}
			}
		}
		else if (alignType == Local) {
			// Perform a front-anchored alignment to extend the alignment to
			// the beginning of the read.
			if (chainAlignment.blocks[0].qPos > 0 and 
					chainAlignment.blocks[0].tPos > 0) {
				qFragment.seq = (Nucleotide*) &query.seq[0];
				qFragment.length = chainAlignment.blocks[0].qPos;

				tFragment.seq = (Nucleotide*) &target.seq[0];
				tFragment.length = chainAlignment.blocks[0].tPos;
				Alignment frontAlignment;
				int frontAlignmentScore;
				// Currently, there might be some space between the beginning
				// of the alignment and the beginning of the read.  Run an
				// EndAnchored alignment that allows free gaps to the start of
				// where the alignment begins, but normal, ungapped alignment
				// otherwise. 
				if (extendFrontByLocalAlignment) {
          if (noRecurseUnder == 0 or qFragment.length * tFragment.length < noRecurseUnder) {
            frontAlignmentScore  = 
              SWAlign(qFragment, tFragment, fragScoreMat, fragPathMat, frontAlignment, scoreFn, EndAnchored);
          }
          else {
            //            cout << "running recursive sdp alignment. " << endl;
            vector<int> recurseFragmentChain;
            SDPAlign(qFragment, tFragment, scoreFn,
                     max(wordSize/2, 5),
                     sdpIns, sdpDel,  indelRate,
                     frontAlignment,
                     fragmentSet,
                     prefixFragmentSet,
                     suffixFragmentSet,
                     targetTupleList,
                     targetPrefixTupleList,
                     targetSuffixTupleList,
                     recurseFragmentChain,
                     alignType, detailedAlignment, extendFrontByLocalAlignment, 0);
          }
					
					int anchorBlock;
					for (anchorBlock = 0; anchorBlock < frontAlignment.blocks.size(); anchorBlock++) {
						//
						// The front alignment needs to be transformed to the
						// coordinate offsets that the chain alignment is in.  This
						// is an alignment starting at position 0 in the target and
						// query.  Currently, the front alignment is offset into the
						// sequences by frontAlignment.[q/t]Pos.
						//
						frontAlignment.blocks[anchorBlock].tPos += frontAlignment.tPos;
						frontAlignment.blocks[anchorBlock].qPos += frontAlignment.qPos;
						alignment.blocks.push_back(frontAlignment.blocks[anchorBlock]);
					}
				}
			}
		}
		
		// 
		// The chain alignment blocks are not complete blocks, so they
		// must be appended to the true alignment and then patched up.
		//

		for (b = 0; b < chainAlignment.size() - 1; b++) {
			alignment.blocks.push_back(chainAlignment.blocks[b]);
			int alignScore;
      
      //
      // Do a detaied smith-waterman alignment between blocks, if this
      // is specified.  
      fragAlignment.Clear();
      qFragment.ReferenceSubstring(query, chainAlignment.blocks[b].qPos + chainAlignment.blocks[b].length);
      qFragment.length = chainAlignment.blocks[b+1].qPos - 
        (chainAlignment.blocks[b].qPos + chainAlignment.blocks[b].length);
                                   
			tFragment.seq    = &(target.seq[chainAlignment.blocks[b].tPos + chainAlignment.blocks[b].length]);
			tFragment.length = (chainAlignment.blocks[b+1].tPos - 
													(chainAlignment.blocks[b].tPos + chainAlignment.blocks[b].length));

			if (qFragment.length > 0 and 
					tFragment.length > 0 and
					detailedAlignment == true) {

        if (noRecurseUnder == 0 or qFragment.length * tFragment.length < noRecurseUnder) {
          alignScore = SWAlign(qFragment, tFragment, fragScoreMat, fragPathMat, fragAlignment, scoreFn, Global);
        }
        else {
          //          cout << "running recursive sdp alignment on " << qFragment.length * tFragment.length << endl;
          vector<int> recurseFragmentChain;
          SDPAlign(qFragment, tFragment, scoreFn,
                   max(wordSize/2, 5),
                   sdpIns, sdpDel,  indelRate,
                   fragAlignment,
                   fragmentSet,
                   prefixFragmentSet,
                   suffixFragmentSet,
                   targetTupleList,
                   targetPrefixTupleList,
                   targetSuffixTupleList,
                   recurseFragmentChain,
                   alignType, detailedAlignment, 0, 0);
        }
        /*
        if (noRecurseUnder and qFragment.length * tFragment.length > noRecurseUnder) {
          cout << "uh oh " << qFragment.length << " " << tFragment.length << " " << (1.0*qFragment.length) / tFragment.length << endl; 
          StickPrintAlignment(fragAlignment, qFragment, tFragment, cout);
        }
        */
				fragAlignment.qPos = 0;
				fragAlignment.tPos = 0;

				int qOffset = chainAlignment.blocks[b].qPos + chainAlignment.blocks[b].length;
				int tOffset = chainAlignment.blocks[b].tPos + chainAlignment.blocks[b].length;

				for (fb = 0; fb < fragAlignment.blocks.size(); fb++) {
					fragAlignment.blocks[fb].qPos += qOffset;
					fragAlignment.blocks[fb].tPos += tOffset;
					alignment.blocks.push_back(fragAlignment.blocks[fb]);
				}
			}
		}
		int lastBlock = chainAlignment.blocks.size() - 1;
		if (alignType == Global or alignType == Local) {
			if (chainAlignment.size() > 0) {
				// Add the last block.
				alignment.blocks.push_back(chainAlignment.blocks[lastBlock]);
				if (alignType == Global and detailedAlignment == true) {
					//
					// When doing a global alignment, the sequence from the end of
					// the last block of the query should be aligned to the end of 
					// the text.
					//
					
          qFragment.ReferenceSubstring(query, chainAlignment.blocks[lastBlock].qPos + chainAlignment.blocks[lastBlock].length,
                                      query.length -  
                                      (chainAlignment.blocks[lastBlock].qPos + chainAlignment.blocks[lastBlock].length));

					tFragment.seq    = &(target.seq[chainAlignment.blocks[lastBlock].tPos + chainAlignment.blocks[lastBlock].length]);
					tFragment.length = (target.length - 
														(chainAlignment.blocks[lastBlock].tPos + chainAlignment.blocks[lastBlock].length));
					if (qFragment.length > 0 and
							tFragment.length > 0 ) {

            if (extendFrontByLocalAlignment) {

              fragAlignment.Clear();
              if (qFragment.length * tFragment.length > 10000) {
                //              cout << "Cautin: slow alignment crossing! " << qFragment.length  << " " << tFragment.length << endl;
              }

              if (noRecurseUnder == 0 or qFragment.length * tFragment.length < noRecurseUnder) {
                SWAlign(qFragment, tFragment, fragScoreMat, fragPathMat, fragAlignment, scoreFn, EndAnchored);
              }
              else {
                vector<int> recurseFragmentChain;
                SDPAlign(qFragment, tFragment, scoreFn,
                         max(wordSize/2, 5),
                         sdpIns, sdpDel,  indelRate,
                         fragAlignment,
                         fragmentSet,
                         prefixFragmentSet,
                         suffixFragmentSet,
                         targetTupleList,
                         targetPrefixTupleList,
                         targetSuffixTupleList,
                         recurseFragmentChain,
                         alignType, detailedAlignment, extendFrontByLocalAlignment, 0);
              }



              int qOffset = chainAlignment.blocks[lastBlock].qPos + chainAlignment.blocks[lastBlock].length;
              int tOffset = chainAlignment.blocks[lastBlock].tPos + chainAlignment.blocks[lastBlock].length;
              unsigned int fb;
              for (fb = 0; fb < fragAlignment.size(); fb++) { 
                fragAlignment.blocks[fb].qPos += qOffset;
                fragAlignment.blocks[fb].tPos += tOffset;
                alignment.blocks.push_back(fragAlignment.blocks[fb]);
              }
            }
					}
				}
			}
		}
	}

	if (alignType == Local) {
		alignment.tPos = alignment.blocks[0].tPos;
		alignment.qPos = alignment.blocks[0].qPos;
		VectorIndex b;
		for (b = 0; b < alignment.blocks.size(); b++) { 
			alignment.blocks[b].qPos -= alignment.qPos;
			alignment.blocks[b].tPos -= alignment.tPos;
		}
	}
	int alignmentScore;
	alignmentScore = ComputeAlignmentScore(alignment, query, target, 
                                         scoreFn.scoreMatrix, scoreFn.ins, scoreFn.del);
	return alignmentScore;
}

#endif
