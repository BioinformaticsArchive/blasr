#ifndef SPARSE_DYNAMIC_PROGRAMMING_H_
#define SPARSE_DYNAMIC_PROGRAMMING_H_

#include <stdlib.h>
#include <assert.h>
#include <set>
#include <limits.h>
#include "SDPSet.h"
#include "SDPFragment.h"
#include "SDPColumn.h"
#include "FragmentSort.h"
#include "../AlignmentUtils.h"
#include "../../../defs.h"


/*******************************************************************************
 *  Sparse dynamic programming implementation of Longest Common Subsequence
 *  
 *  Implementation of method described in Baker and Giancarlo, Journal of
 *  Algorithms 42, 231-254, 2002.
 * 
 *  5/7/09 -- Modified to incorporate different linear cost functions, and 
 *  local alignments.
 *  
 ******************************************************************************/

int IndelPenalty(int x1, int y1, int x2, int y2, int insertion, int deletion) {
  int drift, driftPenalty;
  drift = (x1 - y1) - (x2 - y2);
  if (drift > 0) {
    driftPenalty = drift * insertion;
  }
  else {
    driftPenalty = -drift * deletion;
  }
  return driftPenalty;
}


template<typename T_Fragment>
int SDPLongestCommonSubsequence(DNALength queryLength,
																vector<T_Fragment> &fragmentSet, 
																DNALength fragmentLength,
																int insertion, int deletion, int match,
																vector<int> &maxFragmentChain, AlignmentType alignType=Global) {
	maxFragmentChain.clear();

	if (fragmentSet.size() < 1)
		return 0;
		
	std::sort(fragmentSet.begin(), fragmentSet.end(), LexicographicFragmentSort<T_Fragment>());

	SDPSet<Fragment> sweepSet;
	SDPSet<SDPColumn>   colSet;

	unsigned int sweepRow;
	unsigned int trailRow;
	VectorIndex fSweep, fTrail;
	VectorIndex fi;
	for (fi = 0; fi < fragmentSet.size(); fi++) {
		fragmentSet[fi].index = fi;
	}

	sweepRow = fragmentSet[0].x;

	
	Fragment pred, succ;
	fSweep = 0;
	fTrail = 0;
	unsigned int maxChainLength = 0;
	int maxChainFragment = -1;
	int minFragmentCost, minFragmentIndex;
	minFragmentCost = INF_INT;
	minFragmentIndex = -1;
	for (sweepRow = 0; sweepRow < queryLength + fragmentLength; sweepRow++) {
		//
		// Add all elements on the sweep row to the sweep set.  Note that when
		// fSweep is past query.length.
		int startF = fSweep;
    int fragmentSetSize = fragmentSet.size();
		while (fSweep < fragmentSetSize and 
					 fragmentSet[fSweep].x == sweepRow) {

			//
			// Compute the cost of every fragment in the sweep.
			//
			int cp = INF_INT, cl = INF_INT, ca = INF_INT;
			SDPColumn curCol, predCol;
			curCol.col = fragmentSet[fSweep].y;
			//
			// Search preceeding fragments.
			//

			//
			// Compute the cost of fragment_f
			int foundPrev = 0;
      int drift, driftPenalty;
			if (colSet.Predecessor(curCol, predCol)) {
				//
				// predCol points to the fragment with greatest value less than curCol.
				// 
				/*
					// Baker and Giancarlo LCS cost
				cp = fragmentSet[fSweep].x + fragmentSet[fSweep].y +
					fragmentSet[predCol.optFragment].cost - 
					fragmentSet[predCol.optFragment].x - 
					fragmentSet[predCol.optFragment].y - 2 * Fragment::length;
				*/
        driftPenalty = IndelPenalty(fragmentSet[fSweep].x, fragmentSet[fSweep].y,
                                           fragmentSet[predCol.optFragment].x, fragmentSet[predCol.optFragment].y,
                                           insertion, deletion);
				cp = fragmentSet[predCol.optFragment].cost + driftPenalty;
				foundPrev = 1;
			}

			// Search overlapping fragments.
			if (sweepSet.Predecessor(fragmentSet[fSweep], pred)) {
				/*
					Baker and Giancarlo LCS cost
				cl = pred.cost + 
					(fragmentSet[fSweep].y - fragmentSet[fSweep].x) - 
					(pred.y - pred.x);
				*/
        /*
         * Cost with insertion and deletion penalty.
         */
				cl = (pred.cost + 
              MIN((int)(fragmentLength - (fragmentSet[fSweep].y - pred.y)) * match, 0) + 
              IndelPenalty(fragmentSet[fSweep].x, fragmentSet[fSweep].y,
                           pred.x, pred.y, insertion, deletion));
				foundPrev = 1;
			}
			
			if (sweepSet.Successor(fragmentSet[fSweep], succ) and
					succ.y < fragmentSet[fSweep].y and
					succ.y + fragmentLength >= fragmentSet[fSweep].y) {
				assert(((int)succ.y) - ((int)succ.x) >= ((int)fragmentSet[fSweep].y) - ((int)fragmentSet[fSweep].x));

				/*
					Baker and Giancarlo LCS cost 
				ca = succ.cost +
					(succ.y - succ.x) -
					(fragmentSet[fSweep].y - fragmentSet[fSweep].x);
				*/
				ca = (succ.cost + 
              (fragmentLength - (int)(fragmentSet[fSweep].y - succ.y)) * match + 
              IndelPenalty(fragmentSet[fSweep].x, fragmentSet[fSweep].y, succ.x, succ.y, insertion, deletion));
				foundPrev = 1;
			}
			
			//
			// Now compute the minimum of all these.
			//
			int minCost;
			minCost = MIN(cp, MIN(cl, ca));

			//
			//  If doing a global alignment, chain is always extended.  If local, the chain may not be.
 			// 
			if (foundPrev and 
					(alignType == Global or
					 (alignType == Local and minCost < 0))) {
				fragmentSet[fSweep].cost = minCost - fragmentSet[fSweep].weight;
				if (minCost == cp) {
					fragmentSet[fSweep].chainPrev = predCol.optFragment;
				}
				else if (minCost == cl) {
					fragmentSet[fSweep].chainPrev = pred.index;
				}
				else if (minCost == ca) {
					fragmentSet[fSweep].chainPrev = succ.index;
				}
				assert(fragmentSet[fSweep].chainPrev >= 0 and
							 fragmentSet[fSweep].chainPrev < (int)fragmentSet.size());
				fragmentSet[fSweep].chainLength = fragmentSet[fragmentSet[fSweep].chainPrev].chainLength + 1;
			}
			else {
				if (alignType == Global) {
					fragmentSet[fSweep].chainPrev = (int) -1;
					fragmentSet[fSweep].cost = (fragmentSet[fSweep].x + fragmentSet[fSweep].y) * deletion + fragmentLength * match - fragmentSet[fSweep].weight;
					fragmentSet[fSweep].chainLength = 1;
				}
				else if (alignType == Local) {
					fragmentSet[fSweep].chainPrev = (int) -1;
					fragmentSet[fSweep].cost = fragmentLength * match - fragmentSet[fSweep].weight;
					fragmentSet[fSweep].chainLength = 1;
				}					
			}

			if (minFragmentCost > fragmentSet[fSweep].cost) {
				minFragmentCost = fragmentSet[fSweep].cost;
				minFragmentIndex = fSweep;
				//	maxChainLength = fragmentSet[fSweep].chainLength;
			}

			if (fragmentSet[fSweep].chainLength > maxChainLength) {
				maxChainLength = fragmentSet[fSweep].chainLength;
				maxChainFragment = fSweep;
			}

			// Done computing the optimal score for this fragment.
			fSweep++;
		}
		
		//
		// Insert all fragments in the sweep set 
		//
		fSweep = startF;
		while (fSweep < fragmentSetSize and 
					 fragmentSet[fSweep].x == sweepRow) {
			//			cout << "inserting sweep set with index" << fragmentSet[fSweep].index << endl;
			sweepSet.Insert(fragmentSet[fSweep]);
			++fSweep;
		}
		
		// Remove elements from the sweep set that are too far back.
		if (sweepRow >= fragmentLength + 1) {
			trailRow = sweepRow - fragmentLength - 1;
			while (fTrail < fragmentSetSize and
						 fragmentSet[fTrail].x == trailRow) {
				//
				// These elements are removed from the sweep set since they are done being processed.
				// If they are the lowest cost in the value, update colSet
				//
				SDPColumn col;
				int storeCol = 0;
				col.col = fragmentSet[fTrail].y;

				if (colSet.Member(col)) {
					if (fragmentSet[col.optFragment].cost < fragmentSet[fTrail].cost) {
						storeCol = 1;
					}
				}
				else {
					storeCol = 1;
				}
				if (storeCol) {
					col.col = fragmentSet[fTrail].y;
					col.optFragment = fTrail;
					// 
					// Insert new column or replace col with a more optimal one.
					//
					colSet.Insert(col);

					// 
					// The invariant structure of the colSet is that
					// after inserting a fragment of score S at column col, 
					// the score of all columns greater than 'col' in col set
					// must be less than col. 
					//
					// To preserve this invariant, when an element is inserted
					// at 'col', look to columns greater.  As long as any columns
					// have scores that are greater than col, remove them.
					// Once a column col_next has been found that has a score less than S
					// by the structure of the loop invariant, all columns greater than col_next
					// are guaranteed to have lower score than S, so we can continue searching
					// through this loop.
					//
					// Since fragments are processed at most once, this remains O(M).
				
					SDPColumn successorCol = col;

					while (colSet.Successor(col, successorCol) and
								 fragmentSet[successorCol.optFragment].cost > fragmentSet[fTrail].cost) {
						colSet.Delete(successorCol);
					}
				}

				//
				// Now remove this fragment, it is at the end of the sweep line.
				//
				int deleted;
				deleted = sweepSet.Delete(fragmentSet[fTrail]);
				assert(deleted);

				++fTrail;
			}
		}
	}
	if (alignType == Local) {
		maxChainFragment = minFragmentIndex;
	}
	while (maxChainFragment != -1) {
		maxFragmentChain.push_back(maxChainFragment);
		maxChainFragment = fragmentSet[maxChainFragment].chainPrev;
	}
	std::reverse(maxFragmentChain.begin(), maxFragmentChain.end());
	return maxFragmentChain.size();
}

#endif
