#ifndef ALGORITHMS_ALIGNMENT_BASE_SCORE_FUNCTION_H_
#define ALGORITHMS_ALIGNMENT_BASE_SCORE_FUNCTION_H_

class BaseScoreFunction {
 public:
	int ins;
	int del;
  int substitutionPrior;
  int globalDeletionPrior;
  int affineExtend;
  BaseScoreFunction() {
    ins = del = substitutionPrior = globalDeletionPrior = affineExtend = 0;
  }

  BaseScoreFunction(int insP, int delP, int subPriorP, int delPriorP, int affineExtensionP)  {
    ins = insP;
    del = delP;
    substitutionPrior = subPriorP;
    globalDeletionPrior = delPriorP;
    affineExtend = affineExtensionP;
  }
};

#endif
