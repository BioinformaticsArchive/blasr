#ifndef GUSFIELD_CHAINER_H_
#define GUSFIELD_CHAINER_H_

#include <vector>
#include <string>
#include <iostream>
#include <map>
#include <algorithm>

#include "Rectangle.h"
#include "RectangleTree.h"

using namespace std;
	
class Chainer
{
	
	public:
	
		Chainer();
		Chainer(int diagonalCutoff, int joinCutoff, int minimumChainLength);
		~Chainer();
		
		vector<Rectangle*> rectangles;
		vector<int> maximum_chain;
		// the Chainer owns the rectangle after it's added
		void addRectangle(Rectangle *rectangle);
		void calcMaximumChain(void);
		vector<int> getMaximumChain(void);
		int getHighScore(void);
		int highScore;
		int diagonalCutoff;
		int joinCutoff;
		int minimumChainLength;
        int maxY; 

			
	private:
		void reverseRectangleYcoords(void);
		bool chainable(Rectangle *lowerRectangle, Rectangle *upperRectangle);
};

#endif /*GUSFIELD_CHAINER_H_*/
