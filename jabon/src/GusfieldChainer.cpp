#include "GusfieldChainer.h"

using namespace std;

static bool compare_termini(RectangleTerminus * one, RectangleTerminus * two)
{
  return one->x < two->x;
}

Chainer::Chainer()
	: 	highScore(0),
		diagonalCutoff(3),
		joinCutoff(6),
		minimumChainLength(1),
        maxY(-1)
{
}

Chainer::Chainer(int diagonalCutoff, int joinCutoff, int minimumChainLength)
	: 	highScore(0),
		diagonalCutoff(diagonalCutoff),
		joinCutoff(joinCutoff),
		minimumChainLength(minimumChainLength)	
{
	
}

void Chainer::addRectangle(Rectangle *rectangle)
{
	// cerr << "Adding rectangle: " << rectangle->toString() << endl;
	rectangles.push_back(rectangle);
}

vector<int> Chainer::getMaximumChain(void)
{
	return maximum_chain;
}

bool Chainer::chainable(Rectangle * lowerRectangle, Rectangle * upperRectangle)
{
    // cerr << "trying to chain " << endl << lowerRectangle->toString() << endl <<
    //        " to " << endl << upperRectangle->toString() << endl;
    
    // TODO at some point remove maxY from Chainer
    // it is necessary now to undo what the reverseRectangleYcoords did
    // to conform to Gusfield so that we can 
    // calculate the off diagonal properly here. A little annoying.
	int upperOffDiagonal = upperRectangle->right - (maxY - upperRectangle->bottom);
	int lowerOffDiagonal = lowerRectangle->left - (maxY - lowerRectangle->top);

	int diagonalDifference = abs(lowerOffDiagonal - upperOffDiagonal);

	int joinDistance = 	abs(lowerRectangle->left   - upperRectangle->right) + 
					    abs(upperRectangle->bottom - lowerRectangle->top);
    joinDistance /= 2;

    // cerr << " joinDistance " << joinDistance << " joinCutoff " << joinCutoff << " diagonalDifference " << diagonalDifference << " diagonalCutoff " << diagonalCutoff << endl;
	if (diagonalDifference <= diagonalCutoff && joinDistance <= joinCutoff)
	{
		return true;
	}  
    // cerr << "chainable = false" << endl;
	return false;
}


// the core algorithm
void Chainer::calcMaximumChain(void)
{
	// Construct an array L that will store of good rectangle bottoms. 
  	// Sorted so that rects with higher bottoms in the DP matrix are first.
  	// This is the fundamental data structure for this algorithm.
  	RectangleTree L; // MEMLEAK should be on the stack, hence destroyed at end of method.
	for (int i=0; i < rectangles.size(); i++)
	{
		L.weights.push_back( rectangles[i]->weight );
	}
	
	// traceback structure
	map<int, int> traceback;
	
	// Transform coords so that highest y is now lowest to conform with Gusfield
	reverseRectangleYcoords(); // TODO this confuses things
	
	// Create rectangle x-axis termini on the x-axis ...
	vector<RectangleTerminus*> termini;
	for (int i=0; i < rectangles.size(); i++)
	{
		termini.push_back( new RectangleTerminus(rectangles[i], true) );
		termini.push_back( new RectangleTerminus(rectangles[i], false) );
	}
	std::sort(termini.begin(), termini.end(), compare_termini );
	 
	// ... which we traverse from left to right in the DP matrix
	for (int i=0; i < termini.size(); i++)
	{
		RectangleTerminus * terminus = termini[i];
		// cerr << "Terminus " << terminus.toString() << endl;
		int id = terminus->id;
		if (terminus->left)
		{
			RectangleBottom min_bottom;
			if (L.minimum_higher_bottom(rectangles[id]->top, &min_bottom))
			{
				// cout << "Weight = " << L.weights[ id ] << endl;
				if (chainable(rectangles[id], rectangles[min_bottom.id]))
				{
				    // cout << "traceback id " << id << " -> " << min_bottom.id << endl;
					L.weights[ id ] += L.weights[ min_bottom.id ];
					traceback.insert(pair<int,int>(id, min_bottom.id));
				}
			}
		}
		else
		{	
			int terminus_bottom = rectangles[id]->bottom;
			// The RectangleTree L will own this RectangleBottom
			RectangleBottom max_bottom;
			if (not L.maximum_bottom_not_higher(rectangles[id]->bottom, &max_bottom))
			{
				L.insert_bottom(new RectangleBottom(terminus_bottom, id));
			} 
			else if ((max_bottom.bottom <  terminus_bottom ) 
										||
					((max_bottom.bottom == terminus_bottom) && 
					 (L.weights[ id ] > L.weights[ max_bottom.id ])))
			{
				RectangleBottom * new_bottom = new RectangleBottom(terminus_bottom, id);
				L.insert_bottom( new_bottom );
				L.delete_useless_bottoms( new_bottom );
			} 	
			// cerr << "Max bottom " << max_bottom.bottom << " Terminus bottom " << terminus_bottom << endl; 
		}
	}
	
	// We're done with the termini
	for (int i = 0; i < termini.size(); i++) delete(termini[i]);
	 // 
	 // DEBUG 
	 // for(map<int,int>::const_iterator it = traceback.begin(); it != traceback.end(); ++it)
     //	cout << it->first << " ->  " << it->second << endl; 
        
	// At the end, the last entry in L gives the value of the optimal chain 
    // and specifies the last rectangle in the chain. with a traceback, we 
	// can find the correct path 
	maximum_chain.push_back( L.get_max_id() );
	// cout << "Max id = " << L.get_max_id() << endl;
	while (traceback.find(maximum_chain[ maximum_chain.size() - 1]) != traceback.end())
	{
		maximum_chain.push_back( traceback.find( maximum_chain[ maximum_chain.size() - 1] )->second);
	}
	reverse(maximum_chain.begin(), maximum_chain.end());
	
	// cout << "Path length: " << maximum_chain.size() << endl;
	if (maximum_chain.size() < minimumChainLength)
	{
		maximum_chain.clear();
		return;
	}
	
	highScore = L.weights[maximum_chain[maximum_chain.size()-1]];
	return;
}

int Chainer::getHighScore(void)
{
	return highScore;
}

void Chainer::reverseRectangleYcoords(void)
{
	int maxY = 0;
	for (int i=0; i < rectangles.size(); i++)
	{
		maxY = rectangles[i]->bottom > maxY ? rectangles[i]->bottom : maxY;
	}
	for (int i=0; i < rectangles.size(); i++)
	{
		
		rectangles[i]->bottom = maxY - rectangles[i]->bottom;
		rectangles[i]->top = maxY - rectangles[i]->top;
	}
    this->maxY = maxY;
}

Chainer::~Chainer()
{
	for(int i=0; i < rectangles.size(); i++)
	{
		delete(rectangles[i]);
	}
}




