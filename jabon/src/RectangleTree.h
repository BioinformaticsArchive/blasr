#ifndef TESTAVL_H_
#define TESTAVL_H_

#include <vector>
#include <iostream>
#include <stdlib.h>
#include "Rectangle.h"
#include "AVL.h"

using namespace std; 

class RectangleTree
{
	public:
		// we use two AVL trees, because our library only support GetNext
		// not get previous. Both are kept in sync.
		AVL<int, RectangleBottom*> higher_first_tree;
		AVL<int, RectangleBottom*> lower_first_tree;
		vector<RectangleBottom*> all_bottoms;
		vector<int> weights;
		vector<int> traceback; // this may or may not work

		bool minimum_higher_bottom(int top, RectangleBottom *result);
		bool maximum_bottom_not_higher(int bottom, RectangleBottom *result);
		// RectangleTree will now own this bottom.
		void insert_bottom(RectangleBottom *bottom);
		void delete_useless_bottoms(RectangleBottom *new_bottom);
		int get_max_id();

		RectangleTree();
		~RectangleTree();
		// PORT are these needed?
		string bottomString(RectangleBottom bottom);
		int highScore();		
		void print_tree(void);
};

#endif /*TESTAVL_H_*/
