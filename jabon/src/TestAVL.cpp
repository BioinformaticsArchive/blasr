#include "TestAVL.h"
#include <iostream>
#include <stdlib.h>
#include <vector>

RectangleTree::RectangleTree (vector<int> weights)
{
	this->weights = weights;
}

bool RectangleTree::minimum_higher_bottom(int top, int starting_bottom, RectangleBottom *result)
{
	RectangleBottom *temp_bottom;
	int temp_int;
	AVL<int, RectangleBottom*>::Iterator<int, RectangleBottom*> iter(&this->lower_first_tree);	
	if (!iter.GetFirst(temp_int, temp_bottom))
	{	// empty tree
		return false;
	}

	bool more = iter.Find(starting_bottom, temp_bottom);	
	if (!more)
	{
		// sw error. we should be able to find this bottom!
		cerr << "Strangeness could not find input bottom!" << endl;
		return false;
	}
	more = iter.GetNext(starting_bottom, temp_bottom);
	while (more){
		if (temp_bottom->bottom > top)
		{
			*result = *temp_bottom;
			return true;
		}
		more = iter.GetNext(starting_bottom, temp_bottom);
	}
	return false;
}
		
bool RectangleTree::maximum_bottom_not_higher(int bottom, RectangleBottom *result)
{	
	RectangleBottom *temp_bottom;
	int temp_int;
	AVL<int, RectangleBottom*>::Iterator<int, RectangleBottom*> iter(&this->higher_first_tree);	
	if (!iter.GetFirst(temp_int, temp_bottom))
	{	// empty tree
		return false;
	}

	bool more = iter.Find(-bottom, temp_bottom);	
	if (!more)
	{
		// sw error. we should be able to find this bottom!
		cerr << "Strangeness could not find input bottom!" << endl;
		return false;
	}
	more = iter.GetNext(bottom, temp_bottom);
	if (temp_bottom != NULL)
	{
		*result = *temp_bottom;
	}
	return more;
}

void RectangleTree::insert_bottom(RectangleBottom *bottom)
{
	higher_first_tree.Insert(-bottom->bottom, bottom); 
	// above is negative to sort from top to bottom
	lower_first_tree.Insert(bottom->bottom, bottom);
}

void RectangleTree::delete_useless_bottoms(RectangleBottom new_bottom)
{
	RectangleBottom *temp_bottom;
	AVL<int, RectangleBottom*>::Iterator<int, RectangleBottom*> iter(&this->higher_first_tree);	

	int temp_int = -new_bottom.bottom;
	bool more = iter.Find(temp_int, temp_bottom);	
	while (more){
		if (temp_bottom->bottom <= new_bottom.bottom &&
			weights[ temp_bottom->id ] < weights[ new_bottom.id ])
		{
			// delete this bottom from both trees
			higher_first_tree.Remove(temp_int);
			lower_first_tree.Remove(-temp_int);
		}
		more = iter.GetNext(temp_int, temp_bottom);
	}	
}

void RectangleTree::print_tree(void)
{

	AVL<int, RectangleBottom*>::Iterator<int, RectangleBottom*> iter(&this->higher_first_tree);

	bool more;
	int key = 1;
	RectangleBottom* item;

	cout << "Higher first tree contains:\n";
	
	more = iter.GetFirst(key, item);
	while (more)
	{
		cout << "key - " << key << " item - " << item->toString() << "\n";
		more = iter.GetNext(key, item);
	}	
	
	AVL<int, RectangleBottom*>::Iterator<int, RectangleBottom*> lower_iter(&this->lower_first_tree);

	cout << "Lower first tree contains:\n";	
	more = lower_iter.GetFirst(key, item);
	while (more)
	{
		cout << "key - " << key << " item - " << item->toString() << "\n";
		more = lower_iter.GetNext(key, item);
	}	
	
}

using namespace std; 	//introduces namespace std
int main( void )
{
	int i;
	int myints[] = {1, 2, 1};
  	vector<int> weights (myints, myints + sizeof(myints) / sizeof(int) );
	RectangleTree tree = RectangleTree(weights);

	RectangleBottom bottom0 = RectangleBottom(0, 0);
	RectangleBottom bottom1 = RectangleBottom(1, 1);
	RectangleBottom bottom2 = RectangleBottom(2, 2);	
	// bottom0 is effectively useless, since it is below 1 and has lower weight

	tree.insert_bottom(&bottom0);
	tree.insert_bottom(&bottom1);	
	tree.insert_bottom(&bottom2);
	RectangleBottom bottom_not_higher;
	tree.maximum_bottom_not_higher(2, &bottom_not_higher);
	cout << "Bottom not higher than 2 = " << bottom_not_higher.toString() << endl;
	
	RectangleBottom minimum_higher_bottom;
	tree.minimum_higher_bottom(1, 0, &minimum_higher_bottom);
	cout <<"Minimum higher bottom than 1 starting at 0 = " << minimum_higher_bottom.toString() << endl;
	
	bool exists = tree.minimum_higher_bottom(2, 0, &minimum_higher_bottom);
	cout << "Minimum higher bottom than 2 starting at 0 does";
	if (exists)
		printf(" exist!\n");
	else
		printf(" not exist!\n");
	
	cout << "Tree before deletion:" << endl;	
	tree.print_tree();

	tree.delete_useless_bottoms(bottom1);

	cout << "Tree after deletion:" << endl;	
	tree.print_tree();

	return 0;
}
