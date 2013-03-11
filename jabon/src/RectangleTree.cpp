#include "RectangleTree.h"
#include <iostream>
#include <stdlib.h>
#include <vector>

using namespace std; 	//introduces namespace std

RectangleTree::RectangleTree()
{	
	return;
}

RectangleTree::~RectangleTree()
{
	for (int i = 0; i < all_bottoms.size(); i++)
	{
		delete(all_bottoms[i]);
	}
	AVL<int, RectangleBottom*>::Iterator<int, RectangleBottom*> iter(&this->lower_first_tree);
	int key;
	RectangleBottom* item;
	bool more = iter.GetFirst(key, item);
	while(more){
		lower_first_tree.Remove(key);
		more = iter.GetNext(key, item);
	}
	AVL<int, RectangleBottom*>::Iterator<int, RectangleBottom*> iter2(&this->higher_first_tree);
	more = iter2.GetFirst(key, item);
	while(more){
		higher_first_tree.Remove(key);
		more = iter2.GetNext(key, item);
	}
	// MEMLEAK do we need more?
}

bool RectangleTree::minimum_higher_bottom(int top, RectangleBottom *result)
{
	// cout << "Minimum higher bottom called with top " << top << endl;
	RectangleBottom *temp_bottom;
	int temp_int;
	AVL<int, RectangleBottom*>::Iterator<int, RectangleBottom*> iter(&this->lower_first_tree);	
	if (!iter.GetFirst(temp_int, temp_bottom))
	{
		return false;
	}
	
	bool more = iter.FindGte(top+1, temp_bottom);	 // +1 because do not want equal to.
	
	// cout << "Tree" << endl;
	// this->print_tree();
	if (more)
	{
		// cout << "Found minimum higher bottom " << temp_bottom->id << endl;
		*result = *temp_bottom;
	}
	return more;
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

	bool more = iter.FindGte(-bottom, temp_bottom);	
	if (temp_bottom != NULL)
	{
		*result = *temp_bottom;
	}
	// cout << "Returning from maximum bottom not higher " << more << endl;
	return more;
}

// RectangleTree now owns this bottom
void RectangleTree::insert_bottom(RectangleBottom * bottom)
{
	// cout << "Inserting bottom " << bottom->toString() << endl;
	higher_first_tree.Insert(-bottom->bottom, bottom); 
	// above is negative to sort from top to bottom
	lower_first_tree.Insert(bottom->bottom, bottom);
	all_bottoms.push_back(bottom);
}

void RectangleTree::delete_useless_bottoms(RectangleBottom* new_bottom)
{
	RectangleBottom *temp_bottom;
	AVL<int, RectangleBottom*>::Iterator<int, RectangleBottom*> iter(&this->higher_first_tree);	

	int temp_int = -new_bottom->bottom;
	bool more = iter.Find(temp_int, temp_bottom);	
	while (more){
		if (temp_bottom->bottom <= new_bottom->bottom &&
			weights[ temp_bottom->id ] < weights[ new_bottom->id ])
		{
			// delete this bottom from both trees
			// cout << "Deleting " << temp_int << endl;
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

int RectangleTree::get_max_id()
{
	int max_weight = 0;
	int max_id = 0;
	for (int i=0; i < weights.size(); i++)
	{
		if (weights[i] > max_weight)
		{
			max_weight = weights[i];
			max_id = i;
		}
	}
	return max_id;
}

/*
int main( void )
{
	int i;
	int myints[] = {1, 2, 1};
  	vector<int> weights (myints, myints + sizeof(myints) / sizeof(int) );
	RectangleTree tree = RectangleTree();
	tree.weights = weights;
	RectangleBottom* bottom0 = new RectangleBottom(0, 0);
	RectangleBottom* bottom1 = new RectangleBottom(1, 1);
	RectangleBottom* bottom2 = new RectangleBottom(2, 2);	
	// bottom0 is effectively useless, since it is below 1 and has lower weight

	tree.insert_bottom(bottom0);
	tree.insert_bottom(bottom1);	
	tree.insert_bottom(bottom2);
	RectangleBottom* bottom_not_higher = new RectangleBottom();
	tree.maximum_bottom_not_higher(2, bottom_not_higher);
	cout << "Bottom not higher than 2 = " << bottom_not_higher->toString() << endl;
	
	RectangleBottom* minimum_higher_bottom = new RectangleBottom();
	tree.minimum_higher_bottom(1, minimum_higher_bottom);
	cout <<"Minimum higher bottom than 1 starting at 0 = " << minimum_higher_bottom->toString() << endl;
	
	bool exists = tree.minimum_higher_bottom(2, minimum_higher_bottom);
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
}*/
