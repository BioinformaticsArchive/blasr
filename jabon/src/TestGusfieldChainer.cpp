#include "GusfieldChainer.h"
#include <iostream>
using namespace std;

void test1(void)
{
	cout << endl << "Test1:" << endl << endl;
	Chainer chainer = Chainer();
  	chainer.addRectangle( new Rectangle( 1,2,1,2,3,0) );
  	chainer.addRectangle( new Rectangle( 4,6,3,4,7,1) );
  	chainer.addRectangle( new Rectangle( 3,5,5,6,5,2) );
  	chainer.addRectangle( new Rectangle( 7,8,7,8,9,3) );
  	 
  	chainer.calcMaximumChain();
  	vector<int> path = chainer.getMaximumChain(); 
	for (int i = 0; i < path.size(); i++)
	{
		cout << path[i] << " -> ";
	}
	cout << "score = " << chainer.getHighScore() << endl;
}

void test2(void)
{
	cout << endl << "Test2:" << endl << endl;
	Chainer chainer = Chainer();
  	chainer.addRectangle( new Rectangle( 1,2,1,2,1,0) );
  	chainer.addRectangle( new Rectangle( 3,4,5,6,3,1) );
  	chainer.addRectangle( new Rectangle( 5,6,3,4,1,2) );
  	chainer.addRectangle( new Rectangle( 7,8,5,6,9,3) );

  	chainer.calcMaximumChain();
  	vector<int> path = chainer.getMaximumChain(); 
	for (int i = 0; i < path.size(); i++)
	{
		cout << path[i] << " -> ";
	}
	cout << "score = " << chainer.getHighScore() << endl;
}

void test3(void)
{
	cout << endl << "Test3:" << endl << endl;
	Chainer chainer = Chainer();
  	chainer.addRectangle( new Rectangle( 1,2,1,2,1,0) );
  	chainer.addRectangle( new Rectangle( 3,4,5,6,3,1) );
  	chainer.addRectangle( new Rectangle( 3,4,4,5,3,2) );
  	chainer.addRectangle( new Rectangle( 5,6,3,4,1,3) );
  	chainer.addRectangle( new Rectangle( 7,8,5,6,9,4) );
  	 
  	chainer.calcMaximumChain();
  	vector<int> path = chainer.getMaximumChain(); 
	for (int i = 0; i < path.size(); i++)
	{
		cout << path[i] << " -> ";
	}
	cout << "score = " << chainer.getHighScore() << endl;
}

void test4(void)
{
	cout << endl << "Test4:" << endl << endl;
	Chainer chainer = Chainer();
  	chainer.addRectangle( new Rectangle( 1,3,1,3,1,0) );
  	chainer.addRectangle( new Rectangle( 2,4,2,4,3,1) );
  	chainer.addRectangle( new Rectangle( 3,5,3,5,1,2) );
  	chainer.addRectangle( new Rectangle( 6,7,6,7,4,3) );
  	 
  	chainer.calcMaximumChain();
  	vector<int> path = chainer.getMaximumChain();  
	for (int i = 0; i < path.size(); i++)
	{
		cout << path[i] << " -> ";
	}
	cout << "score = " << chainer.getHighScore() << endl;
}

int main( int argc, char **argv )
{
	test1();
	test2();
	test3();
	test4();
}
