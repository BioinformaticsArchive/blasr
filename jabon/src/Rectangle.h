#ifndef RECTANGLE_H_
#define RECTANGLE_H_

#include <string>

using namespace std;

class Rectangle
{
	public:
		int left;
		int right;
		int top;
		int bottom;
		int weight;
		int id;
		Rectangle(int left, int right, int top, int bottom, int weight, int id);
		Rectangle();
		string toString(void);
};

class RectangleBottom 
{
	public:
		int bottom;
		int id;		
		RectangleBottom() : bottom(0), id(0)
		{  };
		RectangleBottom(int bottom_param, int id_param)
		{
			bottom = bottom_param;
			id = id_param;
		};
		string toString(void);
};

class RectangleTerminus
{
	public:
		int x;
		int id;
		bool left;
		RectangleTerminus(Rectangle *rectangle, bool left);
		string toString(void);
};

#endif /*RECTANGLE_H_*/
