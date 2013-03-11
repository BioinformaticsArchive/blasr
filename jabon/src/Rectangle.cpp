#include "Rectangle.h"

Rectangle::Rectangle() { }

Rectangle::Rectangle(int left, int right, int top, int bottom, int weight, int id)
{
	this->left = left;
	this->right = right;
	this->top = top;
	this->bottom = bottom;
	this->weight = weight;
	this->id = id;
}

RectangleTerminus::RectangleTerminus(Rectangle *rectangle, bool leftboolean)
{	
	x = leftboolean ? rectangle->left : rectangle->right;
	id = rectangle->id;
	left = leftboolean;
};

string Rectangle::toString(void)
{
	char buffer[100];
	int s=sprintf (buffer, "left %i right %i top %i bottom %i weight %i id %i", left, right, top, bottom, weight, id);
	return buffer;
}

string RectangleTerminus::toString(void)
{
	char buffer[100];
	int s=sprintf (buffer, "x=%i id=%i isLeft=%s", x, id, left ? "true" : "false");
	return buffer;
}

string RectangleBottom::toString(void)
{
	char buffer[100];
	int s=sprintf (buffer, "bottom=%i id=%i", bottom, id);
	return buffer;
}


