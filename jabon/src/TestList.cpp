#include <list>
#include <iostream>
#include <time.h>

using namespace::std;
int main( int argc, char **argv){

  list<char> alpha;

  /* for (int i=0; i < 10; i++){
    alpha.push_back( i + 65 );
  } */

  int max = 10000000;
  for (int i=0; i < max; i++)
  {
    alpha.insert( alpha.begin(), 1, 'A' );
  }
  time_t start, end;
  time(&start);

  for (int i=0; i < 10000000; i++)
  {
    alpha.insert( alpha.begin(), 1, 'A' );
  }

  time(&end);
  double dif = difftime(end,start);
  printf("Hey %.2f\n", dif);


  /* for (list<char>::iterator iterator = alpha.begin(); iterator != alpha.end(); iterator++)
  {
    // cout << *iterator;
  } */

  return 0;
}
