#ifndef PELUSAEXCEPTION_H_
#define PELUSAEXCEPTION_H_
#include <string>
using namespace std;

class PelusaException
{
  public:
    PelusaException() : _msg("") {}
    PelusaException( const string& msg ) : _msg(msg) {}
    
    const string& getMessage() const { return _msg; }
    
  private:
    string _msg;
};
#endif /*PELUSAEXCEPTION_H_*/
