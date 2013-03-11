#ifndef JABONEXCEPTION_H_
#define JABONEXCEPTION_H_

#include <string>
using namespace std;

class JabonException
{
  public:
    JabonException() : _msg("") {}
    JabonException( const string& msg ) : _msg(msg) {}
    
    const string& getMessage() const { return _msg; }
    
  private:
    string _msg;
};

#endif /*JABONEXCEPTION_H_*/
