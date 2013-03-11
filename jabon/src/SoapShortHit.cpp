#include "SoapShortHit.h"

SoapShortHit* SoapShortHit::parseLine( const string& line )
{
  typedef boost::tokenizer< boost::char_separator<char> > tokenizer;
  //boost::char_separator<char> sep("\t");
  boost::char_separator<char> sep( " \t" );
  tokenizer tok( line, sep );
  SoapShortHit *hit = new SoapShortHit();
  int j = 0;
  string prefix;
  // cerr << line << endl;
  for( tokenizer::iterator i=tok.begin(); i!=tok.end(); ++i, ++j )
  {
    // cerr << "parsing " << *i << endl;
    switch(j)
    {
      case 0: getSuffix(string(*i), ':', hit->query_id);
              getPrefix(string(*i), ':', prefix );
              if ( prefix.length()>0 )
                hit->query_start = atoi( prefix.c_str() );
              break;
      case 1: hit->num_hits = atoi(i->c_str()); 
      		  break;
      case 2: // skip
              break;
      case 3: hit->query_length = atoi(i->c_str());
              hit->query_end = hit->query_start + hit->query_length; 
              break;
      case 4: hit->target_strand = (*i)[0];
              break;
      case 5: hit->target_id = string(*i); 
      		  break;
      case 6: hit->target_start = atoi(i->c_str()) - 1;
              hit->target_end = hit->target_start + hit->query_length; 
              break;
      case 7: // skip
      		  break;
      case 8:
      		  hit->target_length = atoi(i->c_str());
      		  break;
    }
  }
  
  // Note that coords for strand "-" are stored on the reverse complement
  if (hit->target_strand == '-')
  {
  	hit->reverseTargetCoords();
  }
  return hit;
}

bool SoapShortHit::onPositiveTargetStrand()
{ 
	return this->target_strand == '+' ? true : false;
}

bool SoapShortHit::appendAble(SoapShortHit * otherHit)
{
	// cerr << "Checking hit " << this->toString() << " against " << endl;
	// cerr << " hit " << otherHit->toString() << endl;
	if (this->target_strand != otherHit->target_strand)
	{
		// cerr << "Returning false strand" << endl;
		return false;
	}
	// they have to be on the same diagonal...
	int lowerOffDiagonal = this->query_start - otherHit->query_start;
	int upperOffDiagonal = this->target_start - otherHit->target_start;
	upperOffDiagonal *= this->target_strand == '-' ? -1 : 1;
	int diagonalDifference = abs(lowerOffDiagonal - upperOffDiagonal);
	if (diagonalDifference != 0)
	{	
		// cerr << "Returning false diag " << diagonalDifference << " lower " << lowerOffDiagonal << " upper " << upperOffDiagonal << endl;
		return false;
	}
	// ... and there ends should be adjacent to one another
	int joinCutoff = 2;
	int distance = 	abs(this->query_end - otherHit->query_end) + 
					abs(this->target_end - otherHit->target_end);
	if (distance > joinCutoff)
	{
		// cerr << "Returning false join" << endl;
		return false;
	}  
	// cerr << "Returning true " << endl;
	return true;
}

bool SoapShortHit::append(SoapShortHit * otherHit)
{
	this->query_length += otherHit->query_end - this->query_end;
	this->query_end = otherHit->query_end;
	this->target_end = otherHit->target_end;
}

void SoapShortHit::getPrefix( const string& s, char delim, string& prefix )
{
  char buffer[40];
  buffer[0] = 0;
  int k = 0;
  int sLength = s.length();
  for( int l=0; l<sLength && k<40; l++ )
  {
    char ch = s[l];
    if ( ch==delim )
    {
      buffer[k] = 0;
      break;
    }
    buffer[k++] = ch;
  }
  if ( k==40 ) buffer[0] = 0;
  prefix = string(buffer);
}

void SoapShortHit::getSuffix( const string& s, char delim, string& suffix )
{
  char buffer[40];
  buffer[0] = 0;
  int k = 0;
  int sLength = s.length();
  bool seenDelim = false;
  for( int l=0; l<sLength && k<40; l++ )
  {
    char ch = s[l];
    if ( seenDelim )
 	{
    	buffer[k++] = ch;
 	}
    else if ( ch==delim )
    {
      seenDelim = true;
    }
  }
  buffer[k] = 0;
  if ( k==40 ) buffer[0] = 0;
  suffix = string(buffer);
}


void SoapShortHit::getQueryName( string& queryName ) const
{
  queryName.clear();
  int index = query_id.find( ':' );
  if ( index>=0 )
    queryName = query_id.substr( index+1, query_id.length()-index );
  else
    queryName = query_id;
}

void SoapShortHit::reverseTargetCoords(void)
{
	/*int saved_target_start = this->target_start;
	this->target_start = this->target_length - this->target_end;
	this->target_end = this->target_length - saved_target_start;*/
	
}

const string& SoapShortHit::toString()
{
	char buffer[1000];
	sprintf( buffer, "num_hits:%d  target_strand:%c "
					 "query_id:%s  query_start:%d  query_end:%d  query_length:%d "
					 "target_id:%s target_start:%d target_end:%d",
					 num_hits, target_strand,
					 query_id.c_str(), query_start, query_end, query_length,
					 target_id.c_str(), target_start, target_end);
	rep = string(buffer);
	return rep;
}
