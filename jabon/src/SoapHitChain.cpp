#include"SoapHitChain.h"


SoapHitChain::SoapHitChain(string query_id, string target_id, int score)
	: 	query_id(query_id),
		target_id(target_id),
		query_integer(0),
		target_integer(0),
		score(score),
		positive_strand(true)
{
    this->min_hit.query_start   = INT_MAX; // TODO make a ptr and put in init?
    this->min_hit.query_end     = INT_MAX;
    this->min_hit.target_start  = INT_MAX;
    this->min_hit.target_end    = INT_MAX;
    
    this->max_hit.query_start   = INT_MIN;
    this->max_hit.query_end     = INT_MIN;
    this->max_hit.target_start  = INT_MIN;
    this->max_hit.target_end    = INT_MIN;
   
    std::istringstream query_iss( query_id );
    query_iss >> this->query_integer;
    
    std::istringstream target_iss( target_id ); 
    target_iss >> this->target_integer;
}

void SoapHitChain::addHit( SoapShortHit * hit )
{
	hits.push_back(hit);
	// score += hit->query_length + hit->target_length;
	if (hits.size() == 1)
	{	
		this->positive_strand = hit->onPositiveTargetStrand();
	}
	else if (this->positive_strand != hit->onPositiveTargetStrand())
	{
		cerr << "Tried to add two hits on different strands to a chain!" << endl;
		exit(1);
	}
	if (compare(hit, &min_hit) == -1)
	{		
		min_hit = *hit;
	} 
	if (compare(hit, &max_hit) == 1)
	{
		max_hit = *hit;
	}
}
    
int SoapHitChain::compare( 
	SoapShortHit * hit1, 
	SoapShortHit * hit2)
{
   	return 	hit1->query_start > hit2->query_start  ? 1 :
   			hit1->query_start == hit2->query_start ? 0 :
   												  -1;
}

int SoapHitChain::getLength(void)
{
	return this->hits.size();
}

const string SoapHitChain::toAmosString()
{	
	if (hits.size() == 0)
	{
		return "";
	}
	stringstream stream;
	stream << "{OVL" << endl << "adj:";
	stream << (this->positive_strand ? "N" : "I");
	stream << endl << "rds:" << this->query_integer << "," << this->target_integer << endl;
	stream << "scr:" << this->score << endl; 
	
	int bhg, ahg;
	if (this->positive_strand)
	{	// AAK this logic was checked against hash-overlap.
		ahg = min_hit.query_start - min_hit.target_start;
		bhg = (max_hit.target_length - max_hit.target_end) - (max_hit.full_query_length - max_hit.query_end); 
	} 
	else
	{	// AAK this logic was checked against hash-overlap.
		ahg = min_hit.query_start - (min_hit.target_length - min_hit.target_end);
		bhg = max_hit.target_start - (max_hit.full_query_length - max_hit.query_end);
	}
	
	stream << "ahg:" << ahg << endl;
	stream << "bhg:" << bhg << endl;
	stream << "}" << endl;
	// stream << endl << "Min:" << min_hit.toString() << endl << endl;
	// stream << "Max:" << max_hit.toString() << endl << endl;
	return stream.str();
}
