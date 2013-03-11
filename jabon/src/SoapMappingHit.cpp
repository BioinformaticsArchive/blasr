#include "SoapMappingHit.h"

SoapMappingHit::SoapMappingHit(string query_id, string target_id)
{
    	this->query_id = query_id;
    	this->target_id = target_id;
}



void SoapMappingHit::addHit(SoapShortHit * hit)
{
	if (hits.size() == 0)
	{
		hits.push_back(hit);
		return;	
	}
	
	SoapShortHit * finalHit = hits[hits.size() - 1];
	if (not finalHit->appendAble(hit))
	{
		hits.push_back(hit);
	}
	else
	{
		finalHit->append(hit);
	}
}

SoapHitChain * SoapMappingHit::maximumChain(
								int diagonalCutoff, 
								int joinCutoff, 
								int minChainLength)
{		
	Chainer chainer = Chainer(diagonalCutoff, joinCutoff, minChainLength);
	
	// cerr << "Found " << hits.size() << " hits " << endl;
	
	for (int hit_idx = 0; hit_idx < hits.size(); hit_idx++)
	{
		SoapShortHit * hit = hits[hit_idx];
		int weight = abs(hit->query_end  - hit->query_start) +
					 abs(hit->target_end - hit->target_start);
		Rectangle *new_rect;
		if (hit->onPositiveTargetStrand())
		{
			new_rect = new Rectangle(	
									hit->query_start,  
									hit->query_end,
									hit->target_start, 
									hit->target_end,
									weight,
									hit_idx); // id
		}
		else
		{
			new_rect = new Rectangle(	
									hit->query_start,  
									hit->query_end,
									- hit->target_end,    // so that subsequent hits on the 
									- hit->target_start,  // target are still to the "right"
									weight,
									hit_idx); 
		}
		// cerr << "Rectangle " << hit_idx << " is "  << new_rect->toString() << endl;
		chainer.addRectangle(new_rect);
	}
	
	chainer.calcMaximumChain();
	vector<int> ids = chainer.getMaximumChain();
	int highScore  =  chainer.getHighScore();
	
	// cerr << "Path:" << endl;
	SoapHitChain * hit_chain = new SoapHitChain(this->query_id, this->target_id, highScore);
	for (int i=0; i < ids.size(); i++)
	{
		// cerr << ids[i] << " -> ";
		hit_chain->addHit( hits[ ids[i] ] );
	}
 	// cerr << "Score:" << highScore << endl;
	return hit_chain;
}
