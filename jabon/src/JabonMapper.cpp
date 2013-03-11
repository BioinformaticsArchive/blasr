#include "JabonMapper.h"

void *t_SingleAlign(void *);

extern Param param;

// we have to use global variables until we can switch
// over to a more elegant thread solution
RefSeq          *g_ref;
ReadClass       *g_read_a;
ifstream        *g_fin_a;
pthread_mutex_t g_mutex_fin=PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t g_mutex_fout=PTHREAD_MUTEX_INITIALIZER;
// JMS
ref_id_t        g_longReadUid = 0;
bit32_t         *g_n_aligned;

JabonMapper::JabonMapper( const string& queryFile ) : queryFile(queryFile), ref(),
  n_aligned(0)
{
}

void JabonMapper::initialize( const string& targetFile )
{
  ifstream      fin_db;
  
  Initial_Time();
  
  fin_db.open(targetFile.c_str());
  if(!fin_db) 
  {
    cerr<<"fatal error: failed to open ref file: " << targetFile << endl;
    throw JabonException();
  }
  ref.Run_ConvertBinseq(fin_db);
  
  cerr<<"Load in "<<ref.total_num<<" db seqs, total size "<<ref.sum_length<<" bp. "<<Cal_AllTime()<<" secs passed"<<endl;
}

void JabonMapper::run()
{
  Do_Formatdb();
  RunProcess();
  ref.ReleaseIndex();
  // MEMLEAK valgrind doesn't really like ReleaseIndex (double frees?)
  // but it prevents memleak.
}

void JabonMapper::Do_Formatdb()
{
  ref.CreateIndex();
  cerr << "Create seed table. " <<Cal_AllTime()<<" secs passed"<<endl;
}

JabonMapper::~JabonMapper()
{
}

void JabonMapper::RunProcess()
{
  if(!queryFile.empty()) 
  {
    fin_a.open(queryFile.c_str());
    if(!fin_a) {
      cerr<<"Failed to open file: "<<queryFile<<endl;
      throw JabonException();
    }
  }
  else
  {
    throw JabonException("missing query file(s)");
  }
  
  //cerr<<"Single read alignment:\n";
  //cerr<<"Query: "<<queryFile<<"  Reference: "<<ref_file;
  
  n_aligned=0;
  read_a.InitialIndex();
  
  // set up globals
  g_read_a = &read_a;
  g_fin_a = &fin_a;
  g_n_aligned = &n_aligned;
  g_ref = &ref;
  
  DoSingleAlign();
  fin_a.close();
  cerr<<"Total number of aligned reads: "<<n_aligned<<" ("<<setprecision(2)
    <<100.0*n_aligned/read_a._index<<"%)\n";
}

void JabonMapper::DoSingleAlign()
{
  read_a.CheckFile(fin_a);
  vector<pthread_t> pthread_ids(param.num_procs);
    //create
  for(int i=0; i<param.num_procs; i++)
    pthread_create(&pthread_ids[i], NULL, t_SingleAlign, NULL);
  //join
  for (int i=0; i<param.num_procs; i++)
    pthread_join(pthread_ids[i], NULL);
}

void *t_SingleAlign(void *)
{
  JabonMapperWorker worker = JabonMapperWorker();
  worker.run();
}

JabonMapperWorker::JabonMapperWorker()
{
  
}

void JabonMapperWorker::run()
{
  int n;
  bit32_t cur_at;
  a.ImportFileFormat(g_read_a->_file_format);
  a.SetFlag('a');
  // JMS
  a.short_format = param.short_format;
  while(1)
  {
    pthread_mutex_lock(&g_mutex_fin);
    // JMS
    if(!param.chopReads)
        n = g_read_a->LoadBatchReads(*g_fin_a);
    else
    {
        n = g_read_a->LoadBatchLongReads( *g_fin_a, g_longReadUid );
        g_longReadUid += n;
    }
    cur_at=g_read_a->_index;
    a.ImportBatchReads(g_read_a->num, g_read_a->mreads);
    pthread_mutex_unlock(&g_mutex_fin);
    if(!n)
      break;
        // cerr << "Before do batch\n";
    a.Do_Batch(*g_ref);
        //cerr << "After do batch\n";
    
    postProcessReads();
    
    //pthread_mutex_lock(&g_mutex_fout);
    //cout<<a._str_align;
    //cerr<<cur_at<<" reads finished. "<<Cal_AllTime()<<" secs passed"<<endl;
    //pthread_mutex_unlock(&g_mutex_fout);    
  }
  pthread_mutex_lock(&g_mutex_fout);
  *g_n_aligned+=a.n_aligned;
  pthread_mutex_unlock(&g_mutex_fout);
}

//
// chain the hits for a reads into larger, more significant hits 
//
void JabonMapperWorker::postProcessReads()
{
   typedef boost::tokenizer< boost::char_separator<char> > tokenizer;
   // given how complex the SOAP code is, the best way to get the
   // reads out is to parse the string output
   // really not too much of a performance hit
   // (and the boost tokenizing is streaming, not in memory)
   string& sHits = a._str_align;
   boost::char_separator<char> newlineSep("\n");
   tokenizer tok( sHits, newlineSep );
   
   vector<SoapShortHit*> short_hits;
   
   int number = 0;
   string currentName = "";
   for( tokenizer::iterator i=tok.begin(); i!=tok.end(); ++i )
   {
     // convert to SOAP hit. JabonMapper owns the SoapShortHit.
     SoapShortHit * hit = SoapShortHit::parseLine( *i );
     hit->full_query_length = g_read_a->id2length.find( hit->query_id )->second;
     // group by read name
     // TODO it is not true that all read hits come out of soap simulataneously
     // I'm not sure what to do about this.
     if (currentName != "" && hit->query_id != currentName)
     {
     	postProcessReadHits(short_hits); // MEMLEAK I think this copies the whole vector!
     	for (int i=0; i< short_hits.size(); i++) delete short_hits[i];
     	short_hits.clear();
     }
     short_hits.push_back(hit);
	 currentName = hit->query_id;
   }

   // trigger processing for last read
   postProcessReadHits(short_hits);
   for (int i=0; i < short_hits.size(); i++) delete short_hits[i];
   short_hits.clear();
}

struct strCmp {
	bool operator()( const char* s1, const char* s2 ) const {
	return strcmp( s1, s2 ) < 0;
	}
};

//
// chain the hits for a given read into larger, more significant hits 
//
void JabonMapperWorker::postProcessReadHits(vector<SoapShortHit*> short_hits)
{
   // group by target 
   map<string, SoapMappingHit> target2mapping_hit;
   for(int hit_idx = 0; hit_idx < short_hits.size(); hit_idx++)
   {
		SoapShortHit * hit = short_hits[hit_idx];
		string key = hit->target_id + hit->target_strand;
		if (target2mapping_hit.count(key) == 0)
		{	
			target2mapping_hit.insert(
				pair<string, SoapMappingHit>(
					key, 
					SoapMappingHit(hit->query_id, hit->target_id)));
		}
		
		target2mapping_hit.find(key)->second.addHit(hit);
   }
   
   map<string,SoapMappingHit>::iterator it;
   for ( it=target2mapping_hit.begin(); it != target2mapping_hit.end(); it++ )
   {
   		SoapMappingHit mapping_hit = (*it).second;
    	SoapHitChain * hit_chain = mapping_hit.maximumChain(param.diagonalCutoff, param.joinCutoff, param.minChainLength);
        if (hit_chain->hits.size() > 0)
    	{
    		cout << hit_chain->toAmosString();
    	}
    	delete(hit_chain);
   }	
}

