#include "reads.h"



using namespace std;



extern Param param;



ReadClass::ReadClass()

{

	_index=0;

	mreads.resize(BatchNum);

}



void ReadClass::CheckFile(ifstream &fin)

{

	string s1,s2,s3,s4;

	char ch[1000];

	fin>>s1;

	fin.getline(ch, 1000);

	fin>>s2;

	fin.getline(ch, 1000);

	

	if('>' == s1[0])

		_file_format=1;

	else if('@' == s1[0]) {

		fin>>s3;		

		fin.getline(ch, 1000);

		fin>>s4;

		fin.getline(ch, 1000);

		_file_format=0;

		if(s2.size() != s4.size()) {

			cerr<<"fatal error: fq format, sequence length not equal to quality length\n";

			exit(1);

		}		

	}

	else {

		cerr<<"fatal error: unrecognizable format of reads file.\n";

		exit(1);

	}

	if(s2.size() < param.min_read_size) {

		cerr<<"fatal error: please set smaller seed size, <="<<2*((s2.size()-4+1)/4)<<endl;

		exit(1);

	}

	fin.seekg(0);

}



void ReadClass::InitialIndex()

{

	_index=0;

}

int ReadClass::LoadBatchReads(ifstream &fin)

{

	char ch[1000];

	char c;

	num=0;

	vector<ReadInf>::iterator p=mreads.begin();

	for(; num<BatchNum; p++,num++,_index++) {

		fin>>c;

		if(fin.eof())

			break;

		p->index=_index;

		fin>>p->name;

		fin.getline(ch,1000);

		fin>>p->seq;

		if(!_file_format) {//*.fq

			fin>>ch;

			fin.getline(ch, 1000);

			fin>>p->qual;

		}

		else

			p->qual=string(p->seq.size(), param.zero_qual+param.default_qual);

	}

	return num;

}



//

// 04/08/08 - JMS - allow mapping long reads by splitting into smaller reads

//

int ReadClass::LoadBatchLongReads( ifstream &fin, ref_id_t longReadUid )

{

    int windowSize = param.longReadWindowSize;

    int windowStep = param.longReadWindowStep;

    //cerr << "Using windowSize=" << windowSize << "; windowStep=" << windowStep << endl;

    // FIXME: This is very fragile and causes seg faults when shorter than actual read length

    int MAX_READ_LENGTH = 1500;

    // we need to reduce the effective value of batch num to control memory

    // and not get segfaults when populating mreads

    int newBatchNum = BatchNum / (MAX_READ_LENGTH / windowStep);

    // FIXME: This is a hard-coded buffer size which means that only query sequences

    // with length<=20kb can be handled by this routine

    char buffer[20000];

    // max line size of FASTA file that we can tolerate is 2000

    char ch[2000];

    char ch2[2] = {0, 0};

    char c;

    int nLongReads = 0;

    vector<ReadInf>::iterator p=mreads.begin();

    string name;

    int nShortReads = 0;

    int luid = longReadUid;

    

    for( ; nLongReads<newBatchNum; luid++, nLongReads++ )

    {

        // read '>'

        fin>>c;

        if(fin.eof())

            break;

        // read sequence name

        fin>>name;

        //cerr << name << endl;

        bool inSequence = 1;

        buffer[0] = 0;

        while( inSequence )

        {

            // eat up trailing white-space, \n

            fin.getline(ch,1000);

            fin>>c;

            if ( fin.eof() )

            {

                inSequence = 0;

            }

            else

            {

              if ( c=='>' )

              {

                  fin.putback( c );

                  inSequence = 0;

              }

              else

              {

                  ch2[0] = c;

                  strcat( buffer, ch2 );

                  if ( fin.peek()!='\n' )

                  {

                      fin>>ch;

                      strcat( buffer, ch );

                  }

              }

            }

        }



        // split long read into shorter reads

        int start = 0;

        int end = windowSize;

        string sequence = string(buffer);

        int seqLen = sequence.size();
        
        id2length.insert(pair<string, int>(name, seqLen)); // AAK
		
        // cout << "Name: " << name << " Length: " << seqLen << endl;

        for( ; end<=seqLen; p++, _index++, nShortReads++ )

        {

            p->index = _index;

            sprintf( ch, "%d:%s", start, name.c_str() );

            p->name = string(ch);

            p->seq = sequence.substr( start, windowSize );

            // TODO: There is no support for quality values at this point

            p->qual = string(p->seq.size(), param.zero_qual+param.default_qual);

            p->longReadUid = luid;

            start += windowStep;

            end += windowStep;

            //cerr << p->name << endl << p->seq << endl;

        }

        //cerr << "next long" << endl;

    }

    // cerr << "finished reading batch of long reads" << endl;

    num = nShortReads;

    return nLongReads;

}

