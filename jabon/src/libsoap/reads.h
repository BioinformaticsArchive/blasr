#ifndef _READS_H_

#define _READS_H_



#include<iostream>

#include<fstream>

#include<string>

#include<vector>

#include<map>

#include "param.h"



using namespace std;



const int BatchNum=10000;



struct ReadInf

{

	bit32_t index;

	string name;

	string seq;

	string qual;

    ref_id_t longReadUid;

};



class ReadClass

{

public:

	ReadClass();

	void CheckFile(ifstream &fin);

	void InitialIndex();

	int LoadBatchReads(ifstream &fin);

    int LoadBatchLongReads( ifstream &fin, ref_id_t longReadUid );

public:

	vector<ReadInf> mreads;

	bit32_t num;

	map<string, int> id2length;
	

	int _file_format;  //0: fq; 1: fa;

	bit32_t _index;

};



#endif //_READS_H_

