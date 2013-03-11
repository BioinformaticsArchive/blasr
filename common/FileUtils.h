#ifndef FILE_UTILS_H_
#define FILE_UTILS_H_

#include <iostream>
#include <fstream>
#include <string>

using namespace std;


inline void CriticalOpenRead(string &fileName, ifstream &file, std::ios::openmode mode=std::ios::in) {
	file.open(fileName.c_str(), mode | std::ios::in);
	if (!file.good()) {
		cerr << "Could not open file:"  << fileName << endl;
		exit(1);
	}
}

inline int OpenRead(string &fileName, ifstream &file, std::ios::openmode mode=std::ios::in) {
	file.open(fileName.c_str(), mode | std::ios::in);
	return file.good();
}


inline void CriticalOpenWrite(string &fileName, ofstream &file, std::ios::openmode mode=std::ios::out) {
	file.open(fileName.c_str(), mode | std::ios::out);
	if (!file.good()) { 
		cerr << "Could not open file: " << fileName << endl;
		exit(1);
	}
	
}

inline int OpenWrite(string &fileName, ofstream &file, std::ios::openmode mode=std::ios::out) {
	file.open(fileName.c_str(), mode | std::ios::out);
	return file.good();
}


#endif
