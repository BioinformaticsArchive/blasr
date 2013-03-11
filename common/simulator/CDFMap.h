#ifndef SIMULATOR_CDF_MAP_H_
#define SIMULATOR_CDF_MAP_H_
#include "statistics/statutils.h"
#include <algorithm>

template<typename T_Data>
class CDFMap {
 public:
	vector<int> cdf;
	vector<T_Data> data;
	int SelectRandomValue(T_Data &value) {
		vector<int>::iterator search_it;
		int randomIndex = RandomInt(cdf[cdf.size()-1]);
		search_it = lower_bound(cdf.begin(), cdf.end(), randomIndex);
		assert(search_it != cdf.end());
		int valueIndex = search_it - cdf.begin();
		value = data[valueIndex];
		return valueIndex;
	}
};
		
#endif
