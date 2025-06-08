#ifndef LMC_CE_INCLUDE_GETORBITS_H_
#define LMC_CE_INCLUDE_GETORBITS_H_

#include <set>
#include <string>
#include <vector>
#include <unordered_map>
#include "Config.h"
#include "ClusterExpansion.h"

using namespace std;

unordered_map<string, vector<vector<size_t>>> GetOrbits(
    const Config &config,
    const size_t &maxClusterSize,
    const size_t &maxBondOrder,
    const vector<vector<size_t>> &equivalentSiteEncoding,
    const vector<size_t> &symmetricallSortedVector);

#endif // LMC_CE_INCLUDE_GETORBITS_H_
