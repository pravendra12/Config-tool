#ifndef CONFIGTOOL_GEN_INCLUDE_GENERATESTRUCTURESWITHB2_H_
#define CONFIGTOOL_GEN_INCLUDE_GENERATESTRUCTURESWITHB2_H_

#include "Config.h"
#include <unordered_set>

using namespace std;

void GenerateStructureWithB2(Config &config,
                             size_t numB2Centers,
                             pair<Element, Element> &b2ElementPair);

#endif // CONFIGTOOL_GEN_INCLUDE_GENERATESTRUCTURESWITHB2_H_
