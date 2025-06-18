#ifndef CONFIGTOOL_GEN_INCLUDE_GENERATENEBSTRUCTURE_H_
#define CONFIGTOOL_GEN_INCLUDE_GENERATENEBSTRUCTURE_H_

#include "Config.h"

using namespace std;

// Post-processing required to make the files compatible with LAMMPS
void GenerateNEBStructure(
  const string &filename, 
  Config &config, 
  const pair<size_t, size_t> &latticeIdJumpPair, 
  map<Element, size_t> &elementMap);

#endif // CONFIGTOOL_GEN_INCLUDE_GENERATENEBSTRUCTURE_H_

