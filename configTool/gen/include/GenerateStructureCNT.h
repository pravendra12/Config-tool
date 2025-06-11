#ifndef CONFIGTOOL_GEN_INCLUDE_GENERATESTRUCTURECNT_H_
#define CONFIGTOOL_GEN_INCLUDE_GENERATESTRUCTURECNT_H_

#include <Eigen/Dense>
#include "Config.h"
#include "ConfigEncoding.h"
#include "PotentialEnergyEstimator.h"
#include "B2Ordering.h"
#include <iostream>
#include <fstream>
#include <filesystem>
#include <random>
#include <string>
#include <sstream>

namespace fs = std::filesystem;

using namespace std;
using namespace Eigen;

class GenerateStructureCNT
{
public:
  GenerateStructureCNT(
      const string &predictorFilename,
      const Config &config,
      const Config &trainingConfig,
      const size_t supercellSize,
      set<Element> &elementSet,
      const vector<double> &compositionVector);

  void GenerateRandomStructures(
      const size_t numStructures,
      const bool doComputeEnergy,
      const bool doSaveConfig,
      const string &outputDir);

  void GenerateStructureWithB2(
      Config &config,
      const size_t numB2Centers,
      const size_t numSamplesPerB2Config,
      const pair<Element, Element> &b2OrderedElements,
      const bool doComputeEnergy,
      const string &outputDir);

private:
  void BuildSingleB2(
      Config &config,
      map<Element, int> &elementMap,
      const size_t selectedLatticeId,
      const Element &alphaElement,
      const Element &betaElement,
      unordered_set<size_t> &visitedSites);

  void BuildConfigWithB2(
      Config &config,
      const size_t numB2Center,
      map<Element, int> elementMap,
      const unordered_set<size_t> &visitedSites,
      int randomId,
      const bool doComputeEnergy,
      const string &outputDir);

  const size_t supercellSize_;
  const vector<string> elementVector_;
  const vector<double> compositionVector_;

  const size_t trainingSupercellSize_ = 5;
  const string structureType_ = "BCC";
  const double latticeParam_ = 3.4;
  const vector<double> cutoffs_ = {3.3, 4.7, 5.6};

  const PotentialEnergyEstimator peEstimator_;
};

#endif // CONFIGTOOL_GEN_INCLUDE_GENERATESTRUCTURECNT_H_
