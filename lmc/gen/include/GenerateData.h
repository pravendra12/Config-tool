#ifndef LMC_GEN_INCLUDE_GENERATEDATA_H_
#define LMC_GEN_INCLUDE_GENERATEDATA_H_

#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <filesystem>
#include <sstream>
#include <unordered_map>
#include <eigen3/Eigen/Dense>
#include "Config.h"
#include "Element.hpp"
#include "VacancyMigrationInfo.h"
#include "ClusterExpansion.h"
#include "ClusterExpansionInfo.h"
#include "SymmetryCorrect.h"
#include "B2Ordering.h"
#include <random>

using namespace std;
namespace fs = std::filesystem;

class GenerateData
{
public:
  GenerateData(
      const Config config,
      const size_t supercellSize,
      const vector<string> &elementVector,
      const vector<double> &elementComposition,
      const vector<double> &cutoffs,
      const size_t uniqueConfigId,
      const bool isOrderedStructure);

  GenerateData(const size_t supercellSize,
               const double latticeParam,
               const string structureType,
               const vector<string> &elementVector,
               const vector<double> &elementComposition,
               const vector<double> &cutoffs,
               const size_t uniqueConfigId,
               const bool isOrderedStructure);

  bool saveConfig();

  VacancyMigrationInfo generateNEBStructure() const;

  Eigen::VectorXd getEncodeVector(const Config &config,
                                  const set<Element> &elementSet,
                                  const size_t maxClusterSize,
                                  const size_t maxBondOrder) const;

  // Takes element set such it also includes vacancy if it there in config.
  ClusterExpansionInfo generateCEData(const set<Element> &elementSet,
                                      const size_t maxClusterSize,
                                      const size_t maxBondOrder);

  // Returns a pair of <forwardEncoding, backwardEncoding>
  pair<VectorXd, VectorXd> getMigratingAtomNeighborPairEncodeVector(
      const set<Element> &elementSet);

  void writeVacancyJumpCeConfig(ofstream &outFile,
                                const VacancyMigrationInfo &vacancyMigrationInfo,
                                const ClusterExpansionInfo &ceInfo,
                                const B2OrderingInfo &b2OrderInfo,
                                const set<Element> &elementSet);

  /**
   * @brief Generates a partially B2-ordered structure by assigning elements to a fraction of B2 sites.
   *
   * @param fraction A value between 0 and 1 specifying the fraction of B2 sites to assign elements.
   *                 - A value of 1.0 results in a fully ordered B2 structure.
   *                 - A value of 0.5 means only half of the identified B2 sites will be assigned.
   *
   * @return Config The modified configuration with partial B2 ordering.
   *
   * @warning This function selects a subset of alpha and beta sites at random.
   *          - Since `alphaSites` and `betaSites` are unordered sets, they are converted to vectors for random selection.
   *          - If the fraction is too low, the resulting structure may have significant disorder.
   *          - If strict nearest-neighbor pairing is required for B2 ordering, this method may disrupt local ordering.
   */
  B2OrderingInfo generateB2Structure(vector<double> elementComposition,
                                     const Element elementAtAlphaSite,
                                     const Element elementAtBetaSite);

  B2OrderingInfo AddB2Structure(size_t numB2Center, 
                                vector<Element> &elementSet);


private:
  void writeNeighbourVectorToFile(ofstream &outputFile,
                                  const vector<size_t> &vec) const;

  const string alloySystem_;
  const size_t supercellSize_;

  const vector<double> cutoffs_;

  const std::string configFilename_;

  const size_t uniqueConfigId_;

  const bool isOrderedStructure_;

  mutable Config config_{};

  // Always vacancy will be at central atom
  mutable size_t vacancyId_{};

  // Migrating atom id
  mutable size_t migratingElementLatticeId_{};

  // Migrating atom element
  mutable Element migratingAtomElement_{};

  // Forward Jump Pair
  mutable pair<size_t, size_t> forwardJumpPair_{};

  // Backward Jump Pair
  mutable pair<size_t, size_t> backwardJumpPair_{};
};



#endif // LMC_GEN_INCLUDE_GENERATEDATA_H_

