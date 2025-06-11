#include "GetKRAEncoding.h"

Eigen::RowVectorXd GetKRAEncoding(
    const Config &config,
    const set<Element> &elementSet,
    const string &basisType,
    const size_t maxBondOrder,
    const size_t maxBondOrderOfCluster,
    const size_t maxClusterSize,
    const pair<size_t, size_t> &latticeIdJumpPair)
{
  auto equivalentSitesVector = GetEquivalentSitesUnder3BarSymmetry(config,
                                                                   latticeIdJumpPair,
                                                                   maxBondOrder);

  auto orbitMap = GetOrbits(config,
                            maxClusterSize,
                            maxBondOrderOfCluster,
                            equivalentSitesVector);


  Eigen::RowVectorXd lceEncoding = GetLocalEnvironmentEncoding(config,
                                                               elementSet,
                                                               basisType,
                                                               orbitMap);

  return lceEncoding;
}