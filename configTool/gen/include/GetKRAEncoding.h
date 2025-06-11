#ifndef CONFIGTOOL_GEN_INCLUDE_GETKRAENCODING_H_
#define CONFIGTOOL_GEN_INCLUDE_GETKRAENCODING_H_

#include <Eigen/Dense>
#include "Config.h"
#include "Symmetry.h"
#include "GetOrbits.h"
#include "LocalEnvironmentEncoder.h"

/**
 * @brief Returns Encoding for Kinetically resolved activation barrier fitting
 *
 * @param config
 * @param maxBondOrder
 * @param latticeIdJumpPair
 * @return Eigen::VectorXd
 */

Eigen::RowVectorXd GetKRAEncoding(
    const Config &config,
    const set<Element> &elementSet,
    const string &basisType,
    const size_t maxBondOrder,
    const size_t maxBondOrderOfCluster,
    const size_t maxClusterSize,
    const pair<size_t, size_t> &latticeIdJumpPair);

#endif // CONFIGTOOL_GEN_INCLUDE_GETKRAENCODING_H_
