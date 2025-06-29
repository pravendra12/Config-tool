#ifndef LMC_CE_INCLUDE_SYMMETRYCORRECT_H_
#define LMC_CE_INCLUDE_SYMMETRYCORRECT_H_

#include "Config.h"
#include "eigen3/Eigen/Dense"
#include "Constants.hpp"
#include "EncodingUtility.h"
#include "PrintUtility.h"
#include <boost/functional/hash.hpp>
#include <vector>
#include <cmath>

using namespace std;
using namespace Eigen;

/**
 * Generate all equivalent points under n-fold symmetry around an axis.
 *
 * Args:
 *     start_point: Starting point coordinates
 *     axis: Rotation axis
 *     theta: Angle in degrees between successive rotations in degree
 *     center: Center of rotation (default: nullptr)
 *     tolerance: Distance tolerance for matching points (default: 1e-2)
 *
 * Returns:
 *     Vector of equivalent points
 */

struct Vector3dHash
{
    std::size_t operator()(const Vector3d &v) const
    {
        std::size_t h1 = std::hash<int>{}(static_cast<int>(v.x() * 1e6));
        std::size_t h2 = std::hash<int>{}(static_cast<int>(v.y() * 1e6));
        std::size_t h3 = std::hash<int>{}(static_cast<int>(v.z() * 1e6));
        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};

size_t findClosestMatch(const unordered_map<Vector3d, size_t, Vector3dHash> &positionToIndex, const Vector3d &target, double tolerance);

double toRadians(double degrees);

std::vector<Vector3d> getEquivalentPoints(const Vector3d &start_point,
                                          const Vector3d &axis,
                                          double theta,
                                          const Vector3d &center,
                                          bool applyInversion);

vector<vector<size_t>> GetEquivalentSitesUnderKBarSymmetry(
    const Config &config, 
    size_t maxBondOrder, 
    size_t kFoldRotation);

std::vector<std::vector<size_t>> GetEquivalentSitesUnderKFoldRotation(const Config &config,
                                                                      size_t maxBondOrder,
                                                                      size_t kFoldRotation);

bool PositionCompareState(const pair<size_t, RowVector3d> &lhs,
                          const pair<size_t, RowVector3d> &rhs);

void RotateLatticeVector(
    unordered_map<size_t, RowVector3d> &lattice_id_hashmap,
    const Matrix3d &rotation_matrix);

/*! \brief Computes the center position of a lattice pair while accounting for
             periodic boundary conditions.
   *
   * This function calculates the geometric center of two lattice points specified
   * by their IDs. It adjusts the relative positions of the lattice points to
   * ensure that the computed distance between them is within the range (0, 0.5)
   * in each dimension, considering periodic boundaries.
   *
   *  \param config               A reference to the Config object containing
   *                              lattice configurations and relative positions.
   *  \param lattice_id_jump_pair A pair of lattice IDs representing the two
   *                              lattice points.
   *  \return                     Vector3d The computed center position of
   *                              the lattice pair in three dimensions.
   */
RowVector3d GetLatticePairCenter(
    const Config &config,
    const pair<size_t, size_t> &lattice_id_jump_pair);

Matrix3d GetLatticePairRotationMatrix(
    const Config &config,
    const pair<size_t, size_t> &lattice_id_jump_pair);

vector<size_t> GetSortedLatticeVectorStateOfPair(
    const Config &config,
    const pair<size_t, size_t> &lattice_id_jump_pair,
    const size_t &max_bond_order);

/*! \brief Finds equivalent sites under a 3-fold rotation along the 111 direction.
 *  \param config Reference to the configuration.
 *  \param maxBondOrder Specifies the maximum nearest neighbor shell to consider
 *                     (currently, only maxBondOrder = 2 is supported).
 *  \return A list of equivalent sites encoding under the 3-fold rotation along
 *          the 111 direction.
 *
 * Works only till 2nd NN
 *
 */
vector<vector<size_t>>
GetEquivalentSites3Fold(const Config &config,
                        const size_t maxBondOrder);

// Pairs between Migrating Atom and Nearest Neigbours
VectorXd GetEncodingMigratingAtomPair(
    const Config &config,
    const vector<vector<size_t>> &equivalentSitesEncoding,
    const vector<size_t> &symmetricallySortedVector,
    const unordered_map<string, RowVectorXd> &oneHotEncodingMap,
    const Element &migratingAtom);

#endif // LMC_CE_INCLUDE_SYMMETRY_H_
