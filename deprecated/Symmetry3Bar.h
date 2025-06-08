#ifndef SYMMETRY_3BAR_H
#define SYMMETRY_3BAR_H

#include "Config.h"
#include <Eigen/Dense>
#include <vector>
#include <unordered_map>
#include <utility>

using Eigen::Vector3d;
using Eigen::RowVector3d;
using Eigen::Matrix3d;

// Hash function for Vector3d to use in unordered_map
struct Vector3dHash3Bar {
    std::size_t operator()(const Vector3d& v) const {
        return std::hash<double>()(v.x()) ^ std::hash<double>()(v.y()) ^ std::hash<double>()(v.z());
    }
};

/**
 * @brief Computes the closest <111> direction to a given vector, ensuring consistency.
 * @param pairDirection The direction vector of the jump pair.
 * @return The normalized <111> direction vector in a canonical form.
 */
Vector3d getClosest111Direction(const Vector3d& pairDirection);

/**
 * @brief Computes equivalent points under 3-bar symmetry (3-fold rotation + inversion).
 * @param start_point The starting point to transform.
 * @param axis The rotation axis (typically a <111> direction).
 * @param center The inversion center.
 * @return A vector of equivalent points.
 */
std::vector<Vector3d> getEquivalentPoints3Bar(const Vector3d& start_point,
                                             const Vector3d& axis,
                                             const Vector3d& center);

/**
 * @brief Groups lattice sites equivalent under 3-bar symmetry for a jump pair.
 * @param config The lattice configuration.
 * @param maxBondOrder Maximum bond order for neighboring sites.
 * @param latticeIdPair The jump pair (first and second lattice IDs).
 * @return A vector of vectors, each containing indices of equivalent sites.
 */
std::vector<std::vector<size_t>> GetEquivalentSitesUnder3BarSymmetry(
    const Config& config,
    size_t maxBondOrder,
    const std::pair<size_t, size_t>& latticeIdPair);

/**
 * @brief Computes a symmetrically sorted vector of lattice IDs under 3-bar symmetry.
 * @param config The lattice configuration.
 * @param lattice_id_jump_pair The jump pair (first and second lattice IDs).
 * @param max_bond_order Maximum bond order for neighboring sites.
 * @return A sorted vector of lattice IDs, consistent for forward and backward jumps.
 */
std::vector<size_t> GetSortedLatticeVectorStateOfPair3Bar(
    const Config& config,
    const std::pair<size_t, size_t>& lattice_id_jump_pair,
    const size_t& max_bond_order);

#endif // SYMMETRY_3BAR_H