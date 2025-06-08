#include "Symmetry3Bar.h"
#include "PrintUtility.h"
#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>
#include <iostream>

using namespace std;

// Reused GetLatticePairCenter3Bar (assumed to be in another file, included here for completeness)
RowVector3d GetLatticePairCenter3Bar(const Config& config, const pair<size_t, size_t>& lattice_id_jump_pair) {
  Vector3d first_relative = config.GetRelativePositionOfLattice(lattice_id_jump_pair.first);
  Vector3d second_relative = config.GetRelativePositionOfLattice(lattice_id_jump_pair.second);
  Vector3d center_position;
  for (int kDim = 0; kDim < 3; ++kDim) {
      double distance = first_relative[kDim] - second_relative[kDim];
      int period = static_cast<int>(distance / 0.5);
      if (period != 0) {
          first_relative[kDim] -= period;
      }
      center_position[kDim] = 0.5 * (first_relative[kDim] + second_relative[kDim]);
  }
  return (center_position.transpose() * config.GetBasis()).transpose();
}

// Constants (assumed defined in constants namespace)
// namespace constants {
//     const double kEpsilon = 1e-6;
// }
// 
// Convert degrees to radians
inline double toRadians(double degrees) {
    return degrees * M_PI / 180.0;
}

// Rotate a point around an axis by angle theta (in degrees) using Rodrigues' formula
Vector3d rotatePoints3Bar(const Vector3d& point, const Vector3d& axis, double theta, const Vector3d& center) {
    Vector3d k = axis.normalized();
    double theta_rad = toRadians(theta);
    double cos_theta = cos(theta_rad);
    double sin_theta = sin(theta_rad);

    Vector3d p = point - center;
    Vector3d p_rot = p * cos_theta +
                     k.cross(p) * sin_theta +
                     k * (k.dot(p)) * (1.0 - cos_theta);
    return p_rot + center;
}

// Find the closest match within a tolerance
size_t findClosestMatch(const unordered_map<Vector3d, size_t, Vector3dHash3Bar>& positionToIndex,
                        const Vector3d& target, double tolerance = constants::kEpsilon) {
    for (const auto& [pos, index] : positionToIndex) {
        if ((pos - target).norm() < tolerance) {
            return index;
        }
    }
    return numeric_limits<size_t>::max();
}

Vector3d getClosest111Direction(const Vector3d& pairDirection) {
    Vector3d normDir = pairDirection.normalized();
    vector<Vector3d> directions = {
        Vector3d(1, 1, 1).normalized(),
        Vector3d(-1, 1, 1).normalized(),
        Vector3d(1, -1, 1).normalized(),
        Vector3d(-1, -1, 1).normalized(),
        Vector3d(1, 1, -1).normalized(),
        Vector3d(-1, 1, -1).normalized(),
        Vector3d(1, -1, -1).normalized(),
        Vector3d(-1, -1, -1).normalized()
    };
    
    Vector3d closestDir = directions[0];
    double maxDot = fabs(normDir.dot(closestDir));
    
    for (const auto& dir : directions) {
        double dot = fabs(normDir.dot(dir));
        if (dot > maxDot + constants::kEpsilon) {
            maxDot = dot;
            closestDir = dir;
        } else if (fabs(dot - maxDot) < constants::kEpsilon) {
            for (int i = 0; i < 3; ++i) {
                if (fabs(dir[i]) > constants::kEpsilon && fabs(closestDir[i]) > constants::kEpsilon) {
                    if (dir[i] > closestDir[i]) {
                        closestDir = dir;
                        maxDot = dot;
                    }
                    break;
                }
            }
        }
    }
    
    return closestDir;
}

vector<Vector3d> getEquivalentPoints3Bar(const Vector3d& start_point,
                                         const Vector3d& axis,
                                         const Vector3d& center) {
    vector<Vector3d> equivalent_points;
    equivalent_points.reserve(6);

    vector<double> angles = {0.0, 120.0, 240.0};

    for (double angle : angles) {
        equivalent_points.push_back(rotatePoints3Bar(start_point, axis, angle, center));
    }

    for (double angle : angles) {
        Vector3d rotated = rotatePoints3Bar(start_point, axis, angle, center);
        equivalent_points.push_back(2.0 * center - rotated);
    }

    return equivalent_points;
}

vector<vector<size_t>> GetEquivalentSitesUnder3BarSymmetry(const Config& config,
                                                          size_t maxBondOrder,
                                                          const pair<size_t, size_t>& latticeIdPair) {
    auto ssVector = GetSortedLatticeVectorStateOfPair3Bar(config, latticeIdPair, maxBondOrder);

    Vector3d centralPos = config.GetCartesianPositionOfLattice(latticeIdPair.first);
    Vector3d nnPos = config.GetCartesianPositionOfLattice(latticeIdPair.second);
    Vector3d pairDirection = centralPos - nnPos;
    Vector3d rotationAxis = getClosest111Direction(pairDirection);

    Vector3d center = GetLatticePairCenter3Bar(config, latticeIdPair);

    unordered_map<Vector3d, size_t, Vector3dHash3Bar> positionToIndex;
    vector<Vector3d> cartesianPositionVector;
    cartesianPositionVector.reserve(ssVector.size());
    for (size_t i = 0; i < ssVector.size(); ++i) {
        auto pos = config.GetCartesianPositionOfLattice(ssVector[i]);
        cartesianPositionVector.push_back(pos);
        positionToIndex[pos] = i;
    }

    vector<bool> processed(ssVector.size(), false);
    vector<vector<size_t>> equivalentEncodingVector;
    for (size_t i = 0; i < ssVector.size(); ++i) {
        if (processed[i]) continue;

        vector<size_t> equivalentIndices = {i};
        processed[i] = true;

        auto equivPositions = getEquivalentPoints3Bar(cartesianPositionVector[i], rotationAxis, center);

        for (const auto& equivPos : equivPositions) {
            size_t index = findClosestMatch(positionToIndex, equivPos, constants::kEpsilon);
            if (index != numeric_limits<size_t>::max() && !processed[index]) {
                equivalentIndices.push_back(index);
                processed[index] = true;
            }
        }

        sort(equivalentIndices.begin(), equivalentIndices.end());
        equivalentEncodingVector.push_back(move(equivalentIndices));
    }

    sort(equivalentEncodingVector.begin(), equivalentEncodingVector.end(),
         [](const auto& a, const auto& b) { return a[0] < b[0]; });

    //cout  << "Equivalent Sites for Pair " << latticeIdPair.first << " -> " << latticeIdPair.second << ": ";
    // print2DVector(equivalentEncodingVector);

    return equivalentEncodingVector;
}

vector<size_t> GetSortedLatticeVectorStateOfPair3Bar(const Config& config,
                                                    const pair<size_t, size_t>& lattice_id_jump_pair,
                                                    const size_t& max_bond_order) 
                                                    {
    auto neighboring_lattice_ids = config.GetNeighboringLatticeIdSetOfPair(lattice_id_jump_pair, max_bond_order);
    size_t num_sites = neighboring_lattice_ids.size();

    Vector3d first_pos = config.GetCartesianPositionOfLattice(lattice_id_jump_pair.first);
    Vector3d second_pos = config.GetCartesianPositionOfLattice(lattice_id_jump_pair.second);
    Vector3d pair_direction = first_pos - second_pos;
    Vector3d rotation_axis = getClosest111Direction(pair_direction);

    RowVector3d center = GetLatticePairCenter3Bar(config, lattice_id_jump_pair);

    Vector3d z_axis(0, 0, 1);
    Matrix3d rotation_matrix;
    if (rotation_axis.dot(z_axis) < 1.0 - constants::kEpsilon) {
        Vector3d cross = rotation_axis.cross(z_axis);
        double angle = acos(rotation_axis.dot(z_axis));
        rotation_matrix = Eigen::AngleAxisd(angle, cross.normalized()).toRotationMatrix();
    } else {
        rotation_matrix.setIdentity();
    }

    vector<tuple<size_t, double, RowVector3d>> lattice_id_keys;
    lattice_id_keys.reserve(num_sites);
    for (const auto id : neighboring_lattice_ids) {
        RowVector3d pos = config.GetCartesianPositionOfLattice(id).transpose();
        pos -= center;
        pos = pos * rotation_matrix;
        double key = pos.norm();
        lattice_id_keys.emplace_back(id, key, pos);
    }

    sort(lattice_id_keys.begin(), lattice_id_keys.end(), [](const auto& lhs, const auto& rhs) {
        double key_diff = get<1>(lhs) - get<1>(rhs);
        if (fabs(key_diff) > constants::kEpsilon) {
            return key_diff < 0;
        }
        const auto& pos_lhs = get<2>(lhs);
        const auto& pos_rhs = get<2>(rhs);
        for (int i = 0; i < 3; ++i) {
            double diff = pos_lhs[i] - pos_rhs[i];
            if (abs(diff) > constants::kEpsilon) {
                return diff < 0;
            }
        }
        return false;
    });

    vector<size_t> sorted_lattice_ids;
    sorted_lattice_ids.reserve(num_sites);
    for (const auto& item : lattice_id_keys) {
        sorted_lattice_ids.push_back(get<0>(item));
    }

    //cout  << "Jump Pair: " << lattice_id_jump_pair.first << " -> " << lattice_id_jump_pair.second << endl;
    //cout  << "Rotation Axis: " << rotation_axis.transpose() << endl;
    //cout  << "Center: " << center << endl;
    //cout  << "Sorted Lattice IDs: "; print1DVector(sorted_lattice_ids);

    return sorted_lattice_ids;
}

