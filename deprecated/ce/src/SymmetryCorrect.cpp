#include "SymmetryCorrect.h"

// Function to find closest match within a tolerance
size_t findClosestMatch(const unordered_map<Vector3d, size_t, Vector3dHash> &positionToIndex,
                        const Vector3d &target, double tolerance = 1e-6)
{
  for (const auto &[pos, index] : positionToIndex)
  {
    if ((pos - target).norm() < tolerance) // Check if within tolerance
    {
      return index; // Return the matching index
    }
  }
  return numeric_limits<size_t>::max(); // Return an invalid index if no match is found
}

// Convert degrees to radians
inline double toRadians(double degrees)
{
  return degrees * M_PI / 180.0;
}

Vector3d rotatePoints(const Vector3d &point,
                      const Vector3d &axis,
                      double theta,
                      const Vector3d &center)
{
  /**
   * Rotate a point around an axis by angle theta (in degrees) using Rodrigues' formula.
   *
   * Args:
   *     point: 3D coordinates of the point to rotate
   *     axis: direction vector of the rotation axis
   *     theta: Angle in degrees
   *     center: center of rotation
   *
   * Returns:
   *     Rotated point as Vector3d
   */

  // Normalize axis
  Vector3d k = axis.normalized();
  double theta_rad = toRadians(theta);
  double cos_theta = cos(theta_rad);
  double sin_theta = sin(theta_rad);

  // Translate point relative to center
  Vector3d p = point - center;

  // Rodrigues' rotation formula: p_rot = p cosθ + (k × p) sinθ + k (k · p) (1 - cosθ)
  Vector3d p_rot = p * cos_theta +
                   k.cross(p) * sin_theta +
                   k * (k.dot(p)) * (1.0 - cos_theta);

  // Translate back
  return p_rot + center;
}
/*
vector<Vector3d> getEquivalentPoints(const Vector3d &start_point,
               const Vector3d &axis,
               double theta,
               const Vector3d &center)
{
int n = 360.0 / theta;

// cout << "number of rotations: " << n << endl;

vector<Vector3d> equivalent_points;
equivalent_points.reserve(n); // Pre-allocate space

for (int i = 0; i < n; ++i)
{
double angle = i * theta;
Vector3d rotated = rotatePoints(start_point, axis, angle, center);
equivalent_points.push_back(rotated);
// cout << rotated.transpose() << endl;
}

return equivalent_points;
}
*/

vector<Vector3d> getEquivalentPoints(const Vector3d &start_point,
                                     const Vector3d &axis,
                                     double theta,
                                     const Vector3d &center,
                                     bool applyInversion)
{
  int n = 360.0 / theta;
  vector<Vector3d> equivalent_points;
  equivalent_points.reserve(n);

  for (int i = 0; i < n; ++i)
  {
    double angle = i * theta;
    Vector3d rotated = rotatePoints(start_point, axis, angle, center);

    if (applyInversion)
    {
      rotated = 2.0 * center - rotated; // Inversion through the center
    }

    equivalent_points.push_back(rotated);
  }

  return equivalent_points;
}

vector<Vector3d> getRotoInvertedPoints(const Vector3d &start_point,
                                       const Vector3d &axis,
                                       double theta,
                                       const Vector3d &center)
{
  bool applyInversion = true;
  auto rotated_points = getEquivalentPoints(start_point, axis, theta, center, applyInversion);
  vector<Vector3d> rotoinverted_points;
  rotoinverted_points.reserve(rotated_points.size());

  for (const auto &p : rotated_points)
  {
    rotoinverted_points.push_back(2.0 * center - p); // Apply inversion
  }

  return rotoinverted_points;
}
vector<vector<size_t>> GetEquivalentSitesUnderKBarSymmetry(const Config &config,
                                                           size_t maxBondOrder,
                                                           size_t kFoldRotation)
{
  vector<vector<size_t>> equivalentEncodingVectorKBar;

  // Get lattice pair (central and nearest neighbor)
  const size_t centralLatticeId = config.GetCentralAtomLatticeId();
  const auto nnLatticeIdVector = config.GetNeighborLatticeIdVectorOfLattice(centralLatticeId, 1);

  const size_t nnLatticeId = nnLatticeIdVector[0];
  const pair<size_t, size_t> latticeIdPair = {centralLatticeId, nnLatticeId};

  // Get symmetrically sorted lattice ID vector for the pair
  const auto ssVector = config.GetSortedLatticeVectorStateOfPair(latticeIdPair, maxBondOrder);

  // Transition position
  const Vector3d centralLatticePosition = config.GetCartesianPositionOfLattice(centralLatticeId);
  const Vector3d nnLatticePosition = config.GetCartesianPositionOfLattice(nnLatticeId);
  const Vector3d transitionPosition = 0.5 * (centralLatticePosition + nnLatticePosition);

  // Rotation axis
  const Vector3d rotationAxis = centralLatticePosition - nnLatticePosition;
  // cout << "Transition Position: " << transitionPosition.transpose() << endl;

  // Pre-populate position vector
  vector<Vector3d> cartesianPositionVector;
  cartesianPositionVector.reserve(ssVector.size());

  for (const auto &latticeId : ssVector)
  {
    cartesianPositionVector.emplace_back(config.GetCartesianPositionOfLattice(latticeId));
  }

  // Use a hash map for O(1) lookup of positions to indices
  unordered_map<Vector3d, size_t, Vector3dHash> positionToIndex;
  for (size_t i = 0; i < cartesianPositionVector.size(); ++i)
  {
    positionToIndex[cartesianPositionVector[i]] = i;
    // cout << cartesianPositionVector[i].transpose() << " : " << i << endl;
  }

  auto equivalentEncodingVectorKFold = GetEquivalentSitesUnderKFoldRotation(config,
                                                                            maxBondOrder,
                                                                            kFoldRotation);
  unordered_map<size_t, size_t> invertedSiteMap;

  for (size_t i = 0; i < equivalentEncodingVectorKFold.size(); ++i)
  {

    // set<size_t> allIndices(equivalentEncodingVector[i].begin(), equivalentEncodingVector[i].end());

    // Invert each position in group i
    for (auto id : equivalentEncodingVectorKFold[i])
    {
      Vector3d pos = config.GetCartesianPositionOfLattice(ssVector[id]); // individual site position
      Vector3d inverted = 2.0 * transitionPosition - pos;

      // Find the closest matching representative position
      size_t matchedId = findClosestMatch(positionToIndex, inverted);

      // cout << id << " " << pos.transpose() << " Inverted " << inverted.transpose() << " " << matchedId << endl;

      invertedSiteMap.insert(make_pair(id, matchedId));
    }
  }

  // combine the inversion and 3 bar encodings

  std::unordered_set<size_t> visited;

  for (const auto &equivalentSites : equivalentEncodingVectorKFold)
  {
    std::vector<size_t> combinedEquivalentSites;
    for (auto siteId : equivalentSites)
    {
      if (visited.count(siteId))
        continue;

      combinedEquivalentSites.push_back(siteId);
      visited.insert(siteId);

      if (invertedSiteMap.count(siteId))
      {
        size_t pairId = invertedSiteMap[siteId];
        if (!visited.count(pairId))
        {
          combinedEquivalentSites.push_back(pairId);
          visited.insert(pairId);
        }
      }
    }

    if (!combinedEquivalentSites.empty())
      equivalentEncodingVectorKBar.push_back(combinedEquivalentSites);
  }

  // print2DVector(equivalentEncodingVectorKBar);

  return equivalentEncodingVectorKBar;
}

vector<vector<size_t>> GetEquivalentSitesUnderKFoldRotation(const Config &config,
                                                            size_t maxBondOrder,
                                                            size_t kFoldRotation)
{
  // Get lattice pair (central and nearest neighbor)
  const size_t centralLatticeId = config.GetCentralAtomLatticeId();
  const auto nnLatticeIdVector = config.GetNeighborLatticeIdVectorOfLattice(centralLatticeId, 1);

  const size_t nnLatticeId = nnLatticeIdVector[0];
  const pair<size_t, size_t> latticeIdPair = {centralLatticeId, nnLatticeId};

  // Get symmetrically sorted lattice ID vector for the pair
  const auto ssVector = config.GetSortedLatticeVectorStateOfPair(latticeIdPair, maxBondOrder);

  // Transition position
  const Vector3d centralLatticePosition = config.GetCartesianPositionOfLattice(centralLatticeId);
  const Vector3d nnLatticePosition = config.GetCartesianPositionOfLattice(nnLatticeId);
  const Vector3d transitionPosition = 0.5 * (centralLatticePosition + nnLatticePosition);

  // Rotation axis
  const Vector3d rotationAxis = centralLatticePosition - nnLatticePosition;
  // cout << "Transition Position: " << transitionPosition.transpose() << endl;

  // Pre-populate position vector
  vector<Vector3d> cartesianPositionVector;
  cartesianPositionVector.reserve(ssVector.size());

  for (const auto &latticeId : ssVector)
  {
    cartesianPositionVector.emplace_back(config.GetCartesianPositionOfLattice(latticeId));
  }

  // Store equivalent sites and their encodings
  vector<vector<size_t>> equivalentEncodingVector;
  vector<Vector3d> equivalentEncodingPositionVector;

  // Use a hash map for O(1) lookup of positions to indices
  unordered_map<Vector3d, size_t, Vector3dHash> positionToIndex;
  for (size_t i = 0; i < cartesianPositionVector.size(); ++i)
  {
    positionToIndex[cartesianPositionVector[i]] = i;
    // cout << cartesianPositionVector[i].transpose() << " : " << i << endl;
  }

  double rotationAngle = 360.0 / kFoldRotation; // 120 degrees for 3-fold

  vector<bool> processed(cartesianPositionVector.size(), false);

  for (size_t i = 0; i < cartesianPositionVector.size(); ++i)
  {
    if (processed[i])
      continue;

    const Vector3d &position = cartesianPositionVector[i];
    auto equivalentPositions = getEquivalentPoints(position, rotationAxis, rotationAngle, transitionPosition, false); // without inversion

    vector<size_t> equivalentIndices = {i}; // Include the original position
    processed[i] = true;

    // Check equivalent positions

    for (const auto &equivPos : equivalentPositions)
    {
      size_t index = findClosestMatch(positionToIndex, equivPos);

      if (index != numeric_limits<size_t>::max()) // If a valid index is found
      {
        if (!processed[index])
        {
          equivalentIndices.push_back(index);
          processed[index] = true;
        }
      }
    }

    equivalentEncodingVector.push_back(move(equivalentIndices));
    equivalentEncodingPositionVector.push_back(position);
  }

  // print2DVector(equivalentEncodingVector);

  return equivalentEncodingVector;
}

inline bool PositionCompareState(
    const pair<size_t, RowVector3d> &lhs,
    const pair<size_t, RowVector3d> &rhs)
{
  const auto &relative_position_lhs = lhs.second;
  const auto &relative_position_rhs = rhs.second;

  // Compare individual components (x, y, z)
  const double diff_x = relative_position_lhs[0] - relative_position_rhs[0];
  if (diff_x < -constants::kEpsilon)
  {
    return true;
  }
  if (diff_x > constants::kEpsilon)
  {
    return false;
  }

  const double diff_y = relative_position_lhs[1] - relative_position_rhs[1];
  if (diff_y < -constants::kEpsilon)
  {
    return true;
  }
  if (diff_y > constants::kEpsilon)
  {
    return false;
  }

  const double diff_z = relative_position_lhs[2] - relative_position_rhs[2];
  if (diff_z < -constants::kEpsilon)
  {
    return true;
  }
  if (diff_z > constants::kEpsilon)
  {
    return false;
  }

  return false;
}

void RotateLatticeVector(
    unordered_map<size_t, RowVector3d> &lattice_id_hashmap,
    const Matrix3d &rotation_matrix)
{
  const RowVector3d
      move_distance_after_rotation = RowVector3d(0.5, 0.5, 0.5) -
                                     (RowVector3d(0.5, 0.5, 0.5) * rotation_matrix);

  for (auto &lattice : lattice_id_hashmap)
  {

    auto relative_position = lattice.second;
    // rotate
    relative_position = relative_position * rotation_matrix;
    // move to new center
    relative_position += move_distance_after_rotation;
    relative_position -= relative_position.unaryExpr([](double x)
                                                     { return floor(x); });

    lattice.second = relative_position;
  }
}

vector<size_t> GetSortedLatticeVectorStateOfPair(
    const Config &config,
    const pair<size_t, size_t> &lattice_id_jump_pair,
    const size_t &max_bond_order)
{
  auto neighboring_lattice_ids =
      config.GetNeighboringLatticeIdSetOfPair(lattice_id_jump_pair,
                                              max_bond_order);

  size_t num_sites = neighboring_lattice_ids.size();

  RowVector3d move_distance = RowVector3d(0.5, 0.5, 0.5) -
                              GetLatticePairCenter(config, lattice_id_jump_pair);

  unordered_map<size_t, RowVector3d> lattice_id_hashmap;
  lattice_id_hashmap.reserve(num_sites);

  // Move lattice IDs to center
  for (const auto id : neighboring_lattice_ids)
  {
    RowVector3d relative_position =
        config.GetRelativePositionOfLattice(id).transpose();

    relative_position += move_distance;
    relative_position -= relative_position.unaryExpr([](double x)
                                                     { return floor(x); });
    lattice_id_hashmap.emplace(id, relative_position);
  }

  // Rotate lattice vectors
  auto rotationMatrix = GetLatticePairRotationMatrix(config,
                                                     lattice_id_jump_pair);
  RotateLatticeVector(lattice_id_hashmap, rotationMatrix);

  // Convert unordered_map to vector for sorting
  vector<pair<size_t, RowVector3d>>
      lattice_id_vector(lattice_id_hashmap.begin(), lattice_id_hashmap.end());

  // Sort the lattice vector based on PositionCompareMMM
  sort(lattice_id_vector.begin(), lattice_id_vector.end(), PositionCompareState);

  // Extract and return only the lattice IDs
  vector<size_t> sorted_lattice_ids;
  sorted_lattice_ids.reserve(lattice_id_vector.size());
  for (const auto &pair : lattice_id_vector)
  {
    sorted_lattice_ids.push_back(pair.first);
  }

  return sorted_lattice_ids;
}

RowVector3d GetLatticePairCenter(
    const Config &config,
    const pair<size_t, size_t> &lattice_id_jump_pair)
{
  // Get relative positions of the lattice IDs
  Vector3d first_relative =
      config.GetRelativePositionOfLattice(lattice_id_jump_pair.first);

  Vector3d second_relative =
      config.GetRelativePositionOfLattice(lattice_id_jump_pair.second);

  Vector3d center_position;
  for (int kDim = 0; kDim < 3; ++kDim)
  {
    double distance = first_relative[kDim] - second_relative[kDim];
    int period = static_cast<int>(distance / 0.5);

    // Adjust the positions once based on the period
    if (period != 0)
    {
      first_relative[kDim] -= period;
    }

    center_position[kDim] = 0.5 * (first_relative[kDim] + second_relative[kDim]);
  }

  return center_position.transpose();
}

Matrix3d GetLatticePairRotationMatrix(
    const Config &config,
    const pair<size_t, size_t> &lattice_id_jump_pair)
{
  // Get Pair Direction
  Vector3d
      relative_distance_vector_pair =
          config.GetRelativeDistanceVectorLattice(lattice_id_jump_pair.first,
                                                  lattice_id_jump_pair.second);

  const Matrix3d basis = config.GetBasis();

  RowVector3d
      pair_direction = relative_distance_vector_pair.transpose() * basis;
  pair_direction.normalize();

  RowVector3d vertical_vector;

  // First nearest neighbors
  vector<size_t> nn_list =
      config.GetNeighborLatticeIdVectorOfLattice(lattice_id_jump_pair.first, 1);

  sort(nn_list.begin(), nn_list.end());

  for (const auto nn_id : nn_list)
  {
    // Get Jump Vector
    Vector3d
        relative_distance_vector =
            config.GetRelativeDistanceVectorLattice(lattice_id_jump_pair.first, nn_id);

    RowVector3d jump_vector = relative_distance_vector.transpose() * basis;
    jump_vector.normalize();

    // Check for vertical direction
    const double dot_prod = pair_direction.dot(jump_vector);
    if (abs(dot_prod) < constants::kEpsilon)
    {
      vertical_vector = jump_vector;
      break;
    }
  }

  // Rotation Matrix
  Matrix3d rotation_matrix;
  rotation_matrix.row(0) = pair_direction;
  rotation_matrix.row(1) = vertical_vector;
  rotation_matrix.row(2) = pair_direction.cross(vertical_vector);

  return rotation_matrix.transpose();
}

vector<vector<size_t>> GetEquivalentSites3Fold(const Config &config,
                                               const size_t maxBondOrder)
{
  // Does not depends on jump pairs.
  pair<size_t, size_t> jumpPair = {0, config.GetNeighborLatticeIdVectorOfLattice(0, 1)[0]};

  auto symmetricallySortedVector = GetSortedLatticeVectorStateOfPair(config,
                                                                     jumpPair,
                                                                     maxBondOrder);

  vector<size_t> equivalentSites;
  vector<vector<size_t>> equivalentSiteVector;

  size_t idx = 0;

  while (idx < symmetricallySortedVector.size())
  {
    auto latticeId = symmetricallySortedVector[idx];

    auto directionVectorI =
        config.GetNormalizedDirection(latticeId, jumpPair.first);
    auto directionVectorII =
        config.GetNormalizedDirection(latticeId, jumpPair.second);

    if (directionVectorI == directionVectorII)
    {
      equivalentSites.emplace_back(idx);
      idx++;
    }
    else
    {
      // Sample in threes
      // Works only till 2nd NN
      for (size_t id = idx; id < idx + 3; id++)
      {
        equivalentSites.emplace_back(id);
      }
      idx += 3;
    }

    equivalentSiteVector.emplace_back(equivalentSites);
    equivalentSites = {};
  }

  return equivalentSiteVector;
}

/*
// Pairs between Migrating Atom and Nearest Neigbours
VectorXd GetEncodingMigratingAtomPair (
  const Config &config,
  const vector<vector<size_t>> &equivalentSitesEncoding,
  const vector<size_t> &symmetricallySortedVector,
  const unordered_map<string, RowVectorXd> &oneHotEncodingMap,
  const Element &migratingAtom)
{

// maxBondOrder is used for getting the NN for the jumpPair
// Currently only support till maxBondOrder = 2

// for bondOrder 1, there are 6 equivalent groups
// the two central groups will be the first nn of the migrating atom (6 atoms)
// then the next to each of the central groups will be the second nn atom
// for the migrating atom at transition site
// and third nn for the migrating atom based on the cutoff described
// so in the first nn shell of the jump pair, there will be 14 nn of the migrating Atom

  // Pair of Migrating Atom and the NN
  // Not required storing the element pairs for verification
  // vector<vector<string>> equivalentPairsClusters.reserve();

  size_t sizeEncodeNonSymmPairs = 4;

  VectorXd encodeVector;

  for (auto &equivalentSites : equivalentSitesEncoding)
  {
    vector<string> equivalentPairs;

    // Since non symmetric pairs
    // Migrating Atom will be first site
    // Need to work on this so as to generalize this

    RowVectorXd pairEncodeVector = RowVectorXd::Zero(Index(sizeEncodeNonSymmPairs));

    for (auto &sites : equivalentSites)
    {
      auto latticeId = symmetricallySortedVector[sites];
      auto neighborElement = config.GetElementOfLattice(latticeId);

      string elementPair = migratingAtom.GetElementString() +
                                neighborElement.GetElementString();

      equivalentPairs.emplace_back(elementPair);

      RowVectorXd oneHotVector = oneHotEncodingMap.at(elementPair);

      pairEncodeVector += oneHotVector;
    }

  // equivalentPairsClusters.emplace_back(equivalentPairs);

  encodeVector.conservativeResize(encodeVector.size() + pairEncodeVector.size());

  // Concatenate pairEncodeVector to the resized vector
  encodeVector.tail(pairEncodeVector.size()) = pairEncodeVector;

  }

  // print2DStringVector(equivalentPairsClusters);

  return encodeVector;
}

*/

VectorXd GetEncodingMigratingAtomPair(
    const Config &config,
    const vector<vector<size_t>> &equivalentSitesEncoding,
    const vector<size_t> &symmetricallySortedVector,
    const unordered_map<string, RowVectorXd> &oneHotEncodingMap,
    const Element &migratingAtom)
{
  // Size of non-symmetric pairs
  size_t numElements = 2; // Considering only binary for now
  size_t sizeEncodeNonSymmPairs = pow(numElements, 2);

  // Calculate the total size needed for encodeVector
  size_t totalSize;
  totalSize = equivalentSitesEncoding.size() * sizeEncodeNonSymmPairs;

  // Pre-allocate encodeVector with the total expected size
  VectorXd encodeVector(totalSize);
  size_t offset = 0;

  // Loop through the equivalent sites encoding
  for (const auto &equivalentSites : equivalentSitesEncoding)
  {
    RowVectorXd pairEncodeVector = RowVectorXd::Zero(sizeEncodeNonSymmPairs);

    // Loop through the equivalent sites and build the pair encoding vector
    for (auto &sites : equivalentSites)
    {
      auto latticeId = symmetricallySortedVector[sites];
      auto neighborElement = config.GetElementOfLattice(latticeId);

      string elementPair = migratingAtom.GetElementString() +
                           neighborElement.GetElementString();

      try
      {
        pairEncodeVector += oneHotEncodingMap.at(elementPair);
      }
      catch (const out_of_range &e)
      {
        cout << "Error: Missing Element Pair for " << elementPair << "(" << latticeId << ")" << endl;
        exit(1);
      }
    }
    // Normalized encoding vector
    // Store the result into encodeVector at the correct offset
    encodeVector.segment(offset, pairEncodeVector.size()) = pairEncodeVector / equivalentSites.size();
    offset += pairEncodeVector.size(); // Move the offset for next pair
  }

  // Return the final accumulated encode vector
  return encodeVector;
}
