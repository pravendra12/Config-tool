#ifndef LMC_CFG_CLUSTEREXPANSIONINFO_H_
#define LMC_CFG_CLUSTEREXPANSIONINFO_H_

#include <unordered_set>
#include "eigen3/Eigen/Dense"

/**
 * @brief Stores cluster expansion data for vacancy migration events.
 */
struct ClusterExpansionInfo
{
  /// Unique identifier for the folder storing configuration data.
  size_t uniqueConfigId;

  /// Material system name (e.g., alloy composition).
  std::string materialSystem;

  /// Index of the initial vacancy position.
  size_t initialVacancyIndex;

  /// Coordinates of the initial vacancy in Cartesian space.
  Eigen::Vector3d initialVacancyCoords;

  /// Index of the final vacancy position after migration.
  size_t finalVacancyIndex;

  /// Coordinates of the final vacancy in Cartesian space.
  Eigen::Vector3d finalVacancyCoords;

  /// Chemical element of the migrating atom.
  std::string migratingElement;

  /// Cluster expansion encoding at the initial vacancy position.
  Eigen::VectorXd ceEncodingInitial;

  /// Cluster expansion encoding at the final vacancy position.
  Eigen::VectorXd ceEncodingFinal;
};

#endif //LMC_CFG_CLUSTEREXPANSIONINFO_H_
