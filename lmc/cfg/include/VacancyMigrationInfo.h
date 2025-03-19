#ifndef LMC_CFG_VACANCYMIGRATIONINFO_H_
#define LMC_CFG_VACANCYMIGRATIONINFO_H_

#include <unordered_set>
#include "eigen3/Eigen/Dense"
#include <string>
#include "Element.hpp"

struct VacancyMigrationInfo {
 
    /// Unique identifier for the folder storing configuration data.
    size_t uniqueConfigId;

    /// Material system name (e.g., alloy composition).
    std::string materialSystem;

    /// Initial Vacancy Index
    size_t initialVacancyIndex;

    /// Final Vacancy Index
    size_t finalVacancyIndex;

    /// Initial Vacancy Position
    Eigen::RowVector3d initialVacancyPosition;

    /// Final Vacancy Position
    Eigen::RowVector3d finalVacancyPosition;

    /// Migrating Element
    Element migratingAtom;

    /// First Nearest Neighbours of the Migrating Pair
    std::vector<size_t> firstShellNeighbors ;

    /// Second Nearest Neighbours of the Migrating Pair
    std::vector<size_t> secondShellNeighbors;

    /// Third Nearest Neighbours of the Migrating Pair
    std::vector<size_t> thirdShellNeighbors;
  
};

#endif // LMC_CFG_VACANCYMIGRATIONINFO_H_
