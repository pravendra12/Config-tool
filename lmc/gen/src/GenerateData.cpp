#include "GenerateData.h"

static string getConfigFileName(const size_t supercellSize,
                                const vector<string> &elementSet,
                                const vector<double> &elementComposition)
{
  if (elementSet.size() != elementComposition.size())
  {
    throw runtime_error("Size mismatch between elementSet and elementComposition.");
  }

  string fileName;

  for (size_t i = 0; i < elementSet.size(); i++)
  {
    fileName += elementSet[i];
    fileName += to_string(int(elementComposition[i]));
  }

  ostringstream oss;
  oss << supercellSize << "x"
      << supercellSize << "x"
      << supercellSize;

  fileName += "_" + oss.str();

  return fileName;
}

static string generateFolder(const string &alloySystem,
                             const size_t uniqueConfigID,
                             const bool isOrderedStructure)
{
  // Create the main alloy system directory (e.g., "Ta80W20")
  fs::path alloySystemDir = alloySystem;

  try
  {
    // Check if the alloy system directory exists
    if (!fs::exists(alloySystemDir))
    {
      // Create the alloy system directory
      if (fs::create_directory(alloySystemDir))
      {
        cout << "Directory created for alloy system: " << alloySystemDir << endl;
      }
      else
      {
        cerr << "Error: Failed to create alloy system directory: " << alloySystemDir << endl;
      }
    }
  }
  catch (const fs::filesystem_error &e)
  {
    cerr << "Filesystem error while creating alloy system directory: " << e.what() << endl;
  }

  // Now, create the unique config ID folder inside the alloy system folder (e.g., "Ta80W20/01")
  ostringstream configDirName;
  configDirName << setw(2) << setfill('0') << uniqueConfigID;

  if (isOrderedStructure)
  {
    configDirName << "_O";
  }

  fs::path uniqueConfigDir = alloySystemDir / configDirName.str(); // Combine the paths

  try
  {
    // Check if the unique config directory exists
    if (!fs::exists(uniqueConfigDir))
    {
      // Create the unique config directory
      if (fs::create_directory(uniqueConfigDir))
      {
        cout << "Directory created: " << uniqueConfigDir << endl;
      }
      else
      {
        cerr << "Error: Failed to create unique config directory: " << uniqueConfigDir << endl;
      }
    }
  }
  catch (const fs::filesystem_error &e)
  {
    cerr << "Filesystem error while creating unique config directory: " << e.what() << endl;
  }

  return uniqueConfigDir.string();
}

// Convert cluster set to a map with the number of appearance of each cluster type.
static unordered_map<ClusterType, size_t, boost::hash<ClusterType>> ConvertSetToHashMap(
    const set<ClusterType> &cluster_type_set)
{
  unordered_map<ClusterType, size_t, boost::hash<ClusterType>> cluster_type_count;
  for (const auto &clusterType : cluster_type_set)
  {
    cluster_type_count[clusterType] = 0;
  }
  return cluster_type_count;
}

static string getAlloySystem(const vector<string> &elementSet,
                             const vector<double> &elementComposition)
{
  string alloySystem;

  for (size_t i = 0; i < elementSet.size(); ++i)
  {
    // Convert the composition to an integer and append it to the element symbol
    // A20B80
    alloySystem += elementSet[i] +
                   to_string(static_cast<int>(elementComposition[i]));
  }

  return alloySystem;
}

// Constructor 1
GenerateData::GenerateData(
    const Config config,
    const size_t supercellSize,
    const vector<string> &elementVector,
    const vector<double> &elementComposition,
    const vector<double> &cutoffs,
    const size_t uniqueConfigId,
    const bool isOrderedStructure) : alloySystem_(getAlloySystem(elementVector,
                                                                 elementComposition)),
                                     supercellSize_(supercellSize),
                                     cutoffs_(cutoffs),
                                     configFilename_(getConfigFileName(supercellSize,
                                                                       elementVector,
                                                                       elementComposition)),
                                     uniqueConfigId_(uniqueConfigId),
                                     isOrderedStructure_(isOrderedStructure),
                                     config_(move(config))

{
  vacancyId_ = config_.GetCentralAtomLatticeId();
  migratingElementLatticeId_ = config_.GetNeighborLatticeIdVectorOfLattice(vacancyId_, 1)[0];

  migratingAtomElement_ = config_.GetElementOfLattice(migratingElementLatticeId_);

  // Need to think about this
  forwardJumpPair_ = {migratingElementLatticeId_, vacancyId_};

  backwardJumpPair_ = {vacancyId_, migratingElementLatticeId_};
}

GenerateData::GenerateData(
    const size_t supercellSize,
    const double latticeParam,
    const string structureType,
    const vector<string> &elementVector,
    const vector<double> &elementComposition,
    const vector<double> &cutoffs,
    const size_t uniqueConfigId,
    const bool isOrderedStructure) : alloySystem_(getAlloySystem(elementVector,
                                                                 elementComposition)),
                                     supercellSize_(supercellSize),
                                     cutoffs_(cutoffs),
                                     configFilename_(getConfigFileName(supercellSize,
                                                                       elementVector,
                                                                       elementComposition)),
                                     uniqueConfigId_(uniqueConfigId),
                                     isOrderedStructure_(isOrderedStructure)

{
  config_ = Config::GenerateAlloySupercell(supercellSize,
                                           latticeParam,
                                           structureType,
                                           elementVector,
                                           elementComposition,
                                           uniqueConfigId);
  config_.UpdateNeighborList(cutoffs);

  vacancyId_ = config_.GetCentralAtomLatticeId();
  migratingElementLatticeId_ = config_.GetNeighborLatticeIdVectorOfLattice(vacancyId_, 1)[0];

  migratingAtomElement_ = config_.GetElementOfLattice(migratingElementLatticeId_);

  // Need to think about this
  forwardJumpPair_ = {migratingElementLatticeId_, vacancyId_};

  backwardJumpPair_ = {vacancyId_, migratingElementLatticeId_};
}

bool GenerateData::saveConfig()
{
  bool isConfigSaved;

  string configDir;

  configDir = generateFolder(alloySystem_, uniqueConfigId_, isOrderedStructure_) + "/Config/";

  try
  {
    // Check if the unique config directory exists
    if (!fs::exists(configDir))
    {
      // Create the unique config directory
      if (fs::create_directory(configDir))
      {
        cout << "Directory created: " << configDir << endl;
      }
      else
      {
        cerr << "Error: Failed to create unique config directory: " << configDir << endl;
      }
    }
  }
  catch (const fs::filesystem_error &e)
  {
    cerr << "Filesystem error while creating unique config directory: " << e.what() << endl;
  }

  try
  {
    // Lammps data file
    Config::WriteLAMMPSDataFile(configDir + configFilename_ + ".data", config_);
    // Cfg file
    Config::WriteConfig(configDir + configFilename_ + ".cfg", config_);

    isConfigSaved = true;
  }
  catch (const exception &e)
  {
    cerr << "Error: Not able to save the configurations" << '\n';
    cerr << e.what() << '\n';

    isConfigSaved = false;
  }

  return isConfigSaved;
}

VacancyMigrationInfo GenerateData::generateNEBStructure() const
{

  string uniqueOutputDirectory;

  uniqueOutputDirectory = generateFolder(alloySystem_, uniqueConfigId_, isOrderedStructure_);

  // Temporary configuraiton
  Config config = config_;
  config.UpdateNeighborList(cutoffs_);

  auto configRelativePositionMatrix = config.GetRelativePositionMatrix();
  auto numAtoms = config.GetNumAtoms();
  auto atomVector = config.GetAtomVector();
  auto configBasis = config.GetBasis();

  // Initially vacancy will be at central site in the supercell
  // Final lattice id of the migrating atom
  auto initVacancyId = vacancyId_;
  auto initVacancyRelativePosition = config.GetRelativePositionOfLattice(initVacancyId);

  // Finally vacancy will move to the neighbouring site
  // Initial lattice id of the migrating atom
  auto finalVacancyId = migratingElementLatticeId_;
  auto finalVacancyRelativePosition = config.GetRelativePositionOfLattice(finalVacancyId);

  auto firstNNOfJumpPair = config.GetSortedLatticeVectorStateOfPair({initVacancyId, finalVacancyId}, 1);
  auto secondNNOfJumpPair = config.GetSortedLatticeVectorStateOfPair({initVacancyId, finalVacancyId}, 2);
  auto thirdNNOfJumpPair = config.GetSortedLatticeVectorStateOfPair({initVacancyId, finalVacancyId}, 3);

  // Initial NEB Structure
  Eigen::Matrix3Xd relativePositionMatrix(3, numAtoms - 1);

  relativePositionMatrix << configRelativePositionMatrix.leftCols(initVacancyId),
      configRelativePositionMatrix.rightCols(configRelativePositionMatrix.cols() - initVacancyId - 1);

  auto initAtomVector = atomVector;

  initAtomVector.erase(initAtomVector.begin() + initVacancyId);

  Config initialConfig = Config{configBasis,
                                relativePositionMatrix,
                                initAtomVector};
  // initialConfig.ReassignLattice();
  // initialConfig.Wrap();

  Config::WriteLAMMPSDataFile(uniqueOutputDirectory + "/" + configFilename_ + "_initial.data", initialConfig);

  // Final NEB Structure

  double tolerance = 1e-6;
  for (int i = 0; i < relativePositionMatrix.cols(); ++i)
  {
    if ((relativePositionMatrix.col(i) - finalVacancyRelativePosition).norm() < tolerance)
    {
      relativePositionMatrix.col(i) = initVacancyRelativePosition; // Move the neighbouring atom to the vacancy
      break;
    }
  }

  // Bug : Using same atom vector for initial and final config
  // Possible issues are that the migrating element may not be same
  // Need to work on this
  Config finalConfig = Config{configBasis,
                              relativePositionMatrix,
                              initAtomVector};
  // finalConfig.ReassignLattice();
  // finalConfig.Wrap();

  Config::WriteLAMMPSDataFile(uniqueOutputDirectory + "/" + configFilename_ + "_final.data", finalConfig);

  VacancyMigrationInfo vacancyMigrationInfo(uniqueConfigId_,
                                            alloySystem_,
                                            initVacancyId,
                                            finalVacancyId,
                                            config.GetCartesianPositionOfLattice(initVacancyId),
                                            config.GetCartesianPositionOfLattice(finalVacancyId),
                                            migratingAtomElement_,
                                            firstNNOfJumpPair,
                                            secondNNOfJumpPair,
                                            thirdNNOfJumpPair);

  return vacancyMigrationInfo;
}

Eigen::VectorXd GenerateData::getEncodeVector(const Config &config,
                                              const set<Element> &elementSet,
                                              const size_t maxClusterSize,
                                              const size_t maxBondOrder) const
{
  auto initializedClusterTypeSet(InitializeClusterTypeSet(config,
                                                          elementSet,
                                                          maxClusterSize,
                                                          maxBondOrder));

  // This function is responsible for generating training data.
  // Since the supercell configuration is used for training the CE model,
  // it will remain the same as the configuration to ensure cluster occupancy is preserved.
  auto latticeClusterTypeCount(CountLatticeClusterTypes(config,
                                                        maxClusterSize,
                                                        maxBondOrder));

  auto clusterTypeCountHashMap(ConvertSetToHashMap(initializedClusterTypeSet));

  auto allLatticeClusterHashSet = FindAllLatticeClusters(config,
                                                         maxClusterSize,
                                                         maxBondOrder,
                                                         {});

  for (const auto &latticeCluster : allLatticeClusterHashSet)
  {
    auto atomClusterType = IndentifyAtomClusterType(config,
                                                    latticeCluster.GetLatticeIdVector());
    clusterTypeCountHashMap.at(ClusterType(atomClusterType, latticeCluster.GetClusterType()))++;
  }
  Eigen::VectorXd encodeVector(initializedClusterTypeSet.size());
  int idx = 0;
  for (const auto &clusterType : initializedClusterTypeSet)
  {
    // Count of Cluster Types for given configuration
    auto bondCount = static_cast<double>(clusterTypeCountHashMap.at(clusterType));
    // Count of Cluster Types for normalization
    auto totalBondCount = static_cast<double>(latticeClusterTypeCount.at(clusterType.lattice_cluster_type_));

    encodeVector(idx) = bondCount / totalBondCount;

    // cout << clusterType << " : " << bondCount << " : " << totalBondCount << endl;
    ++idx;
  }

  return encodeVector;
}

ClusterExpansionInfo GenerateData::generateCEData(const set<Element> &elementSet,
                                                  const size_t maxClusterSize,
                                                  const size_t maxBondOrder)
{
  Element vacancy("X");

  // Initial Configuration
  // Central Atom Id : Initial Vacancy Id
  // Neighbor Atom Id : Initial Migraiting Atom Id
  // CE encoding for this has already been saved

  size_t initVacancyId = vacancyId_;
  // size_t initMigratingElementLatticeId = migratingElementLatticeId_;

  // Copy config
  auto startConfig = config_;
  startConfig.UpdateNeighborList(cutoffs_);

  startConfig.SetElementOfLattice(initVacancyId, vacancy);

  Eigen::VectorXd startEncodeVector = getEncodeVector(startConfig,
                                                      elementSet,
                                                      maxClusterSize,
                                                      maxBondOrder);

  // Final Configuration
  // Central Atom Id : Migrating Atom
  // Its Neighbor : Vacancy

  size_t finalVacancyId = migratingElementLatticeId_;
  size_t finalMigratingElementLatticeId = vacancyId_;

  auto finalConfig = config_;
  finalConfig.UpdateNeighborList(cutoffs_);

  // Migrating Atom to the previous vacant site
  finalConfig.SetElementOfLattice(finalMigratingElementLatticeId,
                                  migratingAtomElement_);

  // Making a vacancy at neighbourAtom
  finalConfig.SetElementOfLattice(finalVacancyId, vacancy);

  Eigen::VectorXd endEncodeVector = getEncodeVector(finalConfig,
                                                    elementSet,
                                                    maxClusterSize,
                                                    maxBondOrder);

  ClusterExpansionInfo ceInfo(uniqueConfigId_,
                              alloySystem_,
                              initVacancyId,
                              config_.GetCartesianPositionOfLattice(initVacancyId),
                              finalVacancyId,
                              config_.GetCartesianPositionOfLattice(finalVacancyId),
                              migratingAtomElement_.GetElementString(),
                              startEncodeVector,
                              endEncodeVector);

  return ceInfo;
}

// Other variation of the function can also be there but the current focus is not on that
pair<VectorXd, VectorXd> GenerateData::getMigratingAtomNeighborPairEncodeVector(
    const set<Element> &elementSet)
{
  // Need to work on this
  size_t vacancyMigrationBO = 3;

  unordered_map<string, RowVectorXd> oneHotEncodingHashMap = GetOneHotEncodeHashmap(elementSet);

  // vector<vector<size_t>> equivalentSitesEncoding = GetEquivalentSites3Fold(config_,
  //                                                                         vacancyMigrationBO);

  vector<vector<size_t>> equivalentSitesEncoding = GetEquivalentSitesUnderKFoldRotation(config_, 
                                                                                        vacancyMigrationBO, 
                                                                                        6);

  // Forward Jump
  vector<size_t> forwardSymmSortedVector = GetSortedLatticeVectorStateOfPair(config_,
                                                                             forwardJumpPair_,
                                                                             vacancyMigrationBO);

  VectorXd forwardEncodeVector = GetEncodingMigratingAtomPair(config_,
                                                              equivalentSitesEncoding,
                                                              forwardSymmSortedVector,
                                                              oneHotEncodingHashMap,
                                                              migratingAtomElement_);

  // Backward Jump
  vector<size_t> backwardSymmSortedVector = GetSortedLatticeVectorStateOfPair(config_,
                                                                              backwardJumpPair_,
                                                                              vacancyMigrationBO);

  VectorXd backwardEncodeVector = GetEncodingMigratingAtomPair(config_,
                                                               equivalentSitesEncoding,
                                                               backwardSymmSortedVector,
                                                               oneHotEncodingHashMap,
                                                               migratingAtomElement_);

  return {forwardEncodeVector, backwardEncodeVector};
}


// Helper function to write a list of elements as a comma-separated string
void GenerateData::writeNeighbourVectorToFile(ofstream &outFile,
                                              const vector<size_t> &vec) const
{
  for (size_t i = 0; i < vec.size(); ++i)
  {
    if (i == vec.size() - 1)
    {
      outFile << vec[i] << config_.GetElementOfLattice(vec[i]);
    }
    else
    {
      outFile << vec[i] << config_.GetElementOfLattice(vec[i]) << ", ";
    }
  }
}

void writeEncodingToFile(ofstream &outFile,
                         const VectorXd &vec)
{
  for (size_t i = 0; i < vec.size(); ++i)
  {
    if (i == vec.size() - 1)
    {
      outFile << vec[i];
    }
    else
    {
      outFile << vec[i] << ", ";
    }
  }
}

// Function to write headers to the output file
void writeHeaders(ofstream &outFile)
{
  // Define column headers as a vector for easy modification
  const vector<string> headers = {
      "FolderId",
      "System",
      "Initial Vacancy ID",
      "Initial Vacancy Position",
      "Final Vacancy ID",
      "Final Vacancy Position",
      "Migrating Atom",
      "First Nearest Neighbors",
      "Second Nearest Neighbors",
      "Third Nearest Neighbors",
      "Migration Atom Pairs Forward Encoding",
      "Migration Atom Pairs Backward Encoding",
      "Cluster Expansion Encoding Start",
      "Cluster Expansion Encoding End",
      "isOrdered",
      "Element1",
      "Alpha Site Occupancy Element1",
      "Element2",
      "Alpha Site Occupancy Element2"};

  // Check if file is empty before writing headers
  if (outFile.tellp() == 0)
  {
    // Join headers with tab separator and write in one go
    string headerLine;
    for (const auto &header : headers)
    {
      if (!headerLine.empty())
      {
        headerLine += "\t"; // Add tab separator except for the first header
      }
      headerLine += header;
    }

    // Add newline at the end
    headerLine += "\n";

    outFile << headerLine;
  }
}

void GenerateData::writeVacancyJumpCeConfig(ofstream &outFile,
                                            const VacancyMigrationInfo &vacancyMigrationInfo,
                                            const ClusterExpansionInfo &ceInfo,
                                            const B2OrderingInfo &b2OrderInfo,
                                            const set<Element> &elementSet)
{
  if (outFile.is_open())
  {
    // Check if file is empty to write header only once
    if (outFile.tellp() == 0)
    {
      writeHeaders(outFile);
    }

    // Folder Id
    if (isOrderedStructure_)
    {
      outFile << uniqueConfigId_;
      outFile << "_O" << "\t";
    }
    else
    {
      outFile << uniqueConfigId_ << "\t";
    }

    // system
    outFile << alloySystem_ << "\t";

    // Initial Vacancy ID
    outFile << vacancyMigrationInfo.initialVacancyIndex << "\t";

    // Initial Vacancy Position
    outFile << vacancyMigrationInfo.initialVacancyPosition << "\t";

    // Final Vacancy ID
    outFile << vacancyMigrationInfo.finalVacancyIndex << "\t";

    // Final Vacancy Position
    outFile << vacancyMigrationInfo.finalVacancyPosition << "\t";

    // Migrating Atom
    outFile << vacancyMigrationInfo.migratingAtom.GetElementString() << "\t";

    // First Nearest Neighbors
    writeNeighbourVectorToFile(outFile, vacancyMigrationInfo.firstShellNeighbors);
    outFile << "\t";

    // Second Nearest Neighbors
    writeNeighbourVectorToFile(outFile, vacancyMigrationInfo.secondShellNeighbors);
    outFile << "\t";

    // Third Nearest Neighbors
    writeNeighbourVectorToFile(outFile, vacancyMigrationInfo.thirdShellNeighbors);
    outFile << "\t";

    // Migrating Atom Pairs Encoding
    pair<VectorXd, VectorXd> migratingAtomPairEncoding = getMigratingAtomNeighborPairEncodeVector(elementSet);

    // Migration Atom Pairs Forward Encoding
    writeEncodingToFile(outFile, migratingAtomPairEncoding.first);
    outFile << "\t";

    // Migration Atom Pairs Backward Encoding
    writeEncodingToFile(outFile, migratingAtomPairEncoding.second);
    outFile << "\t";

    // CE Encodings

    // Cluster Expansion Encoding Start
    writeEncodingToFile(outFile, ceInfo.ceEncodingInitial);
    outFile << "\t";

    // Cluster Expansion Encoding End
    writeEncodingToFile(outFile, ceInfo.ceEncodingFinal);
    outFile << "\t";

    // Ordering Info

    if (b2OrderInfo.isOrderedStructure)
    {
      // isOrdered
      outFile << "Yes" << "\t";

      // element 1
      outFile << b2OrderInfo.element1.GetElementString() << "\t";

      // Alpha site fractional occupancy
      outFile << b2OrderInfo.alphaOccupancyFractionElement1 << "\t";

      // element 2
      outFile << b2OrderInfo.element2.GetElementString() << "\t";

      // beta site fractional occupancy
      outFile << b2OrderInfo.alphaOccupancyFractionElement2 << "\t";
    }

    else
    {
      // If not ordered, write "N/A" for all fields
      outFile << "N/A" << "\t";
      outFile << "N/A" << "\t";
      outFile << "N/A" << "\t";
      outFile << "N/A" << "\t";
      outFile << "N/A" << "\t";
    }

    outFile << endl;
  }
}

B2OrderingInfo GenerateData::generateB2Structure(vector<double> elementComposition,
                                                 const Element elementAtAlphaSite,
                                                 const Element elementAtBetaSite)
{
  B2Ordering b2Ordering(config_);

  // Get unordered sets of B2 alpha and beta sites
  unordered_set<size_t> alphaSitesSet = b2Ordering.GetAlphaLatticeSites();
  unordered_set<size_t> betaSitesSet = b2Ordering.GetBetaLatticeSites();

  vector<size_t> alphaSites(alphaSitesSet.begin(), alphaSitesSet.end());
  vector<size_t> betaSites(betaSitesSet.begin(), betaSitesSet.end());

  auto cfg = config_;
  cfg.UpdateNeighborList(cutoffs_);

  double totalComposition = 0;
  for (auto composition : elementComposition)
  {
    totalComposition += composition;
  }
  // Ensure the fraction is within a valid range
  if (totalComposition != 100.0)
  {
    throw invalid_argument("Total composition must be 100%, but got " +
                           to_string(totalComposition) + "%");
  }

  mt19937 mersenne_twister_rng_(uniqueConfigId_);

  // Shuffle sites to ensure randomness in element assignment.
  shuffle(alphaSites.begin(), alphaSites.end(), mersenne_twister_rng_);
  shuffle(betaSites.begin(), betaSites.end(), mersenne_twister_rng_);

  // Determine the number of alpha and beta sites to modify based on the fraction.
  size_t numAlphaToAssign = static_cast<size_t>(elementComposition[0] * alphaSites.size() / 100.0);
  size_t numBetaToAssign = static_cast<size_t>(elementComposition[1] * betaSites.size() / 100.0);

  // Assign elements to a subset of the identified sites.
  for (size_t i = 0; i < numAlphaToAssign; ++i)
  {
    cfg.SetElementOfLattice(alphaSites[i], elementAtAlphaSite);
  }
  for (size_t i = 0; i < numBetaToAssign; ++i)
  {
    cfg.SetElementOfLattice(betaSites[i], elementAtBetaSite);
  }

  B2Ordering b2OrderedStructure(cfg);

  B2OrderingInfo b2OrderInfo(
      true,
      cfg,
      elementAtAlphaSite,
      b2OrderedStructure.GetAlphaSiteOccupancy(elementAtAlphaSite),
      elementAtBetaSite,
      b2OrderedStructure.GetAlphaSiteOccupancy(elementAtBetaSite));

  return b2OrderInfo;
}
