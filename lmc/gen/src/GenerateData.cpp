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

  // Previous
  // forwardJumpPair_ = {migratingElementLatticeId_, vacancyId_};

  // backwardJumpPair_ = {vacancyId_, migratingElementLatticeId_};

  // Now

  forwardJumpPair_ = {vacancyId_, migratingElementLatticeId_};

  backwardJumpPair_ = {migratingElementLatticeId_, vacancyId_};
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

  forwardJumpPair_ = {vacancyId_, migratingElementLatticeId_};

  backwardJumpPair_ = {migratingElementLatticeId_, vacancyId_};
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

  size_t kFoldRotation = 6;
  vector<vector<size_t>> equivalentSitesEncoding = GetEquivalentSitesUnderKFoldRotation(config_,
                                                                                        vacancyMigrationBO,
                                                                                        kFoldRotation);

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
      "Alpha Site Occupancy Element2",
      "numB2Center"};

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

      // number of B2 Centers
      outFile << b2OrderInfo.numB2Centers << "\t";
    }

    else
    {
      // If not ordered, write "N/A" for all fields
      outFile << "N/A" << "\t";
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
/*
void GenerateData::addB2Ordering(double percentB2)
{
  auto atomVector = config_.GetAtomVector();

  set elementSet(atomVector.begin(), atomVector.end());

  if (elementSet.size() == 2)
  {

    // some lattice Id

    unordered_set<size_t> latticeIdAssigned = {};

    size_t b2CenterLatticeId = 0;


    latticeIdAssigned.insert(b2CenterLatticeId);

    auto b2CenterElement = config_.GetElementOfLattice(b2CenterLatticeId);

    Element neighbourElement;

    for (const auto &element : elementSet)
    {
      if (element != b2CenterElement)
      {
        neighbourElement = element;
        break; // Only two elements, so we can stop here
      }
    }

    // Neighbour
    auto firstNNB2Center = config_.GetNeighborLatticeIdVectorOfLattice(b2CenterLatticeId, 1);
    // set all the element of neigbour to the neighbour element and
    // keep track of the b2CenterElement

    // Stores count of of those element which are changed
    // to keep the concentration constant
    int numB2CenterElement = 0;

    for (auto nnId : firstNNB2Center)
    {
      if (config_.GetElementOfLattice(nnId) == b2CenterElement)
      {
        config_.SetElementOfLattice(nnId, neighbourElement);
        // these number of b2 center element removed
        numB2CenterElement += 1;
      }
      latticeIdAssigned.insert(nnId);
    }

    // second nearest nn should be same as the b2 center element

    auto secondNNB2Center = config_.GetNeighborLatticeIdVectorOfLattice(b2CenterLatticeId, 2);

    for (auto nnId : secondNNB2Center)
    {
      if (config_.GetElementOfLattice(nnId) != b2CenterElement)
      {
        config_.SetElementOfLattice(nnId, b2CenterElement);
        numB2CenterElement -= 1;
      }
      latticeIdAssigned.insert(nnId);
    }

    // 15 ids are done

    if (numB2CenterElement != 0)
    {
      if (numB2CenterElement < 0)
      {
        // more number of B2 center element in the config
        for (size_t id = 0; id < config_.GetNumAtoms(); id++)
        {
          if (latticeIdAssigned.find(id) == latticeIdAssigned.end())
          {
            // id not part of the B2 structure
            while (numB2CenterElement)
            {
              if (config_.GetElementOfLattice(id) == b2CenterElement)
              {
                config_.SetElementOfLattice(id, neighbourElement);
                numB2CenterElement+=1;
              }
            }
          }
        }
      }
      else
      {
        // less number of B2 center element in the config
        for (size_t id = 0; id < config_.GetNumAtoms(); id++)
        {
          if (latticeIdAssigned.find(id) == latticeIdAssigned.end())
          {
            // id not part of the B2 structure
            while (numB2CenterElement)
            {
              if (config_.GetElementOfLattice(id) != b2CenterElement)
              {
                config_.SetElementOfLattice(id, neighbourb2CenterElementElement);
                numB2CenterElement-=1;
              }
            }
          }
        }
      }
    }

    // Change the neighbour to opposite
  }
}

B2OrderingInfo GenerateData::AddB2Structure(const int numB2Center,
                                            vector<Element> elementSet)
{
  if (elementSet.size() != 2)
  {
    throw std::runtime_error("B2 ordering requires exactly two elements");
  }

  auto config = config_;

  // Get atom vector and unique elements
  auto atomVector = config.GetAtomVector();

  // Extract the two elements;
  Element centerElement = elementSet[0];
  Element neighborElement = elementSet[1];

  // Calculate number of sites to apply B2 ordering based on percentB2
  size_t totalSites = config.GetNumAtoms();
  size_t numB2Sites = numB2Center * 15;

  // Track assigned lattice sites
  std::unordered_set<size_t> assignedSites;

  // Random number generator for selecting center sites
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<size_t> dist(0, totalSites - 1);

  // Count current elements to maintain concentration
  size_t centerElementCount = std::count(atomVector.begin(), atomVector.end(), centerElement);
  size_t neighborElementCount = totalSites - centerElementCount;

  int deltaCenterElement = 0;   // Net change in centerElement count
  int deltaNeigbourElement = 0; // Net change in neighbour element count

  // Apply B2 ordering to selected sites
  size_t sitesProcessed = 0;
  while (sitesProcessed < numB2Sites && assignedSites.size() < totalSites)
  {
    // Select a random unassigned lattice site as B2 center
    size_t centerId;
    do
    {
      centerId = dist(gen);
    } while (assignedSites.find(centerId) != assignedSites.end());

    assignedSites.insert(centerId);

    // Set center to centerElement
    if (config.GetElementOfLattice(centerId) != centerElement)
    {
      config.SetElementOfLattice(centerId, centerElement);
      deltaCenterElement++;
    }

    // Set first nearest neighbors to neighborElement
    auto firstNN = config.GetNeighborLatticeIdVectorOfLattice(centerId, 1);
    for (size_t nnId : firstNN)
    {
      if (assignedSites.find(nnId) == assignedSites.end())
      {
        if (config.GetElementOfLattice(nnId) != neighborElement)
        {
          config.SetElementOfLattice(nnId, neighborElement);
          deltaNeigbourElement++;
        }
        assignedSites.insert(nnId);
      }
    }

    // Set second nearest neighbors to centerElement
    auto secondNN = config.GetNeighborLatticeIdVectorOfLattice(centerId, 2);
    for (size_t nnId : secondNN)
    {
      if (assignedSites.find(nnId) == assignedSites.end())
      {
        if (config.GetElementOfLattice(nnId) != centerElement)
        {
          config.SetElementOfLattice(nnId, centerElement);
          deltaCenterElement++;
        }
        assignedSites.insert(nnId);
      }
    }

    sitesProcessed += 1 + firstNN.size() + secondNN.size();
  }

  // Correct the composition

  centerElementCount -= deltaCenterElement;
  neighborElementCount -= deltaNeigbourElement;

  int id = 0;
  while ((centerElementCount != neighborElementCount) &&
         id < totalSites)
  {
    auto currentElement = config.GetElementOfLattice(id);

    if (assignedSites.find(id) == assignedSites.end())
    {
      if ((centerElementCount < neighborElementCount) &&
          currentElement == neighborElement)
      {
        config.SetElementOfLattice(id, centerElement);
        centerElementCount++;
      }
      else if ((centerElementCount > neighborElementCount) &&
               currentElement == centerElement)
      {
        config.SetElementOfLattice(id, neighborElement);
        neighborElementCount++;
      }
    }
    id++;
  }

  B2Ordering b2OrderedStructure(config);

  return B2OrderingInfo(
      true,
      config,
      centerElement,
      b2OrderedStructure.GetAlphaSiteOccupancy(centerElement),
      neighborElement,
      b2OrderedStructure.GetAlphaSiteOccupancy(neighborElement),
      numB2Center);
}

*/

B2OrderingInfo GenerateData::AddB2Structure(size_t numB2Center,
                                            vector<Element> &elementSet)
{
  if (elementSet.size() != 2)
  {
    throw std::runtime_error("B2 ordering requires exactly two elements");
  }

  auto config = config_;

  // Get atom vector and unique elements
  // auto atomVector = config.GetAtomVector();

  // Extract the two elements;
  Element alphaElement = elementSet[0];
  Element betaElement = elementSet[1];

  size_t totalSites = config.GetNumAtoms();

  unordered_set<size_t> allAtomIds;

  // Insert 0 to n-1
  for (int id = 0; id < totalSites; ++id)
  {
    allAtomIds.insert(id);
  }

  // random number generator which will select where to insert the B2
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<size_t> selectB2Center(0, totalSites - 1);

  int numAlpha = totalSites / 2;
  int numBeta = totalSites / 2;

  // Simple B2
  if (numB2Center == 1)
  {
    size_t b2CenterAtomId = selectB2Center(gen);

    config.SetElementOfAtom(b2CenterAtomId, alphaElement);

    allAtomIds.erase(b2CenterAtomId);
    numAlpha -= 1;

    auto firstNN = config.GetNeighborAtomIdVectorOfAtom(b2CenterAtomId, 1);

    for (auto id : firstNN)
    {
      config.SetElementOfAtom(id, betaElement);
      allAtomIds.erase(id);
      numBeta -= 1;
    }

    auto secondNN = config.GetNeighborAtomIdVectorOfAtom(b2CenterAtomId, 2);

    for (auto id : secondNN)
    {
      config.SetElementOfAtom(id, alphaElement);
      allAtomIds.erase(id);
      numAlpha -= 1;
    }
  }

  else if (numB2Center == 2)
  {
    // Corner shared B2

    // b2 Center 1
    size_t b2Center1AtomId = selectB2Center(gen);

    config.SetElementOfAtom(b2Center1AtomId, alphaElement);

    allAtomIds.erase(b2Center1AtomId);
    numAlpha -= 1;

    auto firstNN = config.GetNeighborAtomIdVectorOfAtom(b2Center1AtomId, 1);

    for (auto id : firstNN)
    {
      config.SetElementOfAtom(id, betaElement);
      allAtomIds.erase(id);
      numBeta -= 1;
    }

    auto secondNN = config.GetNeighborAtomIdVectorOfAtom(b2Center1AtomId, 2);

    for (auto id : secondNN)
    {
      config.SetElementOfAtom(id, alphaElement);
      allAtomIds.erase(id);
      numAlpha -= 1;
    }

    // Second B2 center will be the first nearest neighbour of the the first B2 center

    size_t b2Center2AtomId = firstNN[0];

    config.SetElementOfAtom(b2Center2AtomId, betaElement);

    if (allAtomIds.find(b2Center2AtomId) != allAtomIds.end())
    {
      allAtomIds.erase(b2Center2AtomId);
      numBeta -= 1;
    }

    firstNN = config.GetNeighborAtomIdVectorOfAtom(b2Center2AtomId, 1);

    for (auto id : firstNN)
    {
      if (allAtomIds.find(id) != allAtomIds.end())
      {

        config.SetElementOfAtom(id, alphaElement);
        allAtomIds.erase(id);
        numAlpha -= 1;
      }
    }

    secondNN = config.GetNeighborAtomIdVectorOfAtom(b2Center2AtomId, 2);

    for (auto id : secondNN)
    {
      if (allAtomIds.find(id) != allAtomIds.end())
      {
        config.SetElementOfAtom(id, betaElement);
        allAtomIds.erase(id);
        numBeta -= 1;
      }
    }
  }

  else if (numB2Center == 3)
  {
    // Corner shared B2

    // b2 Center 1
    size_t b2Center1AtomId = selectB2Center(gen);

    config.SetElementOfAtom(b2Center1AtomId, alphaElement);

    allAtomIds.erase(b2Center1AtomId);
    numAlpha -= 1;

    auto firstNN = config.GetNeighborAtomIdVectorOfAtom(b2Center1AtomId, 1);

    for (auto id : firstNN)
    {
      config.SetElementOfAtom(id, betaElement);
      allAtomIds.erase(id);
      numBeta -= 1;
    }

    auto secondNN = config.GetNeighborAtomIdVectorOfAtom(b2Center1AtomId, 2);

    for (auto id : secondNN)
    {
      config.SetElementOfAtom(id, alphaElement);
      allAtomIds.erase(id);
      numAlpha -= 1;
    }

    // Second B2 center will be the first nearest neighbour of the the first B2 center

    size_t b2Center2AtomId = firstNN[0];
    size_t b2Center3AtomId = firstNN[2];

    config.SetElementOfAtom(b2Center2AtomId, betaElement);

    if (allAtomIds.find(b2Center2AtomId) != allAtomIds.end())
    {
      allAtomIds.erase(b2Center2AtomId);
      numBeta -= 1;
    }

    firstNN = config.GetNeighborAtomIdVectorOfAtom(b2Center2AtomId, 1);

    for (auto id : firstNN)
    {
      if (allAtomIds.find(id) != allAtomIds.end())
      {

        config.SetElementOfAtom(id, alphaElement);
        allAtomIds.erase(id);
        numAlpha -= 1;
      }
    }

    secondNN = config.GetNeighborAtomIdVectorOfAtom(b2Center2AtomId, 2);

    for (auto id : secondNN)
    {
      if (allAtomIds.find(id) != allAtomIds.end())
      {
        config.SetElementOfAtom(id, betaElement);
        allAtomIds.erase(id);
        numBeta -= 1;
      }
    }

    // third B2 center

    if (allAtomIds.find(b2Center3AtomId) != allAtomIds.end())
    {
      allAtomIds.erase(b2Center3AtomId);
      numBeta -= 1;
    }

    firstNN = config.GetNeighborAtomIdVectorOfAtom(b2Center3AtomId, 1);

    for (auto id : firstNN)
    {
      if (allAtomIds.find(id) != allAtomIds.end())
      {

        config.SetElementOfAtom(id, alphaElement);
        allAtomIds.erase(id);
        numAlpha -= 1;
      }
    }

    secondNN = config.GetNeighborAtomIdVectorOfAtom(b2Center3AtomId, 2);

    for (auto id : secondNN)
    {
      if (allAtomIds.find(id) != allAtomIds.end())
      {
        config.SetElementOfAtom(id, betaElement);
        allAtomIds.erase(id);
        numBeta -= 1;
      }
    }
  }

  // Iterate over the rest of atom Ids
  // assign remaining alpha and beta randomly

  std::vector<size_t> remainingIds(allAtomIds.begin(), allAtomIds.end());

  // Shuffle the remaining atom IDs randomly
  std::shuffle(remainingIds.begin(), remainingIds.end(), gen);

  for (auto id : remainingIds)
  {
    if (numAlpha > 0 && numBeta > 0)
    {
      // Randomly decide between alpha and beta
      std::uniform_int_distribution<int> dist(0, 1);
      int choice = dist(gen);
      if (choice == 0)
      {
        config.SetElementOfAtom(id, alphaElement);
        numAlpha -= 1;
      }
      else
      {
        config.SetElementOfAtom(id, betaElement);
        numBeta -= 1;
      }
    }
    else if (numAlpha > 0)
    {
      config.SetElementOfAtom(id, alphaElement);
      numAlpha -= 1;
    }
    else if (numBeta > 0)
    {
      config.SetElementOfAtom(id, betaElement);
      numBeta -= 1;
    }
  }

  B2Ordering b2OrderedStructure(config);

  return B2OrderingInfo(
    true,
    config,
    alphaElement,
    b2OrderedStructure.GetAlphaSiteOccupancy(alphaElement),
    betaElement,
    b2OrderedStructure.GetAlphaSiteOccupancy(betaElement),
    numB2Center);
}
