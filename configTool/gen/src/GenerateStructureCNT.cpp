#include "GenerateStructureCNT.h"

GenerateStructureCNT::GenerateStructureCNT(
    const string &predictorFilename,
    const Config &config,
    const Config &trainingConfig,
    const size_t supercellSize,
    set<Element> &elementSet,
    const vector<double> &compositionVector) : PotentialEnergyEstimator(predictorFilename,
                                                                        config,
                                                                        trainingConfig,
                                                                        elementSet),
                                               supercellSize_(supercellSize),
                                               elementVector_([&]()
                                                              {
                                                  vector<string> vec;
                                                  for (const auto& element : elementSet) 
                                                    if (element != Element("X"))
                                                      vec.emplace_back(element.GetElementString());
                                                  return vec; }()),
                                               compositionVector_(compositionVector),
                                               allLatticeClusterSet_(
                                                   FindAllLatticeClusters(
                                                       config,
                                                       maxClusterSize_,
                                                       maxBondOrder_,
                                                       {}))

{
  if (config.GetNumAtoms() != supercellSize * supercellSize * supercellSize * 2)
  {
    ostringstream oss;
    oss << "Error in GenerateStructureCNT: "
        << "Config contains " << config.GetNumAtoms() << " atoms, "
        << "but expected " << (supercellSize * supercellSize * supercellSize * 2)
        << " atoms for a BCC supercell of size " << supercellSize << "^3.";
    throw runtime_error(oss.str());
  }
}

double GenerateStructureCNT::GetEnergyOfConfig(const Config &config) const
{
  auto clusterTypeCountHashMap = clusterTypeCountHashMap_;

  for (const auto &latticeCluster : allLatticeClusterSet_)
  {
    auto atomClusterType = IdentifyAtomClusterType(config, latticeCluster.GetLatticeIdVector());
    clusterTypeCountHashMap.at(ClusterType(atomClusterType, latticeCluster.GetClusterType()))++;
  }

  Eigen::VectorXd encodeVector(initializedClusterTypeSet_.size());
  int idx = 0;

  for (const auto &clusterType : initializedClusterTypeSet_)
  {
    // Count of Cluster Types for given configuration
    auto count_bond = static_cast<double>(clusterTypeCountHashMap.at(clusterType));
    // Count of Cluster Types for normalization
    auto total_bond = static_cast<double>(latticeClusterTypeCount_.at(clusterType.lattice_cluster_type_));

    encodeVector(idx) = count_bond / total_bond;
    ++idx;
  }

  double energy = betaCE_.dot(encodeVector);

  return energy;
}

void GenerateStructureCNT::GenerateRandomStructures(
    const size_t numStructures,
    const bool doComputeEnergy,
    const bool doSaveConfig,
    const string &outputDir)
{
  for (int i = 0; i < int(numStructures); i++)
  {
    auto randomConfig = Config::GenerateAlloySupercell(supercellSize_,
                                                       latticeParam_,
                                                       structureType_,
                                                       elementVector_,
                                                       compositionVector_,
                                                       i);
    randomConfig.UpdateNeighborList(cutoffs_);

    double energyOfConfig = 0;
    if (doComputeEnergy)
    {
      energyOfConfig = GetEnergyOfConfig(randomConfig);
      cout << i << "\t" << energyOfConfig << endl;
    }

    if (doSaveConfig)
    {
      // Create a directory as randomConfig

      ostringstream folderName;
      folderName << outputDir << "/randomConfig";

      fs::path folderPath = folderName.str();

      if (!fs::exists(folderPath))
      {
        fs::create_directories(folderPath);
      }

      // Create filename
      ostringstream filename;
      filename << "random_" << i << ".xyz.gz";

      fs::path filePath = folderPath / filename.str();

      map<string, Config::ValueVariant> globalList = {{"total_energy", energyOfConfig}};

      Config::WriteXyzExtended(filePath, randomConfig, {}, globalList);
    }
  }
}

void GenerateStructureCNT::GenerateStructureWithB2(
    Config &config,
    const size_t numB2Centers,
    const size_t numSamplesPerB2Config,
    const pair<Element, Element> &b2OrderedElements,
    const bool doComputeEnergy,
    const string &outputDir)
{
  auto atomVector = config.GetAtomVector();

  map<Element, int> elementMap;

  for (const auto &atom : atomVector)
  {
    elementMap[atom]++;
  }

  Element alphaElement = b2OrderedElements.first;
  Element betaElement = b2OrderedElements.second;

  random_device rd;
  mt19937 gen(rd());

  uniform_int_distribution<> dis(0, config.GetNeighborLatticeIdVectorOfLattice(0, 1).size() - 1);

  size_t centralLatticeId = config.GetCentralAtomLatticeId();

  vector<size_t> selectedLatticeIdVector = {centralLatticeId};

  unordered_set<size_t> visitedSites;
  unordered_set<size_t> b2CenterLatticeIdSet;

  for (int numB2 = 1; numB2 < numB2Centers + 1; numB2++)
  {
    size_t selectedLatticeId = selectedLatticeIdVector[numB2 - 1];

    if (config.GetElementOfLattice(selectedLatticeId) == alphaElement)
    {
      BuildSingleB2(config, elementMap, selectedLatticeId, alphaElement, betaElement, visitedSites);
    }
    else
    {
      BuildSingleB2(config, elementMap, selectedLatticeId, betaElement, alphaElement, visitedSites);
    }

    b2CenterLatticeIdSet.insert(selectedLatticeId);

    // How many configuration to save
    for (int numSample = 1; numSample < numSamplesPerB2Config + 1; numSample++)
    {
      // Now we have configuration with B2 embedded in it
      // Next task is to build random config around it such that no B2 is there
      // And also to get random configuration randomize the B2 precipitate
      BuildConfigWithB2(config, numB2, elementMap, visitedSites, numSample, doComputeEnergy, outputDir);
    }

    // possible b2 centers based on the previous selected lattice Id
    vector<size_t> possibleB2Centers = config.GetNeighborLatticeIdVectorOfLattice(selectedLatticeId, 1);
    for (auto id : possibleB2Centers)
    {
      if (b2CenterLatticeIdSet.find(id) == b2CenterLatticeIdSet.end())
      {
        selectedLatticeIdVector.push_back(id);
      }
    }
  }
}

void GenerateStructureCNT::BuildSingleB2(Config &config,
                                         map<Element, int> &elementMap,
                                         const size_t selectedLatticeId,
                                         const Element &alphaElement,
                                         const Element &betaElement,
                                         unordered_set<size_t> &visitedSites)
{
  // assign alpha element to lattice id
  config.SetElementOfLattice(selectedLatticeId, alphaElement);

  if (visitedSites.find(selectedLatticeId) == visitedSites.end())
  {
    elementMap[alphaElement]--;
    visitedSites.insert(selectedLatticeId);
  }

  // assign beta element to first nearest neighbours
  for (const auto &firstNNId : config.GetNeighborLatticeIdVectorOfLattice(selectedLatticeId, 1))
  {
    config.SetElementOfLattice(firstNNId, betaElement);
    if (visitedSites.find(firstNNId) == visitedSites.end())
    {

      elementMap[betaElement]--;
      visitedSites.insert(firstNNId);
    }
  }

  // assign alpha element to second nearest neighbours
  for (const auto &secondNNId : config.GetNeighborLatticeIdVectorOfLattice(selectedLatticeId, 2))
  {
    config.SetElementOfLattice(secondNNId, alphaElement);

    if (visitedSites.find(secondNNId) == visitedSites.end())
    {
      elementMap[alphaElement]--;
      visitedSites.insert(secondNNId);
    }
  }
}

void GenerateStructureCNT::BuildConfigWithB2(Config &config,
                                             const size_t numB2Center,
                                             map<Element, int> elementMap,
                                             const unordered_set<size_t> &visitedSites,
                                             int randomId,
                                             const bool doComputeEnergy,
                                             const string &outputDir)
{

  random_device rd;
  mt19937 gen(rd());

  auto numAtoms = config.GetNumAtoms();

  // Randomly assign the remaining elements to those  lattice Ids which are not visited yet
  vector<size_t> remainingLatticeIds;

  for (size_t id = 0; id < numAtoms; id++)
  {
    if (visitedSites.find(id) == visitedSites.end())
    {
      remainingLatticeIds.emplace_back(id);
    }
  }
  shuffle(remainingLatticeIds.begin(), remainingLatticeIds.end(), gen);

  // Shuffle element pool
  vector<Element> elementPool;
  for (const auto &[el, count] : elementMap)
  {
    for (int i = 0; i < count; ++i)
      elementPool.push_back(el);
  }
  shuffle(elementPool.begin(), elementPool.end(), gen);

  assert(remainingLatticeIds.size() == elementPool.size());

  for (size_t i = 0; i < remainingLatticeIds.size(); ++i)
  {
    Element element = elementPool[i];
    size_t latticeId = remainingLatticeIds[i];

    config.SetElementOfLattice(latticeId, element);

    assert(elementMap[element] > 0);
    elementMap[element]--;
  }

  uniform_int_distribution<> dis(0, remainingLatticeIds.size() - 1);

  // To avoid formation of B2 precipitate in the disordered matrix
  // Swap if a B2 is formed
  // Increase the swap count if supercell is large
  int numSwapTry = 20;
  for (auto id : remainingLatticeIds)
  {
    if (B2Ordering::isB2Structure(config, id)) // Detect local B2-like ordering
    {
      bool swapped = false;

      // Try a limited number of swaps to prevent infinite loops
      for (int attempt = 0; attempt < numSwapTry; ++attempt)
      {
        size_t randomId = remainingLatticeIds[dis(gen)];

        // Don't swap with self
        if (randomId == id)
          continue;

        Element e1 = config.GetElementOfLattice(id);
        Element e2 = config.GetElementOfLattice(randomId);

        // Only swap if elements are different
        if (e1 != e2)
        {
          config.LatticeJump({id, randomId}); // Swaps elements at those lattice IDs

          // Check if swap broke B2 pattern at both sites
          if (!B2Ordering::isB2Structure(config, id) && !B2Ordering::isB2Structure(config, randomId))
          {
            swapped = true;
            break;
          }
          else
          {
            // Undo the swap if it still forms B2
            config.LatticeJump({id, randomId});
          }
        }
      }

      // Optional: report if unable to remove B2 after attempts
      if (!swapped)
      {
        cerr << "Warning: could not break B2 structure at site " << id << endl;
      }
    }
  }

  double energyOfConfig = 0;
  if (doComputeEnergy)
  {
    energyOfConfig = GetEnergyOfConfig(config);
    cout << numB2Center << "\tR" << randomId << "\t" << energyOfConfig << endl;
  }

  // make a new folder numB2_{numB2Center}
  // and inside it save the config file with b2_{numB2Center}_{randomId}.cfg.gz

  // Make folder path

  fs::path folderPath = outputDir + "/orderedConfig";

  if (!fs::exists(folderPath))
  {
    fs::create_directories(folderPath);
  }

  // Create filename
  ostringstream filename;
  filename << "b2_" << numB2Center << "_" << randomId << ".cfg.gz";

  fs::path filePath = folderPath / filename.str();

  map<string, Config::ValueVariant> globalList = {{"total_energy", energyOfConfig}};

  Config::WriteXyzExtended(filePath, config, {}, globalList);
}
