#include "GenerateStructuresWithB2.h"

bool isB2Structure(const Config &config,
                   const size_t latticeId)
{
  // alpha element
  auto centerElement = config.GetElementOfLattice(latticeId);
  // beta element
  Element oppositeElement("X");

  // first NN
  for (const auto firstNNId : config.GetNeighborLatticeIdVectorOfLattice(latticeId, 1))
  {
    if (oppositeElement == Element("X"))
    {
      oppositeElement = config.GetElementOfLattice(firstNNId);

      if (oppositeElement == centerElement)
      {
        return false;
      }
    }
    else
    {
      if (config.GetElementOfLattice(firstNNId) != oppositeElement)
      {
        return false;
      }
    }
  }

  for (const auto secondNNId : config.GetNeighborLatticeIdVectorOfLattice(latticeId, 2))
  {
    if (config.GetElementOfLattice(secondNNId) != centerElement)
    {
      return false;
    }
  }

  // Only B2 if the first NN of the central element are opposite
  // and second NN are same as the central element
  return true;
}

void buildSingleB2(Config &config,
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

void BuildConfigWithB2(Config &config,
                       int i,
                       map<Element, int> elementMap,
                       const unordered_set<size_t> &visitedSites, int randomId)
{

  random_device rd;
  mt19937 gen(rd());

  auto numAtoms = config.GetNumAtoms();
  //
  cout << "After B2" << endl;

  for (auto entry : elementMap)
  {
    cout << entry.first.GetElementString() << " " << entry.second << endl;
  }

  // Randomly assign the remaining elements to those  lattice Ids which are not visited yet

  // Convert unordered_set to vector and shuffle

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

  for (auto id : remainingLatticeIds)
  {
    if (isB2Structure(config, id)) // Detect local B2-like ordering
    {
      bool swapped = false;

      // Try a limited number of swaps to prevent infinite loops
      for (int attempt = 0; attempt < 20; ++attempt)
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
          if (!isB2Structure(config, id) && !isB2Structure(config, randomId))
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
        std::cerr << "Warning: could not break B2 structure at site " << id << std::endl;
      }
    }
  }

  cout << "After assignment" << endl;

  for (auto entry : elementMap)
  {
    cout << entry.first.GetElementString() << " " << entry.second << endl;
  }

  Config::WriteConfig("//media/sf_Phd/CNT/structures/b2_" + to_string(i) + "_" + to_string(randomId) + ".cfg.gz", config);
}

void GenerateStructureWithB2(Config &config,
                             size_t numB2Centers,
                             pair<Element, Element> &b2OrderedElements, int randomId)
{
  auto atomVector = config.GetAtomVector();

  map<Element, int> elementMap;

  for (const auto &atom : atomVector)
  {
    elementMap[atom]++;
  }

  for (auto entry : elementMap)
  {
    cout << entry.first.GetElementString() << " " << entry.second << endl;
  }

  Element alphaElement = b2OrderedElements.first;
  Element betaElement = b2OrderedElements.second;

  size_t numAtoms = config.GetNumAtoms();

  random_device rd;
  mt19937 gen(rd());

  uniform_int_distribution<> dis(0, 7);

  // single B2

  size_t centralLatticeId = config.GetCentralAtomLatticeId();

  vector<size_t> selectedLatticeIdVector = {centralLatticeId};

  unordered_set<size_t> visitedSites;

  unordered_set<size_t> b2CenterLatticeIdSet;

  for (int i = 0; i < numB2Centers; i++)
  {
    size_t selectedLatticeId = selectedLatticeIdVector[i];

    if (config.GetElementOfLattice(selectedLatticeId) == alphaElement)
    {
      buildSingleB2(config, elementMap, selectedLatticeId, alphaElement, betaElement, visitedSites);
    }
    else
    {
      buildSingleB2(config, elementMap, selectedLatticeId, betaElement, alphaElement, visitedSites);
    }

    b2CenterLatticeIdSet.insert(selectedLatticeId);

    BuildConfigWithB2(config, i, elementMap, visitedSites, 0);
    BuildConfigWithB2(config, i, elementMap, visitedSites, 1);
    BuildConfigWithB2(config, i, elementMap, visitedSites, 2);
    BuildConfigWithB2(config, i, elementMap, visitedSites, 3);


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
