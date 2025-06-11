#include "B2Ordering.h"

B2Ordering::B2Ordering(const Config &config) : config_(config)
{
  InitializeAlphaLatticeSites();
  InitializeBetaLatticeSites();
}

double B2Ordering::GetB2OrderParameter(const Element &element)
{
  // fractional occupancy of element at alpha and beta sites
  auto alphaOccupancy = GetAlphaSiteOccupancy(element);
  auto betaOccupancy = GetBetaSiteOccupancy(element);

  double b2OrderParameter = (alphaOccupancy - betaOccupancy) / (alphaOccupancy + betaOccupancy);

  return b2OrderParameter;
}

double B2Ordering::GetAlphaSiteOccupancy(const Element &element)
{
  // Number of occupied alpha sites by element
  size_t numOccupiedASites = 0;
  // size_t numAlphaSites = alphaLatticeSites.size();

  for (const auto &id : alphaLatticeSites)
  {
    auto siteElement = config_.GetElementOfLattice(id);

    if (siteElement == element)
    {
      numOccupiedASites += 1;
    }
  }

  double alphaSiteOccupancy = double(numOccupiedASites) / double(config_.GetNumAtoms());

  return alphaSiteOccupancy;
}

double B2Ordering::GetBetaSiteOccupancy(const Element &element)
{
  // Number of occupied beta sites by element
  size_t numOccupiedBSites = 0;
  // size_t numBetaSites = betaLatticeSites.size();

  for (const auto &id : betaLatticeSites)
  {
    auto siteElement = config_.GetElementOfLattice(id);

    if (siteElement == element)
    {
      numOccupiedBSites += 1;
    }
  }

  // will be between 0 and 0.5
  double betaSiteOccupancy = double(numOccupiedBSites) / double(config_.GetNumAtoms());

  return betaSiteOccupancy;
}

unordered_set<size_t> B2Ordering::GetAlphaLatticeSites()
{
  return alphaLatticeSites;
}

unordered_set<size_t> B2Ordering::GetBetaLatticeSites()
{
  return betaLatticeSites;
}

void B2Ordering::InitializeAlphaLatticeSites()
{

  auto secondNNList = config_.GetNeighborLists()[1];

  for (size_t id1 = 0; id1 < secondNNList.size(); id1++)
  {
    const auto &secondNN = secondNNList[id1];

    for (size_t id2 : secondNN)
    {
      bool id1Valid = true;
      bool id2Valid = true;

      // Check if id1 has any bond order of 1 with sites already in alphaLatticeSites
      for (size_t id3 : alphaLatticeSites)
      {
        if (config_.GetDistanceOrder(id1, id3) == 1)
        {
          id1Valid = false;
          break;
        }
      }

      // Check if id2 has any bond order of 1 with sites already in alphaLatticeSites
      for (size_t id3 : alphaLatticeSites)
      {
        if (config_.GetDistanceOrder(id2, id3) == 1)
        {
          id2Valid = false;
          break;
        }
      }

      // Add valid sites to alphaLatticeSites
      if (id1Valid)
      {
        alphaLatticeSites.emplace(id1);
      }
      if (id2Valid)
      {
        alphaLatticeSites.emplace(id2);
      }
    }
  }
}

void B2Ordering::InitializeBetaLatticeSites()
{
  // betaSites = allSites - alphaSites
  for (size_t id = 0; id < numLattice; id++)
  {
    // O(1) Time Complexity
    if (alphaLatticeSites.find(id) == alphaLatticeSites.end())
    {
      betaLatticeSites.emplace(id);
    }
  }
}

void B2Ordering::AddB2Precipitate(Config &config,
                                  const size_t numB2Centers,
                                  const vector<Element> &elementVector)
{
  if (elementVector.size() != 2)
  {
    throw std::runtime_error("B2 ordering requires exactly two elements");
  }

  // Extract the two elements;
  Element alphaElement = elementVector[0];
  Element betaElement = elementVector[1];

  size_t totalSites = config.GetNumAtoms();

  // random number generator which will select where to insert the B2
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<size_t> selectB2Center(0, totalSites - 1);

  int numAlpha = totalSites / 2;
  int numBeta = totalSites / 2;

  for (size_t i = 0; i < numB2Centers; i++)
  {
    size_t selectedLatticeId = selectB2Center(gen);

    // Assign alpha element to selected lattice Id
    config.SetElementOfLattice(selectedLatticeId, alphaElement);

    // Assign beta element to it neighbours
    for (auto const &firstNNId : config.GetNeighborLatticeIdVectorOfLattice(selectedLatticeId, 1))
    {
      config.SetElementOfLattice(firstNNId, betaElement);
    }

    for (auto const &secondNNId : config.GetNeighborLatticeIdVectorOfLattice(selectedLatticeId, 2))
    {
      config.SetElementOfLattice(secondNNId, alphaElement);
    }
  }
}

void B2Ordering::InitializeB2Structure(Config &config)
{
  auto atomVector = config.GetAtomVector();
  std::map<Element, int> elementCount;

  for (const auto &atom : atomVector)
  {
    elementCount[atom]++;
  }

  // Get available sites
  std::vector<size_t> alphaAvailable(alphaLatticeSites.begin(), alphaLatticeSites.end());
  std::vector<size_t> betaAvailable(betaLatticeSites.begin(), betaLatticeSites.end());

  for (auto &[element, count] : elementCount)
  {
    // Assign to alpha sites first
    while (count > 0 && !alphaAvailable.empty())
    {
      size_t siteId = alphaAvailable.back();
      alphaAvailable.pop_back();
      config.SetElementOfLattice(siteId, element);
      count--;
    }

    // Then assign to beta sites
    while (count > 0 && !betaAvailable.empty())
    {
      size_t siteId = betaAvailable.back();
      betaAvailable.pop_back();
      config.SetElementOfLattice(siteId, element);
      count--;
    }

    // If still atoms left, show warning
    if (count > 0)
    {
      std::cerr << "Warning: Not enough lattice sites to assign all atoms of element "
                << element.GetElementString() << std::endl;
    }
  }
}

bool B2Ordering::isB2Structure(const Config &config,
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
