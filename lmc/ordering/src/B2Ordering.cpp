#include "B2Ordering.h"

B2Ordering::B2Ordering(const Config &config) : config(config)
{
  InitializeAlphaLatticeSites();
  InitializeBetaLatticeSites();
}

double B2Ordering::GetB2OrderParameter(const Element &element)
{
  // fractional occupancy of element at alpha and beta sites 
  auto alphaOccupancy = GetAlphaSiteOccupancy(element);
  auto betaOccupancy = GetBetaSiteOccupancy(element);

  double b2OrderParameter = (alphaOccupancy - betaOccupancy)/(alphaOccupancy + betaOccupancy);

  return b2OrderParameter;
}

double B2Ordering::GetAlphaSiteOccupancy(const Element &element)
{
  // Number of occupied alpha sites by element
  size_t numOccupiedASites = 0;
  // size_t numAlphaSites = alphaLatticeSites.size();

  for (const auto &id : alphaLatticeSites)
  {
    auto siteElement = config.GetElementOfLattice(id);

    if (siteElement == element)
    {
        numOccupiedASites += 1;
    }
  }

  double alphaSiteOccupancy = double(numOccupiedASites)/double(config.GetNumAtoms());

  return alphaSiteOccupancy;
}

double B2Ordering::GetBetaSiteOccupancy(const Element &element)
{
  // Number of occupied beta sites by element
  size_t numOccupiedBSites = 0;
  // size_t numBetaSites = betaLatticeSites.size();

  for (const auto &id : betaLatticeSites)
  {
    auto siteElement = config.GetElementOfLattice(id);

    if (siteElement == element)
    {
        numOccupiedBSites += 1;
    }
  }
  
  // will be between 0 and 0.5
  double betaSiteOccupancy = double(numOccupiedBSites)/double(config.GetNumAtoms());

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
    
  auto secondNNList = config.GetNeighborLists()[1];  
  
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
        if (config.GetDistanceOrder(id1, id3) == 1) 
        {
            id1Valid = false;
            break;
        }
      }

      // Check if id2 has any bond order of 1 with sites already in alphaLatticeSites
      for (size_t id3 : alphaLatticeSites) {
        if (config.GetDistanceOrder(id2, id3) == 1) {
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
