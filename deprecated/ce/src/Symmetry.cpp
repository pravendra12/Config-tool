#include "Symmetry.h"
#include "PrintUtility.h"



// returns the requivalent clusters for 3 fold rotation
std::vector<std::vector<size_t>>
Symmetry::GetEquivalentSites3Fold(const Config &config, 
                                  const size_t maxBondOrder)
{

  // this can be defined one does not depends on the jump pair
  

  std::pair<size_t, size_t> jumpPair = {0, config.GetNeighborLatticeIdVectorOfLattice(0, 1)[0]};

  auto symmetricallySortedVector = 
              config.GetSortedLatticeVectorStateOfPair(jumpPair, maxBondOrder);

  std::vector<size_t> equivalentSites;
  std::vector<std::vector<size_t>> equivalentSiteVector;
  
  size_t idx = 0;

  while (idx < symmetricallySortedVector.size())
  { 
    auto latticeId = symmetricallySortedVector[idx];

    auto directionVectorI = config.GetNormalizedDirection(latticeId, jumpPair.first);
    auto directionVectorII = config.GetNormalizedDirection(latticeId, jumpPair.second);

    if (directionVectorI == directionVectorII)
    {
      equivalentSites.emplace_back(idx);
      idx++;
    }
    else
    {
      // Sample in threes
      // Works only till 2nd NN
      for (size_t id = idx; id < idx+3; id++){
        equivalentSites.emplace_back(id);
      }
      idx += 3;
    }

    equivalentSiteVector.emplace_back(equivalentSites);
    equivalentSites = {};
  }

  return equivalentSiteVector;
}



//std::vector<std::vector<std::vector>>
PairClusters Symmetry::GetEquivalentPairClusters3Fold(const Config &config, 
                                                  size_t maxBondOrder)
{
  auto equivalentSites3F = GetEquivalentSites3Fold(config,
                                                   maxBondOrder);
  
  // Does not depends on the jump pair, it can be anything
  std::pair<size_t, size_t> jumpPair = {74, 54};


  auto symmetricallySortedVector = 
              config.GetSortedLatticeVectorStateOfPair(jumpPair, maxBondOrder);
  
  PairClusters clusters;

  // Symmetric Pairs (Pairs within each subvector)
  for (const auto& group : equivalentSites3F) 
  {
    if (group.size() > 1)
    { // Only form pairs if there are at least two elements

      std::vector<std::vector<size_t>> symPairs;

      for (size_t i = 0; i < group.size(); ++i) 
      {
        for (size_t j = i + 1; j < group.size(); ++j) 
        {
          size_t latticeId_i = symmetricallySortedVector[group[i]];
          size_t latticeId_j = symmetricallySortedVector[group[j]];
              
          auto bondOrder = config.GetDistanceOrder(latticeId_i,
                                                   latticeId_j);
          if (bondOrder <= 3)
          {
            std::vector<size_t> pairCluster = {group[i], group[j]};
            symPairs.push_back(pairCluster);
          }
        }
      }
          if (!symPairs.empty())
          {
            clusters.symmetricPairs.emplace_back(symPairs);
          } 
    }
  }

  // print3DVector(clusters.symmetricPairs);


  // Non-Symmetric Pairs (Pairs between elements from different subvectors)
  for (size_t i = 0; i < equivalentSites3F.size(); ++i) 
  {
    for (size_t j = i + 1; j < equivalentSites3F.size(); ++j) 
    {
      std::vector<std::vector<size_t>> nonSymPairs;
      
      for (size_t idxI : equivalentSites3F[i]) 
      {
        for (size_t idxJ : equivalentSites3F[j]) 
        {
          size_t latticeId_i = symmetricallySortedVector[idxI];
          size_t latticeId_j = symmetricallySortedVector[idxJ];
              
          auto bondOrder = config.GetDistanceOrder(latticeId_i,
                                                   latticeId_j);
                                                
          std::vector<size_t> pairCluster = {idxI, idxJ};

          if (bondOrder <= 3)
          {
            nonSymPairs.push_back(pairCluster);
          }
        }
      }
      if (!nonSymPairs.empty())
      {
        clusters.nonSymmetricPairs.push_back(nonSymPairs);
      }
    }
  }
  // //std::cout << "Non Symmetric Pairs" << std::endl;
  // print3DVector(clusters.nonSymmetricPairs);


  std::vector<std::vector<std::vector<size_t>>> singletClusters;
  singletClusters.emplace_back(equivalentSites3F);
  
  // //std::cout << "Singlet Clusters" << std::endl;
  // print3DVector(singletClusters);       
  
  return clusters;

}



/// One Hot Encoding HashMaps
/*
std::unordered_map<std::string, std::vector<double>> 
GetOneHotEncodeHashmap(const std::set<Element> &elementSet) 
{
  size_t numElements = elementSet.size();
  std::unordered_map<std::string, std::vector<double>> encodeDict;
  size_t ct1 = 0;

  for (const auto &element : elementSet) 
  {
    std::vector<double> elementEncode(numElements, 0);
    elementEncode[ct1] = 1.0;
    encodeDict[element.GetElementString()] = elementEncode;
    ++ct1;
  }

  size_t numPairs = numElements * numElements;
  size_t ct2 = 0;
  for (const auto &element1 : elementSet) 
  {
    for (const auto &element2 : elementSet) 
    {
      std::vector<double> elementEncode(numPairs, 0);
      elementEncode[ct2] = 1.0;
      encodeDict[element1.GetElementString() + 
                 element2.GetElementString()] = elementEncode;
      ++ct2;
    }
  }
  size_t numPairsSymmetry = (numElements + 1) * numElements / 2;
  size_t ct3 = 0;
  for (auto it1 = elementSet.cbegin(); it1 != elementSet.cend(); ++it1) {
    for (auto it2 = it1; it2 != elementSet.cend(); ++it2) {
      std::vector<double> elementEncode(numPairsSymmetry, 0);
      elementEncode[ct3] = 1.0;
      encodeDict[it1->GetElementString() + '-' + 
                 it2->GetElementString()] = elementEncode;
      ++ct3;
    }
  }

  // // printOneHotEncodeHashmap(encodeDict);
  return encodeDict;
}
*/

std::unordered_map<std::string, Eigen::RowVectorXd> 
GetOneHotEncodeHashmap(const std::set<Element> &elementSet) 
{
    size_t numElements = elementSet.size();
    std::unordered_map<std::string, Eigen::RowVectorXd> encodeDict;

    size_t ct1 = 0;
    for (const auto &element : elementSet) 
    {
        Eigen::RowVectorXd elementEncode = Eigen::RowVectorXd::Zero(numElements);
        elementEncode(ct1) = 1.0;
        encodeDict[element.GetElementString()] = elementEncode;
        ++ct1;
    }

    size_t numPairs = numElements * numElements;
    size_t ct2 = 0;
    for (const auto &element1 : elementSet) 
    {
        for (const auto &element2 : elementSet) 
        {
            Eigen::RowVectorXd elementEncode = Eigen::RowVectorXd::Zero(numPairs);
            elementEncode(ct2) = 1.0;
            encodeDict[element1.GetElementString() + element2.GetElementString()] = elementEncode;
            ++ct2;
        }
    }

    size_t numPairsSymmetry = (numElements + 1) * numElements / 2;
    size_t ct3 = 0;
    for (auto it1 = elementSet.cbegin(); it1 != elementSet.cend(); ++it1) 
    {
        for (auto it2 = it1; it2 != elementSet.cend(); ++it2) 
        {
            Eigen::RowVectorXd elementEncode = Eigen::RowVectorXd::Zero(numPairsSymmetry);
            elementEncode(ct3) = 1.0;
            encodeDict[it1->GetElementString() + '-' + it2->GetElementString()] = elementEncode;
            ++ct3;
        }
    }

    return encodeDict;
}

// Given the Singlets, both type of pairs sum up the equivalent one hot encodings

std::vector<double>
GetEncodeVector3F(const Config &config,
                  std::pair<size_t, size_t> jumpPair,
                  std::set<Element> elementSet,
                  size_t maxBondOrder)
{
  auto singletsClusters = Symmetry::GetEquivalentSites3Fold(config, 
                                                            maxBondOrder);


  // print2DVector(singletsClusters);

  auto symmetricallySortedVector = 
              config.GetSortedLatticeVectorStateOfPair(jumpPair, maxBondOrder);

  std::vector<Element> elementVector;
  for (auto id : symmetricallySortedVector)
  {
    auto element = config.GetElementOfLattice(id);
    elementVector.emplace_back(element);
    //std::cout << element.GetElementString() << std::endl;
  }


  auto oneHotEncodeHashMap = GetOneHotEncodeHashmap(elementSet);

  // Combined Vector
  std::vector<double> encodeVector;

  // //std::cout << "--------------Singlets---------------------" << std::endl;
  //std::cout << singletsClusters.size() << std::endl;
  for (auto equivalentSites : singletsClusters)
  {

    // For now since we have two elements assignig it to 2 but need to  
    // generalize it
    size_t sizeEncodingSinglets = 2;


    Eigen::RowVectorXd tempEncodeVector = Eigen::RowVectorXd::Zero(sizeEncodingSinglets);
    size_t equivalentClusterCount = 0;

    //std::cout << equivalentSites.size() << std::endl;

    for (auto site : equivalentSites)
    {
      auto siteElement = elementVector[site];
      auto oneHotEncoding = oneHotEncodeHashMap[siteElement.GetElementString()];

      //std::cout << site << "\t" << siteElement.GetElementString() << "\t";
      
      //std::cout << "{ " << oneHotEncoding << " }" << std::endl;
      //std::cout << tempEncodeVector << std::endl;
      
      tempEncodeVector += oneHotEncoding;
      equivalentClusterCount++;

    }
    // //std::cout << "---------------------------------------" << std::endl;
    
    // //std::cout << tempEncodeVector << std::endl;
    // //std::cout << equivalentClusterCount << std::endl;

    tempEncodeVector /= equivalentClusterCount;

    // //std::cout << "{ " << tempEncodeVector  << " }" << std::endl;

    // Convert Eigen::RowVectorXd to std::vector<double> and append to encodeVector
    encodeVector.insert(encodeVector.end(), 
                        tempEncodeVector.data(), 
                        tempEncodeVector.data() + tempEncodeVector.size());

    // //std::cout << "---------------------------------------" << std::endl;
  }
  

  auto pairClusters = Symmetry::GetEquivalentPairClusters3Fold(config, maxBondOrder);

  // printOneHotEncodeHashmap(oneHotEncodeHashMap);

  auto symmetricPairs = pairClusters.symmetricPairs;

  // //std::cout << encodeVector.size() << std::endl;

  
  // Symmetric Pairs
  for (auto &equivalentPairs : symmetricPairs)
  {
    // Since symmetric pairs
    size_t sizeEncodeSymmPairs = 3;

    Eigen::RowVectorXd tempEncodeVector = Eigen::RowVectorXd::Zero(sizeEncodeSymmPairs);

    size_t equivalentClusterCount = 0;

    for (auto &pair : equivalentPairs)
    {
      if (pair.size() != 2) 
      {
        // //std::cout << "Not a Pair" << std::endl;
        exit(EXIT_FAILURE);  // Properly exits the program
      }
      else
      {
        auto element1 = elementVector[pair[0]].GetElementString();
        auto element2 = elementVector[pair[1]].GetElementString();
        
        std::string elementPair;
        if (element1 < element2)
        {
          elementPair = element1 + "-" + element2;
        }
        else
        {
          elementPair = element2 + "-" + element1;
        }

        auto oneHotEncoding = oneHotEncodeHashMap[elementPair];

        // Print Statement

        // //std::cout << "{ " << pair[0] << element1 << " " << pair[1] << element2 << "} ";

        // //std::cout << "{ ";
        for (auto i : oneHotEncoding)
        {
          // //std::cout << i << " ";
        }
        // //std::cout << "}" << std::endl;

        tempEncodeVector += oneHotEncoding;
        equivalentClusterCount++;        
      }
    }
    // //std::cout << "-------------------" << std::endl;
    
    // //std::cout << tempEncodeVector << std::endl;
    // //std::cout << equivalentClusterCount << std::endl;

    tempEncodeVector /= equivalentClusterCount;

    // //std::cout << tempEncodeVector << std::endl;

    encodeVector.insert(encodeVector.end(), 
                        tempEncodeVector.data(), 
                        tempEncodeVector.data() + tempEncodeVector.size());

    // //std::cout << "-------------------" << std::endl;
  }

  // //std::cout << encodeVector.size() << std::endl;


  auto nonSymmetricPairs = pairClusters.nonSymmetricPairs;

  // //std::cout << "---------------Non Symmetric Pairs----------------" << std::endl;
 
  // Non Symmetric Pairs
  for (auto &equivalentPairs : nonSymmetricPairs)
  {
    // Since non symmetric pairs // Need to work this out
    size_t sizeEncodeNonSymmPairs = 4;

    Eigen::RowVectorXd tempEncodeVector = Eigen::RowVectorXd::Zero(sizeEncodeNonSymmPairs);

    size_t equivalentClusterCount = 0;

    for (auto &pair : equivalentPairs)
    {
      if (pair.size() != 2) 
      {
        // //std::cout << "Not a Pair" << std::endl;
        exit(EXIT_FAILURE);  // Properly exits the program
      }
      else
      {
        auto element1 = elementVector[pair[0]].GetElementString();
        auto element2 = elementVector[pair[1]].GetElementString();
        
        std::string elementPair = element1 + element2;

        auto oneHotEncoding = oneHotEncodeHashMap[elementPair];

        // Print Statement

        // //std::cout << "{ " << pair[0] << element1 << " " << pair[1] << element2 << "} ";

        // //std::cout << "{ ";
        for (auto i : oneHotEncoding)
        {
          // //std::cout << i << " ";
        }
        // //std::cout << "}" << std::endl;

        tempEncodeVector += oneHotEncoding;
        equivalentClusterCount++;        
      }
    }
    // //std::cout << "-------------------" << std::endl;
    
    // //std::cout << tempEncodeVector << std::endl;
    // //std::cout << equivalentClusterCount << std::endl;

    tempEncodeVector /= equivalentClusterCount;

    // //std::cout << tempEncodeVector << std::endl;

    encodeVector.insert(encodeVector.end(), 
                        tempEncodeVector.data(), 
                        tempEncodeVector.data() + tempEncodeVector.size());

    // //std::cout << "-------------------" << std::endl;
  }

  // //std::cout << encodeVector.size() << std::endl;
  return encodeVector;

}


// Pairs between Migrating Atom and Nearest Neigbours

std::vector<double> GetEncodingMigratingAtomPairs (const Config &config,
                                                   size_t migratingAtomId,
                                                   std::pair<size_t, size_t> jumpPair,
                                                   std::set<Element> elementSet, 
                                                   size_t maxBondOrder)
{
  // maxBondOrder is used for getting the NN for the jumpPair
  // And assuming the migrating is at transition site

  auto equivalentSitesClusters = Symmetry::GetEquivalentSites3Fold(config, 
                                                           maxBondOrder);

  auto symmetricallySortedVector = 
              config.GetSortedLatticeVectorStateOfPair(jumpPair, maxBondOrder);

  auto GetOneHotEncoding = GetOneHotEncodeHashmap(elementSet);


  // for bondOrder 1, there are 6 equivalent groups
  // the two central groups will be the first nn of the migrating atom (6 atoms)
  // then the next to each of the central groups will be the second nn atom
  // for the migrating atom at transition site
  // and third nn for the migrating atom based on the cutoff described
  // so in the first nn shell of the jump pair, there will be 14 nn of the migrating Atom 

  // print2DVector(equivalentSitesClusters);
  
  auto migratingAtom = config.GetElementOfLattice(migratingAtomId);
  
  // Pair of Migrating Atom and the NN
  std::vector<std::vector<std::string>> equivalentPairsClusters;

  std::vector<double> encodeVector;

  for (auto &equiivalentSites : equivalentSitesClusters)
  {
    std::vector<std::string> equivalentPairs;
    
    // Since non symmetric pairs 
    // Need to work on this so as to generalize this
    size_t sizeEncodeNonSymmPairs = 4;

    Eigen::RowVectorXd tempEncodeVector = Eigen::RowVectorXd::Zero(sizeEncodeNonSymmPairs);

    size_t equivalentClusterCount = 0;

    for (auto &sites : equiivalentSites)
    {
      auto latticeId = symmetricallySortedVector[sites];
      auto neighborElement = config.GetElementOfLattice(latticeId);

      std::string elementPair = migratingAtom.GetElementString() + 
                  neighborElement.GetElementString();

      equivalentPairs.emplace_back(elementPair);

      tempEncodeVector += GetOneHotEncoding[elementPair];
    }

    equivalentPairsClusters.emplace_back(equivalentPairs);

    encodeVector.insert(encodeVector.end(), 
    tempEncodeVector.data(), 
    tempEncodeVector.data() + tempEncodeVector.size());
  }

  // print2DStringVector(equivalentPairsClusters);

  // print1DVector(encodeVector);


  return encodeVector;
  

}































