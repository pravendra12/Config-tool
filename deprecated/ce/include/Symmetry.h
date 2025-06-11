/**************************************************************************************************
 * Copyright (c) 2023-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 7/18/23 4:13 PM                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 7/18/23 4:21 PM                                                           *
 **************************************************************************************************/

#ifndef LMC_CE_INCLUDE_SYMMETRY_H_
#define LMC_CE_INCLUDE_SYMMETRY_H_

#include "Config.h"

struct PairClusters 
{
      std::vector<std::vector<std::vector<size_t>>> symmetricPairs;
      std::vector<std::vector<std::vector<size_t>>> nonSymmetricPairs;
};

class Symmetry 

{


  public:
    static std::vector<std::vector<size_t>> GetEquivalentSites3Fold(const Config &config, const size_t maxBondOrder);
    // static std::vector<std::vector<std::vector>> GetEquivalentPairClusters3Fold(const Config config, size_t maxBondOrder);
    static PairClusters GetEquivalentPairClusters3Fold(const Config &config, size_t maxBondOrder);

  
};

// std::unordered_map<std::string, Eigen::RowVectorXd>  GetOneHotEncodeHashmap(const std::set<Element> &elementSet);

std::vector<double> GetEncodeVector3F(const Config &config, 
                                      std::pair<size_t, size_t> jumpPair, 
                                      std::set<Element> elementSet, 
                                      size_t maxBondOrder);

std::vector<double> GetEncodingMigratingAtomPairs(const Config &config,
                                                  size_t migratingAtomId, 
                                                  std::pair<size_t, size_t> jumpPair, 
                                                  std::set<Element> elementSet, 
                                                  size_t maxBondOrder);


#endif // LMC_CE_INCLUDE_SYMMETRY_H_

// void GetEncodeVector3F(const Config &config, size_t maxBondOrder);
