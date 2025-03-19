#ifndef LMC_CONSTANT_INCLUDE_PRINTUTILITY_H_
#define LMC_CONSTANT_INCLUDE_PRINTUTILITY_H_

#include <vector>
#include <string>
#include <iostream>
#include <unordered_map>
#include "eigen3/Eigen/Dense"

void print1DVector(const std::vector<double> &vec);
void print2DVector(const std::vector<std::vector<size_t>> &vec);

void print2DStringVector(const std::vector<std::vector<std::string>> &vec);

void print3DVector(const std::vector<std::vector<std::vector<size_t>>> &vec);

void printOneHotEncodeHashmap(const std::unordered_map<std::string, Eigen::RowVectorXd> &hashmap);

#endif // LMC_CONSTANT_INCLUDE_PRINTUTILITY_H_
