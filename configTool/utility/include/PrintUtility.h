#ifndef LMC_CONSTANT_INCLUDE_PRINTUTILITY_H_
#define LMC_CONSTANT_INCLUDE_PRINTUTILITY_H_

#include <vector>
#include <string>
#include <iostream>
#include <unordered_map>
#include "eigen3/Eigen/Dense"

template <typename T>
void print1DVector(const std::vector<T> &vec);

template <typename T>
void print2DVector(const std::vector<std::vector<T>> &vec);

template <typename T>
void print3DVector(const std::vector<std::vector<std::vector<T>>> &vec);

#include "PrintUtility.inl"

#endif // LMC_CONSTANT_INCLUDE_PRINTUTILITY_H_
