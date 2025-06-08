#include "PrintUtility.h"

template <typename T>
void print1DVector(const T &vec) 
{
  std::cout << "{";
  for (size_t i = 0; i < vec.size(); ++i) {
      std::cout << vec[i];
      if (i < vec.size() - 1) {
          std::cout << ", ";
      }
  }
  std::cout << ", }" << std::endl;  // Extra comma before closing brace
};

template <typename T>
void print2DVector(const std::vector<std::vector<T>> &vec)
{
    std::cout << "{\n";
    for (const auto &row : vec)
    {
        std::cout << "  {";
        for (size_t i = 0; i < row.size(); ++i)
        {
            std::cout << row[i] << ", ";
        }
        std::cout << "},\n";
    }
    std::cout << "}\n";
};

template <typename T>
void print3DVector(const std::vector<std::vector<std::vector<T>>> &vec)
{
    std::cout << "{\n";
    for (const auto &matrix : vec)
    { // Iterate over first dimension
        std::cout << "  {\n";
        for (const auto &row : matrix)
        { // Iterate over second dimension
            std::cout << "    {";
            for (size_t i = 0; i < row.size(); ++i)
            { // Iterate over third dimension
                std::cout << row[i] + 1;
                if (i < row.size() - 1)
                    std::cout << ", ";
            }
            std::cout << "},\n";
        }
        std::cout << "  },\n";
    }
    std::cout << "}\n";
}

