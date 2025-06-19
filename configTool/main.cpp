#include "Home.h"

int main(int argc, char *argv[]) {
  if (argc == 1) {
    cout << "No input parameter filename." << endl;
    return 1;
  }
  api::Parameter parameter(argc, argv);
  api::Print(parameter);
  api::Run(parameter);
}

/*
using namespace std;
#include "ConfigEncoding.h"
#include "GenerateStructuresWithB2.h"
#include "PotentialEnergyEstimator.h"
#include "GenerateStructureCNT.h"
#include <chrono>

namespace fs = std::filesystem;

int main()
{
  vector<string> elementVector = {"Ta", "W"};
  vector<double> compositionVector = {50, 50};

  string predictorFilename = "/home/pravendra3/Documents/Config-tool/bin/predictorFileKRA_BO2_WTa.json";

  auto config = Config::GenerateAlloySupercell(10, 3.4, "BCC", elementVector, compositionVector, 1);
  config.UpdateNeighborList({3.3, 4.7, 5.6});

  auto trainingConfig = Config::GeneratePristineSupercell(5, 3.4, "X", "BCC");
  trainingConfig.UpdateNeighborList({3.3, 4.7, 5.6});

  set<Element> elementSet = {Element("Ta"), Element("W"), Element("X")};

  GenerateStructureCNT cnt(
      predictorFilename,
      config,
      trainingConfig,
      10,
      elementSet,
      compositionVector);

  auto start = std::chrono::high_resolution_clock::now();

  double energy = cnt.GetEnergyOfConfig(config);

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;

  std::cout << "Energy: " << energy << std::endl;
  std::cout << "Elapsed time: " << elapsed.count() << " seconds" << std::endl;

  PotentialEnergyEstimator peEstimator(
      predictorFilename,
      config,
      trainingConfig,
      elementSet);

  start = std::chrono::high_resolution_clock::now();

  double energyPE = peEstimator.GetEnergy(config);

  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;

  std::cout << "Energy: " << energyPE << std::endl;
  std::cout << "Elapsed time due to energy from PE: " << elapsed.count() << " seconds" << std::endl;
}
*/

/*
int main() {
    std::vector<size_t> supercellSizes = {5, 7, 10};
    double latticeParam = 3.319;
    std::string baseElement = "Ta";
    std::string soluteElement = "W";
    size_t numSolute = 14;
    std::vector<double> cutoffs = {3.2};

    std::string basePath = "//media/sf_Phd/nebBenchmark/MLIP/";
    std::string systemFolder = baseElement + soluteElement + std::to_string(numSolute);

    // Define element-to-type mapping
    std::map<Element, size_t> elementMap;
    std::vector<std::string> elementSetMLIP = {"Nb", "Mo", "Ta", "W"};
    size_t idx = 1;
    for (const auto& element : elementSetMLIP) {
        elementMap[Element(element)] = idx++;
    }

    for (size_t sSize : supercellSizes) {
        std::string fullPath = basePath + systemFolder + "/ss_" + std::to_string(sSize) + "/Config/";
        fs::create_directories(fullPath);

        // Generate pristine BCC supercell
        auto cfg = Config::GeneratePristineSupercell(sSize, latticeParam, baseElement, "BCC");
        cfg.UpdateNeighborList(cutoffs);

        // Define the jump pair (central atom to 1st neighbor)
        auto centerId = cfg.GetCentralAtomLatticeId();
        auto neighborLattices = cfg.GetNeighborLatticeIdVectorOfLattice(centerId, 1);
        auto latticeJumpPair = std::make_pair(centerId, neighborLattices[0]);

        // Introduce solute
        // cfg.SetElementOfLattice(latticeJumpPair.second, Element(soluteElement));
        for (const auto id : cfg.GetSortedLatticeVectorStateOfPair(latticeJumpPair, 1))
        {
          cfg.SetElementOfLattice(id, Element(soluteElement));
        }

        // Construct filename
        std::string filename = baseElement + soluteElement + std::to_string(numSolute) + "_" +
                               std::to_string(sSize) + "x" +
                               std::to_string(sSize) + "x" +
                               std::to_string(sSize);

        // Write NEB structure
        GenerateNEBStructure(fullPath + filename, cfg, latticeJumpPair, elementMap);

        std::cout << "NEB structure generated at: " << fullPath + filename << std::endl;
    }

    return 0;
}
*/