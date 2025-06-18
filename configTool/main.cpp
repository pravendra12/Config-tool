#include "Home.h"
#include "GenerateNEBStructure.h"
#include "GetKRAEncoding.h"

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

int main()
{
  size_t sSize = 10;
  double latticeParam = 3.2;
  string baseElement = "Mo";
  string soluteElement = "Ta";

  vector<double> cutoffs = {3};

  auto cfg = Config::GeneratePristineSupercell(sSize, latticeParam, baseElement, "BCC");
  cfg.UpdateNeighborList(cutoffs);

  // cout << cfg.GetNeighborLatticeIdVectorOfLattice( 0, 1).size() << endl;

  auto centerId = cfg.GetCentralAtomLatticeId();
  auto latticeJumpPair = make_pair(centerId, cfg.GetNeighborLatticeIdVectorOfLattice(centerId, 1)[0]);

  // for (const auto id : cfg.GetNeighboringLatticeIdSetOfPair(latticeJumpPair, 1))
  // {
  //   cfg.SetElementOfLattice(id, Element(soluteElement));
  // }
  
  string ss = to_string(sSize);

  // cfg.SetElementOfLattice(latticeJumpPair.second, Element("X"));
  // cfg.SetElementOfLattice(latticeJumpPair.first, Element("X"));

  // Config::WriteConfig("//media/sf_Phd/nebBenchmark/TaW14/ss_" + ss + "/Config/TaW14_" +  ss + "x" + ss + "x" + ss + ".cfg.gz", cfg);
  cout << "Here 1" << endl;
  // GenerateNEBStructure("//media/sf_Phd/nebBenchmark/Ta14W_" +  ss + "x" + ss + "x" + ss, cfg, latticeJumpPair);
  //media/sf_Phd/nebBenchmark/TaW14

  map<Element, size_t> elementMap;

  // 1 to 4 are Nb, Mo, Ta and W.
  // vector<string> elementSetMLIP = {"Nb", "Mo", "Ta", "W"};
  // size_t idx = 1;
  // for (const auto element : elementSetMLIP)
  // {
  //   elementMap[Element(element)] = idx;
  //   idx++;
  // }

  Config::WriteLAMMPSDataFileCustom("//media/sf_Phd/nebBenchmark/test.lammps", cfg,elementMap);

}
*/