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
  vector<string> elementVector = {"Ta", "W"};
  vector<double> compositionVector = {50, 50};

  vector<double> cutoffs = {3.3, 4.7, 5.6};

  auto trainingSupercell = Config::GeneratePristineSupercell(5, 3.4, "X", "BCC");
  trainingSupercell.UpdateNeighborList(cutoffs);
  auto elementPair = make_pair(Element("W"), Element("Ta"));
  
  Config cfg;
  int idx = 0;
  do 
  {
    cfg = Config::GenerateAlloySupercell(sSize, 3.4, "BCC", elementVector, compositionVector, idx);
    idx++;
    cout << cfg.GetElementOfLattice(cfg.GetCentralAtomLatticeId()).GetElementString() << endl;
  }
  while (elementPair.first != cfg.GetElementOfLattice(cfg.GetCentralAtomLatticeId()));

  cfg.UpdateNeighborList(cutoffs);



  set<Element> elementSet = {Element("X"), Element("W"), Element("Ta")};

  return 0;

  GenerateStructureCNT cnt(
      "predictorFileKRA_BO2_WTa.json",
      cfg,
      trainingSupercell,
      10,
      elementSet,
      compositionVector);

  cnt.GenerateRandomStructures(2, true, true, "//media/sf_Phd/CNT/structures/WTa");
  cnt.GenerateStructureWithB2(
      cfg,
      1000,
      1,
      elementPair,
      false,
      "//media/sf_Phd/CNT/structures/WTa");



      return 0;
}*/