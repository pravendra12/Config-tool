#include "GenerateNEBStructure.h"

void GenerateNEBStructure(
  const string &filename,
  const Config &config, 
  const pair<size_t, size_t> &latticeIdJumpPair)
{
  auto migratingAtom = config.GetElementOfLattice(latticeIdJumpPair.second);

  // Initial NEB Structure
  auto configInitial = config;
  configInitial.SetElementOfLattice(latticeIdJumpPair.first, Element("X"));

  // Final NEB Structure
  auto configFinal = config;

  configFinal.SetElementOfLattice(latticeIdJumpPair.first, migratingAtom);
  configFinal.SetElementOfLattice(latticeIdJumpPair.second, Element("X"));


  Config::WriteConfig(filename + "_initial.cfg.gz", configInitial);
  Config::WriteLAMMPSDataFile(filename + "_initial.data", configInitial);

  Config::WriteConfig(filename + "_final.cfg.gz", configFinal);
  Config::WriteLAMMPSDataFile(filename + "_final.data", configFinal);

}