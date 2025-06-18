#include "GenerateNEBStructure.h"

void GenerateNEBStructure(
    const string &filename,
    Config &config,
    const pair<size_t, size_t> &latticeIdJumpPair,
    map<Element, size_t> &elementMap)
{
  auto previousAtom = config.GetElementOfLattice(latticeIdJumpPair.first);

  config.SetElementOfLattice(latticeIdJumpPair.first, Element("X"));

  // Initial NEB Structure
  Config::WriteConfig(filename + "_initial.cfg.gz", config);
  Config::WriteLAMMPSDataFileCustom(filename + "_initial.data", config, elementMap);

  // Final NEB Structure
  config.LatticeJump(latticeIdJumpPair);

  Config::WriteConfig(filename + "_final.cfg.gz", config);
  Config::WriteLAMMPSDataFileCustom(filename + "_final.data", config, elementMap);

  config.LatticeJump(latticeIdJumpPair);
  config.SetElementOfLattice(latticeIdJumpPair.first, previousAtom);
}