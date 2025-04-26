#include "Home.h"
#include "Config.h"
#include <iostream>
#include <algorithm>
#include <random>

#include <fstream>
#include <sstream>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <cmath>
#include "Symmetry.h"
#include <SymmetryCustom.h>
#include "PrintUtility.h"
#include "GenerateData.h"
#include "SymmetryCorrect.h"
#include <cmath>
using namespace std;


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
int main()
{

  // vector<Element> elementVector{"Al","W"};
  vector<size_t> ssVector = {5, 7, 10};
  double latticeParam = 3.243;
  const vector<double> cutoffs = {3.22, 4.5, 5.3};
  string structureType = "BCC";
  vector<string> elementVector = {"Ta", "W"};
  vector<double> compositionVector = {50, 50};
  bool isOrdered = true;
  size_t seed = 34432;

  
    auto cfg = Config::GenerateAlloySupercell(ssVector[0],
                                              latticeParam,
                                              structureType,
                                              elementVector,
                                              compositionVector,
                                              seed);

  cfg.UpdateNeighborList(cutoffs);
    
  cout << cfg.GetNeighborLatticeIdVectorOfLattice(0, 1).size() << endl;
  cout << cfg.GetNeighborLatticeIdVectorOfLattice(0, 2).size() << endl;
  cout << cfg.GetNeighborLatticeIdVectorOfLattice(0, 3).size() << endl;


  /*
    GenerateData dataGenerator(supercellSize,
      latticeParam,
      structureType,
      elementVector,
      compositionVector,
      cutoffs,
      seed,
      true);


  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<size_t> selectB2Center(1, 3);

  for (int i = 0; i < 10; i++)
  {
    cout << selectB2Center(gen)<< endl;
  }

  vector<Element>elementSet{Element("Ta"), Element("W")};



    auto b2Info1 = dataGenerator.AddB2Structure(1, elementSet);

    Config::WriteConfig("b2Type1.cfg", b2Info1.orderedConfig);

    auto b2Info2 = dataGenerator.AddB2Structure(2, elementSet);
    Config::WriteConfig("b2Type2.cfg", b2Info2.orderedConfig);


    auto b2Info3 = dataGenerator.AddB2Structure(3, elementSet);
    Config::WriteConfig("b2Type3.cfg", b2Info3.orderedConfig);

  */

  // dataGenerator.AddB2Ordering(6,{Element("Ta"), Element("W")} );


/*
int main(int argc, char *argv[]) {
  if (argc == 1) {
    cout << "No input parameter filename." << endl;
    return 1;
  }
  api::Parameter parameter(argc, argv);
  api::Print(parameter);
  api::Run(parameter);
}

*/
////////////// Testing the barrier encoding till 3rd NN ///////////////////////
/*
void writeEncodingToFileMain(ofstream &outFile,
                         const VectorXd &vec)
{
  for (size_t i = 0; i < vec.size(); ++i)
  {
    if (i == vec.size() - 1)
    {

      outFile << vec[i];
    }
    else
    {
      outFile << vec[i] << ", ";
    }
  }
}

int main()
{
  const vector<double> cutoffs = {3.3, 4.7, 5.6};
  size_t vacancyMigrationBO = 3;
  set<Element> elementSet;
  elementSet.insert(Element("Ta"));
  elementSet.insert(Element("W"));

  unordered_map<string, RowVectorXd> oneHotEncodingHashMap = GetOneHotEncodeHashmap(elementSet);

  bool isEquivalentSitesEncodingDeclared = false;
  vector<vector<size_t>> equivalentSitesEncoding;

  ofstream outFile("outputEncodingWTa_BO_3.txt");

  outFile << "FolderId\t"
          << "MigratingElement\t"
          << "VacancyId\t"
          << "MigratingElementId\t"
          << "ForwardEncoding\t"
          << "BackwardEncoding\t" << endl;

  // read the file

  std::string parentDir = "/media/sf_Phd/ActiveLearning/test1/selectedFolder"; // Set your parent directory path

  // Iterate over all directories inside parentDir
  for (const auto &subdir : fs::directory_iterator(parentDir))
  {
    if (fs::is_directory(subdir))
    {
      fs::path configDir = subdir.path() / "Config";

      if (fs::exists(configDir) && fs::is_directory(configDir))
      {
        for (const auto &entry : fs::directory_iterator(configDir))
        {
          if (entry.is_regular_file() && entry.path().extension() == ".cfg")
          {
            std::string configPath = entry.path().string();
            std::cout << "Reading config: " << configPath << std::endl;

            auto cfg = Config::ReadCfg(configPath);
            cfg.UpdateNeighborList(cutoffs);

            size_t vacancyId = cfg.GetCentralAtomLatticeId();
            size_t migratingElementLatticeId = cfg.GetNeighborLatticeIdVectorOfLattice(vacancyId, 1)[0];

            auto migratingAtomElement = cfg.GetElementOfLattice(migratingElementLatticeId);

            pair<size_t, size_t> forwardJumpPair = {migratingElementLatticeId, vacancyId};
            pair<size_t, size_t> backwardJumpPair = {vacancyId, migratingElementLatticeId};

            if (!isEquivalentSitesEncodingDeclared)
            {

              equivalentSitesEncoding = GetEquivalentSitesUnderKFoldRotation(cfg,
                                                                             vacancyMigrationBO,
                                                                             6);
              isEquivalentSitesEncodingDeclared = true;
            }

            // Forward Jump
            vector<size_t> forwardSymmSortedVector = GetSortedLatticeVectorStateOfPair(cfg,
                                                                                       forwardJumpPair,
                                                                                       vacancyMigrationBO);

            VectorXd forwardEncodeVector = GetEncodingMigratingAtomPair(cfg,
                                                                        equivalentSitesEncoding,
                                                                        forwardSymmSortedVector,
                                                                        oneHotEncodingHashMap,
                                                                        migratingAtomElement);

            // Backward Jump
            vector<size_t> backwardSymmSortedVector = GetSortedLatticeVectorStateOfPair(cfg,
                                                                                        backwardJumpPair,
                                                                                        vacancyMigrationBO);

            VectorXd backwardEncodeVector = GetEncodingMigratingAtomPair(cfg,
                                                                         equivalentSitesEncoding,
                                                                         backwardSymmSortedVector,
                                                                         oneHotEncodingHashMap,
                                                                         migratingAtomElement);

            outFile << subdir.path().filename().string() << "\t";

            outFile << migratingAtomElement.GetElementString() << "\t"
                    << vacancyId << "\t"
                    << migratingElementLatticeId << "\t";

            writeEncodingToFileMain(outFile, forwardEncodeVector);
            outFile << "\t";

            writeEncodingToFileMain(outFile, backwardEncodeVector);
            outFile << "\t";

            outFile << endl;
          }
        }
      }
    }
  }
}

/*
int main()
{
  size_t supercellSize = 10;
  double latticeParam = 3.4;
  string structureType = "BCC";
  vector<string> elementVector = {"Ta", "W"};
  vector<double> compositionVector = {50, 50};
  const vector<double> cutoffs = {3.3, 4.7, 5.6};

  auto cfg = Config::GenerateAlloySupercell(supercellSize,
                                            latticeParam,
                                            structureType,
                                            elementVector,
                                            compositionVector,
                                            1);

  cfg.UpdateNeighborList(cutoffs);

  auto encodingVectorBO2 = GetEquivalentSites3Fold(cfg, 2);

  print2DVector(encodingVectorBO2);

  auto encodingVectorBO3 = GetEquivalentSites3Fold(cfg, 3);

  print2DVector(encodingVectorBO3);

  auto centralAtomId = cfg.GetCentralAtomLatticeId();
  auto neighbourAtomId = cfg.GetNeighborLatticeIdVectorOfLattice(centralAtomId, 1)[0];

  Element vacancy("X");
  cfg.SetElementOfLattice(centralAtomId, vacancy);
  cfg.SetElementOfLattice(neighbourAtomId, vacancy);

  pair<size_t, size_t> latticeIdPair = {centralAtomId, neighbourAtomId};

  auto ssVectorIV = GetSortedLatticeVectorStateOfPair(cfg, latticeIdPair, 3);

  // second nn of the jump pair
  Element fourthNNElement("V");
  for (auto id : ssVectorIV)
  {
    cfg.SetElementOfLattice(id, fourthNNElement);
  }

  auto ssVectorIII = GetSortedLatticeVectorStateOfPair(cfg, latticeIdPair, 3);

  // second nn of the jump pair
  Element thirdNNElement("Mg");
  for (auto id : ssVectorIII)
  {
    cfg.SetElementOfLattice(id, thirdNNElement);
  }

  auto ssVectorII = GetSortedLatticeVectorStateOfPair(cfg, latticeIdPair, 2);

  // second nn of the jump pair
  Element secondNNElement("Al");
  for (auto id : ssVectorII)
  {

    cfg.SetElementOfLattice(id, secondNNElement);
  }

  auto ssVectorI = GetSortedLatticeVectorStateOfPair(cfg, latticeIdPair, 1);

  // first nn of the jump pair
  Element firstNNElement("W");
  for (auto id : ssVectorI)
  {
    cfg.SetElementOfLattice(id, firstNNElement);
  }

  // Config::WriteConfig("vacancyMigration.cfg", cfg);

  // Config::WriteLAMMPSDataFile("vacancyMigration.data", cfg);

  print1DVector(ssVectorII);

  print1DVector(ssVectorIII);

  cout << cfg.GetBasis() << endl;

  Eigen::Vector3d point{18.7, 11.9, 18.7};
  Eigen::Vector3d axis{-1, 1, 1};
  Eigen::Vector3d center{17.85, 16.15, 16.15};

  Eigen::Vector3d point1{18.7, 11.9, 18.7};

  getEquivalentPoints(point, axis, 60, center);

  auto encoding6F = GetEquivalentSitesUnderKFoldRotation(cfg, 3, 6);

  cout << "Size of IV encoding: " << ssVectorIV.size() << endl;

  auto cfgWTa = Config::GenerateAlloySupercell(supercellSize,
                                               latticeParam,
                                               structureType,
                                               elementVector,
                                               compositionVector,
                                               1);

  set<Element> elementSet;

  for (auto ele : elementVector)
  {
    elementSet.insert(Element(ele));
  }

  auto oneHotEncodingMap = GetOneHotEncodeHashmap(elementSet);

  auto encodingMigratingAtomPairs = GetEncodingMigratingAtomPair(cfgWTa,
                                                                 encoding6F,
                                                                 ssVectorIII,
                                                                 oneHotEncodingMap,
                                                                 firstNNElement);


  cout << encodingMigratingAtomPairs.transpose() << endl;
}

/*
int main()
{
  size_t supercellSize = 30;
  double latticeParam = 3.4;
  string structureType = "BCC";
  vector<string> elementVector = {"W", "Ta"};
  vector<double> compositionVector = {50, 50};

  auto cfg = Config::GenerateAlloySupercell(supercellSize,
                                            latticeParam,
                                            structureType,
                                            elementVector,
                                            compositionVector,
                                            1);
  Element vacancy("X");
  cfg.SetElementOfLattice(0, vacancy);

  Config::WriteConfig("start_W50Ta50_30x30x30.cfg", cfg);
}

/*
int main()
{
  size_t supercellSize = 5;
  double latticeParam = 3.4;
  string structureType = "BCC";
  vector<string> elementVector = {"Al", "W"};
  vector<double> elementComposition = {90, 10};
  vector<double> cutoffs = {3.3, 4.7, 5.6};
  size_t maxBondOrder = 3;
  size_t maxClusterSize = 3;

  string outFileName = "outputFile_";

  for (int i = 0; i < elementVector.size(); i++)
  {
    outFileName += elementVector[i];
    outFileName += to_string(int(elementComposition[i]));
  }
  outFileName += ".txt";

  ofstream outFile;
  outFile.open(outFileName);

  size_t unqiueConfigId = 2; // seed

  GenerateData dataGenerator(supercellSize,
                             latticeParam,
                             structureType,
                             elementVector,
                             elementComposition,
                             cutoffs,
                             unqiueConfigId,
                             false);

  auto isSaved = dataGenerator.saveConfig();
  auto vacMigInfo = dataGenerator.generateNEBStructure();

  // elementVector for CE
  vector<string> elementVectorCE = elementVector;
  elementVectorCE.push_back("X");

  set<Element> elementSetCE(elementVectorCE.begin(), elementVectorCE.end());
  auto ceInfo = dataGenerator.generateCEData(elementSetCE,
                                             maxBondOrder,
                                             maxClusterSize);

  B2OrderingInfo disOrderedStructure;
  disOrderedStructure.isOrderedStructure = false;

  dataGenerator.writeVacancyJumpCeConfig(outFile,
                                         vacMigInfo,
                                         ceInfo,
                                         disOrderedStructure,
                                         elementSetCE);

  Element elementAtAlphaSite("Al");
  Element elementAtBetaSite("W");

  auto orderingInfo = dataGenerator.generateB2Structure(elementComposition, elementAtAlphaSite, elementAtBetaSite);

  GenerateData dataGeneratorOrdered(orderingInfo.orderedConfig,
                                    supercellSize,
                                    elementVector,
                                    elementComposition,
                                    cutoffs,
                                    unqiueConfigId,
                                    true);

  auto isSavedOrdered = dataGeneratorOrdered.saveConfig();
  // auto vacMigInfo = dataGenerator.generateNEBStructure();

  auto vacMigInfoOrdered = dataGeneratorOrdered.generateNEBStructure();

  auto ceInfoOrdered = dataGeneratorOrdered.generateCEData(elementSetCE,
                                                           maxBondOrder,
                                                           maxClusterSize);

  dataGeneratorOrdered.writeVacancyJumpCeConfig(outFile,
                                                vacMigInfoOrdered,
                                                ceInfoOrdered,
                                                orderingInfo,
                                                elementSetCE);

  // 10 percent fraction
  // Config::WriteConfig("ordered_Al50W50_2.cfg", orderingInfo.orderedConfig);

  // auto encode =  dataGenerator.getMigratingAtomNeighborPairEncodeVector(elementSetCE);

  // cout << encode.first.transpose() << endl;

  outFile.close();
}
*/
/*
struct ClusterExpansionInfo
{
  string FolderId;
  string system;
  size_t initialVacancyId;
  Eigen::Vector3d initialVacancyPosition;
  size_t finalVacancyId;
  Eigen::Vector3d finalVacancyPosition;
  string migratingAtomElement;
  Eigen::VectorXd ceEncodingStart;
  Eigen::VectorXd ceEncodingEnd;

};

void writeEncodingToFile(string filename,
                         ClusterExpansionInfo &ceInfo)
{
  std::ofstream outputFile(filename, std::ios::app);

  if (outputFile.is_open())
  {
    // Check if file is empty to write header only once
    if (outputFile.tellp() == 0)
    {
        // Writing headers (column names) only if file is empty
        outputFile << "FolderId\t"
                   << "System\t"
                   << "Initial Vacancy ID\t"
                   << "Initial Vacancy Position\t"
                   << "Final Vacancy ID\t"
                   << "Final Vacancy Position\t"
                   << "Migrating Atom\t"
                   << "Cluster Expansion Encoding Start\t"
                   << "Cluster Expansion Encoding End"
                   << std::endl;
    }
    outputFile << ceInfo.FolderId << "\t";
    outputFile << ceInfo.system << "\t";
    outputFile << ceInfo.initialVacancyId << "\t";
    outputFile << ceInfo.initialVacancyPosition.transpose() << "\t";
    outputFile << ceInfo.finalVacancyId << "\t";
    outputFile << ceInfo.finalVacancyPosition.transpose() << "\t";
    outputFile << ceInfo.migratingAtomElement << "\t";

    // Write the ceEncodingStart vector
    for (size_t i = 0; i < ceInfo.ceEncodingStart.size(); ++i)
    {
        if (i == ceInfo.ceEncodingStart.size()-1)
        {
            outputFile << ceInfo.ceEncodingStart[i];
        }
        else
        {
            outputFile << ceInfo.ceEncodingStart[i] << ", ";
        }
    }
    outputFile << "\t";

    // Write the ceEncodingEnd vector
    for (size_t i = 0; i < ceInfo.ceEncodingEnd.size(); ++i)
    {
        if (i == ceInfo.ceEncodingEnd.size()-1)
        {
            outputFile << ceInfo.ceEncodingEnd[i];
        }
        else
        {
            outputFile << ceInfo.ceEncodingEnd[i] << ", ";

        }
    }
    outputFile << std::endl;

    outputFile.close();
  }
  else
  {
      std::cout << "Unable to open file for writing." << std::endl;
  }

};

int main()
{

  const vector<double> cutoffs = {3.3, 4.7, 5.6};

  /// Generating a supercell of size 20 x 20 x 20
  auto supercellWTa = Config::GenerateAlloySupercell(20, 3.4, "BCC",
                                                    {"W", "Ta"},
                                                    {50, 50},
                                                    43);
  supercellWTa.UpdateNeighborList(cutoffs);

  Element vac("X");

  supercellWTa.SetElementOfLattice(123, vac);

  Config::WriteConfig("start_W50Ta50_20x20x20.cfg", supercellWTa);


  size_t supercellSize = 5;

  std::vector<std::string> internal_folders = {};


  string base_path = "/media/sf_Phd/WTaNEB/";


  for (const auto& internal_folder : internal_folders)
  {
      for (int i = 1; i <= 1; ++i) {
          std::string folder_number = (i < 10) ? "0" + std::to_string(i) : std::to_string(i);
          std::string file_path = base_path + internal_folder + "/" + folder_number + "/Config/" + internal_folder + "_5x5x5.cfg";

          if (filesystem::exists(file_path))
          {
              // std::cout << "Reading file: " << file_path << std::endl;
              std::cout << file_path << std::endl;


  // Read the configuration file
  auto cfg = Config::ReadCfg(file_path);

  cfg.UpdateNeighborList(cutoffs);

  auto centralAtomId = cfg.GetCentralAtomLatticeId();
  auto neighoursOfCentralAtom = cfg.GetNeighborLatticeIdVectorOfLattice(centralAtomId, 1);

  auto neighbourAtom = neighoursOfCentralAtom[0];


  // PE estimator initialization

  auto atomVector = cfg.GetAtomVector();
  set<Element> elementSet(atomVector.begin(), atomVector.end());

  size_t maxBondOrderCE = 3;
  size_t maxClusterSizeCE = 3;

  Element vacancy("X");

  elementSet.emplace(vacancy);

  // Migrating Element should be same in initial and final config
  auto migratingElement = cfg.GetElementOfLattice(neighbourAtom);

  PotentialEnergyEstimator peEstimator("predictor_file.json",
                                       cfg,
                                       cfg,
                                       elementSet,
                                       maxClusterSizeCE,
                                       maxBondOrderCE);

  auto encode = peEstimator.GetEncodeVector(cfg);


  // Initial Configuration
  // Central Atom Id : Initial Vacancy Id
  // Its Neighbor : Initial Migraiting Atom Id
  // CE encoding for this has already been saved

  size_t initialVacancyId = centralAtomId;
  size_t initialMigratingAtomId = neighbourAtom;


  auto initialCfg = cfg;
  initialCfg.UpdateNeighborList(cutoffs);

  initialCfg.SetElementOfLattice(initialVacancyId, vacancy);

  Eigen::VectorXd encodeCEInitial = peEstimator.GetEncodeVector(initialCfg);

  /*
  cout << "Initial Configuration" << endl;
  cout << "Vacancy Id: " << initialVacancyId << endl;
  cout << "Migrating Id: " << initialMigratingAtomId << endl;

  cout << "Element of Central Atom Id in Initial Config: "
       << initialCfg.GetElementOfLattice(centralAtomId) << endl;

  cout << "Element of Neighbouring Atom Id Initial Config: "
       << initialCfg.GetElementOfLattice(neighbourAtom) << endl;

  cout << "----------------------------------------------------------" << endl;
  */

// Final Configuration
// Central Atom Id : Migrating Atom
// Its Neighbor : Vacancy
/*
size_t finalVacancyId = neighbourAtom;
size_t finalMigratingAtomId = centralAtomId;
auto finalCfg = cfg;
finalCfg.UpdateNeighborList(cutoffs);

// Migrating Atom to the previous vacant site
finalCfg.SetElementOfLattice(finalMigratingAtomId, migratingElement);

// Making a vacancy at neighbourAtom
finalCfg.SetElementOfLattice(finalVacancyId, vacancy);

Eigen::VectorXd encodeCEFinal = peEstimator.GetEncodeVector(finalCfg);


cout << "Final Configuration" << endl;
cout << "Vacancy Id: " << finalVacancyId << endl;
cout << "Migrating Id: " << finalMigratingAtomId << endl;


cout << "Element of Central Atom Id in Final Config: "
     << finalCfg.GetElementOfLattice(centralAtomId) << endl;

cout << "Element of Neighbouring Atom Id Final Config: "
     << finalCfg.GetElementOfLattice(neighbourAtom) << endl;
*/
/*
ClusterExpansionInfo ceInfo;
ceInfo.FolderId = folder_number;
ceInfo.system = internal_folder;
ceInfo.initialVacancyId = initialVacancyId;
ceInfo.initialVacancyPosition = cfg.GetCartesianPositionOfLattice(initialVacancyId);

ceInfo.finalVacancyId = finalVacancyId;
ceInfo.finalVacancyPosition = cfg.GetCartesianPositionOfLattice(finalVacancyId);

ceInfo.migratingAtomElement = migratingElement.GetElementString();

ceInfo.ceEncodingStart = encodeCEInitial;
ceInfo.ceEncodingEnd = encodeCEFinal;

// writeEncodingToFile("ceEncodingOutput.txt", ceInfo);
}
else
{
  std::cout << "File not found: " << file_path << std::endl;
}
}
}
}
*/

/*
int main(int argc, char* argv[])
{
  const vector<double> cutoffs = {3.3, 4.7, 5.6};
  auto cfg = Config::ReadCfg("01/Config/Nb30Ta70_5x5x5.cfg");
  cfg.UpdateNeighborList(cutoffs);

  std::cout << cfg.GetCentralAtomLatticeId() << std::endl;
  size_t vacancyId = cfg.GetCentralAtomLatticeId();
  size_t migratingAtomId = cfg.GetNeighborLatticeIdVectorOfLattice(vacancyId, 1)[0];

  std::pair<size_t, size_t> forwardJumpPair = {migratingAtomId, vacancyId};

  auto atomVector = cfg.GetAtomVector();
  std::set<Element> setElement(atomVector.begin(), atomVector.end());

  auto encodeVectorForward = GetEncodeVector3F(cfg,
          forwardJumpPair,
          setElement,
          2);







  // finding the staring and ending corners


}
*/
/*
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
int main()
{
  const vector<double> cutoffs = {3.3, 4.7, 5.6};
  const size_t initialVacancyId = 85; // Final Migrating Atom Id
  const size_t finalVacancyId = 74; // Initial Migrating Atom Id
  const size_t maxBondOrder = 2;
  const size_t maxClusterSize = 2;

  auto supercellCfg = Config::GenerateSupercell(4, 3.4, "X", "BCC");
  supercellCfg.UpdateNeighborList(cutoffs);

  string parentDir = "/media/sf_Phd/MEP/NbTa/";


  for (int i = 1; i<=200; i++)
  {
    string folderId = to_string(i);

    if(folderId.size() == 1)
    {
      folderId = "0" + folderId;
    }


    auto cfgMain = Config::ReadCfg(parentDir + folderId + "/NbTa_128.cfg");
    cfgMain.UpdateNeighborList(cutoffs);


    auto atomVector = cfgMain.GetAtomVector();
    set<Element> elementSet(atomVector.begin(), atomVector.end());
    Element vacancy("X");
    elementSet.emplace(vacancy);

    auto clusterTypeSet = InitializeClusterTypeSet(cfgMain, elementSet, maxClusterSize, maxBondOrder);

    PotentialEnergyEstimator peEstimator("predictor_file.json",
                                          cfgMain,
                                          supercellCfg,
                                          elementSet,
                                          maxClusterSize,
                                          maxBondOrder);

    vector<string> files = {"initial", "final"};

    for (auto fileType : files)
    {
      // From initial file
      // Migrating Atom is at its initial position (74)
      Config cfg;
      size_t atomId;
      size_t vacancyId;
      if (fileType == "initial")
      {

        cfg = Config::ReadCfg(parentDir + folderId + "/NbTa_128_initial.cfg");
        cfg.UpdateNeighborList(cutoffs);

        // cout << cfg.GetBasis() << endl;

        atomId = finalVacancyId;
        vacancyId = initialVacancyId;

      }
      else
      {
        cfg = Config::ReadCfg(parentDir + folderId + "/NbTa_128_final.cfg");
        cfg.UpdateNeighborList(cutoffs);

        atomId = initialVacancyId;
        vacancyId = finalVacancyId;

      }

      // cout << cfg.GetCartesianPositionOfLattice(atomId);
      Eigen::VectorXd encodeVector = peEstimator.GetEncodeVectorOfCluster(cfg, {atomId});


      Element vacancyIdElement = cfg.GetElementOfLattice(vacancyId);
      Element migratingElement = cfg.GetElementOfLattice(atomId);

      ofstream outputFile(parentDir + "outputCE.txt", std::ios::app);

      if  (outputFile.is_open())
      {
        if (outputFile.tellp() == 0)
        {
          outputFile << "FolderId\tFile\tVacancyIdElement\tMigratingElement\t";
          for (auto clusterType : clusterTypeSet)
          {
            outputFile << clusterType << "\t" ;
          }

          outputFile << "\n";
        }

        outputFile << folderId << "\t";
        outputFile << fileType << "\t";
        outputFile << vacancyIdElement << "\t";
        outputFile << migratingElement << "\t";

        for (auto encode : encodeVector)
        {
          outputFile << to_string(size_t(encode)) << "\t";
        }
        outputFile << "\n";

        outputFile.close();
      }

    }

    cout << "Done For " << folderId << endl;
  }

}
*/
//  Testing the CE for initial and final structure of Vacancy Migration
//   Manual
//   auto cfg = Config::ReadCfg("testing/NbTa_128.cfg");
//   cfg.UpdateNeighborList(cutoffs);
//
//   size_t initialVacancyId = 85; // Final Migrating Atom Id
//   size_t finalVacancyId = 74; // Initial Migrating Atom Id
//
//   auto migratingAtomElement = cfg.GetElementOfLattice(finalVacancyId);
//   Element vacancy("X");
//
//   auto supercellCfg = Config::GenerateSupercell(4, 3.4, "X", "BCC");
//   supercellCfg.UpdateNeighborList(cutoffs);
//
//
//   // Initial Structure
//
//   cfg.SetElementOfLattice(initialVacancyId, vacancy);
//
//   auto atomVector = cfg.GetAtomVector();
//   set<Element> elementSet(atomVector.begin(), atomVector.end());
//
//   PotentialEnergyEstimator peEstimator("predictor_file.json",
//                                        cfg,
//                                        supercellCfg,
//                                        elementSet,
//                                        2,
//                                        2);
//
//   cout << "Initial Structure" << endl;
//   cout << endl;
//   peEstimator.GetEncodeVectorOfCluster(cfg, {finalVacancyId});
//   cout << endl;
//
//   // Final Structure
//   cfg.SetElementOfLattice(initialVacancyId, migratingAtomElement);
//   cfg.SetElementOfLattice(finalVacancyId, vacancy);
//
//
//   cout << "Final Structure" << endl;
//   cout << endl;
//   peEstimator.GetEncodeVectorOfCluster(cfg, {initialVacancyId});
//   cout << endl;
//
//
//   // From initial file
//   // Migrating Atom is at its initial position (74)
//
//   auto cfgInitial = Config::ReadCfg("testing/NbTa_128_initial.cfg");
//   cfgInitial.UpdateNeighborList(cutoffs);
//
//   cout << endl;
//   cout << "74 :" << cfgInitial.GetElementOfLattice(74) << "\t"
//                  << cfgInitial.GetCartesianPositionOfLattice(74).transpose()
//                  << endl;
//   cout << "85 :" << cfgInitial.GetElementOfLattice(85) << "\t"
//                  << cfgInitial.GetCartesianPositionOfLattice(85).transpose()
//                  << endl;
//
//   cout << "Cluster Probability for the Initial Structure" << endl;
//   cout << endl;
//
//   peEstimator.GetEncodeVectorOfCluster(cfgInitial, {74});

// int main(int argc, char *argv[]) {
//
//   // api::Parameter parameter(argc, argv);
//   // api::Print(parameter);
//   // api::Run(parameter);
//
//   size_t supercellSize = 4;
//   double latticeParam = 3.4;
//   string structureType = "BCC";
//   vector<string> elements = {"W", "Ta", "Nb"};
//   vector<double> elementComposition = {40, 30, 30};
//   vector<double> cutoffs = cutoffs;
//
//   for (int i = 1; i <= 10; ++i)
//   {
//     ostringstream dir_name;
//     dir_name << setw(2) << setfill('0') << i;
//         // Create the directory
//         filesystem::path dir_path = dir_name.str();
//         if (filesystem::create_directory(dir_path)) {
//             cout << "Directory created: " << dir_path << endl;
//         } else {
//             cerr << "Failed to create directory: " << dir_path << endl;
//         }
//
//     auto cfg = Config::GenerateAlloySupercell(supercellSize,
//                                               latticeParam,
//                                               structureType,
//                                               elements,
//                                               elementComposition);
//     cfg.UpdateNeighborList(cutoffs);
//
//     string filename;
//     for (size_t i = 0; i < elements.size(); ++i){
//       filename += elements[i];
//     }
//     filename += "_" + to_string(int(pow(supercellSize, 3)*2));
//
//     cfg.WriteLAMMPSDataFile(dir_name.str() + "/" + filename, cfg);
//
//     auto vacancyMigration = cfg.GenerateNEBStructure(dir_name.str(), filename);
//
//     Config::WriteVacancyMigrationInfo(vacancyMigration,
//                                       cfg,
//                                       elements,
//                                       elementComposition,
//                                       "output.txt",
//                                       dir_name.str());
//
//
//
//   }
//
// }

// int main() {
//
//   size_t supercellSize = 4;
//   double latticeParam = 3.4;
//   string structureType = "BCC";
//   vector<string> elements = {"W", "Ta", "Nb"};
//   vector<double> elementComposition = {40, 30, 30};
//   vector<double> cutoffs = cutoffs;
//
//   unsigned seed = 5;
//
//
//   auto cfg = Config::GenerateAlloySupercell(supercellSize,
//                                             latticeParam,
//                                             structureType,
//                                             elements,
//                                             elementComposition,
//                                             seed);
//   cfg.UpdateNeighborList(cutoffs);
//
//   auto vacancyMigration = cfg.GenerateNEBStructure("0X", "testing.data");
//
//
//
//   Config::WriteVacancyMigrationInfo(vacancyMigration,
//                                       cfg,
//                                       elements,
//                                       elementComposition,
//                                       "output.txt",
//                                       "0" + to_string(seed));
//   cout << endl;
// }
