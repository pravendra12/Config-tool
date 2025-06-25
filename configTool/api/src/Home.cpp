/**************************************************************************************************
 * Copyright (c) 2023-2024. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 3/21/22 3:17 PM                                                                         *
 * @Last Modified by: pravendra12                                                                    *
 * @Last Modified time: 11/28/24 8:30 PM                                                            *
 **************************************************************************************************/

/*! \file  Home.h
 *  \brief File for the Home class implementation.
 */

#include "Home.h"
#include "GenerateNEBStructure.h"

using namespace std;

namespace api
{

  void Print(const Parameter &parameter)
  {
    std::cout << "Parameters" << std::endl;
    std::cout << "simulation_method: " << parameter.method << std::endl;

    if (parameter.method == "GenerateNEBStructures")
    {
      cout << "supercell_size: " << parameter.supercell_size_ << endl;
      std::cout << "lattice_param: " << parameter.lattice_param_ << std::endl;
      std::cout << "structure_type: " << parameter.structure_type_ << std::endl;
      std::cout << "element_set: ";
      std::transform(parameter.element_vector_.begin(), parameter.element_vector_.end(),
                     std::ostream_iterator<std::string>(std::cout, " "),
                     [](auto elements)
                     { return elements; });
      std::cout << std::endl;
      std::cout << "element_composition: ";
      std::transform(parameter.element_composition_.begin(), parameter.element_composition_.end(),
                     std::ostream_iterator<std::string>(std::cout, " "),
                     [](auto element_comp)
                     { return std::to_string(element_comp); });
      cout << endl;
      std::cout << "cutoffs: ";
      std::transform(parameter.cutoffs_.begin(), parameter.cutoffs_.end(),
                     std::ostream_iterator<std::string>(std::cout, " "),
                     [](auto cutoff)
                     { return std::to_string(cutoff); });
      std::cout << std::endl;
      std::cout << "random_seed: " << parameter.random_seed_ << endl;
      cout << "filename: " << parameter.config_filename_ << endl;
    }
    else if (parameter.method == "GenerateAlloySupercell")
    {
      cout << "supercell_size: " << parameter.supercell_size_ << endl;
      std::cout << "lattice_param: " << parameter.lattice_param_ << std::endl;
      std::cout << "structure_type: " << parameter.structure_type_ << std::endl;
      std::cout << "element_set: ";
      std::transform(parameter.element_vector_.begin(), parameter.element_vector_.end(),
                     std::ostream_iterator<std::string>(std::cout, " "),
                     [](auto elements)
                     { return elements; });
      std::cout << std::endl;
      std::cout << "element_composition: ";
      std::transform(parameter.element_composition_.begin(), parameter.element_composition_.end(),
                     std::ostream_iterator<std::string>(std::cout, " "),
                     [](auto element_comp)
                     { return std::to_string(element_comp); });
      cout << endl;
      std::cout << "random_seed: " << parameter.random_seed_ << endl;

      cout << "filename: " << parameter.config_filename_ << endl;
    }

    else if (parameter.method == "GenerateStructuresCNT")
    {
      cout << "structure_type: " << parameter.structure_type_ << endl;
      cout << "supercell_size: " << parameter.supercell_size_ << endl;
      cout << "lattice_param: " << parameter.lattice_param_ << endl;
      std::cout << "element_set: ";
      std::transform(parameter.element_vector_.begin(), parameter.element_vector_.end(),
                     std::ostream_iterator<std::string>(std::cout, " "),
                     [](auto elements)
                     { return elements; });
      std::cout << std::endl;
      std::cout << "element_composition: ";
      std::transform(parameter.element_composition_.begin(), parameter.element_composition_.end(),
                     std::ostream_iterator<std::string>(std::cout, " "),
                     [](auto element_comp)
                     { return std::to_string(element_comp); });
      cout << endl;
      cout << "B2_element_pair: " << parameter.b2_element_pair_.first << " "
           << parameter.b2_element_pair_.second << endl;

      cout << "num_random_structures: " << parameter.num_random_structures_ << endl;
      cout << "max_num_B2_center: " << parameter.max_num_B2_center_ << endl;
      cout << "num_sample_per_B2_config: " << parameter.num_sample_per_B2_config_ << endl;

      cout << "compute_energy: " << parameter.compute_energy_ << endl;
      cout << "output_directory: " << parameter.output_directory_ << endl;

      cout << "json_coefficient_filename: " << parameter.json_coefficients_filename_ << endl;
    }
  }

  void Run(const Parameter &parameter)
  {
    if (parameter.method == "GenerateNEBStructures")
    {
      //  BuildGenerateNEBStructuresFromParameter(parameter);
      auto config = Config::GenerateAlloySupercell(parameter.supercell_size_,
                                                   parameter.lattice_param_,
                                                   parameter.structure_type_,
                                                   parameter.element_vector_,
                                                   parameter.element_composition_,
                                                   parameter.random_seed_);

      config.UpdateNeighborList(parameter.cutoffs_);

      pair<size_t, size_t> latticeIdJumpPair = {config.GetCentralAtomLatticeId(),
                                                config.GetNeighborLatticeIdVectorOfLattice(config.GetCentralAtomLatticeId(), 1)[0]};
      
      // 1 to 4 are Nb, Mo, Ta and W.
      map<Element,size_t> elementMap;
      size_t idx = 1;
      for (const auto element : parameter.element_vector_)
      {
        elementMap[Element(element)] = idx;
        idx++;
      }
      
      GenerateNEBStructure(parameter.config_filename_,
                           config,
                           latticeIdJumpPair, 
                           elementMap);
    }
    else if (parameter.method == "GenerateAlloySupercell")
    {
      auto config = Config::GenerateAlloySupercell(parameter.supercell_size_,
                                                   parameter.lattice_param_,
                                                   parameter.structure_type_,
                                                   parameter.element_vector_,
                                                   parameter.element_composition_,
                                                   parameter.random_seed_);
      Config::WriteConfig(parameter.config_filename_ + ".cfg", config);

      // 1 to 4 are Nb, Mo, Ta and W.
      map<Element,size_t> elementMap;
      size_t idx = 1;
      for (const auto element : parameter.element_vector_)
      {
        elementMap[Element(element)] = idx;
        idx++;
      }

      Config::WriteLAMMPSDataFileCustom(parameter.config_filename_ + ".data", config, elementMap);
    }

    else if (parameter.method == "GenerateStructuresCNT")
    {

      auto trainingSupercell = Config::GeneratePristineSupercell(5,
                                                                 parameter.lattice_param_, "X",
                                                                 parameter.structure_type_);
      trainingSupercell.UpdateNeighborList(parameter.cutoffs_);

      std::pair<Element, Element> b2ElementPair;
      b2ElementPair.first = Element(parameter.b2_element_pair_.first);
      b2ElementPair.second = Element(parameter.b2_element_pair_.second);

      Config config;
      int idx = 0;
      do
      {
        config = Config::GenerateAlloySupercell(parameter.supercell_size_,
                                                parameter.lattice_param_,
                                                parameter.structure_type_,
                                                parameter.element_vector_,
                                                parameter.element_composition_,
                                                idx);
        idx++;
      } while (b2ElementPair.first != config.GetElementOfLattice(config.GetCentralAtomLatticeId()));

      cout << "Central Element of B2: " << config.GetElementOfLattice(config.GetCentralAtomLatticeId()) << endl;

      config.UpdateNeighborList(parameter.cutoffs_);

      std::set<Element> elementSetCE;
      for (const auto &element : parameter.element_vector_)
      {
        elementSetCE.insert(Element(element));
      }
      elementSetCE.insert(Element("X"));

      GenerateStructureCNT cnt(
          parameter.json_coefficients_filename_,
          config,
          trainingSupercell,
          parameter.supercell_size_,
          elementSetCE,
          parameter.element_composition_);

      cout << endl;
      if (parameter.compute_energy_)
      {
        cout << "Energy For Random Configuration" << endl;
        cout << "Id\tEnergy" << endl;
      }
      cnt.GenerateRandomStructures(parameter.num_random_structures_,
                                   parameter.compute_energy_,
                                   true,
                                   parameter.output_directory_);

      cout << endl;
      if (parameter.compute_energy_)
      {
        cout << "Energy For B2 Embedded Configuration" << endl;
        cout << "numB2\trandomId\tEnergy" << endl;
      }
      cnt.GenerateStructureWithB2(
          config,
          parameter.max_num_B2_center_,
          parameter.num_sample_per_B2_config_,
          b2ElementPair,
          parameter.compute_energy_,
          parameter.output_directory_);
    }

    else if (parameter.method == "GenerateNEBStructureWithB2")
    {
      //  BuildGenerateNEBStructuresFromParameter(parameter);
      auto config = Config::GenerateAlloySupercell(parameter.supercell_size_,
                                                   parameter.lattice_param_,
                                                   parameter.structure_type_,
                                                   parameter.element_vector_,
                                                   parameter.element_composition_,
                                                   parameter.random_seed_);

      config.UpdateNeighborList(parameter.cutoffs_);

      pair<Element, Element> b2ElementPair = {Element(parameter.b2_element_pair_.first), 
                                       Element(parameter.b2_element_pair_.second)};

      B2Ordering::AddB2Precipitate(config, 4, b2ElementPair);

      pair<size_t, size_t> latticeIdJumpPair = {config.GetCentralAtomLatticeId(),
                                                config.GetNeighborLatticeIdVectorOfLattice(config.GetCentralAtomLatticeId(), 1)[0]};
      
      // 1 to 4 are Nb, Mo, Ta and W.
      map<Element,size_t> elementMap;
      size_t idx = 1;
      for (const auto element : parameter.element_vector_)
      {
        elementMap[Element(element)] = idx;
        idx++;
      }
      
      GenerateNEBStructure(parameter.config_filename_,
                           config,
                           latticeIdJumpPair, 
                           elementMap);
    }

    else
    {
      std::cout << "No such method: " << parameter.method << std::endl;
    }
  }

  void BuildGenerateNEBStructuresFromParameter(const Parameter &parameter)
  {
    for (size_t i = 1; i <= parameter.num_unique_structure_; ++i)
    {
      std::ostringstream dir_name;
      dir_name << std::setw(2) << std::setfill('0') << i;
      // Create the directory

      std::filesystem::path dir_path = dir_name.str();
      if (std::filesystem::create_directory(dir_path))
      {
        std::cout << "Directory created: " << dir_path << std::endl;
      }
      else
      {
        std::cerr << "Failed to create directory: " << dir_path << std::endl;
      }

      auto cfg = Config::GenerateAlloySupercell(parameter.supercell_size_,
                                                parameter.lattice_param_,
                                                parameter.structure_type_,
                                                parameter.element_vector_,
                                                parameter.element_composition_,
                                                i);

      cfg.UpdateNeighborList(parameter.cutoffs_);

      std::string filename;

      // Ta50W50_4x4x4.data

      for (size_t j = 0; j < parameter.element_vector_.size(); ++j)
      {
        filename += parameter.element_vector_[j];
        filename += std::to_string(int(parameter.element_composition_[j]));
      }

      std::ostringstream oss;
      oss << parameter.supercell_size_ << "x"
          << parameter.supercell_size_ << "x"
          << parameter.supercell_size_;

      // Append the formatted string to `filename`
      filename += "_" + oss.str();

      std::filesystem::path dirConfig = dir_name.str() + "/Config";
      if (std::filesystem::create_directory(dirConfig))
      {
        std::cout << "Directory created: " << dirConfig << std::endl;
      }
      else
      {
        std::cerr << "Failed to create directory: " << dirConfig << std::endl;
      }

      // Lammps Data File
      Config::WriteLAMMPSDataFile(dirConfig.string() + "/" + filename + ".data", cfg);

      // Cfg file
      Config::WriteConfig(dirConfig.string() + "/" + filename + ".cfg", cfg);

      // Central Atom Id
      // Initial Vacant Site
      auto vacancyId = cfg.GetCentralAtomLatticeId();
      auto migratingAtomId =
          cfg.GetNeighborLatticeIdVectorOfLattice(vacancyId, 1)[0];
    }
  }

}