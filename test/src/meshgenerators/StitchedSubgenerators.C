//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "StitchedSubgenerators.h"

#include "CastUniquePointer.h"
#include "MooseUtils.h"
#include "FileMeshGenerator.h"

#include "libmesh/replicated_mesh.h"

registerMooseObject("MooseApp", StitchedSubgenerators);

defineLegacyParams(StitchedSubgenerators);

InputParameters
StitchedSubgenerators::validParams()
{
  InputParameters params = MeshGenerator::validParams();

  MooseEnum algorithm("BINARY EXHAUSTIVE", "BINARY");

  params.addRequiredParam<std::vector<std::string>>("inputs", "The input mesh filenames");
  params.addParam<bool>(
      "clear_stitched_boundary_ids", true, "Whether or not to clear the stitched boundary IDs");
  params.addRequiredParam<std::vector<std::vector<std::string>>>(
      "stitch_boundaries_pairs",
      "Pairs of boundaries to be stitched together between the 1st mesh in inputs and each "
      "consecutive mesh");
  params.addParam<MooseEnum>(
      "algorithm",
      algorithm,
      "Control the use of binary search for the nodes of the stitched surfaces.");
  params.addClassDescription(
      "Allows multiple mesh files to be stiched together to form a single mesh.");

  return params;
}

StitchedSubgenerators::StitchedSubgenerators(const InputParameters & parameters)
  : MeshGenerator(parameters),
    _input_names(getParam<std::vector<std::string>>("inputs")),
    _clear_stitched_boundary_ids(getParam<bool>("clear_stitched_boundary_ids")),
    _stitch_boundaries_pairs(
        getParam<std::vector<std::vector<std::string>>>("stitch_boundaries_pairs")),
    _algorithm(parameters.get<MooseEnum>("algorithm"))
{
  MooseApp &app = this->getMooseApp();
  const std::string sg_name_base = "subgenerator_";

  // InputParameters filemesh_params = FileMeshGenerator::validParams();
  InputParameters filemesh_params = _app.getFactory().getValidParams("FileMeshGenerator");

  // Create and add MeshGenerators for the input meshes
  _mesh_ptrs.reserve(_input_names.size());
  int sg_num = 0;
  for (auto & input_name : _input_names)
    {
      filemesh_params.set<MeshFileName>("file") = input_name;

      const std::string sg_name = sg_name_base + std::to_string(sg_num++);

      app.addMeshGenerator("FileMeshGenerator", sg_name, filemesh_params);

      _mesh_ptrs.push_back(&getMeshByName(sg_name));
    }
}

std::unique_ptr<MeshBase>
StitchedSubgenerators::generate()
{
  // We put the first mesh in a local pointer
  std::unique_ptr<ReplicatedMesh> mesh = dynamic_pointer_cast<ReplicatedMesh>(*_mesh_ptrs[0]);

  // Reserve spaces for the other meshes (no need to store the first one another time)
  _meshes.reserve(_input_names.size() - 1);

  // Read in all of the other meshes
  for (MooseIndex(_input_names) i = 1; i < _input_names.size(); ++i)
    _meshes.push_back(dynamic_pointer_cast<ReplicatedMesh>(*_mesh_ptrs[i]));

  // Stitch all the meshes to the first one
  for (MooseIndex(_meshes) i = 0; i < _meshes.size(); i++)
  {
    auto boundary_pair = _stitch_boundaries_pairs[i];

    boundary_id_type first, second;

    try
    {
      first = MooseUtils::convert<boundary_id_type>(boundary_pair[0], true);
    }
    catch (...)
    {
      first = mesh->get_boundary_info().get_id_by_name(boundary_pair[0]);

      if (first == BoundaryInfo::invalid_id)
      {
        std::stringstream error;

        error << "Boundary " << boundary_pair[0] << " doesn't exist on the first mesh in " << name()
              << "\n";
        error << "Boundary names that do exist: \n";
        error << " ID : Name\n";

        auto & sideset_id_name_map = mesh->get_boundary_info().get_sideset_name_map();

        for (auto & ss_name_map_pair : sideset_id_name_map)
          error << " " << ss_name_map_pair.first << " : " << ss_name_map_pair.second << "\n";

        paramError("stitch_boundaries_pairs", error.str());
      }
    }

    try
    {
      second = MooseUtils::convert<boundary_id_type>(boundary_pair[1], true);
    }
    catch (...)
    {
      second = _meshes[i]->get_boundary_info().get_id_by_name(boundary_pair[1]);

      if (second == BoundaryInfo::invalid_id)
      {
        _meshes[i]->print_info();

        std::stringstream error;

        error << "Boundary " << boundary_pair[1] << " doesn't exist on mesh " << i + 1 << " in "
              << name() << "\n";
        error << "Boundary names that do exist: \n";
        error << " ID : Name\n";

        auto & sideset_id_name_map = _meshes[i]->get_boundary_info().get_sideset_name_map();

        for (auto & ss_name_map_pair : sideset_id_name_map)
          error << " " << ss_name_map_pair.first << " : " << ss_name_map_pair.second << "\n";

        paramError("stitch_boundaries_pairs", error.str());
      }
    }

    const bool use_binary_search = (_algorithm == "BINARY");

    mesh->stitch_meshes(*_meshes[i],
                        first,
                        second,
                        TOLERANCE,
                        _clear_stitched_boundary_ids,
                        /*verbose = */ true,
                        use_binary_search);
  }

  return dynamic_pointer_cast<MeshBase>(mesh);
}
