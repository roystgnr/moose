//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BreakMeshByElementGenerator.h"
#include "CastUniquePointer.h"
#include "MooseMeshUtils.h"

#include "libmesh/partitioner.h"

registerMooseObject("MooseApp", BreakMeshByElementGenerator);
registerMooseObjectRenamed("MooseApp",
                           ExplodeMeshGenerator,
                           "05/18/2024 24:00",
                           BreakMeshByElementGenerator);

InputParameters
BreakMeshByElementGenerator::validParams()
{
  InputParameters params = MeshGenerator::validParams();
  params.addClassDescription("Break all element-element interfaces in the specified subdomains.");
  params.addRequiredParam<MeshGeneratorName>("input", "The mesh we want to modify");
  params.addParam<std::vector<SubdomainID>>(
      "subdomains",
      std::vector<SubdomainID>(),
      "The list of subdomain IDs to explode.  Leave unset to explode all subdomains.");
  params.addParam<BoundaryName>(
      "interface_name",
      "element_boundaries",
      "The boundary name containing all broken element-element interfaces.");
  params.addRangeCheckedParam<unsigned int>(
      "interface_sides",
      1,
      "interface_sides<3",
      "Whether to add no interface boundary, a 1-sided boundary (facing from lower to higher "
      "element id), or a 2-sided boundary");
  return params;
}

BreakMeshByElementGenerator::BreakMeshByElementGenerator(const InputParameters & parameters)
  : MeshGenerator(parameters),
    _input(getMesh("input")),
    _subdomains(getParam<std::vector<SubdomainID>>("subdomains")),
    _interface_name(getParam<BoundaryName>("interface_name")),
    _interface_sides(getParam<unsigned int>("interface_sides"))
{
}

std::unique_ptr<MeshBase>
BreakMeshByElementGenerator::generate()
{
  std::unique_ptr<MeshBase> mesh = std::move(_input);

  if (!mesh->is_prepared())
    mesh->prepare_for_use();

  // check that the subdomain IDs exist in the mesh
  for (const auto & id : _subdomains)
    if (!MooseMeshUtils::hasSubdomainID(*mesh, id))
      paramError("subdomains", "The block ID '", id, "' was not found in the mesh");

  BoundaryInfo & boundary_info = mesh->get_boundary_info();
  if (_interface_sides &&
      boundary_info.get_id_by_name(_interface_name) != Moose::INVALID_BOUNDARY_ID)
    paramError("interface_name", "The specified interface name already exists in the mesh.");

  const auto node_to_elem_map = buildSubdomainRestrictedNodeToElemMap(mesh, _subdomains);

  duplicateNodes(mesh, node_to_elem_map);

  createInterface(*mesh, node_to_elem_map);

  Partitioner::set_node_processor_ids(*mesh);

  // We need to update the global_boundary_ids, and this is faster
  // than a full prepare_for_use()
  boundary_info.regenerate_id_sets();

  return dynamic_pointer_cast<MeshBase>(mesh);
}

BreakMeshByElementGenerator::NodeToElemMapType
BreakMeshByElementGenerator::buildSubdomainRestrictedNodeToElemMap(
    std::unique_ptr<MeshBase> & mesh, const std::vector<SubdomainID> & subdomains) const
{
  NodeToElemMapType node_to_elem_map;
  for (const auto & elem : mesh->active_element_ptr_range())
  {
    // Skip if subdomains are specified and the element is not in them
    if (!subdomains.empty() &&
        std::find(subdomains.begin(), subdomains.end(), elem->subdomain_id()) == subdomains.end())
      continue;

    std::set<const Elem *> neighbors;
    elem->find_point_neighbors(neighbors);

    for (auto n : make_range(elem->n_nodes()))
    {
      // if ANY neighboring element that contains this node is not in specified subdomains,
      // don't add this node to the map, i.e. don't split this node.
      bool should_duplicate = true;
      if (!subdomains.empty())
        for (auto neighbor : neighbors)
          if (neighbor->contains_point(elem->node_ref(n)) &&
              std::find(subdomains.begin(), subdomains.end(), neighbor->subdomain_id()) ==
                  subdomains.end())
          {
            should_duplicate = false;
            break;
          }

      if (should_duplicate)
        node_to_elem_map[elem->node_id(n)].insert(elem->id());
    }
  }

  // If we aren't serial, we may not see every element connected to
  // every node we see.  We don't need to see every such element, but
  // we do need to know they exist so we get our numbering right.
  // We'll push our numbers to them first, then query the unioned
  // numbers second.
  if (!mesh->is_serial())
  {
    std::map<processor_id_type, std::vector<std::pair<dof_id_type, std::vector<dof_id_type>>>>
        submaps_to_push;
    std::map<processor_id_type, std::vector<dof_id_type>> nodes_to_query;
    const processor_id_type my_pid = mesh->processor_id();
    for (const auto & [node_id, connected_elem_ids] : node_to_elem_map)
    {
      const processor_id_type node_pid = mesh->node_ref(node_id).processor_id();
      if (node_pid == my_pid)
        continue;
      submaps_to_push[node_pid].push_back(std::make_pair(
          node_id, std::vector<dof_id_type>(connected_elem_ids.begin(), connected_elem_ids.end())));
      nodes_to_query[node_pid].push_back(node_id);
    }

    auto collect_functor =
        [&node_to_elem_map](
            processor_id_type,
            const std::vector<std::pair<dof_id_type, std::vector<dof_id_type>>> & incoming_submap)
    {
      for (const auto & [node_id, connected_elem_ids] : incoming_submap)
        node_to_elem_map[node_id].insert(connected_elem_ids.begin(), connected_elem_ids.end());
    };

    Parallel::push_parallel_vector_data(mesh->comm(), submaps_to_push, collect_functor);

    auto gather_functor = [&node_to_elem_map](processor_id_type,
                                              const std::vector<dof_id_type> & nodes,
                                              std::vector<std::vector<dof_id_type>> & data)
    {
      const std::size_t query_size = nodes.size();

      data.resize(query_size);
      for (auto i : make_range(query_size))
      {
        auto & elems = libmesh_map_find(node_to_elem_map, nodes[i]);
        data[i].insert(data[i].end(), elems.begin(), elems.end());
      }
    };

    auto action_functor = [&node_to_elem_map](processor_id_type,
                                              const std::vector<dof_id_type> & nodes,
                                              const std::vector<std::vector<dof_id_type>> & data)
    {
      for (auto i : make_range(nodes.size()))
        node_to_elem_map[nodes[i]].insert(data[i].begin(), data[i].end());
    };

    std::vector<dof_id_type> * data_ex = nullptr;
    Parallel::pull_parallel_vector_data(
        mesh->comm(), nodes_to_query, gather_functor, action_functor, data_ex);
  }

  return node_to_elem_map;
}

void
BreakMeshByElementGenerator::duplicateNodes(std::unique_ptr<MeshBase> & mesh,
                                            const NodeToElemMapType & node_to_elem_map) const
{
  // If we're not doing this replicated, we need to manually set new
  // node ids and unique_ids.  Let's figure out what our offsets
  // should be.
  if (!mesh->preparation().has_synched_id_counts)
    mesh->update_parallel_id_counts();
  const dof_id_type max_node_id = mesh->max_node_id();
  const dof_id_type max_unique_id = mesh->parallel_max_unique_id();

  std::size_t max_elems_per_node = 0;

  for (const auto & [node_id, connected_elem_ids] : node_to_elem_map)
  {
    max_elems_per_node = std::max(max_elems_per_node, connected_elem_ids.size());
    unsigned int copy_num = 0;
    for (auto & connected_elem_id : connected_elem_ids)
    {
      Elem * elem = mesh->query_elem_ptr(connected_elem_id);
      if (connected_elem_id != *connected_elem_ids.begin() && elem)
        duplicateNode(mesh, elem, mesh->node_ptr(node_id), copy_num, max_node_id, max_unique_id);
      ++copy_num;
    }
  }

  mesh->set_next_unique_id(max_unique_id * max_elems_per_node);

  // We'll want to renumber and we'll need to sync id counts later.
  mesh->unset_has_synched_id_counts();
}

void
BreakMeshByElementGenerator::duplicateNode(std::unique_ptr<MeshBase> & mesh,
                                           Elem * elem,
                                           const Node * node,
                                           unsigned int copy_num,
                                           dof_id_type max_node_id,
                                           dof_id_type max_unique_id) const
{
  std::unique_ptr<Node> new_node = Node::build(*node, Node::invalid_id);
  new_node->processor_id() = elem->processor_id();
  new_node->set_id(max_node_id * copy_num + node->id());
  new_node->set_unique_id(max_unique_id * copy_num + node->unique_id());
  Node * added_node = mesh->add_node(std::move(new_node));
  elem->set_node(elem->get_node_index(node), added_node);

  // Add boundary info to the new node
  BoundaryInfo & boundary_info = mesh->get_boundary_info();
  std::vector<boundary_id_type> node_boundary_ids;
  boundary_info.boundary_ids(node, node_boundary_ids);
  boundary_info.add_node(added_node, node_boundary_ids);
}

void
BreakMeshByElementGenerator::createInterface(MeshBase & mesh,
                                             const NodeToElemMapType & node_to_elem_map) const
{
  std::set<std::pair<dof_id_type, unsigned int>> sides_breaking;

  for (const auto & node_to_elems : node_to_elem_map)
    for (const auto & elem_id_i : node_to_elems.second)
    {
      Elem * elem_i = mesh.elem_ptr(elem_id_i);
      for (const auto & elem_id_j : node_to_elems.second)
      {
        Elem * elem_j = mesh.elem_ptr(elem_id_j);
        if (elem_i != elem_j && elem_i->has_neighbor(elem_j))
          sides_breaking.insert(std::make_pair(elem_id_i, elem_i->which_neighbor_am_i(elem_j)));
      }
    }

  if (_interface_sides)
  {
    BoundaryInfo & boundary_info = mesh.get_boundary_info();

    // libMesh should only need the boundary_info prepared to call
    // get_global_boundary_ids(), but still asserts too much there, so
    // just prepare anything that isn't.
    mesh.complete_preparation();

    const auto & existing_boundary_ids = boundary_info.get_global_boundary_ids();
    const boundary_id_type interface_id =
        existing_boundary_ids.empty() ? 0 : *existing_boundary_ids.rbegin() + 1;
    boundary_info.sideset_name(interface_id) = _interface_name;

    for (const auto & [elem_id, side] : sides_breaking)
      if (_interface_sides > 1 || elem_id > mesh.elem_ptr(elem_id)->neighbor_ptr(side)->id())
        boundary_info.add_side(elem_id, side, interface_id);
  }

  // Remove element neighbor connections on the broken sides
  for (const auto & [elem_id, side] : sides_breaking)
    mesh.elem_ref(elem_id).set_neighbor(side, nullptr);
}
