//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ActivateElementTemp.h"
#include "libmesh/quadrature.h"
#include "libmesh/parallel_algebra.h"
#include "libmesh/parallel.h"
#include "libmesh/point.h"

registerMooseObject("MooseApp", ActivateElementTemp);

template <>
InputParameters
validParams<ActivateElementTemp>()
{
  InputParameters params = validParams<ElementUserObject>();
  params.addClassDescription("Determine activated elements.");
  params.addRequiredParam<int>("active_subdomain_id", "The active subdomain ID.");
  params.addParam<Real>("activate_tol", 1e-4, "The spatial tolerance for activating an element.");
  params.addParam<FunctionName>("function_x", "The x component heating spot travel path");
  params.addParam<FunctionName>("function_y", "The y component heating spot travel path");
  params.addParam<FunctionName>("function_z", "The z component heating spot travel path");

  return params;
}

ActivateElementTemp::ActivateElementTemp(const InputParameters & parameters)
  : ElementUserObject(parameters),
    _active_subdomain_id(getParam<int>("active_subdomain_id")),
    _tol(getParam<Real>("activate_tol")),
    _function_x(getFunction("function_x")),
    _function_y(getFunction("function_y")),
    _function_z(getFunction("function_z"))
{
}

void
ActivateElementTemp::execute()
{
  /*
    Check if current element is activated
  */

  // activate center (assume position of the activate center is only time dependent)
  Real x_t = _function_x.value(_t, _q_point[0]);
  Real y_t = _function_y.value(_t, _q_point[0]);
  Real z_t = _function_z.value(_t, _q_point[0]);

  if(_current_elem->contains_point( Point (x_t, y_t, z_t) ) && _current_elem->subdomain_id()!=_active_subdomain_id)
  {
    /*
      _current_elem subdomain id is not assignable
      create a copy of this element from MooseMesh
    */
    dof_id_type ele_id= _current_elem->id();
    Elem * ele = _mesh.elemPtr(ele_id);
    std::cout<<"====>  Current element info:\n";
    ele->print_info();
    /*
      Add element to the activate subdomain
    */
    ele->subdomain_id()=_active_subdomain_id;
    /*
      Reinit equation systems
    */
    _mesh.getMesh().prepare_for_use();
    _mesh.meshChanged();
    _fe_problem.es().reinit_solutions();
    _fe_problem.es().reinit();
    _fe_problem.es().reinit_systems();
    _mesh.getMesh().prepare_for_use();

    // std::cout<<"====>  Neighbor element info:\n";
    // for (auto s : ele->side_index_range())
    // {
    //   if (ele->neighbor_ptr(s))
    //   {
    //     dof_id_type neighbor_ele_id=ele->neighbor_ptr(s)->id();
    //     Elem * neighbor_ele = _mesh.elemPtr(neighbor_ele_id);
    //     neighbor_ele->print_info();
    //   }
    // }
    // std::cout<<"====>  Current element info:\n";
    // ele->print_info();

  }


}

void
ActivateElementTemp::finalize()
{
  // _communicator.set_union(_activated_elem_map);
}
