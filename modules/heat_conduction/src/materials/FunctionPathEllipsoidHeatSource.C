//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FunctionPathEllipsoidHeatSource.h"

registerMooseObject("HeatConductionApp", FunctionPathEllipsoidHeatSource);

template <>
InputParameters
validParams<FunctionPathEllipsoidHeatSource>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredParam<Real>("power", "laser power");
  params.addParam<Real>("efficienty", 1, "process efficienty");
  params.addRequiredParam<Real>("a", "transverse ellipsoid axe");
  params.addRequiredParam<Real>("b",  "depth ellipsoid axe");
  params.addRequiredParam<Real>("c", "longitudinal ellipsoid axe");
  params.addRequiredParam<Real>("velocity",  "heating spot travel speed");
  params.addRequiredParam<Real>("factor",  "scaling factor");
  params.addParam<FunctionName>("function_x", "The x component heating spot travel path");
  params.addParam<FunctionName>("function_y", "The y component heating spot travel path");
  params.addParam<FunctionName>("function_z", "The z component heating spot travel path");
  params.addClassDescription("Double ellipsoid volumetric source heat with function path.");

  return params;
}

FunctionPathEllipsoidHeatSource::FunctionPathEllipsoidHeatSource(const InputParameters & parameters)
  : Material(parameters),
    _P(getParam<Real>("power")),
    _eta(getParam<Real>("efficienty")),
    _a(getParam<Real>("a")),
    _b(getParam<Real>("b")),
    _c(getParam<Real>("c")),
    _v(getParam<Real>("velocity")),
    _f(getParam<Real>("factor")),
    _function_x(getFunction("function_x")),
    _function_y(getFunction("function_y")),
    _function_z(getFunction("function_z")),
    _volumetric_heat(declareADProperty<Real>("volumetric_heat"))
{
}

void
FunctionPathEllipsoidHeatSource::computeQpProperties()
{
  const Real & x = _q_point[_qp](0);
  const Real & y = _q_point[_qp](1);
  const Real & z = _q_point[_qp](2);

  // center of the heat source
  Real x_t = _function_x.value(_t, _q_point[_qp]);
  Real y_t = _function_y.value(_t, _q_point[_qp]);
  Real z_t = _function_z.value(_t, _q_point[_qp]);

  _volumetric_heat[_qp] = 6.0 * std::sqrt(3.0) * _P * _eta * _f /
                          (_a * _b * _c * std::pow(libMesh::pi, 1.5)) *
                          std::exp(-(3.0 * std::pow(x - x_t, 2.0) / std::pow(_a, 2.0) +
                                     3.0 * std::pow(y - y_t, 2.0) / std::pow(_b, 2.0) +
                                     3.0 * std::pow(z - z_t, 2.0) / std::pow(_c, 2.0)));
}
