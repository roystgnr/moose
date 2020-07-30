//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"

// Forward Declarations
class FunctionPathEllipsoidHeatSource;

template <>
InputParameters validParams<FunctionPathEllipsoidHeatSource>();

/**
 * Double ellipsoid heat source distribution.
 */
class FunctionPathEllipsoidHeatSource : public Material
{
public:
  FunctionPathEllipsoidHeatSource(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// power
  const Real _P;
  /// process efficienty
  const Real _eta;
  /// transverse ellipsoid axe
  const Real _a;
  /// depth ellipsoid axe
  const Real _b;
  /// longitudinal ellipsoid axe
  const Real _c;
  /// heating spot travel speed
  const Real _v;
  /// scaling factor
  const Real _f;
  /// path of the heat source, x, y, z components
  const Function & _function_x;
  const Function & _function_y;
  const Function & _function_z;

  ADMaterialProperty<Real> & _volumetric_heat;
};
