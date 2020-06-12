//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADRadialReturnStressUpdate.h"

/**
 * This class uses the stress update material in a radial return isotropic creep
 * model.  This class is one of the basic radial return constitutive models; more complex
 * constitutive models combine creep and plasticity.
 *
 * This class inherits from RadialReturnCreepStressUpdateBase and must be used
 * in conjunction with ComputeMultipleInelasticStress.  This class calculates
 * creep based on stress, temperature, and time effects.  This class also
 * computes the creep strain as a stateful material property.
 */
class ADRateTempDependentStressUpdate : public ADRadialReturnStressUpdate
{
public:
  static InputParameters validParams();

  ADRateTempDependentStressUpdate(const InputParameters & parameters);

protected:
  virtual void computeStressInitialize(const ADReal & effective_trial_stress,
                                       const ADRankFourTensor & elasticity_tensor) override;
  virtual ADReal computeResidual(const ADReal & effective_trial_stress,
                                 const ADReal & scalar) override;
  virtual ADReal computeDerivative(const ADReal & effective_trial_stress,
                                   const ADReal & scalar) override;

  void computeFlowRule(const ADReal & effective_trial_stress,
                         const ADReal & scalar = 0.0);

  virtual void updateInternalStateVariables(const ADReal & effective_trial_stress,
                                            const ADReal & scalar=0.0,
                                            const ADReal & scalar_increment=0.0) override;

  virtual void initQpStatefulProperties() override;
  virtual void propagateQpStatefulProperties() override;
  virtual void computeStressFinalize(const ADRankTwoTensor & plastic_strain_increment) override;

  // Real computeReferenceResidual(const ADReal & effective_trial_stress, const ADReal & scalar_effective_inelastic_strain) override;

  /// Temperature variable value
  const ADVariableValue * _temperature;

  /// Simulation start time
  const Real _start_time;

  ///
  /// Rate and temperature dependent plasticity model parameters
  ///

  /// Rate independent yield constant [Pa]
  const Real _Y0;

  /// Rate independent yield temperature dependencies [K], [1/K], [K], [-]
  const Real _Y1, _Y2, _Y3, _Y4;

  /// Isotropic hardening shear coefficient [-]
  const Real _Hmu;

  /// Flow rule coefficient constants [1/s], [K]
  const Real _f1, _f2;

  /// Flow rule exponent constant [-]
  const Real _n1;

  /// Flow rule exponent temperature dependence [K]
  const Real _n2;

  /// Isotropic dynamic recovery constant [Pa]
  const Real _Rd1;

  /// Isotropic dynamic recovery temperature dependence [K]
  const Real _Rd2;

  /// Misorientation variable hardening constant [m/(s Pa)]
  const Real _hxi;


  /// Components for computing the derivatives
  ADReal _C1, _C2;

  /// Flow rule function
  ADReal _phi;

  /// Isotropic harderning internal state variable
  ADReal _r, _r_old, _dr;

  /// Temperature dependent yield stress
  ADReal _Y;

  /// Temperature dependent shear modulus
  ADReal _G;

  /// Derivative of temperature dependent shear modulus
  ADReal _dG;


  /// Plastic strain material property
  ADMaterialProperty<RankTwoTensor> & _plastic_strain;
  const MaterialProperty<RankTwoTensor> & _plastic_strain_old;
};
