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
#include "LinearInterpolation.h"

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

  void computePlasticStrainRate(const ADReal & effective_trial_stress,
                         const ADReal & scalar = 0.0);

  void computeMisorientationVariable();
  void computeShearStressDerivative(const ADRankFourTensor & elasticity_tensor);

  virtual void updateInternalStateVariables(const ADReal & effective_trial_stress,
                                            const ADReal & scalar=0.0,
                                            const ADReal & scalar_increment=0.0) override;

  virtual void initQpStatefulProperties() override;
  virtual void propagateQpStatefulProperties() override;
  virtual void computeStressFinalize(const ADRankTwoTensor & plastic_strain_increment) override;

  /**
   * Toggle between solid/fluid stress updates based on element temperature
   * if solid: compute stress from J2 radial return mapping method
   * if fluid: compute stress from the Newtonian fluid constitutive eqn
   * @param strain_increment              Sum of elastic and inelastic strain increments
   * @param inelastic_strain_increment    Inelastic strain increment calculated by this class
   * @param rotation increment            Not used by this class
   * @param stress_new                    New trial stress from pure elastic calculation
   * @param stress_old                    Old state of stress
   * @param elasticity_tensor             Rank 4 C_{ijkl}, must be isotropic
   * @param elastic_strain_old            Old state of total elastic strain
   */
  virtual void updateState(ADRankTwoTensor & strain_increment,
                           ADRankTwoTensor & inelastic_strain_increment,
                           const ADRankTwoTensor & rotation_increment,
                           ADRankTwoTensor & stress_new,
                           const RankTwoTensor & stress_old,
                           const ADRankFourTensor & elasticity_tensor,
                           const RankTwoTensor & elastic_strain_old) override;

  virtual Real computeReferenceResidual(const ADReal & effective_trial_stress, const ADReal & scalar_effective_inelastic_strain) override;

  /// Temperature variable value
  const ADVariableValue * _temperature;

  /// Simulation start time
  const Real _start_time;

  ///
  /// Model parameters (fluid part)
  ///

  /// Melt temperature [K]
  const Real _theta_melt;

  /// Melt viscosity [Pa*s]
  const Real _mu_melt;

  /// Bulk modulus melt [Pa]
  const Real _K_melt;

  ///
  /// Model parameters (solid part)
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

  /// Misorientation variable hardening exponent [-]
  const Real _r;

  /// Components for computing the derivatives
  ADReal _C1, _C2;

  /// Plastic strain rate (flow rule function)
  ADReal _plastic_strain_rate;

  /// Temperature dependent yield stress
  ADReal _yield_stress;

  /// Temperature dependent shear modulus
  ADReal _shear_modulus;

  /// Derivative of shear modulus w.r.t. temperature
  ADReal _shear_modulus_derivative;

  /// LinearInterpolation objects for Young's modulus and Poisson's ratio
  std::unique_ptr<LinearInterpolation> _data_youngs_modulus, _data_poissons_ratio;

  /// Isotropic harderning internal state variable
  ADMaterialProperty<Real> & _hardening_variable;
  const MaterialProperty<Real> & _hardening_variable_old;

  /// Misorientation Variable
  ADMaterialProperty<Real> & _misorientation_variable;
  const MaterialProperty<Real> & _misorientation_variable_old;

  /// Plastic strain material property
  ADMaterialProperty<RankTwoTensor> & _plastic_strain;
  const MaterialProperty<RankTwoTensor> & _plastic_strain_old;

  /// Pressure
  ADMaterialProperty<Real> & _pressure;
  const MaterialProperty<Real> & _pressure_old;

  /// Strain from fluid deformation
  ADMaterialProperty<RankTwoTensor> & _strain_fluid;
  const MaterialProperty<RankTwoTensor> & _strain_fluid_old;
};
