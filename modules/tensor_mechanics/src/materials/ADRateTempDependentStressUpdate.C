//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADRateTempDependentStressUpdate.h"

registerMooseObject("TensorMechanicsApp", ADRateTempDependentStressUpdate);

InputParameters
ADRateTempDependentStressUpdate::validParams()
{
  InputParameters params = ADRadialReturnStressUpdate::validParams();
  params.addClassDescription(
      "This class uses the stress update material in a radial return isotropic power law creep "
      "model. This class can be used in conjunction with other creep and plasticity materials "
      "for more complex simulations.");

  params.set<std::string>("effective_inelastic_strain_name") = "effective_plastic_strain";

  // Linear strain hardening parameters
  params.addCoupledVar("temperature", "Coupled temperature");
  params.addParam<Real>("start_time", 0.0, "Start time (if not zero)");

  // temperature dependent Young's modulus and Poisson's ratio
  params.addParam<std::vector<Real>>("Ex", "The temperature values");
  params.addParam<std::vector<Real>>("Ey", "The Young's modulus values");
  params.addParam<std::vector<Real>>("nux", "The temperature values");
  params.addParam<std::vector<Real>>("nuy", "The Poisson's ratio values");

  // Rate dependent plasticity parameters
  // Default values are from the paper by Michael E. Stender, et. al, 2018.
  // fluid parameters
  params.addParam<Real>("theta_melt", 1700, "Melt temperature [K]");
  params.addParam<Real>("mu_melt", 1.0e-6, "Melt viscosity [Pa*s]");
  params.addParam<Real>("K_melt", 2.2e9, "Bulk modulus melt [Pa]");

  // solid parameters
  params.addParam<Real>("Y0", 5.264e09, "Rate independent yield constant [Pa]");
  params.addParam<Real>("Y1", 2.688e05, "First rate independent yield temperature dependency [K]");
  params.addParam<Real>("Y2", 1.87e-03, "Second rate independent yield temperature dependency [1/K]");
  params.addParam<Real>("Y3", 8.683e02, "Third rate independent yield temperature dependency [K]");
  params.addParam<Real>("Y4", 3.316e01, "Fourth rate independent yield temperature dependency [-]");
  params.addParam<Real>("Hmu", 0.01, "Isotropic hardening shear coefficient [-]");
  params.addParam<Real>("f1", 9.178e-02, "First flow rule coefficient constant [1/s]");
  params.addParam<Real>("f2", 0.0, "Second flow rule coefficient constant [K]");
  params.addParam<Real>("n1", 0.0, "Flow rule exponent constant [-]");
  params.addParam<Real>("n2", 5.699e03, "Flow rule exponent temperature dependence [K]");
  params.addParam<Real>("Rd1", 8.565e02, "Isotropic dynamic recovery constant [Pa]");
  params.addParam<Real>("Rd2", 5.419e03, "Isotropic dynamic recovery temperature dependence [K]");
  params.addParam<Real>("hxi", 1.670e-03, " Misorientation variable hardening constant [m/(s Pa)]");
  params.addParam<Real>("r", 1.0, " Misorientation variable hardening exponent, 0.5<=r<=1 [-]");

  return params;
}

ADRateTempDependentStressUpdate::ADRateTempDependentStressUpdate(const InputParameters & parameters)
  : ADRadialReturnStressUpdate(parameters),
    _temperature(isParamValid("temperature") ? &adCoupledValue("temperature") : nullptr),
    _start_time(getParam<Real>("start_time")),
    _theta_melt(getParam<Real>("theta_melt")),
    _mu_melt(getParam<Real>("mu_melt")),
    _K_melt(getParam<Real>("K_melt")),
    _Y0(getParam<Real>("Y0")),
    _Y1(getParam<Real>("Y1")),
    _Y2(getParam<Real>("Y2")),
    _Y3(getParam<Real>("Y3")),
    _Y4(getParam<Real>("Y4")),
    _Hmu(getParam<Real>("Hmu")),
    _f1(getParam<Real>("f1")),
    _f2(getParam<Real>("f2")),
    _n1(getParam<Real>("n1")),
    _n2(getParam<Real>("n2")),
    _Rd1(getParam<Real>("Rd1")),
    _Rd2(getParam<Real>("Rd2")),
    _hxi(getParam<Real>("hxi")),
    _r(getParam<Real>("r")),
    _hardening_variable(declareADProperty<Real>(_base_name + "hardening_variable")),
    _hardening_variable_old(getMaterialPropertyOld<Real>(_base_name + "hardening_variable")),
    _misorientation_variable(declareADProperty<Real>(_base_name + "misorientation_variable")),
    _misorientation_variable_old(getMaterialPropertyOld<Real>(_base_name + "misorientation_variable")),
    _plastic_strain(declareADProperty<RankTwoTensor>(_base_name + "plastic_strain")),
    _plastic_strain_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "plastic_strain")),
    _pressure(declareADProperty<Real>(_base_name + "pressure")),
    _pressure_old(getMaterialPropertyOld<Real>(_base_name + "pressure")),
    _strain_fluid(declareADProperty<RankTwoTensor>(_base_name + "strain_fluid")),
    _strain_fluid_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "strain_fluid"))
{
  // Get linear interpolation of Young's modulus and Poisson'ratio
  // The goal is to compute _shear_modulus_derivative
  // Todo: how to get ADMaterialProperty derivative w.r.t coupled variable
  std::vector<Real> Ex, Ey, nux, nuy;
  if (!(parameters.isParamValid("Ex")&& parameters.isParamValid("Ey")&& parameters.isParamValid("nux") && parameters.isParamValid("nuy")))
    mooseError("Both 'x' and 'y' data must be specified for the Young's modulus and the Poisson's ratio. ");

  Ex = getParam<std::vector<Real>>("Ex");
  Ey = getParam<std::vector<Real>>("Ey");
  nux = getParam<std::vector<Real>>("nux");
  nuy = getParam<std::vector<Real>>("nuy");

  _data_youngs_modulus=libmesh_make_unique<LinearInterpolation>(Ex, Ey, false);
  _data_poissons_ratio=libmesh_make_unique<LinearInterpolation>(nux, nuy, false);
}

void
ADRateTempDependentStressUpdate::computeStressInitialize(const ADReal & effective_trial_stress,
                                                     const ADRankFourTensor & elasticity_tensor)
{
  _yield_stress = 0.5*_Y0*(1.0 + std::tanh(_Y2*(_Y3 - (*_temperature)[_qp])))/(_Y4 + std::exp(-_Y1/(*_temperature)[_qp]));

  _shear_modulus = ElasticityTensorTools::getIsotropicShearModulus(elasticity_tensor);

  computeShearStressDerivative(elasticity_tensor);

  _hardening_variable[_qp] = _hardening_variable_old[_qp];
  _misorientation_variable[_qp] = _misorientation_variable_old[_qp];
  _pressure[_qp] = _pressure_old[_qp];
  _strain_fluid[_qp]=_strain_fluid_old[_qp];

  updateInternalStateVariables(effective_trial_stress);
}

ADReal
ADRateTempDependentStressUpdate::computeResidual(const ADReal & effective_trial_stress,
                                             const ADReal & scalar)
{
  computePlasticStrainRate(effective_trial_stress, scalar);
  return _plastic_strain_rate * _dt - scalar - _C1*_C2*_hardening_variable[_qp]*_shear_modulus_derivative/_shear_modulus*_dt;
}

ADReal
ADRateTempDependentStressUpdate::computeDerivative(const ADReal & effective_trial_stress,
                                               const ADReal &  scalar )
{
  computePlasticStrainRate(effective_trial_stress, scalar);

  const ADReal theta = (*_temperature)[_qp];
  const ADReal creep_rate_derivative = _C1*(-3.0*_shear_modulus/(_hardening_variable[_qp] + _yield_stress)) - _C1*_C2*(_Hmu*_shear_modulus*(1.0+_misorientation_variable[_qp]/_hardening_variable[_qp]) -_Rd1*std::exp(-_Rd2/theta)*_hardening_variable[_qp]);

  ADReal tmp = creep_rate_derivative * _dt - 1.0;

  return creep_rate_derivative * _dt - 1.0;
}

void
ADRateTempDependentStressUpdate::computePlasticStrainRate(const ADReal & effective_trial_stress,
                                                      const ADReal & scalar)
{
  const ADReal theta = (*_temperature)[_qp];
  const ADReal stress_delta = effective_trial_stress - 3.0 * _shear_modulus * scalar;
  const ADReal ratio = stress_delta/(_hardening_variable[_qp] + _yield_stress);

  if (ratio>1.0)
  {
    _plastic_strain_rate = _f1*std::exp(-_f2*theta)*std::pow(std::sinh(ratio-1.0), _n1+_n2/theta);
    _C1 = (_n1+_n2/theta)*_plastic_strain_rate*std::cosh(ratio-1.0)/std::sinh(ratio-1.0);
    _C2 = (3.0*_shear_modulus*scalar-effective_trial_stress)/(_hardening_variable[_qp] + _yield_stress)/(_hardening_variable[_qp] + _yield_stress);
  }
  else
  {
    _plastic_strain_rate= 0.0;
    _C1=0.0;
    _C2=0.0;
  }

  // check value
  if(std::isinf(_plastic_strain_rate.value()) || std::isinf(_C1.value()) || std::isinf(_C2.value()))
    mooseError("Plastic strain variable out of bound.. check trial stress");
}

void
ADRateTempDependentStressUpdate::computeShearStressDerivative(const ADRankFourTensor & elasticity_tensor)
{
  Real dE = _data_youngs_modulus->sampleDerivative((*_temperature)[_qp].value());
  Real dnu = _data_poissons_ratio->sampleDerivative((*_temperature)[_qp].value());

  ADReal poissons_ratio = ElasticityTensorTools::getIsotropicPoissonsRatio(elasticity_tensor);
  ADReal youngs_modulus = ElasticityTensorTools::getIsotropicYoungsModulus(elasticity_tensor);

  _shear_modulus_derivative = (2.0*dE*(1.0+poissons_ratio) - 2.0*dnu*youngs_modulus)/4.0/(1.0+poissons_ratio)/(1.0+poissons_ratio);

  // std::cout<<"E: "<<youngs_modulus.value()<<", dE: "<<dE<<"; nu: "<<poissons_ratio.value()<<"; dnu: "<<dnu<<"; dG: "<<_shear_modulus_derivative.value()<<std::endl;
}

void
ADRateTempDependentStressUpdate::initQpStatefulProperties()
{
  _plastic_strain[_qp].zero();

  /// initilize _hardening_variable
  // @ t=0, _hardening_variable=exp(_shear_modulus_derivative/_shear_modulus*t)->_hardening_variable=1.0
  // similarly for the misorientation_variable

  if(_hardening_variable[_qp]<1e-10)
    _hardening_variable[_qp]=1.0;

  if(_misorientation_variable[_qp]<1e-10)
    _misorientation_variable[_qp]=1.0;

  ADRadialReturnStressUpdate::initQpStatefulProperties();
}

void
ADRateTempDependentStressUpdate::propagateQpStatefulProperties()
{
  _plastic_strain[_qp] = _plastic_strain_old[_qp];
  _hardening_variable[_qp]=_hardening_variable_old[_qp];
  _misorientation_variable[_qp]=_misorientation_variable_old[_qp];
  _pressure[_qp]= _pressure_old[_qp];

  propagateQpStatefulPropertiesRadialReturn();
}

void
ADRateTempDependentStressUpdate::computeStressFinalize(
    const ADRankTwoTensor & plastic_strain_increment)
{
  _plastic_strain[_qp] = _plastic_strain_old[_qp] + plastic_strain_increment;

  // if (_qp==1)
  //   std::cout<<"\t\t[qp= "<< _qp<<"], After Finalize: r="<<_hardening_variable[_qp].value()<<std::endl;
}

void
ADRateTempDependentStressUpdate::updateInternalStateVariables(
                                          const ADReal & effective_trial_stress,
                                          const ADReal & scalar,
                                          const ADReal & /*scalar_increment*/)
{
  const ADReal theta = (*_temperature)[_qp];

  /// Compute increment of isotropic harderning internal state variable
  ADReal hardening_variable_increment= _hardening_variable[_qp]*(_shear_modulus_derivative/_shear_modulus)+(_Hmu*_shear_modulus*(1.0+_misorientation_variable[_qp]/_hardening_variable[_qp])-_Rd1*std::exp(-_Rd2/theta)* _hardening_variable[_qp])*scalar;
  _hardening_variable[_qp]=_hardening_variable_old[_qp]+hardening_variable_increment;

  /// Compute increment of misorientation variable
  ADReal misorientation_variable_increment;
  const Real n_power = 1.0 - 1.0/_r;
  if (n_power<1e-10)
    misorientation_variable_increment = _misorientation_variable[_qp]*(_shear_modulus_derivative/_shear_modulus) + _hxi*_shear_modulus*std::abs(scalar);
  else
    misorientation_variable_increment = _misorientation_variable[_qp]*(_shear_modulus_derivative/_shear_modulus) + _hxi*_shear_modulus*std::pow(_misorientation_variable[_qp]/_shear_modulus , n_power)*std::abs(scalar);
  _misorientation_variable[_qp] = _misorientation_variable_old[_qp] + misorientation_variable_increment;

  // if (_qp==0)
  //   std::cout<<"\t\t[qp= "<< _qp<<"], Update: r="<<_hardening_variable[_qp].value()<<", dr= "<<hardening_variable_increment.value()<<", xi= "<<_misorientation_variable[_qp].value()<<", dxi= "<<misorientation_variable_increment.value()<<", plas_strain_rate: "<<std::abs(scalar).value()<<std::endl;

  computePlasticStrainRate(effective_trial_stress, scalar);
}

void
ADRateTempDependentStressUpdate::updateState(ADRankTwoTensor & strain_increment,
                                        ADRankTwoTensor & inelastic_strain_increment,
                                        const ADRankTwoTensor & rotation_increment,
                                        ADRankTwoTensor & stress_new,
                                        const RankTwoTensor & stress_old,
                                        const ADRankFourTensor & elasticity_tensor,
                                        const RankTwoTensor & elastic_strain_old)
{

  // accumulate pressure in preparation for calculations after melting
  ADRankTwoTensor strain_increment_total = strain_increment; //+inelastic_strain_increment;

  ADReal p_increment=_K_melt * strain_increment_total.trace();
  _pressure[_qp]=_pressure_old[_qp] + p_increment;

  // get average temperature
  ADReal temp = 0;
  for (unsigned int qp = 0; qp < _qrule->n_points(); ++qp)
    temp += (*_temperature)[qp];

  temp /= _qrule->n_points();

  if (temp >= _theta_melt)
  {
    ADRankTwoTensor strain_increment_rate = 1.0/_dt * strain_increment_total;
    RankTwoTensor I; I.setToIdentity();
    stress_new = _pressure[_qp]*I + 2.0*_mu_melt*strain_increment_rate.deviatoric();
    _strain_fluid[_qp] = elastic_strain_old + strain_increment;
  }
  else
  {
    // compute trial stress using the strain caused only by solid deformation
    stress_new = stress_new - elasticity_tensor*(_strain_fluid[_qp]);

    ADRadialReturnStressUpdate::updateState(strain_increment,
                                            inelastic_strain_increment,
                                            rotation_increment,
                                            stress_new,
                                            stress_old,
                                            elasticity_tensor,
                                            elastic_strain_old);
  }
}

Real
ADRateTempDependentStressUpdate::computeReferenceResidual(
  const ADReal & /*effective_trial_stress*/, const ADReal & scalar_effective_inelastic_strain)
{
  return  scalar_effective_inelastic_strain.value();
}
