#include "SwitchingFunctionConstraintEta.h"

template<>
InputParameters validParams<SwitchingFunctionConstraintEta>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Lagrange multiplier kernel to constrain the sum of all switching functions in a multiphase system. This kernel acts on a non-conserved order parameter eta_i.");
  params.addParam<std::string>("h_name", "Switching Function Materials that provides h(eta_i)");
  params.addRequiredCoupledVar("lambda", "Lagrange multiplier");
  return params;
}

SwitchingFunctionConstraintEta::SwitchingFunctionConstraintEta(const std::string & name, InputParameters parameters) :
    DerivativeMaterialInterface<Kernel>(name, parameters),
    _h_name(getParam<std::string>("h_name")),
    _eta_name(_var.name()),
    _dh(getMaterialPropertyDerivative<Real>(_h_name, _eta_name)),
    _d2h(getMaterialPropertyDerivative<Real>(_h_name, _eta_name, _eta_name)),
    _lambda(coupledValue("lambda")),
    _lambda_var(coupled("lambda"))
{
}

Real
SwitchingFunctionConstraintEta::computeQpResidual()
{
  return _lambda[_qp] * _dh[_qp] * _test[_i][_qp];
}

Real
SwitchingFunctionConstraintEta::computeQpJacobian()
{
  return _lambda[_qp] * _d2h[_qp] * _phi[_j][_qp] * _test[_i][_qp];
}

Real
SwitchingFunctionConstraintEta::computeQpOffDiagJacobian(unsigned int j_var)
{
  if (j_var == _lambda_var)
    return _phi[_j][_qp] * _dh[_qp] * _test[_i][_qp];
  else
    return 0.0;
}
