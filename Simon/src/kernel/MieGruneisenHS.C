#include "MieGruneisenHS.h"

registerMooseObject("beaverApp", MieGruneisenHS);

InputParameters
MieGruneisenHS::validParams()
{
  InputParameters params = HeatSource::validParams(); //we use AD to avoid complicated Jacobian computation
  params.addClassDescription("compute elastic thermomechanical coupling heat source term");
  params.addRequiredParam<Real>("Gamma", "Grunseisen Parameter");
  params.addRequiredParam<Real>("T_ref", "reference temperature for thermal expansion");
  params.addRequiredParam<Real>("C0", "artificial viscosity C0 parameter");
  params.addRequiredParam<Real>("C1", "artificial viscosity C1 parameter");
  params.addRequiredParam<Real>("beta_av", "disipation coefficient for artificial viscosity");
  params.addRequiredCoupledVar("temperature", "temperature");
  params.addRequiredParam<Real>("viscosity_type", "viscosity_type");
  params.addRequiredParam<Real>("element_size", "element_size");
  return params;
}

MieGruneisenHS::MieGruneisenHS(const InputParameters & parameters)
  : HeatSource(parameters),
    _Gamma(getParam<Real>("Gamma")),
    _T_ref(getParam<Real>("T_ref")),
    _rho(getMaterialProperty<Real>("density")),
    _Cv(getMaterialProperty<Real>("specific_heat")),
    _T(coupledValue("temperature")),
    _mechanical_strain(getMaterialProperty<RankTwoTensor>("mechanical_strain")),
    _mechanical_strain_old(getMaterialPropertyOld<RankTwoTensor>("mechanical_strain")),
    _C0(getParam<Real>("C0")),
    _C1(getParam<Real>("C1")),
    _current_elem(_assembly.elem()),
    _beta_av(getParam<Real>("beta_av")),
    _elasticity_tensor(getMaterialProperty<RankFourTensor>("elasticity_tensor")),
    _viscosity_type(getParam<Real>("viscosity_type")),
    _Le(getParam<Real>("element_size"))
{}

Real
MieGruneisenHS::computeQpResidual()
{
  //initialize the symmetric identity tensors
  RankTwoTensor I2(RankTwoTensor::initIdentity);

  //compute sound speed and bulk modulus from elasticity tensor
  Real K0 = (1.0 / 9.0) * I2.doubleContraction(_elasticity_tensor[_qp] * I2);
  Real ss = std::sqrt(K0 / _rho[_qp]);

  //Compute Mie Gruneisen pressure
  Real P_eos;
  Real res_eos;
  Real J = 1.0 + _mechanical_strain[_qp].trace();
  //Real J_prime = 1.0 / (1.0 + J); //this is equivalent to v0/v
  RankTwoTensor epsilon_dot = (_mechanical_strain[_qp] - _mechanical_strain_old[_qp]) / _dt;
  Real J_dot = epsilon_dot.trace();
  res_eos = - _Gamma * _rho[_qp] * _Cv[_qp] * _T[_qp] * epsilon_dot.trace();

  //Compute artificial viscosity term
  Real P_av;
  Real res_av;
  Real Le;
  //assign element size for viscosity
  if(_viscosity_type == 1){ //representing constant
    Le = _Le;
  }
  else {
    Le = _current_elem->hmax();
  }

  P_av = _C0 * _rho[_qp] * (J_dot * std::abs(J_dot) / std::pow(J, 2.0)) * std::pow(Le, 2.0);
  P_av += _C1 * _rho[_qp] * ss * (J_dot / J) * Le;
  RankTwoTensor P_av_tensor = P_av * I2; //create a rank two tensor for the volumetric av stress
  res_av = _beta_av * P_av_tensor.doubleContraction(epsilon_dot);

  Real res_tot = res_eos + res_av;

  return res_tot * _test[_i][_qp];

}

Real
MieGruneisenHS::computeQpJacobian()
{
  Real Jac_eos;
  Real Jac_av;
  Real Jac_tot;

  RankTwoTensor epsilon_dot;
  epsilon_dot = (_mechanical_strain[_qp] - _mechanical_strain_old[_qp]) / _dt;

  Jac_eos = - _Gamma * _rho[_qp] * _Cv[_qp] * epsilon_dot.trace();
  Jac_av = 0.0;
  Jac_tot = Jac_eos + Jac_av;
  return Jac_tot * _test[_i][_qp] * _phi[_j][_qp];
}

Real
MieGruneisenHS::computeQpOffDiagonalJacobian(unsigned int jvar)
{
  Real OffJac;
  OffJac = 0.0;
  return OffJac * _test[_i][_qp] * _phi[_j][_qp];
}