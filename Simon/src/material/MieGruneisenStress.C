#include "MieGruneisenStress.h"

registerMooseObject("TensorMechanicsApp", MieGruneisenStress);

InputParameters
MieGruneisenStress::validParams()
{
  InputParameters params = ComputeStressBase::validParams();
  params.addClassDescription("Compute Elastic MieGruneisen EOS Pressure with Thermomechanical Coupling");
  params.addRequiredParam<Real>("Gamma", "Grunseisen Parameter");
  params.addRequiredParam<Real>("T_ref", "reference temperature for thermal expansion");
  params.addRequiredParam<Real>("C0", "artificial viscosity C0 parameter");
  params.addRequiredParam<Real>("C1", "artificial viscosity C1 parameter");
  params.addRequiredCoupledVar("temperature", "temperature");
  params.addRequiredParam<Real>("viscosity_type", "viscosity_type");
  params.addRequiredParam<Real>("element_size", "element_size");
  params.addRequiredParam<Real>("slope", "slope");
  params.addRequiredParam<Real>("A", "A");
  params.addRequiredParam<Real>("B", "B");
  params.addRequiredParam<Real>("R1", "R1");
  params.addRequiredParam<Real>("R2", "R2");
  params.addRequiredParam<Real>("omega", "omega");
  params.addRequiredCoupledVar("Y1", "Y1");
  params.addRequiredCoupledVar("Y2", "Y2");
  params.addRequiredCoupledVar("Y3", "Y3");
  params.addRequiredCoupledVar("Y4", "Y4");
  params.addRequiredParam<Real>("use_JWL", "use_JWL");
  return params;
}

MieGruneisenStress::MieGruneisenStress(
    const InputParameters & parameters)
  : ComputeStressBase(parameters),
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
    _elasticity_tensor(getMaterialProperty<RankFourTensor>("elasticity_tensor")),
    _viscosity_type(getParam<Real>("viscosity_type")),
    _Le(getParam<Real>("element_size")),
    _pressure_eos(declareProperty<Real>("pressure_eos")),
    _pressure_reacted(declareProperty<Real>("pressure_reacted")),
    _s(getParam<Real>("slope")),
    _J(declareProperty<Real>("J")),
    _A(getParam<Real>("A")),
    _B(getParam<Real>("B")),
    _R1(getParam<Real>("R1")),
    _R2(getParam<Real>("R2")),
    _omega(getParam<Real>("omega")),
    _Y1(coupledValue("Y1")),
    _Y2(coupledValue("Y2")),
    _Y3(coupledValue("Y3")),
    _Y4(coupledValue("Y4")),
    _use_JWL(getParam<Real>("use_JWL"))
{
}

void
MieGruneisenStress::initQpStatefulProperties()
{
  //initialize stress
  //_stress[_qp] = _stress_old[_qp]; //keep stress
  //OPTIONAL
}

void
MieGruneisenStress::computeQpStress()
{
  //initialize the symmetric identity tensors
  RankTwoTensor I2(RankTwoTensor::initIdentity);

  //compute sound speed and bulk modulus from elasticity tensors
  Real K0 = (1.0 / 9.0) * I2.doubleContraction(_elasticity_tensor[_qp] * I2);
  Real ss = std::sqrt(K0 / _rho[_qp]);

  //Compute Mie Gruneisen pressure
  Real P_eos;
  Real delta = _mechanical_strain[_qp].trace();
  Real eta = - delta;
  Real J = 1.0 + _mechanical_strain[_qp].trace();
  RankTwoTensor epsilon_dot = (_mechanical_strain[_qp] - _mechanical_strain_old[_qp]) / _dt;
  Real J_dot = epsilon_dot.trace();
  //P_eos = - _Gamma * _rho[_qp] * _Cv[_qp] * (_T[_qp] - _T_ref) * (1.0 / J); //thermal expansion
  //P_eos -= ((K0 * eta) / std::pow((1.0 - _s * eta), 2.0)) * ((_Gamma / 2.0) * ((1.0 / J) - 1.0) - 1.0);
  P_eos = (- K0 * eta * (1.0 - (_Gamma * eta / 2.0)) / std::pow((1.0 - _s * eta), 2.0)) - (_Gamma * _rho[_qp] * _Cv[_qp] * (_T[_qp] - _T_ref) * (1.0 / J));
  _pressure_eos[_qp] = P_eos;
  
  //compute reacted EOS using JWL formulation
  Real P_reacted;
  Real CV = _rho[_qp] * _Cv[_qp]; //extensive heat capacity
  P_reacted = _A * std::exp(- _R1 / (J)) + _B * std::exp(- _R2 / (J)) + _omega * CV * _T[_qp]; //temperature dependent JWL EOS
  _pressure_reacted[_qp] = P_reacted;
  
  Real P_HMX;
  
  //write partial pressures into stress tensor
  if (_use_JWL == 1)
  {
  	P_HMX = ((_Y1[_qp] + _Y2[_qp] + _Y3[_qp]) * P_eos) + (_Y4[_qp] * P_reacted);

  }
  else
  {
  	P_HMX = P_eos;
  }
  //write into stress tensor
  _stress[_qp] = P_HMX * I2;
	
  //Compute artificial viscosity term
  Real P_av;
  Real Le;

  //assign element size for viscosity
  if(_viscosity_type == 1){ //representing constant
    Le = _Le;
  }
  else {
    Le = _current_elem->hmax();
  }
  
  //construct the total pressure as the sum of partial pressures times mass fractions		
  P_av = _C0 * _rho[_qp] * (J_dot * std::abs(J_dot) / std::pow(J, 2.0)) * std::pow(Le, 2.0);
  P_av += _C1 * _rho[_qp] * ss * (J_dot / J) * Le;
  _stress[_qp] += P_av * I2;
}
