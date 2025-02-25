// Computes heat generation from EOS pressure change, and artificial viscosity 
// RKM 2025

#include "Heat_EOS.h"

registerMooseObject("SolidMechanicsApp",Heat_EOS);

InputParameters
Heat_EOS::validParams()
{
  //InputParameters params = validParams<HeatSource>();
  InputParameters params = HeatSource::validParams();
  params.addClassDescription("Thermal expansion heat source kernel generic kernel for finite strain for Mie Gruneisen equation of state (Menon, 2014) (Zhang, 2011)");

  params.addRequiredParam<std::string>("dPress_dT_name", "The property name for the temperature contribution dP/dT");

  params.addRequiredParam<Real>("beta_av", "artificial viscosity beta parameter");
  return params;
}

Heat_EOS::Heat_EOS(const InputParameters & parameters)
  : HeatSource(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),

    _mechanical_strain(    getMaterialPropertyByName   <RankTwoTensor>(_base_name + "mechanical_strain")),
    _mechanical_strain_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "mechanical_strain")),

    _total_strain(    getMaterialPropertyByName<RankTwoTensor>   (_base_name + "total_strain")),
    _total_strain_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "total_strain")),

    _dPress_dT_name(getParam<std::string>("dPress_dT_name")),
    _dPress_dT(getMaterialPropertyByName<Real>(_dPress_dT_name)),

    _beta_av(getParam<Real>("beta_av")),

    _stress_artificial(     getMaterialPropertyByName         <RankTwoTensor>(_base_name + "stress_artificial"))
{}

Real
Heat_EOS::computeQpResidual()
{
  Real mechanical_strain_trace;
  Real total_strain_trace;

  RankTwoTensor mechanical_strain_dot;
  RankTwoTensor total_strain_dot;

  Real mechanical_strain_trace_dot;
  Real total_strain_trace_dot;

  // compute total, elastic strain with time derivatives
  mechanical_strain_trace = _mechanical_strain[_qp].trace();
  total_strain_trace      =      _total_strain[_qp].trace();

  mechanical_strain_dot = (_mechanical_strain[_qp] - _mechanical_strain_old[_qp]) / _dt;
  total_strain_dot      = (_total_strain[_qp]      -      _total_strain_old[_qp]) / _dt;

  mechanical_strain_trace_dot = (mechanical_strain_trace - _mechanical_strain_old[_qp].trace()) / _dt;
  total_strain_trace_dot      = (total_strain_trace      -      _total_strain_old[_qp].trace()) / _dt;

  Real q_tot;
  q_tot = 0.0;

  //Thermomechanical Coupling
    Real q_eos;
    q_eos = -_u[_qp] * _dPress_dT[_qp] * mechanical_strain_trace_dot;
    q_tot += q_eos;    


  //Artificial Visocosity:Heating
    Real J;
    J = 1 + mechanical_strain_trace;
    if (J < 1.0) {
       q_tot += std::abs(_beta_av * _stress_artificial[_qp].doubleContraction(mechanical_strain_dot));

    }

  return -_test[_i][_qp] * q_tot;
}

Real
Heat_EOS::computeQpJacobian()
{
  Real mechanical_strain_trace;
  Real total_strain_trace;

  RankTwoTensor mechanical_strain_dot;
  RankTwoTensor total_strain_dot;

  Real mechanical_strain_trace_dot;
  Real total_strain_trace_dot;

  // compute total, elastic strain with time derivatives
  mechanical_strain_trace = _mechanical_strain[_qp].trace();
  total_strain_trace      =      _total_strain[_qp].trace();

  mechanical_strain_dot = (_mechanical_strain[_qp] - _mechanical_strain_old[_qp]) / _dt;
  total_strain_dot      = (_total_strain[_qp]      -      _total_strain_old[_qp]) / _dt;

  mechanical_strain_trace_dot = (mechanical_strain_trace - _mechanical_strain_old[_qp].trace()) / _dt;
  total_strain_trace_dot      = (total_strain_trace      -      _total_strain_old[_qp].trace()) / _dt;


  Real q_tot;
  q_tot = 0.0;

  // compute thermomechanical coupled term
  Real q_eos;
    q_eos = -_u[_qp] * _dPress_dT[_qp] * mechanical_strain_trace_dot;
  q_tot += q_eos; 

  // compute artificial viscosity stress and heat source terms
  Real J;
  J = 1 + _mechanical_strain[_qp].trace();  
  if (J < 1.0) {
      q_tot += std::abs(_beta_av * _stress_artificial[_qp].doubleContraction(mechanical_strain_dot));  // to ensure that only when in compression the av term contributes to heat generation
  }


  
  return - q_tot * _phi[_j][_qp] * _test[_i][_qp];
}
