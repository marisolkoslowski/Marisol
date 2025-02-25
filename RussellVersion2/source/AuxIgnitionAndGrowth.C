#include "AuxIgnitionAndGrowth.h"
#include "MathUtils.h"
registerMooseObject("SolidMechanicsApp", AuxIgnitionAndGrowth);
//https://pubs.aip.org/aip/acp/article/1793/1/040015/581730/Shock-initiation-experiments-with-ignition-and

InputParameters
AuxIgnitionAndGrowth::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Pressure Dependent IG Reaction Model");
  params.addRequiredParam<Real>("F_igMin", "Parameter for Pressure Dependent Reaction Rate");
  params.addRequiredParam<Real>("F_igMax", "Parameter for Pressure Dependent Reaction Rate");
  params.addRequiredParam<Real>( "I_chem", "Parameter for Pressure Dependent Reaction Rate");
  params.addRequiredParam<Real>( "a_chem", "Parameter for Pressure Dependent Reaction Rate");
  params.addRequiredParam<Real>( "b_chem", "Parameter for Pressure Dependent Reaction Rate");
  params.addRequiredParam<Real>( "x_chem", "Parameter for Pressure Dependent Reaction Rate");
  params.addRequiredParam<Real>("F_G1Min", "Parameter for Pressure Dependent Reaction Rate");
  params.addRequiredParam<Real>("F_G1Max", "Parameter for Pressure Dependent Reaction Rate");
  params.addRequiredParam<Real>("G1_chem", "Parameter for Pressure Dependent Reaction Rate");
  params.addRequiredParam<Real>("c1_chem", "Parameter for Pressure Dependent Reaction Rate");
  params.addRequiredParam<Real>("d1_chem", "Parameter for Pressure Dependent Reaction Rate");
  params.addRequiredParam<Real>("y1_chem", "Parameter for Pressure Dependent Reaction Rate");
  params.addRequiredParam<Real>("F_G2Min", "Parameter for Pressure Dependent Reaction Rate");
  params.addRequiredParam<Real>("F_G2Max", "Parameter for Pressure Dependent Reaction Rate");
  params.addRequiredParam<Real>("G2_chem", "Parameter for Pressure Dependent Reaction Rate");
  params.addRequiredParam<Real>("c2_chem", "Parameter for Pressure Dependent Reaction Rate");
  params.addRequiredParam<Real>("d2_chem", "Parameter for Pressure Dependent Reaction Rate");
  params.addRequiredParam<Real>("y2_chem", "Parameter for Pressure Dependent Reaction Rate");
  return params;
}

AuxIgnitionAndGrowth::AuxIgnitionAndGrowth(const InputParameters & parameters) : AuxKernel(parameters),
  _F_igMin(getParam<Real>("F_igMin")),
  _F_igMax(getParam<Real>("F_igMax")),
  _I_chem( getParam<Real>( "I_chem")),
  _a_chem( getParam<Real>( "a_chem")),
  _b_chem( getParam<Real>( "b_chem")),
  _x_chem( getParam<Real>( "x_chem")),
  _F_G1Min(getParam<Real>("F_G1Min")),
  _F_G1Max(getParam<Real>("F_G1Max")),
  _G1_chem(getParam<Real>("G1_chem")),
  _c1_chem(getParam<Real>("c1_chem")),
  _d1_chem(getParam<Real>("d1_chem")),
  _y1_chem(getParam<Real>("y1_chem")),
  _F_G2Min(getParam<Real>("F_G2Min")),
  _F_G2Max(getParam<Real>("F_G2Max")),
  _G2_chem(getParam<Real>("G2_chem")),
  _c2_chem(getParam<Real>("c2_chem")),
  _d2_chem(getParam<Real>("d2_chem")),
  _y2_chem(getParam<Real>("y2_chem")),
 // Chemical:Pressure Dependent IG Reaction Model
  _mu(          getMaterialProperty<Real>(          "mu")),
  _pressure_eos(getMaterialProperty<Real>("pressure_eos")),
  _u_old(uOld())
{}


//Real
Real
AuxIgnitionAndGrowth::computeValue() {
  Real rate_piece1 = 0; Real rate_piece2 = 0; Real rate_piece3 = 0;
  Real lambda = _u_old[_qp];
  if((_u[_qp] >= _F_igMin) && (_u[_qp] < _F_igMax)){
    rate_piece1 =  _I_chem * std::pow(1.0-_u[_qp], _b_chem) * std::pow(_mu[_qp]-_a_chem, _x_chem);
  }
  if((_u[_qp] >= _F_G1Min) && (_u[_qp] < _F_G1Max)){
    rate_piece2 = _G1_chem * std::pow(1.0-_u[_qp],_c1_chem) * std::pow(_u[_qp],_d1_chem) * std::pow(_pressure_eos[_qp],_y1_chem);
  }
  if((_u[_qp] >= _F_G2Min) && (_u[_qp] < _F_G2Max)){
    rate_piece3 = _G2_chem * std::pow(1.0-_u[_qp],_c2_chem) * std::pow(_u[_qp],_d2_chem) * std::pow(_pressure_eos[_qp],_y2_chem);
  }
  
  if (lambda > 1){ lambda = 1;}
  if(_pressure_eos[_qp] > 0) {lambda += _dt * (rate_piece1 + rate_piece2 + rate_piece3);}

  return lambda;
}
















//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//