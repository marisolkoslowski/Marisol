/// Compute Melting Temperature using the Simon-Glatzel Equation

#include "SimonGlatzelMelt.h"

registerMooseObject("SolidMechanicsApp", SimonGlatzelMelt);

InputParameters
SimonGlatzelMelt::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Compute Melting Temperature using the Simon-Glatzel Equation");

  params.addRequiredParam<Real>("temp_melt_ref", "Reference Melting Temperature");
  params.addRequiredParam<Real>("pres_melt_ref", "Reference Melting Pressure");
  params.addRequiredParam<Real>(   "a_melt"    , "Melting Coefficient a");
  params.addRequiredParam<Real>(   "c_melt"    , "Melting Exponent c");

  return params;
}

SimonGlatzelMelt::SimonGlatzelMelt(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters), //HeatSource(parameters),
    _temp_melt_ref(getParam<Real>( "temp_melt_ref")),
    _pres_melt_ref(getParam<Real>( "pres_melt_ref")),
    _a_melt(       getParam<Real>( "a_melt")),
    _c_melt(       getParam<Real>( "c_melt")),
    _pressure_eos( getMaterialProperty<Real>( "pressure_eos")),
    _temp_melt(        declareProperty   <Real>("temp_melt"))
{}

void
SimonGlatzelMelt::initQpStatefulProperties()
{
  _temp_melt[_qp] = _temp_melt_ref * std::pow(1,1/_c_melt);
}

void
SimonGlatzelMelt::computeQpProperties()//computeQpStress()
{
  _temp_melt[_qp] = _temp_melt_ref * std::pow(1 + (_pressure_eos[_qp]-_pres_melt_ref)/_a_melt,1/_c_melt);
}























//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//