/// Computes pressure from a jones-wilkens-lee equation of state.

#include "ComputeJWL_Grun.h"

registerMooseObject("SolidMechanicsApp", ComputeJWL_Grun);

InputParameters
ComputeJWL_Grun::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Temperature Dependent Johnson Wilkins Lee");
  params.addParam<std::string>("base_name", "Optional parameter that allows the user to define multiple mechanics material systems on the same block, i.e. for multiple phases");
  params.addRequiredParam<std::string>("property_name", "The property name to declare");
  params.addRequiredParam<std::string>("property_dT_name", "The property name for the temperature contribution dP/dT");
  params.addParam<std::string>("specific_heat", "specific_heat", "Name of Material Property that provides the specific heat");
  params.addRequiredCoupledVar("temperature","temperature");
  params.addRequiredParam<Real>(              "A",                   "A parameter for JWL eos");
  params.addRequiredParam<Real>(              "B",                   "B parameter for JWL eos");
  params.addRequiredParam<Real>(             "R1",                  "R1 parameter for JWL eos");
  params.addRequiredParam<Real>(             "R2",                  "R2 parameter for JWL eos");
  params.addRequiredParam<Real>(          "omega",               "omega parameter for JWL eos");
  params.addRequiredParam<Real>("ref_temperature",                            "O energy State");
  return params;
}

ComputeJWL_Grun::ComputeJWL_Grun(const InputParameters & parameters)
  :  Material(parameters), //DerivativeMaterialInterface<Material>(parameters), //HeatSource(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _mechanical_strain(getMaterialProperty   <RankTwoTensor>( "mechanical_strain")),
    _property_name(getParam<std::string>("property_name")),
    _property(declareProperty<Real>(_property_name)),
    _property_dT_name(getParam<std::string>("property_dT_name")),
    _property_dT(declareProperty<Real>(_property_dT_name)),
    _specific_heat_name(getParam<std::string>("specific_heat")),
    _Cv(getMaterialPropertyByName<Real>(_specific_heat_name)),
    _temperature(coupledValue("temperature")),
    _A(              getParam<Real>(              "A")),
    _B(              getParam<Real>(              "B")),
    _R1(             getParam<Real>(             "R1")),
    _R2(             getParam<Real>(             "R2")),
    _omega(          getParam<Real>(          "omega")),
    _ref_temperature(getParam<Real>("ref_temperature")),
    _density(      getMaterialProperty<Real>(     "density"))
{}

void
ComputeJWL_Grun::initQpStatefulProperties()
{
}

void
ComputeJWL_Grun::computeQpProperties()//computeQpStress()
{
  Real V = _mechanical_strain[_qp].trace() + 1.0; //relative volume change;
  _property[_qp] = _A*(1-_omega/(_R1*V))*std::exp(-_R1*V) + _B*(1-_omega/(_R2*V))*std::exp(-_R2*V) + (_omega*_density[_qp]) * (_Cv[_qp] * (_temperature[_qp]-_ref_temperature))/V; //Density may have a /V term to get to current density
  _property_dT[_qp] = _omega*_density[_qp]* _Cv[_qp]/V; //Density may have a /V term to get to current density
  if(_property[_qp]<0){_property[_qp]=0;_property_dT[_qp] = 0;}
}










//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//