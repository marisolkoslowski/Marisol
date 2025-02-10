/// Computes pressure from a jones-wilkens-lee equation of state.

#include "ComputeJWLdependentT.h"

registerMooseObject("SolidMechanicsApp", ComputeJWLdependentT);

InputParameters
ComputeJWLdependentT::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Temperature Dependent Johnson Wilkins Lee");
  params.addParam<std::string>("base_name", "Optional parameter that allows the user to define multiple mechanics material systems on the same block, i.e. for multiple phases");
  params.addRequiredParam<std::string>("property_name", "The property name to declare");
  params.addRequiredCoupledVar("temperature","temperature");
  params.addRequiredParam<Real>(              "A",                   "A parameter for JWL eos");
  params.addRequiredParam<Real>(              "B",                   "B parameter for JWL eos");
  params.addRequiredParam<Real>(             "R1",                  "R1 parameter for JWL eos");
  params.addRequiredParam<Real>(             "R2",                  "R2 parameter for JWL eos");
  params.addRequiredParam<Real>(          "omega",               "omega parameter for JWL eos");
  params.addRequiredParam<Real>(             "Cv", "heat capacity in units energy/temperature");
  params.addRequiredParam<Real>("ref_temperature",                            "O energy State");
  return params;
}

ComputeJWLdependentT::ComputeJWLdependentT(const InputParameters & parameters)
  :  Material(parameters), //DerivativeMaterialInterface<Material>(parameters), //HeatSource(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _property_name(getParam<std::string>("property_name")),
    _property(declareProperty<Real>(_property_name)),
    _temperature(coupledValue("temperature")),
    _A(              getParam<Real>(              "A")),
    _B(              getParam<Real>(              "B")),
    _R1(             getParam<Real>(             "R1")),
    _R2(             getParam<Real>(             "R2")),
    _omega(          getParam<Real>(          "omega")),
    _Cv(             getParam<Real>(             "Cv")),
    _ref_temperature(getParam<Real>("ref_temperature")),
    _V(            getMaterialProperty<Real>(           "V")),
    _density(      getMaterialProperty<Real>(     "density"))
{}

void
ComputeJWLdependentT::initQpStatefulProperties()
{
  Real pres_piece1, pres_piece2, pres_piece3;
  pres_piece1 = _A*std::exp(-_R1*1);
  pres_piece2 = _B*std::exp(-_R2*1);
  pres_piece3 = _omega * _Cv * (_temperature[_qp]-_ref_temperature) / 1; // ref
  _property[_qp] = (pres_piece1+pres_piece2+pres_piece3);
  if(_property[_qp]<0){_property[_qp]=0;}
}

void
ComputeJWLdependentT::computeQpProperties()//computeQpStress()
{
  Real pres_piece1, pres_piece2, pres_piece3; Real press;
  pres_piece1 = _A*std::exp(-_R1*_V[_qp]);
  pres_piece2 = _B*std::exp(-_R2*_V[_qp]);
  pres_piece3 = _omega *_Cv * (_temperature[_qp]-_ref_temperature) / _V[_qp]; // ref
  _property[_qp] = (pres_piece1+pres_piece2+pres_piece3);
  if(_property[_qp]<0){_property[_qp]=0;}
}










//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//