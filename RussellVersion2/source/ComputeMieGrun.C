// Computes Mie Grunessien, and Hugoniot Pressure 
// RKM 2025

#include "ComputeMieGrun.h"

registerMooseObject("SolidMechanicsApp", ComputeMieGrun);

InputParameters
ComputeMieGrun::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Temperature Dependent Johnson Wilkins Lee");
  params.addParam<std::string>("base_name", "Optional parameter that allows the user to define multiple mechanics material systems on the same block, i.e. for multiple phases");

  params.addRequiredParam<std::string>("property_name",     "The property name to declare");
  params.addRequiredParam<std::string>("property_dT_name", "The property name for the temperature contribution dP/dT");
  params.addParam<std::string>("hugo_press_name","hugo_press",     "The property name to declare");

  params.addRequiredCoupledVar("temperature","Temperature");

  params.addParam<std::string>(      "density","density", "Name of Material Property that provides the initial density");
  params.addParam<std::string>( "bulk_modulus",  "bulk_modulus", "Name of Material Property that provides the bulk_modulus");
  params.addParam<std::string>("specific_heat", "specific_heat", "Name of Material Property that provides the specific heat");

  params.addRequiredParam<Real>(               "Gamma0", "Preshock Gruneisen coefficient G (or Gamma) in Mie-Gruneisen EOS");
  params.addRequiredParam<Real>(           "slope_UsUp", "Us-Up Hugonot :Us=s*Up+c0");
  params.addParam<Real>(                   "C0", "Unshocked Bulk Speed of Sound :Us=s*Up+c0");
  params.addParam<Real>("reference_temperature", 300, "Temperature reference");

  return params;
}

ComputeMieGrun::ComputeMieGrun(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters), //HeatSource(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),

    _property_name(getParam<std::string>("property_name")),
    _property(declareProperty<Real>(_property_name)),
    _property_dT_name(getParam<std::string>("property_dT_name")),
    _property_dT(declareProperty<Real>(_property_dT_name)),

    _hugo_press_name(getParam<std::string>("hugo_press_name")),
    _hugo_press(declareProperty<Real>(_hugo_press_name)),

    _elasticity_tensor_name(_base_name + "elasticity_tensor"),
    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_elasticity_tensor_name)),

    _temperature(coupledValue("temperature")),

    _density_name(      getParam<std::string>(      "density")),
    _specific_heat_name(getParam<std::string>("specific_heat")),

    _density(      getMaterialPropertyByName<Real>(      _density_name)),
    _specific_heat(getMaterialPropertyByName<Real>(_specific_heat_name)),

    _Gamma(                getParam<Real>(               "Gamma0")),
    _slope_UsUp(           getParam<Real>(           "slope_UsUp")),
    _C0(                   getParam<Real>(                   "C0")),
    _reference_temperature(getParam<Real>("reference_temperature")),

    _mechanical_strain(getMaterialPropertyByName<RankTwoTensor>(_base_name + "mechanical_strain"))
{}

void
ComputeMieGrun::initQpStatefulProperties()
{  
}

void
ComputeMieGrun::computeQpProperties()//computeQpStress()
{
  Real   J = _mechanical_strain[_qp].trace() + 1.0; //relative volume change;
  Real eta = 1-J;
  Real bulk_modulus;

  if(_C0) { 
    bulk_modulus = _density[_qp]*std::pow(_C0,2);
    }
  else {
    RankTwoTensor I2(RankTwoTensor::initIdentity);
    bulk_modulus = (1.0 / 9.0) * I2.doubleContraction(_elasticity_tensor[_qp] * I2);  //Voight Average Bulk Modulus
    } 

  _hugo_press[_qp] = bulk_modulus*eta / std::pow((1-_slope_UsUp*eta),2);

  _property[_qp] = _hugo_press[_qp] * (1 - (_Gamma*eta)/2)+_Gamma*(_density[_qp]*_specific_heat[_qp])*(_temperature[_qp] - _reference_temperature);
  _property_dT[_qp] = _Gamma*(_density[_qp]*_specific_heat[_qp]); 

  if(_property[_qp]<0){_property[_qp]=0;}
}










//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//