#include "ComputeFiniteStress_EOS_ART.h"

registerMooseObject("SolidMechanicsApp", ComputeFiniteStress_EOS_ART);

InputParameters ComputeFiniteStress_EOS_ART::validParams() {
  InputParameters params = ComputeStressBase::validParams();
  params.addClassDescription("Computes the stress and free energy derivatives for the phase field fracture model, with small strain considers Mie Gruneisen EOS and artificial viscosity damping");

  params.addRequiredParam<std::string>("pressure_eos_name",     "The property name to declare");
  params.addRequiredParam<Real>("C0", "C0 Artificial Visosity Parameter");
  params.addRequiredParam<Real>("C1", "C1 Artificial Visosity Parameter");

  return params;
}

ComputeFiniteStress_EOS_ART::ComputeFiniteStress_EOS_ART(const InputParameters & parameters) : ComputeStressBase(parameters),
  //Name Finder
    _elasticity_tensor_name("elasticity_tensor"),
    _elasticity_tensor(    getMaterialPropertyByName   <RankFourTensor>(      _elasticity_tensor_name)),
    _rotation_total(           declareProperty         <RankTwoTensor> (_base_name + "rotation_total")),
    _rotation_total_old(   getMaterialPropertyOldByName<RankTwoTensor> (_base_name + "rotation_total")),
    _strain_increment(     getMaterialPropertyByName   <RankTwoTensor> (_base_name + "strain_increment")),
    _rotation_increment(   getMaterialPropertyByName   <RankTwoTensor> (_base_name + "rotation_increment")),
    _stress_old(           getMaterialPropertyOldByName<RankTwoTensor> (_base_name + "stress")),
    _elastic_strain_old(   getMaterialPropertyOldByName<RankTwoTensor> (_base_name + "elastic_strain")),
    _mechanical_strain_old(getMaterialPropertyOldByName<RankTwoTensor> (_base_name + "mechanical_strain")),

  //Specific Volume
    _V(                 getMaterialProperty                  <Real>(                          "V")),
    _mu(                getMaterialProperty                  <Real>(                         "mu")),

  //Pressure EOS
    _pressure_eos_name(getParam<std::string>("pressure_eos_name")),
    _pressure_eos(getMaterialPropertyOldByName<Real>(_pressure_eos_name)),
    //_pressure_eos(      getMaterialPropertyOld               <Real>(               "pressure_eos")),

  //Artificial Viscosity
    _density(           getMaterialProperty                  <Real>(               "density")),
    _stress_artificial(     declareProperty         <RankTwoTensor>(_base_name + "stress_artificial")),
    _C0(getParam<Real>("C0")),
    _C1(getParam<Real>("C1")),
    _current_elem(_assembly.elem()),

  //Debugging
    _report_output(             declareProperty         <Real> ("report_output"))
{}

void ComputeFiniteStress_EOS_ART::computeQpStress() {
  RankTwoTensor intermediate_stress;
  RankFourTensor elasticity_tensor_rotated = _elasticity_tensor[_qp];
  //elasticity_tensor_rotated.rotate(_rotation_total_old[_qp]);

  intermediate_stress = elasticity_tensor_rotated * (_elastic_strain_old[_qp] + _strain_increment[_qp]);
  _stress[_qp] = intermediate_stress;
  //_rotation_total[_qp] = _rotation_increment[_qp] * _rotation_total_old[_qp];
  //_Jacobian_mult[_qp] = elasticity_tensor_rotated; // This is NOT the exact jacobian
  //_stress[_qp] = _rotation_increment[_qp] * intermediate_stress * _rotation_increment[_qp].transpose();
  //_elastic_strain[_qp] = _mechanical_strain[_qp];

  EOSUpdate(); //Pressure Equation of State: Finding stress from Pressure EOS 
  ArtificialViscosity(); //Curve Smoothing

  _stress[_qp] = stress_hyd + stress_dev;
}

void ComputeFiniteStress_EOS_ART::initQpStatefulProperties() {
  ComputeStressBase::initQpStatefulProperties();
  RankTwoTensor identity_rotation(RankTwoTensor::initIdentity);
  _rotation_total[_qp] = identity_rotation;
}

void ComputeFiniteStress_EOS_ART::EOSUpdate  () {//Combining partial pressure, and calculating resulting stress
  //EOS Update
  stress_dev = _stress[_qp] - (_stress[_qp].trace()/3) * I2;
  stress_hyd = -1 * _pressure_eos[_qp] * I2; //pressure is opposite stress

  //Real bulk_modulus = (1.0 / 9.0) * I2.doubleContraction(_elasticity_tensor[_qp] * I2); //Luscher 2017
  //RankTwoTensor _stress_eos_elastic = -_pressure_eos[_qp] * I2;
  //RankTwoTensor _stress_cpl_elastic = _elasticity_tensor[_qp] * (_elastic_strain_old[_qp]+_strain_increment[_qp]) - bulk_modulus * (_mechanical_strain[_qp].trace()) * I2;
  //_stress_previous_calc_method[_qp] = _stress_eos_elastic + _stress_cpl_elastic;
}

void ComputeFiniteStress_EOS_ART::ArtificialViscosity() { //ArtificialViscosity: 
  //Initialize
  Real jacob, jacob_dot, stress_artificial, Le, K0, ss_bulk;
  Le = _current_elem->hmax(); //"Compute the element size using Elem::hmin() or Elem::hmax() from libMesh."

  jacob = 1.0 + _mechanical_strain[_qp].trace();
  jacob_dot = (_mechanical_strain[_qp].trace() - _mechanical_strain_old[_qp].trace()) / _dt;

  K0 = (1.0 / 9.0) * I2.doubleContraction(_elasticity_tensor[_qp] * I2);
  ss_bulk = std::sqrt(K0 / _density[_qp]); //store sound speed as a function of actual bulk modulus from elasticity tensor

  _stress_artificial[_qp] = 0.0*I2;
  if (jacob < 1.0) {
    _stress_artificial[_qp] = I2 * ((_C0*_density[_qp]*jacob_dot*std::abs(jacob_dot)/std::pow(jacob,2.0))*std::pow(Le,2.0) + (_C1*_density[_qp]*ss_bulk*jacob_dot/jacob)*Le);
  }
  stress_hyd += _stress_artificial[_qp]; //Only want smoothing in the volumetric section

  _report_output[_qp] = Le; 
}






















//--Have-A-Nice-Day------------------------------------------------------------------------------------------------------------------------