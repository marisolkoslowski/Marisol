
#include "SpecificVolume.h"


registerMooseObject("SolidMechanicsApp", SpecificVolume);

InputParameters
SpecificVolume::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Computes the stress and free energy derivatives for the phase field fracture model, with small strain considers Mie Gruneisen EOS and artificial viscosity damping");
  return params;
}

SpecificVolume::SpecificVolume(const InputParameters & parameters)
  : Material(parameters), //DerivativeMaterialInterface<Material>(parameters),   ///ComputeGeneralStressBase

    _bulk_modulus(declareProperty<Real>("bulk_modulus")),
    _ss_prop(     declareProperty<Real>(     "ss_prop")),
    _V(           declareProperty<Real>(           "V")),
    _density_corr(declareProperty<Real>("density_corr")),
    _mu(          declareProperty<Real>(          "mu")),

    _density(getMaterialProperty<Real>("density")),
    _mechanical_strain(getMaterialProperty   <RankTwoTensor>( "mechanical_strain")),
    _elasticity_tensor(getMaterialPropertyByName <RankFourTensor>( "elasticity_tensor"))
{
}


void
SpecificVolume::initQpStatefulProperties()
{
  _V[_qp] = 1.0;
}


void
SpecificVolume::computeQpProperties()
{
  Real delta,  eta, peos, trD, q_bv;
  RankTwoTensor stress_eos, stress, stress_cpl;
  RankTwoTensor I2(RankTwoTensor::initIdentity);
  RankFourTensor I4sym(RankFourTensor::initIdentitySymmetricFour);

  //Bulk Modulus from Elasticity Tensor
   _bulk_modulus[_qp] = (1.0 / 9.0) * I2.doubleContraction(_elasticity_tensor[_qp] * I2);

  //Speed of Sound from Bulk Modulus
   _ss_prop[_qp] = std::sqrt(_bulk_modulus[_qp] / _density[_qp]);

  //Specific volume and related values
   _V[_qp] = _mechanical_strain[_qp].trace() + 1.0; //relative volume change;
   _density_corr[_qp] = _density[_qp] / _V[_qp];
   _mu[_qp] = (1/_V[_qp]) - 1;
}








