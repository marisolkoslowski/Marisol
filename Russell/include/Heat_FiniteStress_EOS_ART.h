
#include "HeatSource.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"


// Forward Declarations
class Heat_FiniteStress_EOS_ART;

//template <>
//InputParameters validParams<ThermalExpansionHeatSourceFiniteStrainMieGruneisenNew>();

/**
 * This kernel calculates the heat source term corresponding to thermoelasticity
 * Mie Gruneisen equation of state (Menon, 2014) (Zhang, 2011)
 */
class Heat_FiniteStress_EOS_ART : public HeatSource
{
public:
  Heat_FiniteStress_EOS_ART(const InputParameters & parameters);
  static InputParameters validParams();
  
protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

private:
  
  std::string _base_name;

  //const MaterialProperty<RankTwoTensor> & _deformation_gradient; // deformation gradient
  //const MaterialProperty<RankTwoTensor> & _deformation_gradient_old; // deformation gradient, previous timestep

  //const MaterialProperty<RankFourTensor> & _elasticity_tensor; //elasticity tensor
  const MaterialProperty<RankTwoTensor> & _mechanical_strain;
  const MaterialProperty<RankTwoTensor> & _mechanical_strain_old;
  //const MaterialProperty<RankTwoTensor> & _elastic_strain;
  //const MaterialProperty<RankTwoTensor> & _elastic_strain_old;
  const MaterialProperty<RankTwoTensor> & _total_strain;
  const MaterialProperty<RankTwoTensor> & _total_strain_old;

  const Real _gamma;

  const MaterialProperty<Real> & _specific_heat;
  const MaterialProperty<Real> & _density;

  const Real _beta_av;
  const MaterialProperty<RankTwoTensor> & _stress_artificial;

};
























//###########################################################################################################
