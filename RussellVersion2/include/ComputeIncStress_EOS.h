/// Calculates stress for anisortopic crack propagation
/// Includes artificial viscosity and Mie Gruneisen Equation of State

#pragma once

#include "ComputeStressBase.h"
#include "MathUtils.h"
#include "RankTwoTensor.h"
#include "PiecewiseBilinear.h"
#include "RankFourTensor.h"

class ComputeIncStress_EOS : public ComputeStressBase//, public Material
{
public:
  static InputParameters validParams();
  ComputeIncStress_EOS(const InputParameters & parameters);

protected:
  virtual void          computeQpStress() override;
  virtual void initQpStatefulProperties() override;
  void  EOSUpdate();
  void  ArtificialViscosity();

 //Moose Properties/Variables
  //Name Finder
   const std::string _elasticity_tensor_name;
  //Stress Update:
   const MaterialProperty<RankFourTensor> & _elasticity_tensor;
   MaterialProperty       <RankTwoTensor> & _rotation_total;
   const MaterialProperty <RankTwoTensor> & _rotation_total_old;
   const MaterialProperty <RankTwoTensor> & _strain_increment;
   const MaterialProperty <RankTwoTensor> & _rotation_increment;
   const MaterialProperty <RankTwoTensor> & _stress_old;
   const MaterialProperty <RankTwoTensor> & _elastic_strain_old;
   const MaterialProperty <RankTwoTensor> & _mechanical_strain_old;

  //Specific Volume
   const MaterialProperty          <Real> & _V;
   const MaterialProperty          <Real> & _mu;

  //Pressure EOS
   const MaterialPropertyName               _pressure_eos_name;
   const MaterialProperty          <Real> & _pressure_eos;

  //Artificial Viscosity
   const MaterialProperty          <Real> & _density;
   MaterialProperty       <RankTwoTensor> & _stress_artificial;
   MaterialProperty                <Real> & _Le;
   const Real _C0;
   const Real _C1;
   const Elem * const & _current_elem;

  //Internal Class Update
  RankTwoTensor stress_dev; //Deviatoric Stress
  RankTwoTensor stress_hyd; //Hydrostatic Stress
  RankTwoTensor I2 = RankTwoTensor::Identity();

  //Debugging
  usingTensorIndices(i_, j_, k_, l_);
};












//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
