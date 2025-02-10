/// Calculates the change in phase 1 
#pragma once
#include "AuxKernel.h"

class AuxIgnitionAndGrowth : public AuxKernel {
public:
  static InputParameters validParams();
  AuxIgnitionAndGrowth(const InputParameters & parameters);
protected:
  virtual Real computeValue() override;
private:  
  const Real _F_igMin;
  const Real _F_igMax;
  const Real  _I_chem;
  const Real  _a_chem;
  const Real  _b_chem;
  const Real  _x_chem;

  const Real _F_G1Min;
  const Real _F_G1Max;
  const Real _G1_chem;
  const Real _c1_chem;
  const Real _d1_chem;
  const Real _y1_chem;

  const Real _F_G2Min;
  const Real _F_G2Max;
  const Real _G2_chem;
  const Real _c2_chem;
  const Real _d2_chem;
  const Real _y2_chem;

  const MaterialProperty<Real> &            _mu;
  const MaterialProperty<Real> &  _pressure_eos;
  const VariableValue &  _u_old;
};


















