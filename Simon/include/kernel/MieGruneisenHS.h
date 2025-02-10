#include "HeatSource.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include <vector>

//forward declarte the class object to be acted upon

class MieGruneisenHS : public HeatSource
{
public:
  MieGruneisenHS(const InputParameters & parameters);
  static InputParameters validParams();
  
protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagonalJacobian(unsigned int jvar);

private:
  const Real _Gamma;
  const Real _T_ref;
  const MaterialProperty<Real> &_rho;
  const MaterialProperty<Real> &_Cv;
  const VariableValue &_T;
  const MaterialProperty<RankTwoTensor> &_mechanical_strain;
  const MaterialProperty<RankTwoTensor> &_mechanical_strain_old;
  const Real _C0;
  const Real _C1;

  //variable element length computation
  const Elem * const &_current_elem;
  const Real _beta_av;
  const MaterialProperty<RankFourTensor> &_elasticity_tensor;
  const Real _viscosity_type;
  const Real _Le;
};
