#include "TimeDerivative.h"

// Forward Declarations
class dY4_dt;

//this kernel computes the term $\frac{\partial (rho Y_1)}{\partial t}$
//term for the species conservation equation

class dY4_dt : public TimeDerivative
{
public:
  dY4_dt(const InputParameters & parameters);
  static InputParameters validParams();
  
protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  //virtual Real computeQpOffDiagonalJacobian(unsigned int jvar);

private:
  std::string _base_name;
};

