#include "dY4_dt.h"

registerMooseObject("beaverApp", dY4_dt);

InputParameters
dY4_dt::validParams()
{
  InputParameters params = TimeDerivative::validParams();
  params.addClassDescription("Y1 species time derivate term of the conservation equation"
                              "this kernel computes only the time derivative contribution"
                              "to the residual");
  return params;
}

dY4_dt::dY4_dt(const InputParameters & parameters)
  : TimeDerivative(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : "")
{}

Real
dY4_dt::computeQpResidual()
{
  //compute the weak form contribution of the term $\frac{\partial (rho Y_1)}{\partial t}$
  //the chain rule yields -> d(rho Y1) / dt = rho dY1/dt + Y1 drho/dt
  Real res;
  res = _u_dot[_qp] * _test[_i][_qp]; //rho dY_dt
  return res;
}

Real dY4_dt::computeQpJacobian()
{
  //compute the Jacobian term dR_dY where R is the residual of this kernel
  Real Jac = 0.0;
  Jac += _du_dot_du[_qp] * _phi[_j][_qp] * _test[_i][_qp];
  return Jac;
}

