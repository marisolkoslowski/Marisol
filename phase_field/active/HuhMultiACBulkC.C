//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HuhMultiACBulkC.h"

registerMooseObject("PhaseFieldApp", HuhMultiACBulkC);

InputParameters
HuhMultiACBulkC::validParams()
{
  InputParameters params = HuhMultiACBulkBase::validParams();
  params.addClassDescription("Multi-phase KKS model kernel (part 2 of 2) for the Bulk Allen-Cahn. "
                             "This includes all terms dependent on chemical potential.");
  params.addRequiredCoupledVar(
      "ci_name", "Phase concentration ci corresponding to the variable.");
  params.addRequiredCoupledVar(    
      "ck_names", "Array of phase concentrations ck. Place in same order as Fk_names!");

  return params;
}

HuhMultiACBulkC::HuhMultiACBulkC(const InputParameters & parameters)
  : KKSMultiACBulkBase(parameters),
    _ci(coupledValue("ci_name"))
    _ci_var(coupled("ci_name",0)),
    _c1_name(getVar("ck_name", 0)
                 ->name()), // Can use any dFk/dck since they are equal so pick first ck in the list
    _cks(coupledValues("ck_names")),
    _cks_var(coupledIndices("ck_names")),
    _prop_dF1dc1(getMaterialPropertyDerivative<Real>(_Fk_names[0],
                                                     _c1_name)) // Use first Fj in list for dFj/dcj
{
  if (_num_k != coupledComponents("ck_names"))
    paramError("ck_names", "Need to pass in as many ck_names as Fk_names");
}

Real
HuhMultiACBulkC::computeDFDOP(PFFunctionType type)
{
  Real sum = 0.0;

  switch (type)
  {
    case Residual:
      for (unsigned int n = 0; n < _num_k; ++n){
        // sum += (*_prop_dhdetai[n])[_qp] * (*_cks[n])[_qp];

        // Residual = dF1dc1 * si*sk*Mik*(-ci+ck)
        sum += _si[_qp] * (*_sk[n])[_qp] * (*_Mik[n])[_qp] * ( -_ci[_qp] + (*_cks[n])[_qp]); 
      }
      return  _prop_dF1dc1[_qp] * sum/_num_phases[_qp];

    case Jacobian:
      // For when this kernel is used in the Lagrange multiplier equation
      // In that case the Lagrange multiplier is the nonlinear variable
      // if (_etai_var != _var.number())
      //   return 0.0;

      // For when eta_i is the nonlinear variable
      for (unsigned int n = 0; n < _num_j; ++n){
        // sum += (*_prop_d2hjdetai2[n])[_qp] * (*_cjs[n])[_qp];
        
        // There is no Jacobian for this term.
        sum += 0;
      }
      return sum/_num_phases[_qp];
  }

  mooseError("Invalid type passed in");
}

Real
HuhMultiACBulkC::computeQpOffDiagJacobian(unsigned int jvar)
{
  // // first get dependence of mobility _L on other variables using parent class
  // // member function
  // Real res = ACBulk<Real>::computeQpOffDiagJacobian(jvar);

  // Real sum = 0.0;
  // // Then add dependence of KKSACBulkC on other variables
  // // Treat cj variables specially, as they appear in the residual
  // if (jvar == _cjs_var[0])
  // {
  //   for (unsigned int n = 0; n < _num_j; ++n)
  //     sum += (*_prop_dhjdetai[n])[_qp] * (*_cjs[n])[_qp];

  //   res -= _L[_qp] * (sum * _prop_d2F1dc12[_qp] + _prop_dF1dc1[_qp] * (*_prop_dhjdetai[0])[_qp]) *
  //          _phi[_j][_qp] * _test[_i][_qp];
  //   return res;
  // }

  // for (unsigned int i = 1; i < _num_j; ++i)
  // {
  //   if (jvar == _cjs_var[i])
  //   {
  //     res -=
  //         _L[_qp] * _prop_dF1dc1[_qp] * (*_prop_dhjdetai[i])[_qp] * _phi[_j][_qp] * _test[_i][_qp];
  //     return res;
  //   }
  // }

  // //  for all other vars get the coupled variable jvar is referring to
  // const unsigned int cvar = mapJvarToCvar(jvar);

  // for (unsigned int n = 0; n < _num_j; ++n)
  //   sum += _prop_dF1dc1[_qp] * (*_prop_d2hjdetaidarg[n][cvar])[_qp] * (*_cjs[n])[_qp] +
  //          (*_prop_d2F1dc1darg[cvar])[_qp] * (*_prop_dhjdetai[n])[_qp] * (*_cjs[n])[_qp];

  // res -= _L[_qp] * sum * _phi[_j][_qp] * _test[_i][_qp];

  Real res = 0;

  return res;
}
