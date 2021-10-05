//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "HuhMultiACBulkF.h"

registerMooseObject("PhaseFieldApp", HuhMultiACBulkF);

InputParameters
HuhMultiACBulkF::validParams()
{
  InputParameters params = HuhMultiACBulkBase::validParams();
  params.addClassDescription("KKS model kernel (part 1 of 2) for the Bulk Allen-Cahn. This "
                             "includes all terms NOT dependent on chemical potential.");
  params.addRequiredParam<Real>("wik_names", "Double well height parameters");
  return params;
}

HuhMultiACBulkF::HuhMultiACBulkF(const InputParameters & parameters)
  : HuhMultiACBulkBase(parameters),
    _wik_names(getParam<std::vector<MaterialPropertyName>>("wik_names")),
    _prop_wik(_num_k))
{
  for (unsigned int n = 0; n<_num_k; ++n)
  {
    _prop_wik[n] = &getMaterialPropertyByName<Real>(_wik_names[n]);
  }
}

Real
HuhMultiACBulkF::computeDFDOP(PFFunctionType type)
{
  Real sum = 0.0;

  switch (type)
  {
    case Residual:
      // for (unsigned int n = 0; n < _num_j; ++n)
      //   sum += (*_prop_dhjdetai[n])[_qp] * (*_prop_Fj[n])[_qp];

      // return sum + _wi * _prop_dgi[_qp];
      
      // Get the value of Fi
      Real Fi = 
      for (unsigned int n = 0; n<_num_k; ++n)
      {
        // Residual += si*sk*Mik*(Fi-Fk)
        sum += _prop_si[_qp] * (*_prop_sk[n])[_qp]* (*prop_Mik[n])[_qp]* ((*prop_Fi[n])[_qp] - (*prop_Fk[n])[_qp]));

        // Residual += si*sk*Mik*wik*(phi_k-phi_i) 
        sum += _prop_si[_qp] * (*_prop_sk[n])[_qp]* (*prop_Mik[n])[_qp] * 
                  (*prop_wik[n])[_qp] * ((*_prop_hk[n])[_qp] - _u[_qp]);
      }
      return sum/_num_phases[_qp];


    case Jacobian:
      // For when this kernel is used in the Lagrange multiplier equation
      // In that case the Lagrange multiplier is the nonlinear variable
      // if (_etai_var != _var.number())
      //   return 0.0;

      // For when eta_i is the nonlinear variable
      // for (unsigned int n = 0; n < _num_j; ++n)
      //   sum += (*_prop_d2hjdetai2[n])[_qp] * (*_prop_Fj[n])[_qp];
      for (unsigned int n = 0; n < _num_k; ++n)
      {
        // Jacobian += -si*sk*Mik*wik*1
        sum += -_prop_si[_qp] * (*_prop_sk[n])[_qp]* (*_prop_Mik[n])[_qp]*(*_prop_wik[n])[_qp];
      }

      return _phi[_j][_qp]/ _num_phases *  sum;
  }

  mooseError("Invalid type passed in");
}

Real
HuhMultiACBulkF::computeQpOffDiagJacobian(unsigned int jvar)
{
  // // get the coupled variable jvar is referring to
  // const unsigned int cvar = mapJvarToCvar(jvar);

  // // first get dependence of mobility _L on other variables using parent class
  // // member function
  // Real res = ACBulk<Real>::computeQpOffDiagJacobian(jvar);

  // // Then add dependence of KKSMultiACBulkF on other variables
  // Real sum = 0.0;
  // for (unsigned int n = 0; n < _num_k; ++n)
  //   sum += (*_prop_d2hjdetaidarg[n][cvar])[_qp] * (*_prop_Fj[n])[_qp] +
  //          (*_prop_dhjdetai[n])[_qp] * (*_prop_dFjdarg[n][cvar])[_qp];

  // // Handle the case when this kernel is used in the Lagrange multiplier equation
  // // In this case the second derivative of the barrier function contributes
  // // to the off-diagonal Jacobian
  // if (jvar == _etai_var)
  //   sum += _wi * _prop_d2gi[_qp];

  // res += _L[_qp] * sum * _phi[_j][_qp] * _test[_i][_qp];
  Real res = 0; // bruh i don't know how to do this

  return res;
}
