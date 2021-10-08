#
# KKS simple example in the split form
#

[Debug]
  show_var_residual_norms = true
[]


[Mesh]
  type = GeneratedMesh
  dim = 2
  elem_type = QUAD4
  nx = 2000
  ny = 2
  nz = 0
  xmin = -10
  xmax = 10
  ymin = 0
  ymax = 0.4
  zmin = 0
  zmax = 0
[]

[AuxVariables]
  [./Energy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./eta1_dot]
  [../]
  [./eta2_dot]
  [../]
  [./c_dot]
  [../]
  [./diffusion_c1]
  [../]
  [./diffusion_c2]
  [../]
[]

[Variables]
  # order parameter
  [./eta1]
    order = FIRST
    family = LAGRANGE
    scaling = 1e4
  [../]

  [./eta2]
    order = FIRST
    family = LAGRANGE
    scaling = 1e4
  [../]

  # hydrogen concentration
  [./c]
    order = FIRST
    family = LAGRANGE
  [../]

  # Cu6Sn5 phase solute concentration
  [./c2]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.4105
    scaling = 1e10
  [../]

  # Sn phase solute concentration
  [./c1]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.9769
    scaling = 1e5
  [../]

[]



[Functions]
  [./ic_func_eta1]
    type = ParsedFunction
    value = 'if(x<=0, 1, 0)'
  [../]
  [./ic_func_eta2]
    type = ParsedFunction
    value = 'if(x<=0, 0, 1)'
  [../]

  [./ic_func_c]
    type = ParsedFunction
    value = 'if(x<=0, 0.9769, 0.4105)'
  [../]
[]

[ICs]
  [./eta1]
    variable = eta1
    type = FunctionIC
    function = ic_func_eta1
  [../]
  [./eta2]
    variable = eta2
    type = FunctionIC
    function = ic_func_eta2
  [../]
  [./c]
    variable = c
    type = FunctionIC
    function = ic_func_c
  [../]
[]

[BCs]
  [./Periodic]
    [./all]
      auto_direction = 'y'
    [../]
  [../]

[]

[Materials]

  [./consts]
      type = GenericConstantMaterial
      prop_names  = ' Vm         sigma12  sigma23  R      T    w        scale '
      prop_values = ' 16.29e-6   0.1      0.3     8.314  523   0.1e-6   0.1e-6 '
  [../]

  [./FreeEnergyConstants]
    type = GenericConstantMaterial
    prop_names = ' GAlphaCu   GSERSn    GLCu      GLSn       LL0         LL1        LL2'
    prop_values = '-1.9073e4  -2.8716e4  -1.1083e4 -2.8963e4  -1.0487e4  -1.8198e4  -10528.3'
  [../]

  [./FreeEnergyConstants2]
    type = GenericConstantMaterial
    prop_names = 'GAlphaSn   Lalpha0     Lalpha1 '
    prop_values = '-2.7281e4 -1.1448e4   -1.1694e4 '
  [../]


  # Free energy of the liquid Sn
  [./f1]
    type = DerivativeParsedMaterial
    f_name = f1
    args = 'c1'
    #function = '((1-c1)*GLCu + c1*GLSn +R*T*((1-c1)*log(1-c1)+c1*log(c1))
    #          +c1*(1-c1)*(LL0+LL1*(1-2*c1)+LL2*(1-4*c1+4*c1^2)))/Vm*scale^3' 
    function = '(1.5e4*(0.9769-c1)^2/Vm-1.8e9)*scale^2'
    material_property_names = 'GLCu GLSn LL0 LL1 LL2 R T Vm scale'
    outputs = exodus
  [../]


  # Free energy of the Cu6Sn5
  [./f2]
    type = DerivativeParsedMaterial
    f_name = f2
    args = 'c2'
    function = '(2.0e5*(c2-0.435)^2 + 0.545*GAlphaCu + 0.455*GSERSn-6689.5-0.1589*T)/Vm*scale^2'
    material_property_names = 'GAlphaCu GSERSn T Vm scale'
    outputs = exodus
  [../]


#  [./f3]
#    type = DerivativeParsedMaterial
#    f_name = f3
#    args = 'c3'
#    function = '((1-c3)*GAlphaCu + c3*GAlphaSn + R*T*((1-c3)*log(1-c3)+c3*log(c3))
#              + c3*(1-c3)*(Lalpha0+Lalpha1*(1-2*c3)))/Vm*scale^2'
#    material_property_names = 'GAlphaCu GAlphaSn Lalpha0 Lalpha1 R T Vm scale'
#    outputs = exodus
#  [../]


  [./h1]
    type = DerivativeParsedMaterial
    f_name = h1
    function = 'eta1'
    args = 'eta1'
    derivative_order = 2
  [../]
  [./h2]
    type = DerivativeParsedMaterial
    f_name = h2
    function = 'eta2'
    args = 'eta2'
    derivative_order = 2
  [../]



  [./diffusionConsts]
      # Dimensionless
      type = GenericConstantMaterial
      prop_names =  'DSn DCu6Sn5 DCu      D12      D23'
      prop_values = '1   2.0e-5  1.0e-6   2.0e-1   2.0e-3 '
  [../]



  [./Dfunc]
      type = ParsedMaterial
      f_name = Dfunc
      function = 'if(eta1>0.8, DSn,
                      if(eta2>0.8, DCu6Sn5,
                        if( (eta2>0.2) & (eta1>0.2), D12, D12))))'

      #function = 1
      args = 'eta1 eta2'
      material_property_names = 'DSn D12 DCu6Sn5 D23 DCu'
      derivative_order = 2
      outputs = exodus
  [../]


  [./Dh1]
    type = DerivativeParsedMaterial
    material_property_names = 'Dfunc(eta1,eta2) h1 w'
    function = 'Dfunc*h1'
    #function = 1
    args = 'eta1 eta2'
    f_name = Dh1
    outputs = exodus
  [../]
  [./Dh2]
    type = ParsedMaterial
    material_property_names = 'Dfunc(eta1,eta2) h2 w'
    function = 'Dfunc*h2'
    #function = 1
    args = 'eta1 eta2'
    f_name = Dh2
    outputs = exodus
  [../]


  [./mobilityConsts]
    type = GenericConstantMaterial
    prop_names = 'M12 M23    Mgb'
    prop_values = '1e6 7.0e4 7.0e4'
  [../]
  

  [./eps12_sq]
    type = ParsedMaterial
    f_name = eps12_sq
    function = '(sigma12*w/3.14)*scale'
    material_property_names = 'sigma12 w scale'
    outputs = exodus
  [../]

  [./w12]
    type = ParsedMaterial
    f_name = w12
    function = '4*sigma12/w *scale^2'
    material_property_names = 'sigma12 w scale'
  [../]

  ####
  # Define step functions
  ####

  [./s1]  
    type = ParsedMaterial
    f_name = s1
    function = 'if(eta1>1e-3, 1, 0)'
    args = 'eta1'
  [../]

  [./s2]  
    type = ParsedMaterial
    f_name = s2
    function = 'if(eta2>1e-3, 1,0)'
    args = 'eta2'
  [../]



  ## Number of phases

  [./num_phases]
    type = ParsedMaterial
    f_name = num_phases
    function = '(if(eta1>1e-3, 1, 0) + if(eta2>1e-3, 1, 0))/2'
    args = 'eta1 eta2' 
  [../]





[]

[Kernels]
  # enforce c = (1-h(eta))*cl + h(eta)*cs
  [./phaseconcentration]
    type = KKSMultiPhaseConcentration
    variable = c1
    cj = 'c1 c2'
    hj_names = 'h1 h2'
    etas = 'eta1 eta2'
    c = c
  [../]

  # enforce pointwise equality of chemical potentials
  [./ChemPotSolute]
    type = KKSPhaseChemicalPotential
    variable = c2
    cb       = c1
    fa_name  = f2
    fb_name  = f1
  [../]


  #
  # Cahn-Hilliard Equation
  #
  #Kernels for diffusion equation
  [./diff_time]
    type = TimeDerivative
    variable = c
    save_in = c_dot
  [../]
  [./diff_c1]
    type = MatDiffusion
    variable = c
    diffusivity = Dh1
    v = c1
    args = 'eta1 eta2'
    save_in = diffusion_c1
  [../]
  [./diff_c2]
    type = MatDiffusion
    variable = c
    diffusivity = Dh2
    v = c2
    args = 'eta1 eta2'
    save_in = diffusion_c2
  [../]




  #
  # Allen-Cahn Equation
  #
  [./deta1dt]
    type = TimeDerivative
    variable = eta1
    save_in = eta1_dot
  [../]
  [./ACBulkF1]
    type = HuhMultACBulkF
    variable  = eta1
    Fi_names  = 'f1'
    si_name = 's1'
    Fk_names = 'f2'
    sk_names = 's2'
    hk_names = 'h2'
    Mik_names = 'M12'
    wik_names = 'w12'
    num_phases = 'num_phases'
  [../]
  [./ACBulkC1]
    type = KKSMultiACBulkC
    variable  = eta1
    Fj_names  = 'f1 f2 f3'
    hj_names  = 'h1 h2 h3'
    cj_names  = 'c1 c2 c3'
    eta_i     = eta1
    args      = 'eta2 eta3'
    mob_name = Madd1
  [../]
  [./ACInterface1_1]
    type = ACInterface
    variable = eta1 
    args = 'eta2 eta3' # inputs of mobility Madd
    mob_name = Madd1
    kappa_name = eps_sq
  [../]

  [./ACBulkF1_2]
    type = KKSMultiACBulkF
    variable  = eta1
    Fj_names  = 'f1 f2 f3'
     hj_names  = 'h1 h2 h3'
     gi_name   = g2
     eta_i     = eta2
     wi        = 4e-8
     #wi        = 1
     args      = 'c1 c2 c3 eta1 eta3'
     mob_name = Msubtract1
   [../]
   [./ACBulkC1_2]
     type = KKSMultiACBulkC
     variable  = eta1
     Fj_names  = 'f1 f2 f3'
     hj_names  = 'h1 h2 h3'
     cj_names  = 'c1 c2 c3'
     eta_i     = eta2
     args      = 'eta3'
     mob_name = Msubtract1
   [../]
   [./ACInterface1_2]
     type = SimpleCoupledACInterface
     variable = eta1
     v = eta2
     mob_name = Msubtract1
     kappa_name = eps_sq
   [../]



  [./deta2dt]
    type = TimeDerivative
    variable = eta2
    save_in = eta2_dot
  [../]
   [./ACBulkF2_2]
     type = KKSMultiACBulkF
     variable  = eta2
     Fj_names  = 'f1 f2 f3'
     hj_names  = 'h1 h2 h3'
     gi_name   = g2
     eta_i     = eta2
     wi        = 4e-8
     #wi        = 1
     args      = 'c1 c2 c3 eta1 eta3'
     mob_name = Madd2
   [../]
   [./ACBulkC2_2]
     type = KKSMultiACBulkC
     variable  = eta2
     Fj_names  = 'f1 f2 f3'
     hj_names  = 'h1 h2 h3'
     cj_names  = 'c1 c2 c3'
     eta_i     = eta2
     args      = 'eta1 eta3'
     mob_name = Madd2
   [../]
   [./ACInterface2_2]
     type = ACInterface
     variable = eta2
     args = 'eta1 eta3'
     mob_name = Madd2
     kappa_name = eps_sq
   [../]

   [./ACBulkF2_1]
     type = KKSMultiACBulkF
     variable  = eta2
     Fj_names  = 'f1 f2 f3'
     hj_names  = 'h1 h2 h3'
     gi_name   = g1
     eta_i     = eta1
     wi        = 4e-8
     #wi        = 1
     args      = 'c1 c2 c3 eta2 eta3'
     mob_name = Msubtract1
   [../]
   [./ACBulkC2_1]
     type = KKSMultiACBulkC
     variable  = eta2
     Fj_names  = 'f1 f2 f3'
     hj_names  = 'h1 h2 h3'
     cj_names  = 'c1 c2 c3'
     eta_i     = eta1
     args      = 'eta2 eta3'
     mob_name = Msubtract1 
   [../]
   [./ACInterface2_1]
     type = SimpleCoupledACInterface
     variable = eta2
     v = eta1
     mob_name = Msubtract1
     kappa_name = eps_sq
    [../]

   [./ACBulkF2_3]
     type = KKSMultiACBulkF
     variable  = eta2
     Fj_names  = 'f1 f2 f3'
     hj_names  = 'h1 h2 h3'
     gi_name   = g3
     eta_i     = eta3
     wi        = 4e-8
     #wi        = 1
     args      = 'c1 c2 c3 eta1'
     mob_name = Msubtract3
   [../]
   [./ACBulkC2_3]
     type = KKSMultiACBulkC
     variable  = eta2
     Fj_names  = 'f1 f2 f3'
     hj_names  = 'h1 h2 h3'
     cj_names  = 'c1 c2 c3'
     eta_i     = eta3
     args      = 'eta1 eta2'
     mob_name = Msubtract3
   [../]
   [./ACInterface2_3]
     type = SimpleCoupledACInterface
     variable = eta2
     v = eta3
     mob_name = Msubtract3
     kappa_name = eps_sq
   [../]


  [./deta3dt]
    type = TimeDerivative
    variable = eta3
    save_in = eta3_dot
  [../]
  [./ACBulkF3_3]
    type = KKSMultiACBulkF
    variable  = eta3
    Fj_names  = 'f1 f2 f3'
    hj_names  = 'h1 h2 h3'
    gi_name   = g3
    eta_i     = eta3
    wi        = 4e-8
    #wi        = 1
    args      = 'c1 c2 c3 eta1 eta2'
    mob_name = Madd3
  [../]
  [./ACBulkC3_3]
    type = KKSMultiACBulkC
    variable  = eta3
    Fj_names  = 'f1 f2 f3'
    hj_names  = 'h1 h2 h3'
    cj_names  = 'c1 c2 c3'
    eta_i     = eta3
    args      = 'eta1 eta2'
    mob_name = Madd3
  [../]
  [./ACInterface3_3]
    type = ACInterface
    variable = eta3
    args = 'eta1 eta2'
    mob_name = Madd3
    kappa_name = eps_sq
  [../]

  [./ACBulkF3_2]
    type = KKSMultiACBulkF
    variable  = eta3
    Fj_names  = 'f1 f2 f3'
    hj_names  = 'h1 h2 h3'
    gi_name   = g2
    eta_i     = eta2
    wi        = 4e-8
    #wi        = 1
    args      = 'c1 c2 c3 eta1'
    mob_name = Msubtract3
  [../]
  [./ACBulkC3_2]
    type = KKSMultiACBulkC
    variable  = eta3
    Fj_names  = 'f1 f2 f3'
    hj_names  = 'h1 h2 h3'
    cj_names  = 'c1 c2 c3'
    eta_i     = eta2
    args      = 'eta1'
    mob_name = Msubtract3
  [../]
  [./ACInterface3_2]
    type = SimpleCoupledACInterface
    variable = eta3
    v = eta2
    mob_name = Msubtract3
    kappa_name = eps_sq
  [../]



[]

[AuxKernels]
  [./Energy_total]
    type = KKSMultiFreeEnergy
    Fj_names = 'f1 f2 f3'
    hj_names = 'h1 h2 h3'
    gj_names = 'g1 g2 g3'
    variable = Energy
    w        = 4e-8
    #w        = 1
    interfacial_vars =  'eta1  eta2     eta3'
    kappa_names =       'eps_sq eps_sq  eps_sq'
  [../]
[]
[Executioner]
  type = Transient
  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type'
  petsc_options_value = 'asm      ilu          nonzero'

  l_max_its = 100
  nl_max_its = 100
  nl_abs_tol = 5e-4
  nl_rel_tol = 1e-6


  num_steps = 2000000
  dt = 1e-4

  #num_steps = 30000
  #dt = 2e-3
[]

#
# Precondition using handcoded off-diagonal terms
#
[Preconditioning]
  [./full]
    type = SMP
    full = true
  [../]
[]



[Outputs]
  #interval = 600
  interval = 1
  exodus = true
  checkpoint = true
[]
