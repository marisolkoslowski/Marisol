# Units:
#+00:  K:Temperature
#-09:  s:Time 
#-09:  m:Length
#-18: kg:Mass 
#+00:  A:Current 
#+00:mol:Number 
#+00: cd:Luminosity

#https://mooseframework.inl.gov/modules/heat_transfer/tutorials/introduction/therm_step03.html
################################################################################
[GlobalParams] displacements = 'disp_x disp_y' []
[Mesh]
	#type = GeneratedMesh
	#dim  =  2
	#nx   =  21
	#ny   =  21
	#xmax =  2000
	#xmin = -2000
	#ymax =  2000
	#ymin = -2000
  [./ccmg]
    type = ConcentricCircleMeshGenerator
    num_sectors = 6
    radii = '600 3000'
    rings = '2   4'
    has_outer_square = false
    preserve_volumes = false
  []
[]

[Adaptivity]
  marker = errorfrac #specify the reference marker
  cycles_per_step = 0 #number of iterations to run the refinement scheme
  max_h_level = 1 #2 #amount of times a single element can be refined/coarsened
  initial_steps = 0 #8
  recompute_markers_during_cycles = false#

  #define indicators 
  [./Indicators]
    #indicator to be extracted from the marker at each element
    [./error]
      type = GradientJumpIndicator #ValueJumpIndicator #GradientJumpIndicator
      variable = temperature
      outputs = Exodus
    [../]
  [../]
  #Create a marker that determines which elements to refine/coarsen based on error estimates
  [./Markers]
    [./errorfrac]
      type = ErrorFractionMarker
      indicator = error
      refine  = 0.25
      coarsen = 0.25
      clear_extremes = true
      outputs = Exodus
    [../]
    #[./errorfrac]
    #  type = ErrorToleranceMarker
    #  indicator = error
    #  coarsen = 0.5 # 800
    #  refine  = 0.5 #1000
    #  outputs = Exodus
    #[../]
  [../]
[]
###############################################################################
[Variables]
  [temperature] order = FIRST family = LAGRANGE []
  [Y_RDX]          order = FIRST family = LAGRANGE []
  [Y_INT]          order = FIRST family = LAGRANGE []
  #[c     ]     order = FIRST family = LAGRANGE [] #Damage 0-1
  [disp_x]     order = FIRST family = LAGRANGE []
  [disp_y]     order = FIRST family = LAGRANGE []
[]
[Kernels]
  #Temperature
    [dTdt]    type  = HeatConductionTimeDerivative  variable = temperature  []
    [T_cond]  type  = HeatConduction                variable = temperature  []
    [T_RDX]
      type  = MatCoupledForce
      variable = temperature
      v = Y_RDX #Y_RDX* w1 is Y_RDX
      material_properties = Heat_RDX
    []
    [T_INT]
      type  = MatCoupledForce
      variable = temperature
      v = Y_INT #Y_INT*WINt is Y_PRD_dot
      material_properties = Heat_INT
    []
    [Heat_FiniteStress_EOS_ART]
      type = Heat_EOS
      dPress_dT_name = pressure_dT_eos
      variable = temperature
      beta_av = 0.4
    []

  #Massfraction 1
    [Y_RDX_Consv]   type = TimeDerivative  variable = Y_RDX  []
    [Y_RDX_Loss]   
      type = MatReaction
      variable = Y_RDX
      reaction_rate = w_RDX_loss #w1 needs to be -
      v = Y_RDX
    []
  
  #Massfraction 2
    [Y_INT_Consv]   type = TimeDerivative  variable = Y_INT  []
    [Y_INT_Make]   
      type = MatReaction
      variable = Y_INT
      reaction_rate = w_INT_make #w1 needs to be +
      v = Y_RDX
    []
    [Y_INT_Loss] 
      type = MatReaction
      variable = Y_INT
      reaction_rate = w_INT_loss #w2 needs to be -
      v = Y_INT
    []

  ##Mechanics##
    [./inertia_x]
      type = InertialForce
      variable     = disp_x
      velocity     = velo_x
      acceleration = acce_x
      beta  = 0.3025 ###from dandekar 2019 sec 2.1
      gamma = 0.6000 ###from dandekar 2019 sec 2.1
      use_displaced_mesh = true
    [../]
    [./inertia_y]
      type = InertialForce
      variable     = disp_y
      velocity     = velo_y
      acceleration = acce_y
      beta  = 0.3025 ###from dandekar 2019 sec 2.1
      gamma = 0.6000 ###from dandekar 2019 sec 2.1
      use_displaced_mesh = true
    [../]

    [TM_all0]
      type = StressDivergenceTensors
      block = ''
      component = 0
      displacements = 'disp_x disp_y '
      use_displaced_mesh = true
      variable = disp_x
    []
    [TM_all1]
      type = StressDivergenceTensors
      block = ''
      component = 1
      displacements = 'disp_x disp_y '
      use_displaced_mesh = true
      variable = disp_y
    []
[]
################################################################################
[AuxVariables]
  [Y_PRD] order = FIRST family = LAGRANGE  []

 #Implementation
  [velo_x]          order = FIRST    family = LAGRANGE []
  [acce_x]          order = FIRST    family = LAGRANGE []
  [velo_y]          order = FIRST    family = LAGRANGE []
  [acce_y]          order = FIRST    family = LAGRANGE []

 #Reporting
  [w_RDX_loss]      order = CONSTANT family = MONOMIAL []
  [w_INT_make]      order = CONSTANT family = MONOMIAL []
  [w_INT_loss]      order = CONSTANT family = MONOMIAL []
  [Q_RDX]           order = CONSTANT family = MONOMIAL []
  [Q_INT]           order = CONSTANT family = MONOMIAL []
  [Heat_RDX]        order = CONSTANT family = MONOMIAL []
  [Heat_INT]        order = CONSTANT family = MONOMIAL []


  [ElemLength]               order = CONSTANT family = MONOMIAL []
  [V]               order = CONSTANT family = MONOMIAL []
  [pressure_eos]    order = CONSTANT family = MONOMIAL []
  #[pressure_Y_RDX] order = CONSTANT family = MONOMIAL []

  [stress_xx]       order = CONSTANT family = MONOMIAL []
  [stress_yy]       order = CONSTANT family = MONOMIAL []
  [stress_xy]       order = CONSTANT family = MONOMIAL []
  [stress_vol]      order = CONSTANT family = MONOMIAL []
  [stress_vonmises] order = CONSTANT family = MONOMIAL []

  [strain_mech_xx]  order = CONSTANT family = MONOMIAL []
  [strain_mech_yy]  order = CONSTANT family = MONOMIAL []
  [strain_mech_xy]  order = CONSTANT family = MONOMIAL []
  [strain_mech_vol] order = CONSTANT family = MONOMIAL []

  [report_output]   order = CONSTANT family = MONOMIAL []
[]
[AuxKernels]
  [Y_PRD]
    type = ParsedAux
    variable = Y_PRD
    coupled_variables = 'Y_RDX Y_INT'
    expression = '1 - Y_RDX - Y_INT'
  []

 #Implementation
  [velo_x]          type = NewmarkVelAux     variable = velo_x  acceleration = acce_x                     gamma = 0.6    []
  [velo_y]          type = NewmarkVelAux     variable = velo_y  acceleration = acce_y                     gamma = 0.6    []
  [acce_x]          type = NewmarkAccelAux   variable = acce_x  displacement = disp_x  velocity = velo_x  beta = 0.3025  []
  [acce_y]          type = NewmarkAccelAux   variable = acce_y  displacement = disp_y  velocity = velo_y  beta = 0.3025  []

 #Reporting 
  [w_RDX_loss]      type = MaterialRealAux   variable = w_RDX_loss       property = w_RDX_loss           []
  [w_INT_make]      type = MaterialRealAux   variable = w_INT_make       property = w_INT_make            []
  [w_INT_loss]      type = MaterialRealAux   variable = w_INT_loss       property = w_INT_loss            []
  [Q_RDX]           type = MaterialRealAux   variable = Q_RDX            property = Q_RDX            []
  [Q_INT]           type = MaterialRealAux   variable = Q_INT            property = Q_INT            []
  [Heat_RDX]        type = MaterialRealAux   variable = Heat_RDX         property = Heat_RDX            []
  [Heat_INT]        type = MaterialRealAux   variable = Heat_INT         property = Heat_RDX            []

  [ElemLength]      type = MaterialRealAux   variable = ElemLength       property = Le            []
  [V]               type = MaterialRealAux   variable = V                property = V            []
  [pressure_eos]    type = MaterialRealAux   variable = pressure_eos     property = pressure_eos []
  #[pressure_Y_RDX ]type = MaterialRealAux   variable = pressure_Y_RDX   property = pressure_Y_RDX  []


  [stress_xx]       type = RankTwoAux        variable = stress_xx        rank_two_tensor = stress            index_j = 0    index_i = 0     []
  [stress_yy]       type = RankTwoAux        variable = stress_yy        rank_two_tensor = stress            index_j = 1    index_i = 1     []
  [stress_xy]       type = RankTwoAux        variable = stress_xy        rank_two_tensor = stress            index_j = 1    index_i = 0     []
  [stress_vol]      type = RankTwoScalarAux  variable = stress_vol       rank_two_tensor = stress            scalar_type = Hydrostatic      []
  [stress_vonmises] type = RankTwoScalarAux  variable = stress_vonmises  rank_two_tensor = stress            scalar_type = VonMisesStress   []

  [strain_mech_xx]  type = RankTwoAux        variable = strain_mech_xx   rank_two_tensor = mechanical_strain            index_j = 0    index_i = 0     []
  [strain_mech_yy]  type = RankTwoAux        variable = strain_mech_yy   rank_two_tensor = mechanical_strain            index_j = 1    index_i = 1     []
  [strain_mech_xy]  type = RankTwoAux        variable = strain_mech_xy   rank_two_tensor = mechanical_strain            index_j = 1    index_i = 0     []
  [strain_mech_vol] type = RankTwoScalarAux  variable = strain_mech_vol  rank_two_tensor = mechanical_strain            scalar_type = VolumetricStrain []

  #[report_output]    type = MaterialRealAux variable = report_output  property = report_output []

[]
################################################################################
#[AuxVariables]
#  [bounds_dummy]  order = FIRST  family = LAGRANGE  []
#[]
#[Bounds]
#  [c_upper_bound]  type = ConstantBounds          variable = bounds_dummy  bounded_variable = c  bound_type = upper  bound_value = 1  []
#  [c_lower_bound]  type = VariableOldValueBounds  variable = bounds_dummy  bounded_variable = c  bound_type = lower []
#[]
[ICs]
  [IC_Y_RDX]  type = ConstantIC  variable = Y_RDX  value = 1.0  []
  [IC_Y_INT]  type = ConstantIC  variable = Y_INT  value = 0.0  []

  [T_IC]
	type = SmoothCircleIC
	variable = temperature 
	invalue  =  301 #1200 #1920 
	outvalue =  300 
	radius   = '500'
	#int_width = 2.5
	x1 = 0
	y1 = 0
  []
[]
[BCs]
  #[Tleft]	type = NeumannBC  variable = temperature  value = 0  boundary = left  []
  #[Tright]	type = NeumannBC  variable = temperature  value = 0  boundary = right  []
  #[Ttop]	type = NeumannBC  variable = temperature  value = 0  boundary = top  []
  #[Tbottom]	type = NeumannBC  variable = temperature  value = 0  boundary = bottom  []
  #[./Pressure]
  #  [./Side2]
  #    boundary = 1
  #    displacements = 'disp_x disp_y'
  #    factor = 5.0
  #  [../]
  #[../]
  [fix_x]	type = DirichletBC  variable = disp_x  value = 0  boundary = 1  []
  [fix_y]	type = DirichletBC  variable = disp_y  value = 0  boundary = 1  []
[]
################################################################################
[Materials]
  #General Constants
    [General_Constant]
      type = GenericConstantMaterial
      prop_names = 'density'
      prop_values = 1.807e-06
    []
    [General_Heat]
      type = HeatConductionMaterial
      specific_heat        = 2.357e03
      thermal_conductivity = 3.610e-01
    []

  #Chemical Properties
    [3_Phase_Chemistry_Props]
      type = GenericConstantMaterial
      prop_names  = '    Z_RDX     Z_INT   E_RDX_R   E_INT_R      a_RDX     a_INT     b_RDX      b_INT T_ref_3phase'
      prop_values = '3.210E+04 1.181E+04 1.228E+04 1.183E+04 -5.783E-02 1.415E+01 1.068E-03 -1.305E-03      1.736e3'
    []

  #Elastic Tensor
    [plate_stiffne]
      type = ComputeElasticityTensorCP
      fill_method = symmetric21
      ##C_ijkl = 'C11 C12 C13 C14 C15 C16 
      ##              C22 C23 C24 C25 C26 
      ##                  C33 C34 C35 C36
      ##                      C44 C45 C46 
      ##                          C55 C56 
      ##                              C66'
      C_ijkl = '25.1  9.7 12.8 0.0 -1.300 0.00 
                     22.3 11.8 0.0  4.600 0.00 
                          21.8 0.0  1.400 0.00 
                               9.7  0.000 3.18 
                                   11.036 0.00 
                                          8.66'
    [] #HMX Modulus in GPA

  #Fracture
    [Fracture]
      type = GenericConstantMaterial
      prop_names = 'Gc'
      prop_values = 0.03
    []

  #Mie Grunessien Equation Of State
    [pressure_Y_RDX]
      type = ComputeMieGrun
      property_name = pressure_Y_RDX
      property_dT_name = pressure_dT_Y_RDX
      density = density
      temperature = temperature
      specific_heat = specific_heat
      Gamma0 = 0.667
      slope_UsUp = 1.905  #c0 is , will need to check against the value calculated in code, or rewrite to avoid using bulk modulus
      C0 = 2961.49
      reference_temperature = 300
    []
    #[pressure_Y_RDX]
    #  type = ComputeJWL_Grun
    #  property_name = pressure_Y_RDX
    #  property_dT_name = pressure_dT_Y_RDX
    #  specific_heat = specific_heat
    #  temperature = temperature
    #  A     =   2.01e5
    #  B     =  -5.20
    #  R1    =  12.4
    #  R2    =   1.24
    #  omega =   0.890
    #  ref_temperature = 300
    #[]

    [pressure_Y_PRD]
      type = ComputeJWL_Grun
      property_name = pressure_Y_PRD
      property_dT_name = pressure_dT_Y_PRD
      specific_heat = specific_heat
      temperature = temperature
      A     = 1242
      B     =  65.8
      R1    =   5
      R2    =   2.1
      omega =   0.34
      ref_temperature = 300
    []

  #Johnson Cook EOS


  #Chemistry Calculations 
    [w_RDX_loss]         
      type = ParsedMaterial
      property_name = w_RDX_loss
      coupled_variables = "temperature"
      material_property_names = 'Z_RDX E_RDX_R'
      expression = '-1.0*Z_RDX*exp(-1.0*E_RDX_R/temperature)'
    []
    [w_INT_make]         
      type = ParsedMaterial
      property_name = w_INT_make
      material_property_names = 'w_RDX_loss'
      expression = '-w_RDX_loss'
    []
    [w_INT_loss]         
      type = ParsedMaterial
      property_name = w_INT_loss
      coupled_variables = "temperature"
      material_property_names = 'Z_INT E_INT_R'
      expression = '-1.0*Z_INT*exp(-1.0*E_INT_R/temperature)'
    []
    [w_PRD_make]         
      type = ParsedMaterial
      property_name = w_PRD_make
      material_property_names = 'w_INT_loss'
      expression = '-w_INT_loss'
    []

    [Q_RDX]
      type = ParsedMaterial
      property_name = Q_RDX
      coupled_variables = 'temperature'
      material_property_names = 'a_RDX b_RDX T_ref_3phase'
      expression = 'if(temperature<=T_ref_3phase, a_RDX, a_RDX - b_RDX *(temperature-T_ref_3phase)) '
    []
    [Q_INT]
      type = ParsedMaterial
      property_name = Q_INT
      coupled_variables = 'temperature'
      material_property_names = 'a_INT b_INT T_ref_3phase'
      expression = 'if(temperature<=T_ref_3phase, a_INT, a_INT - b_INT *(temperature-T_ref_3phase)) '
    []

    [Heat_RDX]         
      type = ParsedMaterial
      property_name = Heat_RDX
      coupled_variables = 'temperature'
      material_property_names = 'Q_RDX Q_INT w_RDX_loss'
      expression = '1.0*(-Q_RDX)*w_RDX_loss' #*-1.0 as MatReaction is -LY
    []
    [Heat_INT]         
      type = ParsedMaterial
      property_name = Heat_INT
      coupled_variables = 'temperature'
      material_property_names = 'Q_INT w_PRD_make'
      expression = '1.0*(Q_INT)*w_PRD_make' #*-1.0 as MatReaction is -LY
    [../]

  [all_strain]
    type = ComputeIncrementalStrain #ComputeFiniteStrain
    
    displacements = 'disp_x disp_y'
    eigenstrain_names = ''
  []

  #[plate_stress]   type = ComputeFiniteStrainElasticStress  []
  [plate_stress]    
   type = ComputeIncStress_EOS 
   C0 = 0.1  
   C1 = 1.0
   pressure_eos_name = 'pressure_eos'
  []

 ###Eos###
  [SpecificVolume]   type = SpecificVolume  v_start = 0.85 []


  [pressure_eos] #Implements phase*pressure
    type = PhaseAddUp
    
    property_name = pressure_eos
    sum_materials = 'pressure_Y_RDX pressure_Y_PRD pressure_Y_PRD'
    values        = '         Y_RDX          Y_INT          Y_PRD'
  []

  [pressure_dT_eos] #Implements phase*pressure
    type = PhaseAddUp
    
    property_name = pressure_dT_eos
    sum_materials = 'pressure_dT_Y_RDX pressure_dT_Y_PRD pressure_dT_Y_PRD'
    values        = '         Y_RDX          Y_INT          Y_PRD'
  []

[]################################################################################

[Preconditioning][smp]  type = SMP  full = true  [][]
#[Postprocessors]
#  [ParticleVelo]    type = AverageNodalVariableValue  variable = velo_x  block = '200'  []
#  [ParticleAcce]    type = AverageNodalVariableValue  variable = acce_x  block = '200'  []
#  [initial_residual]  type = Residual  []
#[]
[Functions][./dts] 
  type = PiecewiseLinear 
  x = '0           3'
  y = '0.00001 0.00001' # 0.05
[../][]
[Executioner]
  type = Transient
  start_time = 0.0
  end_time = 99.00 #99 #200
  [./TimeStepper]  type = FunctionDT  function = dts  min_dt = 1e-6  [../]
  nl_rel_tol = 1E-6
  nl_abs_tol = 1E-6
  solve_type = PJFNK

  petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type -pc_hypre_boomeramg_max_iter   -snes_type'
  petsc_options_value = '               201    hypre      boomeramg                           20 vinewtonrsls'
  line_search = 'none'
[]

[Outputs]
  [Checkpoint]  type = Checkpoint  wall_time_interval           = 86400  num_files = 1[]
  [    Exodus]  type = Exodus      min_simulation_time_interval = 0.01                []
  #[       Csv]  type = CSV         time_step_interval           = 0.1                  []
[]
################################################################################################################################################################################################################################################

