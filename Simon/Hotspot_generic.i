#unit system is standar
#mass = 10⁻¹⁵ kg
#time = 10⁻⁹ s
#length = 10⁻⁶ m
#pressure = 10⁹ Pa


###NOTE###
#this model considers only chemical and thermal behavior

#the mechanics, plasticity, viscosity are neglected

#density is assumed to be constant for simplification purposes, as well as the k and Cp of the mixture

#all transport properties are assumed to be constant ---> This has to be updated as heat capacity and
#conductivity depend strongly on temperature and phase 

#the thermo-mechanical coupling is neglected. In future versions of this model, thermal expansion and
#linear elastic behavior will be included


[GlobalParams]
	#number of displacements and directions
	displacements = 'disp_x disp_y'
[]

[Mesh]
	type = GeneratedMesh
	dim = 2
	nx = 100
	ny = 100
	xmax = 3000
	xmin = 0
	ymax = 3000
	ymin = 0
	boundary_id = '0 1 2 3'
	boundary_name = 'bottom right top left'
[]

#[Adaptivity]
#  marker = errorfrac # this specifies which marker from 'Markers' subsection to use
#  steps = 1 # run adaptivity 1 times, recomputing solution, indicators, and markers each time
#  max_h_level = 2
#
#  # Use an indicator to compute an error-estimate for each element:
#  [./Indicators]
#    # create an indicator computing an error metric for the convected variable
#    [./error]
#      # arbitrary, use-chosen name
#      type = GradientJumpIndicator
#      variable = temperature
#      outputs = none
#    [../]
#  [../]
#
#  # Create a marker that determines which elements to refine/coarsen based on error estimates
#  # from an indicator:
#  [./Markers]
#    [./errorfrac]
#      # arbitrary, use-chosen name (must match 'marker=...' name above
#     type = ErrorFractionMarker
#      indicator = error # use the 'error' indicator specified above
#      refine = 0.5 # split/refine elements in the upper half of the indicator error range
#     coarsen = 0 # don't do any coarsening
#      outputs = none
#    [../]
#  [../]
#[]

[Variables]
	###variables are often used to compute kernel contributions
	###to the weak form or residual of a specific equation
	#mechanical vars
	[./disp_x]
	    	order = FIRST
	    	family = LAGRANGE
	[../]
	[./disp_y]
		order = FIRST
		family = LAGRANGE
	[../]

	#thermal vars
	[./temperature]
		order = FIRST
		family = LAGRANGE
	[../]

	#chemical vars
	[./Y1]
		order = FIRST
		family = LAGRANGE
	[../]
	[./Y2]
		order = FIRST
		family = LAGRANGE
	[../]
	[./Y3]
		order = FIRST
		family = LAGRANGE
	[../]
	[./Y4]
		order = FIRST
		family = LAGRANGE
	[../]
		

[]

[AuxVariables]
	#[./density]
	#    order = CONSTANT
	#    family = MONOMIAL
	#[../]

    [saved_x]
    []
    [saved_y]
    []
	
	##dynamics variables
	[./vx]
	[../]
	[./ax]
	[../]
	[./vy]
	[../]
	[./ay]
	[../]
	
	[./Pressure]
		order = CONSTANT
		family = MONOMIAL
	[../]
	[./J]
		order = CONSTANT
		family = MONOMIAL
	[../]

[]

[Modules]
    [./TensorMechanics/DynamicMaster]
        [all]
            add_variables = true 
            strain = SMALL
            use_automatic_differentiation = false
            generate_output = 'stress_xx stress_xy stress_xz stress_yy stress_yz stress_zz
                                strain_xx strain_xy strain_xz strain_yy strain_yz strain_zz
                                elastic_strain_xx elastic_strain_xy elastic_strain_xz elastic_strain_yy
                                elastic_strain_yz elastic_strain_zz mechanical_strain_xx mechanical_strain_xy
                                mechanical_strain_xz mechanical_strain_yy mechanical_strain_yz mechanical_strain_zz
                                vonmises_stress'
            save_in = 'saved_x saved_y'
        []
    []
[]

[Kernels]

	##### temperature dependent kernels
	#these kernels build the terms for the heat equation
    
	[./dT_dX] #spatial derivatives (nabla x)
        type = HeatConduction
        variable = temperature
        diffusion_coefficient = 'thermal_conductivity'
	[../]
	[./dT_dt] #temperature time derivative
        type = HeatConductionTimeDerivative
        variable = temperature
        density_name = 'density'
        specific_heat = 'specific_heat'
	[../]
	
	###chemical decomposition heat source
	##the chemical decomposition kernels compute the diffusion equation for each species
	[./chemHeatSource]
		type = QYdot
		variable = temperature
		Y1 = Y1
		Y2 = Y2
		Y3 = Y3
		Y4 = Y4
		#heats of reaction
		Q1 = -412.4 #check units and sign
		Q2 = 1255.2 #check units and sign
		Q3 = 5020.8 #check units and sign
	[../]
	
	[./EOSHeatSource]
		type = MieGruneisenHS
		variable = temperature
		temperature = temperature
		Gamma = 0.7
		T_ref = 300
		C0 = 1.0
		C1 = 0.1
		beta_av = 0.4
		viscosity_type = '0' #1 is constant, otherwise uses elem.size
		element_size = 3
	[../]

	###Chemical species diffusion equations
	##each specie is modeled using a poisson's diffusion equation
	##of the type dY_dt + nabla^2 Y = Y_dot

	##this computes the first term of the conservation equation \frac{\partial Y_i}{\partial t}, different from \dot{Y}_i
	[./dY1_dt]
		type = TimeDerivative
		variable = Y1
	[../]
	[./dY2_dt]
		type = TimeDerivative
		variable = Y2
	[../]
	[./dY3_dt]
		type = TimeDerivative
		variable = Y3
	[../]
	[./dY4_dt]
		type = TimeDerivative
		variable = Y4
	[../]
	####

	
	##this computes the right hand side of the conservation equation, rate of change of Y_i species using the Tarver 3-step model
	[./Y1dot]
		type = Y1_dot
		variable = Y1
		temperature = temperature
	[../]
	[./Y2dot]
		type = Y2_dot
		variable = Y2
		Y1 = Y1
		temperature = temperature
	[../]
	[./Y3dot]
		type = Y3_dot
		variable = Y3
		Y2 = Y2
		temperature = temperature
	[../]
	[./Y4dot]
		type = Y4_dot
		variable = Y4
		Y3 = Y3
		temperature = temperature
	[../]

	#velocity and acceleration kernels
	#using newmark beta integration scheme
	[./inertia_x]
    	type = InertialForce
    	variable = disp_x
    	velocity = vx
    	acceleration = ax
    	beta = 0.3025 ###from dandekar 2019 sec 2.1
    	gamma = 0.6 ###from dandekar 2019 sec 2.1
	[../]
	[./inertia_y]
    	type = InertialForce
    	variable = disp_y
    	velocity = vy
    	acceleration = ay
    	beta = 0.3025 ###from dandekar 2019 sec 2.1
    	gamma = 0.6 ###from dandekar 2019 sec 2.1
	[../]
[]

[AuxKernels]
	[./vx]
    	type = NewmarkVelAux
    	variable = vx
    	acceleration = ax
    	gamma = 0.6
	[../]
	[./vy]
    	type = NewmarkVelAux
    	variable = vy
    	acceleration = ay
    	gamma = 0.6
	[../]
	[./ax]
    	type = NewmarkAccelAux
    	variable = ax
    	displacement = disp_x
    	velocity = vx
    	beta = 0.3025
	[../]
	[./ay]
    	type = NewmarkAccelAux
    	variable = ay
    	displacement = disp_y
    	velocity = vy
    	beta = 0.3025
	[../]
	
	[./Pressure]
		type = MaterialRealAux
		variable = Pressure
		property = pressure_eos
	[../]
	[./J]
		type = MaterialRealAux
		variable = J
		property = J
	[../]
[]

[BCs]
	[./Periodic]
		[./x_temp]
			variable = temperature
			primary = left
			secondary = right
			translation = '3000 0 0'
		[../]
		[./y_temp]
			variable = temperature
			primary = bottom
			secondary = top
			translation = '0 3000 0'
		[../]
		[./x_disp]
			variable = disp_x
			primary = left
			secondary = right
			translation = '3000 0 0'
		[../]
		[./y_disp]
			variable = disp_y
			primary = bottom
			secondary = top
			translation = '0 3000 0'
		[../]
	[../]
[]

[ICs]
	[./chemY1_IC]
		type = ConstantIC
		variable = Y1
		value = 1.0
	[../]
	[./chemY2_IC]
		type = ConstantIC
		variable = Y2
		value = 0.0
	[../]
	[./chemY3_IC]
		type = ConstantIC
		variable = Y3
		value = 0.0
	[../]
	[./chemY4_IC]
		type = ConstantIC
		variable = Y4
		value = 0.0
	[../]
	[./temp_hotspot]
		#this is done to replicate the chemY1 initial condition
		type = SmoothCircleIC
		variable = temperature
		invalue = 600
		outvalue = 300
		int_width = 10
		x1 = 1500
		y1 = 1500
		radius = 250
	[../]
	#[./single_element_temp]
	#	type = ConstantIC
	#	variable = temperature
	#	value = 2000
	#[../]
[]

[Functions]
[]

[Materials] 
    ### Mechanics ###
    
	[./const_mat_HMX]
		type = GenericConstantMaterial
		prop_names = 'density specific_heat thermal_conductivity'
		prop_values = '1905e-9 2357 0.31e-3' #assumed constant
	[../]
	#[./aniso_HMX]
	#	type = ComputeElasticityTensor
	#	fill_method = symmetric21
	#	C_ijkl = '25.1 9.7 12.8 0 -1.3 0 22.3 11.8 0 4.6 0 21.8 0 1.4 0 9.7 0 3.18 11.036 0 8.66'
	#[../]
	[./Isotropic_HMX]
		type = ComputeIsotropicElasticityTensor
		youngs_modulus = 25.12
		poissons_ratio = 0.24
	[../]

	[./ComputeStress]
	    	type = MieGruneisenStress
	    	temperature = temperature
	    	viscosity_type = '0' #1 is constant, otherwise uses elem.size
	    	element_size = 3
	    	Gamma = 0.7
	    	T_ref = 300
	    	C0 = 1.0
	    	C1 = 0.1
	    	slope = 2.29
	    	A = 1668.9 #GPa
	    	B = 59.69 #GPa
	    	R1 = 5.9
	    	R2 = 2.1
	    	omega = 0.45
	    	Y1 = Y1
	    	Y2 = Y2
	    	Y3 = Y3
	    	Y4 = Y4
	    	use_JWL = 0
	[../]
	
	[./ComputeYdots]
		type = ComputeYdots
		temperature = temperature
		Y1 = Y1
		Y2 = Y2
		Y3 = Y3
		Y4 = Y4

		Z1 = 1.41E+12 #check units, these units are 1/ns and is Log(Z)
		Z2 = 1.58E+07 #check units, these units are 1/ns and is Log(Z)
		Z3 = 1.60E+03 #check units, these units are 1/ns and is Log(Z)

		E1 = 734.2 #check units, these units are kcal/mol
		E2 = 614.4 #check units, these units are kcal/mol
		E3 = 475.1 #check units, these units are kcal/mol

		Rg = 0.028 #check units, units are kJ/(kg HMX) - K
		[../]
[]

[Preconditioning]
	[./smp]
		type = SMP
		full = true
	[../]
	#[./vcp]
	#type = VCP
	#full = true
	#[../]
[]

[Executioner]
	type = Transient
	start_time = 0.0
	dt = 1e-2
	dtmin = 1e-50
	end_time = 2
	nl_rel_tol = 1E-6
	nl_abs_tol = 1E-6
	solve_type = PJFNK
	petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type -pc_hypre_boomeramg_max_iter'
	petsc_options_value = '201 hypre boomeramg 50'
	line_search = 'none'
[]

[VectorPostprocessors]
[]

[Outputs]
	exodus = true
	interval = 1
	checkpoint = false
[]

[Debug]
	#show_material_props = true
	show_var_residual_norms = true
	#show_execution_order = ALWAYS
	check_jacobian = true
[]