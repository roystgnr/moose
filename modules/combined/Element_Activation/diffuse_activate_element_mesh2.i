T_room = 300

[Mesh]
  [./gen]
    type = GeneratedMeshGenerator
    dim = 3
    xmin =0
    xmax =2
    ymin =0
    ymax =2
    zmin =0
    zmax =1
    nx=2
    ny=2
    nz=1
  [../]
  # [./rename_bnd]
  #   type = RenameBoundaryGenerator
  #   input = gen
  #   # bottom = 0, front = 1, right = 2, back = 3, left = 4, top = 5
  #   old_boundary_id='0 1 2 3 4 5'
  #   new_boundary_name='bottom front right back left top'
  # [../]
  [./subdomain_id]
    input = gen
    type = ElementSubdomainIDGenerator
    subdomain_ids = '1 2
                     1 2'
  [../]
[]

[Variables]
  [./temp]
    initial_condition = ${T_room}
    block = '1'
  [../]
  [./dummy]
    block = '2'
  [../]
[]

[Kernels]
  [./time]
    type = ADHeatConductionTimeDerivative
    variable = temp
  [../]
  [./heat_conduct]
    type = ADHeatConduction
    variable = temp
    use_displaced_mesh = true
    thermal_conductivity = thermal_conductivity
  [../]
  [./heatsource]
    type = ADMatHeatSource
    material_property = volumetric_heat
    variable = temp
    scalar = 1
    use_displaced_mesh = true
  [../]

  [./dummy_diffusion]
    type = Diffusion
    variable = dummy
  [../]
  [./dummy_bforce]
    type = BodyForce
    variable = dummy
    function = forcing_fn_dummy
  [../]
[]

[Functions]
  [./heat_source_x]
    type = ParsedFunction
    value= '1.5'
  [../]
  [./heat_source_y]
    type = ParsedFunction
    value= '2.5*t'
  [../]
  [./heat_source_z]
    type = ParsedFunction
    value= '0.25'
  [../]
  [./forcing_fn_dummy]
    type = ConstantFunction
    value = 0.0
  [../]
[]

[BCs]
  [./temp_bottom_fix]
    type = ADDirichletBC
    variable = temp
    boundary = 'back'
    value = ${T_room}
  [../]

  [./dummy_bottom_fix]
    type = ADDirichletBC
    variable = dummy
    boundary = 'back'
    value = 0.0
  [../]
[]

[Materials]
  [./volumetric_heat]
    type = FunctionPathEllipsoidHeatSource
    a = 2
    b = 2
    c = 2
    power = 1000
    efficienty = 1.0
    factor = 2
    velocity = 5
    function_x= heat_source_x
    function_y= heat_source_y
    function_z= heat_source_z
  [../]
  [./density]
    type = ADDensity
    density = 4.43e-6
  [../]
  [./heat]
    type = ADHeatConductionMaterial
    specific_heat = 603
    thermal_conductivity = 10e-2
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]


[Executioner]
  type = Transient

  automatic_scaling = true

  #Preconditioned JFNK (default)
  solve_type = 'NEWTON'

  # petsc_options = '-snes_ksp'
  petsc_options_iname = '-pc_type -ksp_type'
  petsc_options_value = 'lu  preonly'

  line_search = 'none'

  l_max_its = 10
  nl_max_its = 20
  nl_rel_tol = 1e-4
  nl_abs_tol = 1e-12
  # l_tol = 1e-5

  start_time = 0.0
  end_time = 0.1
  dt = 1e-1
  dtmin = 1e-4
[]

# [AuxVariables]
#   [./activated_elem]
#     order = CONSTANT
#     family = MONOMIAL
#   [../]
# []

# [AuxKernels]
#   [./activated_elem]
#     type = ActivatedElementsMarker
#     melt_temperature = 480
#     temp_aux = temp
#     variable = activated_elem
#     execute_on = timestep_begin
#   [../]
# []

[UserObjects]
  [./activated_elem_uo]
    type = ActivateElementTemp
    execute_on = timestep_begin
    function_x= heat_source_x
    function_y= heat_source_y
    function_z= heat_source_z
    # activate_tol=1e-2
    active_subdomain_id = 1
  [../]
[]

# [Postprocessors]
#   [./temperature]
#     type = ElementAverageValue
#     variable = temp
#   [../]
# []

[Outputs]
  file_base = 'temp_out'
  [./exodus]
    type = Exodus
  [../]
[]
