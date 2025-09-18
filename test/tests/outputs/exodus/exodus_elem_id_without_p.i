!include exodus_elem_id_base.i

[Mesh]
  [order_conversion]
    type = ElementOrderConversionGenerator
    input = tempid_4
    conversion_type = SECOND_ORDER
  []
[]

[Adaptivity]
  switch_h_to_p_refinement = true
  initial_marker = uniform
  initial_steps = 1
  [Markers/uniform]
    type = UniformMarker
    mark = REFINE
    block = 1
  []
[]

[Variables]
  [u]
    family = HIERARCHIC
    order = FIRST
  []
[]

[Kernels]
  [src]
    type = BodyForce
    variable = u
  []
  [l2]
    type = Reaction
    variable = u
  []
[]
