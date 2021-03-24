[Mesh]
  [./smg]
    type = StitchedSubgenerators
    inputs = 'left.e center.e right.e'
    clear_stitched_boundary_ids = true
    stitch_boundaries_pairs = 'right left;
                               right left'
    parallel_type = 'replicated'
  []
[]

[Outputs]
  exodus = true
[]
