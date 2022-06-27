# Expected read counts and proportion of reads can be calculated

    Code
      rc_out
    Output
      class: PhIPData 
      dim: 10 10 
      metadata(8): seed a_pi ... b_c fc
      assays(9): counts logfc ... true_theta expected_rc
      rownames(10): 1 2 ... 9 10
      rowData names(2): a_0 b_0
      colnames(10): 1 2 ... 9 10
      colData names(5): group n_init n true_c true_pi
      beads-only name(4): beads

# Bayes factors can be calculated

    Code
      bf_out
    Output
      class: PhIPData 
      dim: 10 10 
      metadata(8): seed a_pi ... b_c fc
      assays(9): counts logfc ... true_theta bayes_factors
      rownames(10): 1 2 ... 9 10
      rowData names(2): a_0 b_0
      colnames(10): 1 2 ... 9 10
      colData names(7): group n_init ... c pi
      beads-only name(4): beads

