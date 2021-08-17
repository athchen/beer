# warns when overwriting matrices [plain]

    Code
      edgeR(sim_data, assay.names = c("logfc", "counts"))
    Message <cliMessage>
      ! The following assays are not empty and will be overwritten: counts
    Output
      class: PhIPData 
      dim: 50 10 
      metadata(8): seed a_pi ... b_c fc
      assays(8): counts logfc ... true_b true_theta
      rownames(50): 1 2 ... 49 50
      rowData names(2): a_0 b_0
      colnames(10): 1 2 ... 9 10
      colData names(5): group n_init n c pi
      beads-only name(4): beads

# warns when overwriting matrices [ansi]

    Code
      edgeR(sim_data, assay.names = c("logfc", "counts"))
    Message <cliMessage>
      [33m![39m The following assays are not empty and will be overwritten: counts
    Output
      class: PhIPData 
      dim: 50 10 
      metadata(8): seed a_pi ... b_c fc
      assays(8): counts logfc ... true_b true_theta
      rownames(50): 1 2 ... 49 50
      rowData names(2): a_0 b_0
      colnames(10): 1 2 ... 9 10
      colData names(5): group n_init n c pi
      beads-only name(4): beads

# warns when overwriting matrices [unicode]

    Code
      edgeR(sim_data, assay.names = c("logfc", "counts"))
    Message <cliMessage>
      ! The following assays are not empty and will be overwritten: counts
    Output
      class: PhIPData 
      dim: 50 10 
      metadata(8): seed a_pi ... b_c fc
      assays(8): counts logfc ... true_b true_theta
      rownames(50): 1 2 ... 49 50
      rowData names(2): a_0 b_0
      colnames(10): 1 2 ... 9 10
      colData names(5): group n_init n c pi
      beads-only name(4): beads

# warns when overwriting matrices [fancy]

    Code
      edgeR(sim_data, assay.names = c("logfc", "counts"))
    Message <cliMessage>
      [33m![39m The following assays are not empty and will be overwritten: counts
    Output
      class: PhIPData 
      dim: 50 10 
      metadata(8): seed a_pi ... b_c fc
      assays(8): counts logfc ... true_b true_theta
      rownames(50): 1 2 ... 49 50
      rowData names(2): a_0 b_0
      colnames(10): 1 2 ... 9 10
      colData names(5): group n_init n c pi
      beads-only name(4): beads

