# warns when overwriting sampleInfo [plain]

    Code
      brew(sim_data)
    Message <cliMessage>
      ! Values in the following assays will be overwritten: sampleInfo
      
      -- Running JAGS ----------------------------------------------------------------
      Sample runs
      
      -- Summarizing results ---------------------------------------------------------
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

# warns when overwriting sampleInfo [ansi]

    Code
      brew(sim_data)
    Message <cliMessage>
      [33m![39m Values in the following assays will be overwritten: sampleInfo
      
      [36m--[39m [1m[1mRunning JAGS[1m[22m [36m----------------------------------------------------------------[39m
      Sample runs
      
      [36m--[39m [1m[1mSummarizing results[1m[22m [36m---------------------------------------------------------[39m
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

# warns when overwriting sampleInfo [unicode]

    Code
      brew(sim_data)
    Message <cliMessage>
      ! Values in the following assays will be overwritten: sampleInfo
      
      ── Running JAGS ────────────────────────────────────────────────────────────────
      Sample runs
      
      ── Summarizing results ─────────────────────────────────────────────────────────
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

# warns when overwriting sampleInfo [fancy]

    Code
      brew(sim_data)
    Message <cliMessage>
      [33m![39m Values in the following assays will be overwritten: sampleInfo
      
      [36m──[39m [1m[1mRunning JAGS[1m[22m [36m────────────────────────────────────────────────────────────────[39m
      Sample runs
      
      [36m──[39m [1m[1mSummarizing results[1m[22m [36m─────────────────────────────────────────────────────────[39m
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

