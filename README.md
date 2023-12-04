
# Sofa.jl
code to compute improved upper bounds on the moving sofa problem

### Instructions

To compute upper bounds on the sofa problem, run:

```
julia --threads auto src/run.jl
```

Intermediate results will be saved to the folder `state` every 1000 iterations.

### Continuing with Precomputed Boxes

To continue lowering the bound using precomputed boxes, follow these steps:

1. Locate the compressed file in the `boxes/sup2.37` folder.
2. Extract this file into a `state` folder at the root directory of the project. If the `state` folder does not exist, create it first.

   Example:
   ```bash
   mkdir -p state
   zstd -d boxes/sup2.37/4209000_2023-12-03T22:55:32.431.dat.zst -o state/extracted_file
   ```

   Replace `extracted_file` with the desired name for the extracted file.

This will set up the environment with the precomputed boxes, allowing the computation to resume from a specific state.
