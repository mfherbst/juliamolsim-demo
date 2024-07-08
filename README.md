# JuliaMolSim demo for JuliaCon 2024

This repository contains the Julia source code used
to generate the aluminium simulation shown in the JuliaCon 2024 keynote
[Materials Modeling: Bonding across atoms, code, and people](https://michael-herbst.com/slides/juliacon2024)
[(abstract)](https://pretalx.com/juliacon2024/talk/RDGSFV/).

## Reproducing the aluminium simulation

1. Add the `ACEsuit` Julia registry and instantiate the environment:
```julia
using Pkg
pkg"registry add https://github.com/ACEsuit/ACEregistry"
pkg"instantiate"
```

2. Execute `1_run_dft.jl`. This script runs the DFT calculations
   and takes a few days. You can skip this if you want to use the precomputed
   data stored in the repository.

3. Execute `2_fit_potential.jl`. This script fits the ACE model from the
   DFT data.

4. Execute `3_run_md.jl`. This script uses the ACE ML model generated in 3.
   and runs a molecular dynamics simulation on aluminium for 1000 time steps.
   The resulting video is stored as `simulation.mp4`.

Note: Due to various random initialisations, the precise potential data
and the precise video will look different each time you run steps 3 or 4.

## Further reading
- [AtomsBase documentation](https://juliamolsim.github.io/AtomsBase.jl/stable)
- [DFTK documentation](https://docs.dftk.org)
- [ACEpotentials documentation](https://acesuit.github.io/ACEpotentials.jl/)
  (in particular the part on the [AtomsBase interface](https://acesuit.github.io/ACEpotentials.jl/dev/tutorials/AtomsBase_interface/)).
- [Molly documentation](https://juliamolsim.github.io/Molly.jl/stable/)
