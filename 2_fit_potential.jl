#!/bin/bash
#=
BN=$(basename "$0" .jl)
julia --project -t1 $BN.jl
exit $?
=#
using ACEpotentials
using ExtXYZ
all_data = ExtXYZ.Atoms.(ExtXYZ.read_frames("Al_training_data.extxyz"))
model = acemodel(
    elements=[:Al, ],
    order=3,
    totaldegree=[20, 16, 12],
    Eref = [:Al => -2.180568046, ],
)
acefit!(model, all_data; solver=ACEfit.BLR(), prior=smoothness_prior(model; p=4))
save_potential("Al_potential.json", model)
