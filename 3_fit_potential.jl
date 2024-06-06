#!/bin/bash
#=
BN=$(basename "$0" .jl)
julia --project -t1 $BN.jl |& tee $BN.log
exit $?
=#

using ACEpotentials
using ExtXYZ
using Random
using Unitful

# Load data and split into training and testing
all_data = ExtXYZ.Atoms.(ExtXYZ.read_frames("Al_DFTK_dataset.extxyz"))
println("Read $(length(all_data)) systems for building the potential.")
perm = randperm(length(all_data))
data_train = all_data[perm[1:floor(Int, 0.8length(all_data))]]
data_test  = all_data[perm[floor(Int, 0.8length(all_data))+1:end]]
@assert length(data_train) + length(data_test) == length(all_data)
@show length(data_train)

# Set up ACEmodel for fitting
model = acemodel(
    elements=[:Al, ],
    order=3,
    totaldegree=[20, 16, 12],
    Eref = [:Al => -2.180568046, ],
)
@show length(model.basis);

# Fit the model
acefit!(model, data_train; solver=ACEfit.BLR(), prior=smoothness_prior(model; p=4))

@info("Training Error Table") 
ACEpotentials.linear_errors(data_train, model);
@info("Test Error Table")  
ACEpotentials.linear_errors(data_test, model);

# Save result
save_potential("Al_DFTK_dataset_potential.json", model)
