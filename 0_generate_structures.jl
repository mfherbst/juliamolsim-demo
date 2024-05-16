#!/bin/sh
#=
BN=$(basename "$0" .jl)
julia --project -t1 $BN.jl
exit $?
=#
using AtomsBuilder
using AtomsIO
using AtomsBase

systems = AbstractSystem[]
for repeat in 1:3
    system = bulk(:Al, cubic=true) * (repeat, repeat, repeat)
    push!(systems, rattle!(FlexibleSystem(system), 1e-3))
end
save_trajectory("Al_bulk.extxyz", systems)
