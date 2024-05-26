#!/bin/sh
#=
BN=$(basename "$0" .jl)
julia --project -t1 $BN.jl
exit $?
=#

using AtomsIO
using JSON3

dftk_results = filter!(endswith(".json"), readdir("dftk_output"))
systems = map(dftk_results) do dftk
    prefix, _ = splitext(dftk)
    structure, index = split(prefix, "-")
    system = load_system(structure * ".extxyz", parse(Int, index))
    data = open(JSON3.read, joinpath("dftk_output", dftk))

    properties = (; energy=data[:energies][:total],
                    forces=reduce(hcat, data[:forces]))  # 3xn_atom array
    if data[:stresses] != "nothing"
        properties = merge(properties,
                           (; stresses=Array(reshape(data[:stresses], 3, 3))))
    end

    FlexibleSystem(system; properties...)
end

# I think there is a bug in extxyz this is why this does not work.
# elseif v isa AbstractArray{<:ExtxyzType} should be added as an option
# to the type check
save_trajectory("Al_DFTK_dataset.extxyz", systems)
