#!/bin/sh
#=
BN=$(basename "$0" .jl)
julia --project -t1 $BN.jl
exit $?
=#

using AtomsIO
using JSON3
using ExtXYZ
using Unitful
using UnitfulAtomic
using Logging

dftk_results = filter!(endswith(".json"), readdir("dftk_output"))
println("Found $(length(dftk_results)) result files.")

# Concatenate systems
systems = map(dftk_results) do dftk
    prefix, _ = splitext(dftk)
    structure, index = split(prefix, "-")
    load_system(structure * ".extxyz", parse(Int, index))
end
save_trajectory("Al_structures.extxyz", systems)

# Concatenate fitting results
systems = map(dftk_results) do dftk
    prefix, _ = splitext(dftk)
    structure, index = split(prefix, "-")
    system = load_system(structure * ".extxyz", parse(Int, index))
    data = open(JSON3.read, joinpath("dftk_output", dftk))

    properties = (; energy=ustrip(u"eV", data[:energies][:total] * u"hartree"))
    if data[:stresses] != "nothing"
        virial = Array(reshape(data[:stresses], 3, 3))
        properties = merge(properties, (; virial=ustrip.(u"eV", virial * u"eV")))
    end
    atoms = map(system, data[:forces]) do at, at_force
        Atom(; pairs(at)..., force=ustrip.(u"eV/Å", at_force * u"hartree/bohr"))
    end
    return FlexibleSystem(system; atoms, properties...)
end
println("Extracted $(length(systems)) systems for building the potential.")
save_trajectory("Al_DFTK_dataset.extxyz", systems)
