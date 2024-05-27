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
systems = map(dftk_results) do dftk
    prefix, _ = splitext(dftk)
    structure, index = split(prefix, "-")
    system = load_system(structure * ".extxyz", parse(Int, index))
    data = open(JSON3.read, joinpath("dftk_output", dftk))

    energy = data[:energies][:total]     * u"hartree"
    force  = reduce(hcat, data[:forces]) * u"hartree/bohr"  # 3 x n_atom array
    properties = (; energy=ustrip(u"eV",   energy),
                    force=ustrip.(u"eV/Å", force))
    if data[:stresses] != "nothing"
        virial = Array(reshape(data[:stresses], 3, 3))
        properties = merge(properties, (; virial=ustrip.(u"eV", virial * u"eV")))
    end
    return FlexibleSystem(system; properties...)
end
save_trajectory("Al_DFTK_dataset.extxyz", systems)
