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

    # TODO Use the right units here
    energy = data[:energies][:total]     # in Hartree
    forces = reduce(hcat, data[:forces]) # in Hartree / Bohr, 3xn_atom array
    properties = (; energy=ustrip(auconvert(u"eV", energy)),
                    forces=ustrip.(auconvert.(u"eV/Å", forces)))
    if data[:stresses] != "nothing"
        virial = Array(reshape(data[:stresses], 3, 3))
        properties = merge(properties,
                           (; virial=ustrip.(auconvert.(u"eV", virial))))
    end

    # Work around a bug in ExtXYZ where just using FlexibleSystem does not work
    #
    # elseif v isa AbstractArray{<:ExtxyzType} should be added as an option
    # to the type check in write_dict
    #
    # return FlexibleSystem(system; properties...)

    sys_with_props = FlexibleSystem(system; properties...)

    extxyz = with_logger(NullLogger()) do
        ExtXYZ.write_dict(sys_with_props)
    end

    extxyz["info"]["forces"] = properties.forces
    if haskey(properties, :virial)
        extxyz["info"]["virial"] = properties.virial
    end
    extxyz
end
ExtXYZ.write_frames("Al_DFTK_dataset.extxyz", systems)
# save_trajectory("Al_DFTK_dataset.extxyz", systems)
