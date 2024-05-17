#!/bin/sh
#=
BN=$(basename "$0" .jl)
julia --project -t1 $BN.jl
exit $?
=#
using AtomsBuilder
using AtomsIO
using AtomsBase
using LinearAlgebra


function rattle_cell(system, r_cell; update_position=true)
    F = I + r_cell * (rand(3, 3) .- 0.5)
    bbox = Ref(F) .* bounding_box(system)

    if update_position
        # Update atomic positions and bounding box
        atoms = [Atom(at; position=F * position(at)) for at in system]
        FlexibleSystem(system; bounding_box=bbox, atoms)
    else
        # Just update bounding box
        FlexibleSystem(system; bounding_box=bbox)
    end
end
function rattle_system(system, r_pos, r_cell)
    system = rattle_cell(system, r_cell; update_position=true)
    rattle!(system, r_pos)  # Rattle positions (inplace)
    system
end

# Some parameters
maxrattle_pos  = 0.1
maxrattle_cell = 0.1
max_supercell  = 3

# Category 1 structures: MD stability
let
    n_structures = 100
    file = "Al_supercells_1.extxyz"

    if !isfile(file)
        systems = map(1:n_structures) do i
            n = rand(1:max_supercell)
            system = bulk(:Al, cubic=true) * (n, n, n)
            rattle_system(system,
                          maxrattle_pos  * rand(),
                          maxrattle_cell * rand())
        end
        save_trajectory(file, systems)
    end
end

# Category 2 structures: Virials, elastic constants, long-range elastic fields
let
    n_structures = 100
    file = "Al_bulk_1.extxyz"

    if !isfile(file)
        systems = map(1:n_structures) do i
            rattle_cell(bulk(:Al), maxrattle_cell * rand(); update_position=false)
        end
        save_trajectory(file, systems)
    end
end

# Category 3 structures: Defects
# TODO
