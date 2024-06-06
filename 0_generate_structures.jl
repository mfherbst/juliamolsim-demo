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
using Unitful

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
maxrattle_pos  = 0.1u"Å"
maxrattle_cell = 0.1
max_supercell  = 3

# Category 1 structures: MD stability
let file = "Al_supercells_1.extxyz", n_structures = 100
    if !isfile(file)
        systems = map(1:n_structures) do i
            nx = rand(1:max_supercell)
            ny = rand(1:max_supercell)
            nz = rand(1:max_supercell)
            system = bulk(:Al, cubic=true) * (nx, ny, nz)
            system = rattle_system(system,
                                   maxrattle_pos  * rand(),
                                   maxrattle_cell * rand())
            FlexibleSystem(system; repeat=[nx, ny, nz])
        end
        save_trajectory(file, systems)
    end
end

let file = "Al_supercells_2.extxyz", n_structures = 100
    if !isfile(file)
        systems = map(1:n_structures) do i
            nx = ny = nz = 2
            system = bulk(:Al, cubic=true) * (nx, ny, nz)
            system = rattle_system(system,
                                   maxrattle_pos  * rand(),
                                   maxrattle_cell * rand())
            FlexibleSystem(system; repeat=[nx, ny, nz])
        end
        save_trajectory(file, systems)
    end
end

# Category 2 structures: Virials, elastic constants, long-range elastic fields
let file = "Al_bulk_1.extxyz", n_structures = 100
    if !isfile(file)
        systems = map(1:n_structures) do i
            system = rattle_cell(bulk(:Al), maxrattle_cell * rand(); update_position=false)
            FlexibleSystem(system; repeat=[1, 1, 1])
        end
        save_trajectory(file, systems)
    end
end

# Category 3 structures: Defects
let file = "Al_defect_1.extxyz", n_structures = 100
    if !isfile(file)
        systems = map(1:n_structures) do i
            nx = rand(1:max_supercell)
            ny = rand(1:max_supercell)
            nz = rand(1:max_supercell)
            system = bulk(:Al, cubic=true) * (nx, ny, nz)

            idx = rand(1:length(system))
            system = deleteat!(system, idx)
            system = rattle_system(system,
                                   maxrattle_pos  * rand(),
                                   maxrattle_cell * rand())
            FlexibleSystem(system; repeat=[nx, ny, nz])
        end
        save_trajectory(file, systems)
    end
end

# Category 4 structures: Defect test systems:
let file = "Al_test_1.extxyz", n_structures = 2
    if !isfile(file)
        systems = map(1:n_structures) do i
            nx = rand(1:max_supercell)
            ny = rand(1:max_supercell)
            nz = rand(1:max_supercell)
            system = bulk(:Al, cubic=true) * (nx, ny, nz)

            if i == 1
                system = deleteat!(system, ceil(Int, length(system)/2))
            end
            system = rattle!(system, maxrattle_pos * rand())  # Rattle positions (inplace)
            FlexibleSystem(system; repeat=[nx, ny, nz])
        end
        save_trajectory(file, systems)
    end
end

