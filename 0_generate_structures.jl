using AtomsBuilder
using AtomsIO
using AtomsBase

systems = []
for repeat in 4:4
    system = bulk(:Al, cubic=true) * (repeat, repeat, repeat)
    push!(systems, system)
    # push!(systems, rattle!(FlexibleSystem(system)), 1e-4)
end
save_trajectory("Al_bulk.extxyz", systems)
