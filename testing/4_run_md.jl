#!/bin/bash
#=
BN=$(basename "$0" .jl)
xvfb-run -s '-screen 0 1024x768x24' julia --project -t1 $BN.jl |& tee $BN.log
exit $?
=#

using ACEpotentials
using AtomsBase
using AtomsBuilder
using Molly
using AtomsIO
using GLMakie

# Load potential and system
potential = load_potential("Al_DFTK_dataset_potential.json"; new_format=true)

defect = true
if defect
    println("Defect system")
    structure = load_system("Al_test_1.extxyz", 1)
else
    println("No defect")
    structure = load_system("Al_test_1.extxyz", 2)
end

## Molly simulation
# Pack data to Molly compatible format, note this is a custom 
# Molly.System builder, that is not part of Molly (ACEmd-Molly extension).
# Could make a PR for it...
sys = Molly.System(structure, potential)

# Set up temperature and velocities
temp = 298.0u"K"
vel = random_velocities!(sys, temp)

# Add loggers
sys = Molly.System(sys; loggers=(
    temp=TemperatureLogger(100),
    coords=CoordinateLogger(1),
))

# Set up simulator
simulator = VelocityVerlet(; dt=0.5u"fs", 
                             coupling=AndersenThermostat(temp, 0.5u"ps"), )

println("Performing 1000 MD steps")
simulate!(sys, simulator, 1000)
@show sys.loggers.temp.history

println("Visualising results")
visualize(sys.loggers.coords, sys.boundary, "simulation.mp4")
