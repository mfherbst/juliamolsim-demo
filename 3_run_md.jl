#!/bin/bash
#=
BN=$(basename "$0" .jl)
xvfb-run -s '-screen 0 1024x768x24' julia --project -t1 $BN.jl
exit $?
=#
using ACEpotentials
using AtomsBase
using AtomsBuilder
using Molly
using AtomsIO
using GLMakie

potential = load_potential("Al_potential.json"; new_format=true)
structure = load_system("Al_test.extxyz", 2)
sys = Molly.System(structure, potential)

temp = 298.0u"K"
random_velocities!(sys, temp)
sys = Molly.System(sys; loggers=(; coords=CoordinateLogger(1)))
simulator = VelocityVerlet(; dt=0.5u"fs", coupling=AndersenThermostat(temp, 0.5u"ps"), )
simulate!(sys, simulator, 1000)  # Simulate 1000 time steps
visualize(sys.loggers.coords, sys.boundary, "simulation.mp4";
          color=:grey, markersize=0.2, show_boundary=false)
