using ACEpotentials
using AtomsBase
using AtomsBuilder
using Molly
using AtomsIO

# Load potential and system
potential = load_potential("Al_DFTK_dataset_potential.json"; new_format=true)
structure = load_system("Al_test_1.extxyz", 1)

## Molly simulation
# Pack data to Molly compatible format, note this is a custom 
# Molly.System builder, that is not part of Molly (ACEmd-Molly extension).
# Could make a PR for it...
sys = Molly.System(structure, potential)

# Set up temperature and velocities
temp = 298.0u"K"
vel = random_velocities!(sys, temp)

# Add loggers
sys = Molly.System(sys;
                   loggers=(temp=TemperatureLogger(100), )  # add more loggers here
)

# Set up simulator
simulator = VelocityVerlet(; dt=0.5u"fs", 
                             coupling=AndersenThermostat(temp, 0.5u"ps"), )

# Perform MD for 1000 step
simulate!(sys, simulator, 1000)
@show sys.loggers.temp.history
