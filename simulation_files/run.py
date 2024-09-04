import openmm as mm
import openmm.app as app
import openmm.unit as unit
from openmmplumed import PlumedForce

from sys import stdout

import time
import random

random.seed(int(time.time()))

coordinates = "coord.pdb"
forcefield = "argon.xml"
plumed_file = "plumed.inp"
temperature = 300 * unit.kelvin
pressure = 1 * unit.atmosphere
timestep = 1e-3 * unit.picosecond
trel = 1 / unit.picosecond
taut = 1 * unit.picosecond
thermo_file = "output.out"
nthermo = 1000
traj_file = "trajectory.dcd"
ntraj = 10000
nsteps = 8000000


pdb = app.PDBFile(coordinates)
forcefield = app.ForceField(forcefield)
platform = mm.Platform.getPlatformByName("CPU")
# properties = {"Precision": "mixed"}
properties = {}

# General system settings
system = forcefield.createSystem(
    pdb.topology,
    nonbondedMethod=app.CutoffNonPeriodic,
    nonbondedCutoff=0.9 * unit.nanometer,
    useDispersionCorrection=False,
    constraints=None,
    rigidWater=False,
)


for n, force in enumerate(system.getForces()):
    if force.getName() == "CMMotionRemover":
        system.removeForce(n)
        break


with open(plumed_file, "r") as f:
    script = f.read()
system.addForce(PlumedForce(script))

# Ensemble
integrator = mm.LangevinMiddleIntegrator(temperature, trel, timestep)

# Create simulation object
simulation = app.Simulation(pdb.topology, system, integrator, platform, properties)

simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(temperature, random.randrange(99999))

simulation.reporters.append(
    app.StateDataReporter(
        thermo_file,
        nthermo,
        separator="\t",
        step=False,
        time=True,
        potentialEnergy=True,
        kineticEnergy=False,
        totalEnergy=False,
        temperature=True,
        volume=True,
        density=True,
        progress=False,
        remainingTime=False,
        speed=True,
        elapsedTime=False,
    )
)

simulation.reporters.append(app.DCDReporter(traj_file, ntraj, enforcePeriodicBox=False))

# Screen output
simulation.reporters.append(
    app.StateDataReporter(
        stdout,
        int(nsteps / 50),
        separator="\t",
        step=False,
        time=False,
        potentialEnergy=True,
        kineticEnergy=False,
        totalEnergy=False,
        temperature=True,
        volume=True,
        density=False,
        progress=True,
        remainingTime=True,
        speed=True,
        elapsedTime=True,
        totalSteps=nsteps,
    )
)


simulation.step(nsteps)
