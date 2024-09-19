<h3>Input files</h3>

All files and data required to perform this tutorial are on the [GitHub repository](https://github.com/blake-armstrong/binding-tutorial-PLUMED). They can either be accessed by cloning the repository to your machine;

```
git clone https://github.com/blake-armstrong/binding-tutorial-PLUMED
```

or by following the download instructions in the `README.md` page on the repository to download the basic ZIP file.

Once downloaded, the input files required to run a simulation with OpenMM and PLUMED are provided in the `simulation_files` directory. Included in the `scripts` directory are bash scripts that will automate setting up the multiple-walker simulation files and will analyse the output for you. These scripts are not strictly necessary to perform this tutorial, as the actions they perform can be carried out manually by the user. All the data has already been produced and analysed, and is presented in the 'data' directory. We will discuss how that data was produced and analysed so you can do it all yourself.

<h3>OpennMM & PLUMED Installation</h3>

The first thing we need to do is get a working installation of OpenMM and openmmplumed. [Comprehensive installation instructions for OpenMM](http://docs.openmm.org/latest/userguide/application/01_getting_started.html). We recommend using `conda`, as it is the most straight forward way to install the software without needing to worry about dependencies. See [instructions for installing conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html). After testing the installation with:

```
python -m openmm.testInstallation
```

All platforms available to you will be shown. In this tutorial we use the CPU platform, and the input files provided set up a simulation to use the CPU platform, which every user should have available by default. The [instructions to install the OpenMM PLUMED plugin](https://github.com/openmm/openmm-plumed) are also very simple when using conda. This will also install the command line plumed tool which will be used for analysis later on.

<h3>Running a simulation</h3>

Running the bash script `setup_walkers.sh`

```
bash scripts/setup_walkers.sh
```

or

```
bash scripts/setup_walkers.sh NUMWALKERS
```

will create a set of directory trees within the `data` directory labelled with cylinder radii ranging from 0.0 nm to 0.5 nm. `NUMWALKERS` is an integer for the number of walkers desired with the default being 4 (NB this can be changed but must be greater than 1 for multiple walker metadynamics). These directories are where the multiple walker simulations will be run for systems with various radii for the flat-bottom cylindrical restraining potential. This setup script has been included so that the user can play around with different numbers of walkers in an easy way, if they so desire. By default, it will create directories and setups for four walkers. The details on how to run a single walker are discussed further later.

If we navigate to the directory for a radius of 0.0 nm,

```
cd data/radius_0.0
```

here we will find a `run.sh` script here that will handle submitting the run.py job for each walker. The `run.py` file (which is located in the `simulation_files` directory) will allow us to run our simulation using OpenMM through its Python API. The `plumed.inp` file which contains the PLUMED specifics. There is an 'analysis.sh' script here as well, which will handle the analysis of our simulation output, which is discussed [here](analysis.md).

There is [extensive documentation for using the OpenMM Python API](http://docs.openmm.org/latest/userguide/application/02_running_sims.html) which we direct you to if you cannot understand some component of the `run.py` file while we briefly discuss it here.

<h4>run.py</h4>

```
import openmm as mm
import openmm.app as app
import openmm.unit as unit
from openmmplumed import PlumedForce
```

These lines give us access to the OpenMM code, and the OpenMM PLUMED plugin that can parse a PLUMED input file, interpret it, and provide the instructions to OpenMM.

```
coordinates = "coord.pdb"           # Name of our coordinates file.    
forcefield = "argon.xml"            # Name of our force field file.
plumed_file = "plumed.inp"          # Name of our PLUMED input file.
temperature = 300 * unit.kelvin     # Simulation temperature.
timestep = 1e-3 * unit.picosecond   # Simulation timestep.
trel = 1 / unit.picosecond          # Thermostat relaxation constant for the Langevin thermostat.
thermo_file = "output.out"          # Name of the file to write our periodic simulation output to.
nthermo = 10000                     # How often to periodicly write simulation output (in steps).
traj_file = "trajectory.dcd"        # Name of the file to periodicly write our system coordinates to.
ntraj = 10000                       # How often to periodicly write our system coordinates (in steps).
nsteps = 2000000                    # Total number of steps to simulate for (this corresponds to 2 ns).
```

Simulation parameter definitions that will be used further in the file. When providing a unit with a dimension, we make use of the OpenMM unit object. By default, each walker will run for 2 ns which should take 5 - 10 minutes to complete. This should be enough to produce a PMF that is sufficiently converged to analyse for our purposes due to the simplicity of the model system.


```
pdb = app.PDBFile(coordinates)
forcefield = app.ForceField(forcefield)
platform = mm.Platform.getPlatformByName("CPU")
properties = {}

system = forcefield.createSystem(
    pdb.topology,
    nonbondedMethod=app.CutoffNonPeriodic,      # Run a non-periodic simulation.
    nonbondedCutoff=0.9 * unit.nanometer,
    useDispersionCorrection=False,
    constraints=None,
    rigidWater=False,
)
```

The details of this section are better explained in the [OpenMM docs](http://docs.openmm.org/latest/userguide/application/02_running_sims.html).


```
for n, force in enumerate(system.getForces()):
    if force.getName() == "CMMotionRemover":
        system.removeForce(n)
        break
```
This section is needed to prevent the thermal energy from being removed from the system as we have only one particle. It is not important to understand, as it would not generally be necessary for a real simulation.

```
with open(plumed_file, "r") as f:
    script = f.read()
system.addForce(PlumedForce(script))
```
This is where a PLUMED input script is read and our forces will be created and added to the system.

```
integrator = mm.LangevinMiddleIntegrator(temperature, trel, timestep)

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
```

The details of this section are better explained in the [OpenMM docs](http://docs.openmm.org/latest/userguide/application/02_running_sims.html).

<h4>plumed.inp</h4>

Here we will discuss all of the sections of the PLUMED input file (plumed.inp).


```plumed
RESTART
```
Required for multiple walkers. Ensures each walker reads HILLS from all walkers when initialised.

```plumed
UNITS ENERGY=kj/mol LENGTH=nm TIME=ps
```
Simulation units. Same as OpenMM units for consistency. Not a necessity.

```plumed
fixed1: FIXEDATOM AT=2.50,2.50,2.50
d1: DISTANCE ATOMS=fixed1,1 NOPBC COMPONENTS
```
Here we define our z-distance that we will use for the collective variable. We define a fixed point at the middle of our cubic simulation box (50 Å in each direction) and calculate the distance between the fixed point and our one particle, and save it to the `d1` variable.

```plumed
#HIDDEN
fixed1: FIXEDATOM AT=2.50,2.50,2.50
d1: DISTANCE ATOMS=fixed1,1 NOPBC COMPONENTS
#ENDHIDDEN
rd: CUSTOM ...

   ARG=d1.x,d1.y
   VAR=x,y
   FUNC=0.5*14473*(max(0,sqrt(x^2+y^2)-0.0))^2
   PERIODIC=NO
...

rdb: BIASVALUE ARG=rd
```
Here we use the x and y components of the d1 variable to create `rd`, the flat-bottomed cylindrical restraining potential. We then tell PLUMED to actually use the potential to apply forces throughout the simulation, and not just record its value.

```plumed
#HIDDEN
fixed1: FIXEDATOM AT=2.50,2.50,2.50
d1: DISTANCE ATOMS=fixed1,1 NOPBC COMPONENTS
#ENDHIDDEN
uwall: UPPER_WALLS ARG=d1.z  AT=1.500 KAPPA=14473
lwall: LOWER_WALLS ARG=d1.z AT=-0.125 KAPPA=14473
```

Here we define upper and lower harmonic walls in the z-dimension. The lower wall emulates the repulsion of a surface, and starts acting on the particle at -0.125 nm in the z direction. The upper wall prevents the binding atom from going too far from the `surface` in the z-direction, and starts acting at 1.5 nm.


```plumed
#HIDDEN
fixed1: FIXEDATOM AT=2.50,2.50,2.50
d1: DISTANCE ATOMS=fixed1,1 NOPBC COMPONENTS
#ENDHIDDEN
fixed2: FIXEDATOM AT=2.50,2.50,2.80
d2: DISTANCE ATOMS=fixed2,1 NOPBC COMPONENTS

rd1: CUSTOM ...

   ARG=d1.x,d1.y,d1.z
   VAR=ax,ay,az
   FUNC=sqrt(ax^2+ay^2+az^2)
   PERIODIC=NO
...

rd2: CUSTOM ...

   ARG=d2.x,d2.y,d2.z
   VAR=bx,by,bz
   FUNC=sqrt(bx^2+by^2+bz^2)
   PERIODIC=NO
...

g: CUSTOM ...

   ARG=rd1,rd2
   VAR=x,y
   FUNC=-30*exp(-((x)^2/(2*(0.1^2))))+15*exp(-((y)^2/(2*(0.1^2))))
   PERIODIC=NO
...

pes: BIASVALUE ARG=g
```

The following is used to create an arbitrary potential energy surface that mimics the binding profile of an atom normal to a surface. There will be a minimum at [2.5, 2.5, 2.5] nm to represent a stable surface site and a maximum at [2.5, 2.5, 2.8] nm to represent a barrier for leaving the surface.


```plumed
PRINT FILE=colvar STRIDE=1000 ARG=*
```
This will output a collective variable file to monitor that everything in the simulation is proceeding as expected.

```plumed
FLUSH STRIDE=1000
```
This will flush the output every 1000 steps, allowing us to monitor that everything within PLUMED is working as it should while the simulation is running.

```plumed
#HIDDEN
fixed1: FIXEDATOM AT=2.50,2.50,2.50
d1: DISTANCE ATOMS=fixed1,1 NOPBC COMPONENTS
#ENDHIDDEN
#SOLUTIONFILE=simulation_files/plumed_example.inp
METAD ...
  ARG=d1.z                      # Active CV.
  HEIGHT=2.5                    # Gaussian height. Approximately kBT in this case.
  SIGMA=0.01                    # gaussian width.
  PACE=1000                     # Gaussian deposition pace.
  BIASFACTOR=18 TEMP=300        # Well-tempered parameters.
  GRID_MIN=-0.5                 # Grid parameters for d1.z CV.
  GRID_MAX=1.8                  #
  WALKERS_DIR=../BIAS           # Multiple walkers directory for reading/writing BIAS.
  WALKERS_RSTRIDE=1000          # Multiple walkers reading pace.
  WALKERS_ID=REPLACEWALKERID    # ID of current walker.
  WALKERS_N=REPLACETOTALWALKERS # Total number of walkers.
... METAD
```

This tells PLUMED to run multiple walker metadynamics with a collective variable defined as the z-component of the `d1` distance. WALKERS_ID and WALKERS_N are set to placeholder names which will get replaced when submitting the job through the `run.sh` script. While the multiple walkers run, the HILLS files will be written to and read from the 'BIAS' directory. If you wish to run a single walker simulation, the `RESTART` keyword should be removed, alongside the last four lines inside the `METAD` section. The `HILLS` file will then be written inside the current directory instead of a dedicated `BIAS` directory. The analysis scripts will then need to be modified accordingly.

<h4>run.sh</h4>

The `run.sh` file will handle submitting all of the walkers for a given simulation.

```
export OPENMM_CPU_THREADS=1
```
The `OPENMM_CPU_THREADS` environmental variable is used for OpenMM to assign the number of threads per job on the CPU. If this isn't set OpenMM will use all of your available CPU cores. In our case where our system only has one atom, it is the most efficient to run with only one thread per walker.

```
pids=""
for w in Walker_*; do
  cd ${w}
  pid=$(python3 ../../scripts/run.py > screen.out & echo $!)
  echo "Running ${w} with PID ${pid}"
  pids="${pids}${pid} "
  cd ..
done
echo "To kill walkers call ' kill -9 ${pids} '"
```

This for loop iterates through every walker directory and runs `python3 ../../simulation_files/run.py` to begin the metadynamics simulation, and submits the job as a background process with `&`. The script saves the process ID and outputs it to the screen afterwards to make it easy to kill the background processes if something goes wrong. On UNIX-based operating systems you can run either the `top` or `htop` command to monitor the background processes.

<h4>Beginning the simulation</h4>

To begin your simulations for a given radius run the command `bash run.sh` from the `data/radius_0.x` directory. If you do not modify the `run.py` file, each walker will run for 100 ns to give a combined time of 400 ns. This amount of time is enough for the PMFs to adequately converge. Feel free to run for longer to get an even smoother profile! Once your simulations have finished running, it is time to [analyse](analysis.md).

---

[Back to binding tutorial main page.](../NAVIGATION.md)
