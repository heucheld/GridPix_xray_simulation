# Simulation for IAXO- and CAST-like GridPix X-ray detectors
This [Garfield++](https://garfieldpp.web.cern.ch/garfieldpp/) and
[Degrad](https://degrad.web.cern.ch/degrad/) based tool can simulate events for
X-ray detectors that are based on GridPix. It simulates the absorption of the
X-rays in the gas volume, the track of the photoelectron including its
scattering, the drift of the secondary electron towards the grid if the GridPix
and the gas amplification. The events are stored in the same format as TOS
stores frame data, so that the same analysis tools (like
[TimepixAnalysis](https://github.com/Vindaar/TimepixAnalysis)) can be used.
Additionally, the photon conversion points and the starting points of the
secondary electrons are stored.

## Requirements
To use this tool Garfield++ and Degrad are required:

### Garfield++
Garfield++ can be installed by following the
[installation guide](https://garfieldpp.web.cern.ch/garfieldpp/getting-started).
If you already have an installation of Garfield++ make sure that you use commit
[57ffa3f9](https://gitlab.cern.ch/garfield/garfieldpp/-/commit/57ffa3f9814dd8325de6951320e36cf21d9b232a)
or newer, as there was a bug in the photon transport before.

### Degrad
For Degrad you need to put an executable of it into the folder where you run the
simulation. You can download the source code
[here](https://degrad.web.cern.ch/degrad/) and you should use version 3.15 or
newer. Degrad can be compiled with (for version 3.15)
```
gfortran degrad-3.15.f
```
Depending on the version of the compiler and degrad it can run into an error like this:
```
degrad-3.15.f:2683:55:

 2683 |      /QION,PEQION,EION,NION,QATT,NATT,QNULL,NNULL,SCLN,NC0,EC0,WK,EFL,
      |                                                       1
Error: Rank mismatch in argument 'nco' at (1) (scalar and rank-1)
degrad-3.15.f:2683:59:

 2683 |      /QION,PEQION,EION,NION,QATT,NATT,QNULL,NNULL,SCLN,NC0,EC0,WK,EFL,
      |                                                           1
Error: Rank mismatch in argument 'eco' at (1) (scalar and rank-1)
```
As a workaround add `-fallow-argument-mismatch` for the compiler.

## Usage
To use the simulation tool it must be compiled in a first step. Therefore create
a new folder "build" and change into it:
```
mkdir build
cd build
```
Then cmake is used to generate the makefile (make sure that the compiled Degrad
executable (`degrad.run`) is in this folder):
```
cmake <path to simulation folder>
```
Then the tool need to be compiled:
```
make
```

Now the tool is ready to be used. It has several arguments that define the detector
parameters. The tool can be stated with:
```
./simulation <path> <job> <energy> <gas1> <gas2> <percentage1> <percentage2> <temperature> <pressure> <field>
```
These are the parameters:
- `path` defines the path to a gas file which is used for the drift of the
   secondary electrons. If no path is provided automatically a new gasfile is
   generated at the start. This can then be used in the future for simulations
   with the same gas mixture, the same temperature and pressure and the same
   drift field.
- `job` defines a number for this simulation. It is used as run number in the
   output files.
- `energy` defines the energy of the incoming photons.
- `gas1` defines the name of the first gas component. For possible gases see
   below.
- `gas2` defines the name of the second gas component. For possible gases see
   below.
- `percentage1` defines the percentage of gas1 in the gas mixture.
- `percentage2` defines the percentage of gas2 in the gas mixture.
- `temperature` defines the temperature of the gas in Celsius.
- `pressure` defines the pressure of the gas in Torr.
- `field` defines the driftfield in V/cm.

## Gases
The following gases are supported:

| Name as argument | Secondary name |
| ---------------- | -------------- |
| CF4              |                |
| Ar               | Argon          |
| He               | Helium         |
| He3              | Helium-3       |
| Ne               | Neon           |
| Kr               | Krypton        |
| Xe               | Xenon          |
| CH4              | Methane        |
| Ethane           |                |
| Propane          |                |
| Isobutane        |                |
| CO2              |                |
| Neo-Pentane      |                |
| H2O              | Water          |
| O                | Oxygen         |
| N                | Nitrogen       |
| Nitric Oxide     |                |
| Nitrous Oxide    |                |
| Ethene           |                |
| Acetylene        |                |
| H2               | Hydrogen       |
| D2               | Deuterium      |
| CO               | Carbon Monoxide|
| Methylal         |                |
| DME              |                |
| Reid Step Model  |                |
| Maxwell Model    |                |
| Reid Ramp Model  |                |
| C2F6             |                |
| SF6              |                |
| NH3              | Ammonia        |
| C3H6             | Propene        |
| Cylopropane      |                |
| CH3OH            | Methanol       |
| C2H5OH           | Ethanol        |
| C3H7OH           | Iso-Propanol   |
| Cs               | Cesium         |
| F                | Flourine       |
| CS2              |                |
| COS              |                |
| CD4              |                |
| BF3              | Boron          |
| C2HF5            | C2H2F4         |
| TMA              |                |
| C3H7OH           | N-Propanol     |
| CHF3             |                |
| CF3BR            |                |
| C3F8             |                |
| O3               | Ozone          |
| Hg               | Mercury        |
| H2S              |                |
| N-Butane         |                |
| N-Pentane        |                |
| CCL4             |                |

## Planned for the future
- Possibility to choose between AvalancheMC and AvalancheMicroscopic
- Timepix3/GridPix3 support
- Mone carlo of angular distribution of photoelectrons
- Consideration of z-boost of photoelectrons
- Selectable parameters for the polya of the gas gain monte carlo
- Optional gas gain simulation based on tracking
- Monte carlo of the Grid collection efficiency
- Consideration of absorptions in detector windows
- Support for field maps
