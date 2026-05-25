# Analysis modules

## Generating grids of model results (module GRI)[sec:gri-module]

Module GRI provides a facility for the construction of grids (e.g. for contour plotting) in a format compatible with SURFER, Matlab, or other software packages. Within the grid module, the window to be gridded must be specified, along with the number of points on the long axis of the grid.

5.7.1 Selecting the output file type (opt directive) The opt directive instructs the grid module which type of output file to create. Usage:

opt grid-type The opt directive expects one parameter, as follows. grid-type Choose surfer for an ASCII SURFER-compatible grid (with the extension .grd) and matlab for

an ASCII MATLABTM-compatible grid (with the extension .m). If the OPT directive is omitted, the grid-type will default to surfer.

Example To select a MatlabTM-compatible output file, issue the directive

opt matlab

5.7.2 Defining the grid window (win directive) Defines the window to be gridded. Usage :

win (x1,y1) (x2,y2) The win directive expects the following parameters. (x1,y1) The lower-left corner of the (rectangular) region to be gridded. (x2,y2) The upper-right corner of the (rectangular) region to be gridded.

Example

win (-100.0,-100.0) (100.0,100.0)

sets the lower-left and upper-right corners of the window for the GRI module at the coordinates \Gamma \Delta \Gamma 100\Delta \Delta \Gamma 100\Theta and\Gamma

100\Delta 100\Theta , respectively.

5.7.3 Choosing the grid resolution (dim directive) Sets the number of grid points along the long axis of the window. Usage:

dim npts The dim directive expects one parameter as follows. npts The number of evenly-spaced points to compute along the long axis of the rectangular grid region. Module GRI

will compute an appropriate number of points along the short axis to ensure that the grid has regular spacing in both directions.

Example To make grid(s) with a resolution of 50 points along the long axis, issue the directive

gri 50

5.7.4 Computing a grid and writing it to a file (directives hea, pot, psi, q_x, and q_y) Once the grid type, grid region, and grid resolution are specified (see the directives opt, dim, and win above), grid files may be computed for a variety of model output values. Currently the following directives are available:

HEA - Create a grid of heads Usage:

hea base-filename generates a grid of the potentiometric head f on the file base-filename_head.grd (or base-filename_head.m)

POT - Create a grid of potentials Usage:

pot base-filename generates a grid of the discharge potential \Phi on the file base-filename_potential.grd (or base-filename_potential.m)

PSI - Create a grid of stream functions Usage:

psi base-filename generates a grid of the streamfunction \Psi on the file base-filename_psi.grd (or base-filename_psi.m)

Q_X - Create a grid of discharges in the x-direction Usage:

q_x base-filename generates a grid of the total aquifer discharge in the x-direction potentiometric heads on the file base-filename_qx.grd (or base-filename_qx.m)

Q_Y - Create a grid of discharges in the y-direction Usage:

q_y base-filename generates a grid of the total aquifer discharge in the y-direction on the file base-filename_qy.grd (or base-filename_qy.m)

## Data extraction (module INQ)

Module INQ provides a facility for extracting computed values from a ModAEM solution and writing them to output files. Two output formats are supported:

- **Traditional (`opt inq`)** — all output is written to a single file opened by the `fil` directive. This is the default.
- **CSV (`opt csv`)** — each directive type writes to its own comma-separated value file, suitable for direct import into Python (`pandas.read_csv()`), R, Excel, or other tools.

The INQ module is entered with the `inq` directive and terminated with `end`.

### Selecting the output format (`opt` directive)

```
opt inq | csv
```

- `opt inq` — write all output to the single file named by `fil` (default behaviour).
- `opt csv` — write each directive type to a separate CSV file. The `fil` directive sets a file-name prefix; if omitted, files are named `<type>_inq.csv`.

The `opt` directive should appear before `fil` and before any extraction directives.

### Naming the output file (`fil` directive)

```
fil filename
```

In **traditional mode** (`opt inq`), `fil` opens `filename` for writing and all subsequent extraction output is written to it. If a file was previously open, it is closed first.

In **CSV mode** (`opt csv`), `fil` sets a prefix for all CSV file names. For example, `fil mymodel` causes extraction directives to write to files named `mymodel_hea.csv`, `mymodel_ls1.csv`, etc. If `fil` is omitted in CSV mode, files are named `hea_inq.csv`, `ls1_inq.csv`, etc.

### Spatial extraction directives

These directives evaluate the solution at one or more points or paths. In traditional mode all output goes to the file opened by `fil`. In CSV mode each directive type appends to its own CSV file; a header row is written the first time a given type is used.

#### `hea` — extract head at a point

```
hea (x,y) label
```

Computes the potentiometric head at `(x,y)`. `label` is an integer that is echoed in the output to identify the record.

**CSV file:** `hea.csv` (columns: `label, x, y, head`)

#### `dis` — extract discharge at a point

```
dis (x,y) label
```

Computes the depth-integrated discharge vector at `(x,y)`.

**CSV file:** `dis.csv` (columns: `label, x, y, discharge_x, discharge_y`)

#### `vel` — extract velocity at a point

```
vel (x,y) label
```

Computes the seepage velocity vector at `(x,y)`.

**CSV file:** `vel.csv` (columns: `label, x, y, velocity_x, velocity_y`)

#### `pot` — extract discharge potential at a point

```
pot (x,y) label
```

Computes the discharge potential at `(x,y)`.

**CSV file:** `pot.csv` (columns: `label, x, y, potential_r, potential_i`)

#### `flo` — extract integrated flow across a path

```
flo (x1,y1) (x2,y2) label
```

Computes the total volumetric flow rate crossing the line segment from `(x1,y1)` to `(x2,y2)`.

**CSV file:** `flo.csv` (columns: `label, x1, y1, x2, y2, flow`)

#### `bdy` — extract head and flux across a boundary segment

```
bdy (x1,y1) (x2,y2) label
```

Computes the head at the midpoint of the segment and the integrated flow crossing it.

**CSV file:** `bdy.csv` (columns: `label, x1, y1, x2, y2, head, flux`)

#### `rch` — extract net recharge rate at a point

```
rch (x,y) label
```

**CSV file:** `rch.csv` (columns: `label, x, y, recharge`)

#### `sat` — extract saturated thickness at a point

```
sat (x,y) label
```

**CSV file:** `sat.csv` (columns: `label, x, y, sat_thickness`)

#### `gag` — extract routed flow at a stream gage

```
gag gage-id label
```

Returns the routed flow at the downstream end of the LS2 string identified by `gage-id`.

**CSV file:** `gag.csv` (columns: `gage_id, label, flow`)

### Module inquiry directives

These directives write a table of element-level check values for a given module. In traditional mode, output goes to the file opened by `fil`. In CSV mode, each directive writes a separate file; a header row listing `tag` and all column names is written at the top.

The first column in every module CSV row is `tag`, which carries the record-type identifier (e.g. `LS1`, `WL0`). This allows modules that write multiple record types (such as AQU and CW0) to be distinguished within a single file.

#### `wl0` / `wl1` — specified-discharge / specified-head wells

**CSV files:** `wl0.csv`, `wl1.csv`

#### `ls0` / `ls1` / `ls2` — line-sink packages

**CSV files:** `ls0.csv`, `ls1.csv`, `ls2.csv`

#### `hb0` — no-flow barrier strings

**CSV file:** `hb0.csv`

#### `as0` — area-sink package

**CSV file:** `as0.csv`

#### `aqu` — aquifer reference and boundary elements

Writes up to three sections: reference point (`AQU` records), boundary elements (`BDY` records), and inhomogeneity interface elements (`IN0` records). In CSV mode the sections are separated by a blank line; each section begins with its own header row.

**CSV file:** `aqu.csv`

#### `cw0` — collector wells

Writes two sections: well summary (`CW0` records) and radial arm elements (`CWR` records). In CSV mode the sections are separated by a blank line.

**CSV file:** `cw0.csv`

### Example — traditional mode

```
inq
  fil mymodel.inq
  hea (0.0,0.0) 1
  hea (100.0,0.0) 2
  flo (-50.0,-100.0) (-50.0,100.0) 10
  ls1
  aqu
end
```

### Example — CSV mode

```
inq
  opt csv
  fil mymodel
  hea (0.0,0.0) 1
  hea (100.0,0.0) 2
  flo (-50.0,-100.0) (-50.0,100.0) 10
  ls1
  aqu
end
```

This writes `mymodel_hea.csv`, `mymodel_flo.csv`, `mymodel_ls1.csv`, and `mymodel_aqu.csv`. Reading results in Python:

```python
import pandas as pd
heads  = pd.read_csv("mymodel_hea.csv", skipinitialspace=True)
ls1    = pd.read_csv("mymodel_ls1.csv", skipinitialspace=True)
```

## Streamline tracing (module TR0)

5.8 Module TR0 - Trace The TR0 directive instructs ModAEM to enter the trace module, which is used to trace 2-D streamlines. The TR0 directive must have a matching END directive. Within the TR0 module, the following directives are valid:

WIN - Set the tracing window. Default tuning parameters are derived from the window size.

TUN - Set tuning parameters Sets tuning parameters for the tracing algorithm. Usage:

TUN step prox frac small step The base step size prox The proximity (in terms of the current step size) to boundary conditions for reducing the step size frac The factor for step size reductions small Smallest allowable step size

TIM - Specify maximum time allowed for particle tracing POI - Release a single particle at the specified location LIN - Release particles along a line N particles along a line

GRI - Release a grid of particles in the sub-window WL0 - Release N particles in reverse from the well bore of a WL0 (discharge-specified well) element. WL1 - Release N particles in reverse from the well bore of a WL1 (head-specified well) element.
