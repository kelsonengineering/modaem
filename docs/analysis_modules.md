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

## Streamline tracing (module TR0)

5.8 Module TR0 - Trace The TR0 directive instructs ModAEM to enter the trace module, which is used to trace 2-D streamlines. The TR0 directive must have a matching END directive. Within the TR0 module, the following directives are valid:

WIN - Set the tracing window. Default tuning parameters are derived from the window size.

TUN - Set tuning parameters Sets tuning parameters for the tracing algorithm. Usage:

TUN step prox frac small step The base step size prox The proximity (in terms of the current step size) to boundary conditions for reducing the step size frac The factor for step size reductions small Smallest allowable step size

TIM - Specify maximum time allowed for particle tracing POI - Release a single particle at the specified location LIN - Release particles along a line N particles along a line

GRI - Release a grid of particles in the sub-window WL0 - Release N particles in reverse from the well bore of a WL0 (discharge-specified well) element. WL1 - Release N particles in reverse from the well bore of a WL1 (head-specified well) element.
