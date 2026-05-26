# Processing Directives 

This chapter describes processing directives that control the solution process and general imspection of values for ModAEM. Many of the directives described here are most commonly used for program debugging, experimentation and development of example problems. When computing results that are to be used for analysis purposes (e.g. within GUI–based modeling environments), the modules GRI (grid/contour generation), EXT (data extraction), INQ (element results inquiry), and OBS (observation data) [^1]. See Chapter [cha:Analysis-modules] for details about the analysis modules.

Each of these directives is given outside the “AEM section” of the ModAEM script. The general layout of the script is as follows, with the processing directives marked in bold face.

        aem
          aqu...
            ref...
            bdy...
            in0...
              ...
            end
          end
          # Other elements go here
        end
        \textbf{\#processingdirectivesgohere}
        \textbf{sol...}
        \textbf{hea...}
        eod

In nearly every case, the first processing directive issued should be **`sol x`**, where **`x`** is the number of iterations (if a solution is to be loaded, use **`sol 0`** to simply load the saved results file). After the solution is complete, the other directives will be useful.

##### Processing directives for solving and reporting solution results 

`SOL`  
Solves the model, based on the available input and results from a previous solution loaded from disk [^2].

`RPT`  
Generates a report of all solution information in HTML format.

##### Processing directives that retrieve data from the model at specific points

`HEA`  
Reports the head at some location in the model to the error log file.

`POT`  
Reports the complex potential at some location in the model to the error log file.

`DIS`  
Reports the discharge at some location in the model to the error log file.

`VEL`  
Reports the velocity at some location in the model to the error log file.

`RCH`  
Reports the net exfiltration rate at some location in the model to the error log file.

`THK`  
Reports the saturated thickness at some location in the model to the error log file.

##### Processing directives that compute numerical approximations for model testing

`GRA`  
Reports the numerical gradient at some location in the model to the error log file.

`LAP`  
Reports the numerical laplacian at some location in the model to the error log file.

##### Processing directives that compute values for a line segment

`FLO`  
Reports the integrated flow across a line segment in the model to the error log file.

The remainder of this chapter describes the general–purpose processing directives in detail.

## Directive SOL – solve the model 

After a model has been defined (using the AEM module input section), it must be solved prior to performing any analyses. ModAEM uses an iterative solution scheme – at each iteration, the solution is improved based on the previous iteration, including the incorporation of non–linear elements such as resistance line sinks or inhomogeneity boundaries in regions where the flow is unconfined.

##### Usage:

        sol niter relaxation

##### Parameters for the sol directive:

niter  
The number of iterations to be performed. $\{integer\}$\
The number of iterations to be used depends strongly on the problem to be solved. The following list describes the issues that determine the number of iterations needed for model convergence. Note that this list is a rule–of–thumb; the modeler should look closely at the solution to ensure that it has fully converged.

1–2iterations  
Simple models that are linear at all points in the domain. Typically, this means that the flow is confined everywhere and that no resistance line sinks that can be removed from solution (that is, “river” and “drain” line sinks) are present. These problems should give accurate solutions with only 1–2 iterations.

3–8iterations  
More complex models that make use of “river” and “drain” line sinks, or have large areas in the model domain in which the flow is unconfined, but where the aquifer base elevation and thickness are constant everywhere will typically converge in 8 iterations or less.

moreiterations  
Very complex models in which the aquifer base elevation and/or thickness varies and the flow is unconfined, or problems where baseflow routing is in use. These problems may require 10 or more i terations to achieve convergence.

relaxation  
The relaxation factor to be used. $\{integer\}$\
This parameter instructs ModAEM to relax the solution during iterations. If $relaxation=1$, then ModAEM applies all of the computed adjustments in the strength coefficients on each iteration. For $0<relaxation<1$, the value of $relaxation$ is multiplied by the strenght adjustments. In some very complex models, this reduces the “stiffness” of the solution algorithm and may reduce the possibility of oscillations during iterations. In nearly all cases, this parameter should be set to 1.0.

### Loading a previous solution

In ModAEM-1.8 and later versions, the modeler has the option of loading a previous solution. This is an advantage for large complex models; the previous solution may be loaded prior to post-processing and analysis, e.g. particle tracking or computing grids of heads. The problem definition, i.e. the aquifer definition and creation of all elements must be complete prior to loading the results.

Because the previous solution is loaded in response to the `sol` directive a special version of the directive, **`sol 0`**, is provided. By issuing the `sol` directive with no iterations, the previous solution will be loaded, and all of ModAEM’s internal data structures will be restored, but no addition solution step will be performed. Whenever a previous solution is to be used, the **`sol 0`** directive should be the first directive in the ModAEM script after the AEM section is complete.

### Saving a solution for future re–use

If the modeler provides a name for a “solution save file” on line 3 of the ModAEM name file, ModAEM stores the solution there as the final step in the “post–solve” procedure. No directives are required in the ModAEM script file.

## Directive RPT – report the solution in HTML format

After the solution is complete (see directive `SOL` above),

## Directives that compute analytic values at a single point

The following procesing directives compute a value and report it in a human–readable format to the run–time message file.

### HEA - report the modeled head

Reports the potentiometric head at a specified point to the message file and to the console. A solution must be present (see the `sol` directive) prior to issuing this directive.

##### Usage:

    hea(x,y)

##### Parameters for the `hea` directive 

`(x,y)`  
The desired location in the complex plane $x+iy$. $\{complex\}$

##### Example:

The head at the location $(100,100)$ is reported in response to the following directive:

    HEA(100.0,100.0)

### POT - report the modeled complex potential

Reports the complex potential at a specified point to the message file and to the console. The complex potential is defined as $\Omega=\Phi+i\Psi$, where $\Phi$ is the discharge potential and $\Psi$ is the streamfunction [^3]. A solution must be present (see the `sol` directive) prior to issuing this directive.

##### Usage:

    pot(x,y)

##### Parameters for the `pot` directive 

`(x,y)`  
The desired location in the complex plane $x+iy$. $\{complex\}$

##### Example

The complex potential at the location $(100,100)$ is reported in response to the following directive:

    pot(100.0,100.0)

### DIS - report the total aquifer discharge 

Reports the total aquifer discharge $Q_{x}+iQ_{y}$ at a specified point to the message file and to the console. The total discharge is a two-dimensional analogue for the specific discharge:
``` math
Q_{i}=Q_{x}+iQ_{y}=\int_{z_{bot}}^{z_{top}}(q_{x}+iq_{y})dz
```
where $q_{x}(z)$ and $q_{y}(z)$ are the horizontal components of the specific discharge $q_{i}=-k\partial_{i}H$ at the elevation $z$, and $z_{top}$ and $z_{bot}$ are the elevations of the aquifer top and bottom, respectively [^4]. A solution must be present (see the `sol` directive) prior to issing this directive.

##### Usage:

    dis(x,y)

##### Parameters for the `dis` directive 

`(x,y)`  
The desired location in the complex plane $x+iy$. $\{complex\}$

##### Example

The complex discharge at the location $(100,100)$ is reported in response to the following directive:

    dis(100.0,100.0)

### VEL - report the horizontal groundwater velocity 

Reports the horizontal aquifer velocity $v_{x}+iv_{y}$ at a specified point to the message file and to the console. The velocity should be interpreted as the *vertically–averaged* velocity at the point in question:
``` math
\bar{v}_{i}=\bar{v}_{x}+i\bar{v}_{y}=\frac{Q_{i}}{hn_{e}}
```
where $\bar{v}_{x}(z)$ and $\bar{v}_{y}$ are the reported horizontal components of the average velocity, $h$ is the saturated thickness, and $n_{e}$ is the effective porosity. A solution must be present (see the `sol` directive) prior to issing this directive.

##### Usage:

    vel(x,y)

##### Parameters for the `vel` directive 

`(x,y)`  
The desired location in the complex plane $x+iy$. $\{complex\}$

##### Example

The average groundwater velocity at the location $(100,100)$ is reported in response to the following directive:

    vel(100.0,100.0)

### RCH - report the net recharge rate

Reports the net rate of areal recharge at a specified point to the message file and to the console. The reported recharge rate has a sign consistent with the ModAEM conventions as described in e.g., Sections[sec:as0-module] and [sec:pd0-module]. The reported rate is negative if the net recharge rate injects water into the aquifer, positive if water is removed. A solution must be present (see the `sol` directive) prior to issuing this directive.

##### Usage:

    rch(x,y)

##### Parameters for the `rch` directive 

`(x,y)`  
The desired location in the complex plane $x+iy$. $\{complex\}$

##### Example

The net recharge rate (negative for recharge) at the location $(100,100)$ is reported in response to the following directive:

    rch(100.0,100.0)

## Directives that compute numerical approximations for testing

Two directives are provided that approximately compute the gradient in the potential and laplacian of the potential, for the purpose of testing the analytic functions that underlie all ModAEM computations.

### GRA - report the modeled numerical gradient in potential

Reports the *numerical* gradient in the discharge potential) at a specified point to the message file and console. This directive is commonly used in program debugging; the numerical gradient should have approximately the same value as the total discharge (see the `dis` directive in Section [DIS_directive]). A solution must be present (see the `sol` directive) prior to issing this directive.

The aproximate gradient $(\hat{q}_{x},\hat{q}_{y})$ is computed at $(x,y)$, a point in the 2D plane, using a spacing $\delta$ provided by the user. The calculations are performed according to the following expression.

``` math
\begin{eqnarray*}
\hat{q}_{x} & = & \frac{1}{\delta}\left[\Phi(x+\delta,y)-\Phi(x-\delta,y)\right]\\
\hat{q}_{y} & = & \frac{1}{\delta}\left[\Phi(x,y+\delta)-\Phi(x,y-\delta)\right]
\end{eqnarray*}
```
in general, it is expected that a more accurate estimate of the gradient is obtained by using a small value of $\delta$.

##### Usage:

    gra z delta

Reports the numerical gradient in the potential at the complex coordinate, $z=(x,y)$, using the spacing $\delta$. Note that in Fortran free-format reads, the two parts of the complex coordinate are provided as $(x,y)$ pairs.

##### Example:

The numerical gradient at the coordinate $(100,100)$, using a spacing of $1.0$ is reported in response to the following directive:

    gra (100.0,100.0) 1.0

## Directives that compute a net analytic value for a line segment

ModAEM provides a directive that is used to report model results along a line segment. Future versions of ModAEM may include additional directives for line segments.

### FLO - report the total flow across a path

Directs ModAEM to report the total integrated groundwater flow across a linear path connecting two specified points to the message file (.err file) and console. The value is computed as

``` math
Q(z_{1},z_{2})=\int_{z_{1}}^{z_{2}}Q_{n}ds
```
where

$z_{1}$ and $z_{2}$  
are the points at the end-points of the line segment, $z_{j}=x_{j}+iy_{j}$

$Q_{n}$  
is the discharge normal to the line segment $z_{1}z_{2}$

$s$  
is the direction along the line segment $z_{1}z_{2}$

A solution must be present (see the `sol` directive) prior to issing this directive.

##### Usage:

    flo (x1,y1) (x2,y2)

Reports the integrated groundwater flux across the line-segment connecting $(x_{1},y_{1})$ and $(x_{2,}y_{2})$ in units of $L^{3}/T$. The sign of the result is determined by the right-hand rule (Section [sec:The-right--hand-rule]).

##### Example:

The total integrated flux across the line segment containing $(50,50)$ and $(100,100)$ is reported in response to the directive

    flo (50.0,50.0) (100.0,100.0)

## Extracting model results in machine-readable format (module INQ)

Module INQ provides the ability to extract results from the model in a machine-readable format that is convenient for post-processing tools and graphical user interface (GUI) programs.

This section needs to be completed

[^1]: Available in ModAEM-1.4.1 and later versions.

[^2]: Available in ModAEM-1.4.1 and later versions.

[^3]: Note that the streamfunction does not exist for problems that include recharge; the reported streamfunction value is not useful in regions where a recharge element (see module AS0 or PD0) is present.

[^4]: Since ModAEM is a 2–D code, it does not explicitly compute the integral. The value is computed analytically by differentiating the discharge potential $\Phi$ at the coordinate $x+iy$. See e.g. Strack (1989) or Haitjema (1995) for details. A detailed description of the formulation of ModAEM will be included in a future version of this book.
