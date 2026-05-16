Once the aquifer is defined using module `AQU` (see Chapter [cha:aqu-module]), additional boundary conditions are added to the model by superposition. This chapter describes the elements that are available in ModAEM. Currently, the following element modules are provided:

Discharge–specified wells (`WL0`)  
These are wells for which the pumping rate is known, for example, water supply wells or irrigation wells. This is the most commonly–used type of well; both pumping wells and injection wells are supported. These are analougous to wells created with the MODFLOW WEL package.

Head–specified wells (`WL1`)  
These are wells for which the pumping rate is not *a priori* known to the modeler. The best example is a dewatering well that turns on and off according to a water level measurement. For WL1 elements, the modeler provides the location and radius of the well, plus a location in the aquifer and a head value; ModAEM computes the pumping rate of the well[^1]. These are somewhat analogous to MODFLOW constant–head cells.

Discharge–specified line sinks (`LS0`)  
These are line segments that add or remove a specific amount of water along their length, where the modeler provides the amount of water to be added. Some examples of discharge–specified line sinks are infiltration galleries or rivers in arid climates that are typically dry, but infiltrate water at certain times of the year. These are analogous to a group of wells created with the MODFLOW WEL package.

Head–specified line sinks (`LS1`)  
These are line segments that add or remove water along their lengths, but for which the pumping rate is not *a priori* known to the modeler, and where there is no “entry resistance”, e.g a silty stream bed, between the line sink and the aa group ofquifer. `LS1` line sinks are often used to represent rivers in the far field when the modeler wishes to use an unbounded aquifer with a modeled far field (see Section [sub:reference-flow-field]). These are somewhat analogous to a group of MODFLOW constant–head cells; they differ because in MODFLOW, constant head cells are specified as part of the BAS package, not as a separate component.

Resistance line sinks (`LS2`)  
These are line segments that add or remove water along their lengths, but for which the pumping rate is not *a priori* known to the modeler, and where an “entry resistance”, e.g a silty stream bed, between the line sink and the aquifer is present. `LS2` line sinks are often used to represent surface waters or drains in the near field. `LS2` line sinks may be created as rivers (analogous to the MODFLOW RIV package), drains (analogous to the MODFLOW DRN package), or general–head boundaries (analogous to the MODFLOW GHB package). Streamflow routing may be performed for `LS2` elements (in a manner similar to the MODFLOW STR package), using the analysis module `RT0` (see Section [sec:rt0-module]).

Horizontal no–flow boundaries (`HB0)`  
These are elements that create a linear no–flow condition within the active area of a ModAEM model. These may be used to model sheet pilings, slurry walls, faults, and other linear no–flow conditions that require an active aquifer domain on both sides of the line. This should not be used for bounded models or for “islands” in an aquifer domain; use the `BDY` elements included in module `AQU` (Section [sub:bounded-aquifers]) for details. These are analougous to the MODFLOW HFB package.

Circular area–sinks (`PD0`)  
These are elements that provide an areal infiltration or exfiltration rate over a circular sub–domain, using a circular “pond” function described by Strack (1989). The “sink density”, or rate of infitration per unit of surface area, is specified by the modeler. These are typically used in example problems and for conceptual models, although they are convenient for representing circular irrigators. For most practical modeling applications, module `PD0` is superseded by the polygonal area–sink module `AS0` (Section [sec:as0-module]).

Polygonal area–sinks (AS0)  
These are elements that provide an areal infiltration or exfiltration rate over a polygonal sub–domain. The “sink density”, or rate of infitration per unit of surface area, is specified by the modeler. These are typically used as sources of areal recharge, e.g. from rainfall, or for infiltration galleries.

# Discharge-specified wells (module `WL0`)

Module `WL0` creates discharge-specified wells. ModAEM processing enters module `WL0` in response to the `wl0` directive.

#### Usage:

        wl0 nwells
          (xw,yw) pumping-rate radius id 
             # Optional features follow if desired
             ppw ...
             eff ...
             dhd ...
             ddn ...
          (xw, yw) ...
          ...
        end

#### Parameters for the `wl0` directive

  
The maximum number of wells in the problem

#### Specifying well elements

The `wl0` directive is followed by one record for every well in the model. If more than nwells well elements is provided, an ModAEM will terminate and report the error. Each well record has the following form:

        (xw,yw) pumping-rate radius id

#### Parameters

  
The coordinates of the center of the well.

  
The pumping rate of the well. The value is positive if the element removes water from the aquifer, negative if it adds water to the aquifer [^2].

  
The radius of the well.

  
A unique identification number for the well.

#### Optional `wl0` features

ModAEM offers several optional features for modeling discharge-specified wells. Two

# Head-specified wells (module WL1)

Module WL1 creates head-specified wells. Head-specified wells should be used in cases where a well is to be used to maintain a particular water level at some location in the aquifer. The model computes the pumping rate for each head–specified well element that allows the model to match the specified head; these elements *must not* be used for the purpose of calibrating the model to known heads. ModAEM processing enters module WL1 in response to the wl1 directive.

#### Usage:

        wl1 nwells
          (xw,yw) head (xc,yc) radius id
          ... 
        end

#### Parameters

The maximum number of wells in the problem

## Specifying head–specified well elements

The `wl1` directive is followed by one record for every well in the model. If more than nwells well elements is provided, an ModAEM will terminate and report the error. Each well record has the following form:

        (xw,yw) head (xc,yc) radius id

#### Parameters

  
The coordinates of the center of the well.

  
The specified head at the control point.

  
The coordinates of the point where the head condition is to be met.

  
The radius of the well.

  
A unique (integer) identification number for the well.

# Discharge-specified line-sinks (module LS0) 

Module LS0 creates discharge-specified line-sinks. These are line segments that add or remove a specific amount of water along their length, where the modeler provides the amount of water to be added. Some examples of discharge–specified line sinks are infiltration galleries or rivers in arid climates that are typically dry, but infiltrate water at certain times of the year. These are analogous to a group of wells created with the MODFLOW WEL package.

#### Usage:

        ls0 nstrings
          str npts id
            (x,y) strength 
            ...
          str ...
            ...
        end

#### Parameters for the `ls0` directive

  
The maximum number of discharge-specified line-sink strings in the problem.

Following the `ls0` directive, the user provides strings of line sinks using the `str` directive. A string of line sinks is composed of a list of vertices. One line sink element is created for each pair of consecutive vertices.

#### Creating a string of elements

        str npts id
          (x,y) strength
          ...

#### Parameters

  
The maximum number of vertices in the line-sink string.

  
A unique identification number for the string.

#### Specifying vertices for line–sink strings

Following the `str` directive, two or more data records define the vertices of the line–sink string. The parameters provided for each vertex are as follows.

        (x,y) strength



#### Parameters

Coordinates of the vertex.

The sink density $\sigma$ of the line–sink string at the vertex. The total volumetric infiltration or abstraction along the line will be the value $(\sigma_{1}+\sigma_{2})\frac{L}{2}$ $[L^{3}/T]$ where $\sigma_{1}$ and $\sigma_{2}$ are the strengths at consecutive vertices and $L$ is the distance between the vertices. The sink density is defined as the total extraction rate of the line-sink per unit length. The value is positive if the element removes water from the aquifer, negative if it adds water to the aquifer [^3].



# Head-specified line-sinks (module `LS1`)

Module LS1 creates head–specified line–sink elements. These are line segments that add or remove water along their lengths, but for which the pumping rate is not *a priori* known to the modeler, and where there is no “entry resistance”, e.g a silty stream bed, between the line sink and the aquifer. LS1 line sinks are often used to represent rivers in the far field when the modeler wishes to use an unbounded aquifer with a modeled far field (see Section [sub:reference-flow-field]). These elements are somewhat analogous to a group of MODFLOW constant–head cells; they differ because in MODFLOW, constant head cells are specified as part of the aquifer, not as a separate element.

#### Usage:

        ls1 nstrings
          str npts id
            (x,y) head 
              ...
          str ...
              ...
        end

#### Parameters for the `ls1` directive

  
The maximum number of discharge-specified line–sink strings in the problem.

Following the `ls1` directive, the user provides strings of line sinks using the `str` directive. A string of line sinks is composed of a list of vertices. One line sink element is created for each pair of consecutive vertices.

#### Creating a string of elements

        str npts id
          (x,y) head
          ...

#### Parameters for the `str` directive 

  
The maximum number of vertices in the line-sink string.

  
A unique identification number for the string.

#### Specifying vertices for line–sink strings

Following the `str` directive, two or more data records define the vertices of the line–sink string. The parameters provided for each vertex are as follows.

  
Coordinates of the vertex.

  
The specified head at the center of the line-sink string.

# Line-sinks with entry resistance (module LS2)

Module LS2 creates head-specified line-sinks. These are line segments that add or remove water along their lengths, but for which the pumping rate is not *a priori* known to the modeler, and where an “entry resistance”, e.g a silty stream bed, between the line sink and the aquifer is present. LS2 line sinks are often used to represent surface waters or drains in the near field. LS2 line sinks may be created as rivers (analogous to the MODFLOW RIV package), drains (analogous to the MODFLOW DRN package), or general–head boundaries (analogous to the MODFLOW GHB package). Streamflow routing may be performed for LS2 elements (in a manner similar to the MODFLOW STR package), using the analysis module RT0 (see Section [sec:rt0-module]).

#### Usage:

        ls2 nstrings
          str npts mode conductance id
            (x,y) head bottom
            ...
          str ...
            ...
        end

#### Parameters for the `ls2` directive

  
The maximum number of discharge-specified line–sink strings in the problem.

Following the `ls2` directive, the user provides strings of line sinks using the `str` directive. A string of line sinks is composed of a list of vertices. One line sink element is created for each pair of consecutive vertices.

#### Creating a string of elements

        str npts mode conductance id
          (x,y) head bottom
          ...

#### Parameters for the `str` directive 

  
The maximum number of vertices in the line-sink string.

  
A unique identification number for the string.

  
Defines the behavior of the line sink. The value provided is an integer from the list below:

`general-head boundary (0)`  
The boundary is always active, in a manner analogous with the MODFLOW GHB package.

`river (1)`  
The line sink becomes a discharge–specified feature if the head in the aquifer falls below the bottom of the resistance layer (a “percolating” condition). In this case, the infiltration density is computed as $\sigma=c\times(h_{r}-h_{b})$ where $c$ is the conductance of the line sink (see below), $h_{r}$ is the specified stage in the river, and $h_{b}$is the elevation of the bottom of the resistance layer for this line sink. This is analogous to the MODFLOW RIV package.

`drain (2)`  
The line sink will be removed from the solution with a sink density of zero when the head in the aquifer falls below the bottom of the drain. This is analogous to the MODFLOW DRN package.

  
The conductance for the line sink. The conductance is defined in a manner consistent with MODFLOW; for a river, the conductance is defined to be $c=w\times\frac{k_{c}}{t_{c}}\,[L/T]$, where $c$ is the conductance, $w$ is the width of the stream, $k_{c}$ is the vertical hydraulic conductivity of the resistance layer, and $t_{c}$ is the thickness of the resistance layer.

  
A unique identification number for this string.

#### Specifying vertices for line–sink strings

Following the `str` directive, two or more data records define the vertices of the line–sink string. The parameters provided for each vertex are as follows.

  
Coordinates of the vertex.

  
The specified head at this vertex. The model interpolates the head along the line segment. Since the head is specified at the center of the line segment, the average of the bottom elevation between adjacent vertices is used.

  
real The elevation of the bottom of the resistance layer at this vertex. The model interpolates the bottom elevation along the line segment. Since the head is specified at the center of the line segment, the average of the bottom elevation between adjacent vertices is used.

# No–flow boundary walls (module HB0)

Module HB0 creates no–flow boundary walls (e.g. sheet pilings or slurry walls). These are elements that create a linear no–flow condition within the active area of a ModAEM model. These may be used to model sheet pilings, slurry walls, faults, and other linear no–flow conditions that require an active aquifer domain on both sides of the line. This should not be used for bounded models or for “islands” in an aquifer domain; use the BDY elements included in module AQU (Section [sub:bounded-aquifers]) for details. These are analogous to the MODFLOW HFB package.

#### Usage:

        hb0 nstrings
          str npts id
            (x,y)
            ...
          str ...
            ...
        end

#### Parameters for the `hb0` directive

  
The maximum number of no–flow strings in the problem.

Following the `hb0` directive, the user provides strings of line sinks using the `str` directive. A no–flow string is composed of a list of vertices. One line doublet element is created for each pair of consecutive vertices.

#### Creating a string of elements

        str npts id
          (x,y) strength
          ...

#### Parameters for the str directive 

  
The maximum number of vertices in the line-sink string.

  
A unique identification number for the string.

#### Specifying vertices for line–sink strings

Following the `str` directive, two or more data records define the vertices of the no–flow string. The parameters provided for each vertex are as follows.

  
Coordinates of the vertex.

# Discharge-specified circular area-sinks (module PD0)

Module PD0 creates discharge-specified circular area-sinks (colloquially known as ponds).

#### Usage:

        pd0 nponds
          (xc,yc) strength radius id 
          ... 
        end

#### Parameters for the `pd0` directive

`nwells`  
The maximum number of wells in the problem

#### Specifying pond elements

The `wl0` directive is followed by one record for every pond in the model. If more than `nponds` elements are provided, ModAEM will terminate and report the error. Each pond record has the following parameters:

        (xc,yc) sink-density radius id

  
The coordinates of the center of the pond.

  
The sink density of the pond. This is the value $\gamma=Q_{p}/A_{p}$, where $\gamma\,[L/T]$ is the sink density, $Q_{p}\,[L^{3}/T]$ is the total amount of water infiltrated or abstracted by the pond, and $A_{p}\,[L^{2}]$ is the area of the pond. The value is positive if the element removes water from the aquifer, negative if it adds water to the aquifer [^4]

  
The radius of the pond.

  
A unique identification number for the pond.

# Polygonal area-sinks (module AS0)

Module AS0 creates discharge-specified polygonal area–sink elements. These are elements that provide an areal infiltration or exfiltration rate over a polygonal sub–domain. The “sink density”, or rate of infiltration per unit of surface area, is specified by the modeler. These are typically used as sources of areal recharge, e.g. from rainfall, or for infiltration galleries.

#### Usage:

        as0 top/bottom nareas
          str npts strength id
            (x,y) 
            ...
          str ...
            ...
        end

Please note: In versions of ModAEM prior to version 1.8, the “bottom” and “nareas” parameters were reversed in the AS0 input.

#### Parameters for the `as0` directive

  
An integer that specifies whether the area sinks are to be placed at the aquifer top (**`0`**) or bottom (**`1`**).

  
The maximum number of discharge–specified area–sink polygons in the problem.

Following the `as0` directive, the user provides polygons using the `str` directive. A polygon is composed of a list of vertices. It is not necessary to duplicate the first vertex to close the polygon; ModAEM automatically closes the polygon.

#### Creating an element

        str npts sink-density id
          (x,y) 
          ...

#### Parameters for the str directive 

  
The maximum number of vertices in the line-sink string.

`sink-density`  
The sink density of the area sink. This is the value $\gamma=Q_{p}/A_{p}$, where $\gamma\,[L/T]$ is the sink density, $Q_{p}\,[L^{3}/T]$ is the total amount of water infiltrated or abstracted by the element, and $A_{p}\,[L^{2}]$ is the area of the polygon. The value is positive if the element removes water from the aquifer, negative if it adds water to the aquifer. [^5]

  
A unique identification number for the string.

#### Specifying vertices for polygons

Following the `str` directive, three or more data records define the vertices of the area–sink perimeter. The parameters provided for each vertex are as follows.

  
Coordinates of the vertex.

[^1]: As will be discussed in Section [sec:wl1-module], `WL1` wells should not be confused with calibration targets or inverse models.

[^2]: Users who make use of ModAEM using the GMS preprocessor will note that GMS makes use of the MODFLOW convention that abstraction of water from the aquifer is negative, while injection is positive. Conversion to ModAEM’s convention is handled transparently by GMS.

[^3]: Users who make use of ModAEM using the GMS preprocessor will note that GMS makes use of the MODFLOW convention that abstraction of water from the aquifer is negative, while injection is positive. Conversion to ModAEM’s convention is handled transparently by GMS.

[^4]: Users who make use of ModAEM using the GMS preprocessor will note that GMS makes use of the MODFLOW convention that abstraction of water from the aquifer is negative, while injection is positive. Conversion to ModAEM’s convention is handled transparently by GMS.

[^5]: Users who make use of ModAEM using the GMS preprocessor will note that GMS makes use of the MODFLOW convention that abstraction of water from the aquifer is negative, while injection is positive. Conversion to ModAEM’s convention is handled transparently by GMS.
