"""modaem.py - Functions for building modaem models. This module provides the reference 
   implementation for simply building the features of a modaem modelgroundwater simulation.
   This code should be the basis for constructing more-sophisticated tools for preprocessing
   and postprocessing. 
   
   This code may also utilized to build validation and test simulations for modaem testing,
   and conceptual models for new modaem users. 

   Modelers are encouraged to utilize Marimo for a testable, notebook-oriented, modeling 
   experience. An example Marimo notebook is included in the distribution.

   Copyright (c) 2026 Vic Kelson, Kelson Environmental LLC
   All Rights Reserved
"""

from math import sqrt
from typing import Iterable, Callable, Generator

class AquReferencePoint:
    ...


class AquDomain:
    ...


class AquString:
    ...


class AquBoundary:
    ...


class Aquifer:
    """Contains the configuration of a modaem2 aquifer."""
    def __init__(self, 
                 bottom: float = 0.0,
                 top: float = 100.0,
                 conductivity: float = 100.0,
                 porosity: float = 0.2,
                 average_head: float = None,):
        assert top > bottom, ValueError("Top elevation must be greater than bottom elevation.")
        assert conductivity > 0, ValueError("Conductivity must be greater than zero.")
        assert porosity > 0, ValueError("Porosity must be greater than zero.")
        assert average_head > bottom, ValueError("Average head must be above the aquifer bottom.")

        self.bottom = bottom
        self.top = top
        self.conductivity = conductivity
        self.porosity = porosity
        # Default the average head to the aquifer top, if it's missing
        if average_head is not None:
            self.average_head = average_head
        else:
            self.average_head = self.bottom + self.top

        self.reference: AquReferencePoint | None = None

        self.domains = list[AquDomain]
        self.strings = list[AquString]
        self.boundary = list[AquBoundary]

    def add_reference(self, ref: AquReferencePoint) -> None:
        self.reference = ref

    def add_domain(self, domain: AquDomain) -> None:
        self.domains.append(domain)

    def add_string(self, string: AquString) -> None:
        self.strings.append(string)

    def to_modaem(self) -> Generator[str]:
        yield (f"aqu {len(self.domains)} {len(self.strings)} {self.bottom} "
               f"     {self.top - self.bottom} {self.conductivity} {self.porosity} {self.average_head}"
               )

        # Use the domain entries to write out the polygons of varying properties
        if len(self.domains) > 0:
            yield "dom"
            for dom in iter(self.domains):
                yield from dom.to_modaem()
            yield "end"

        # Write the aquifer perimeter
        if len(self.perimeter) > 0:
            yield "prm"
            for bdy in iter(self.perimeter):
                yield from bdy.to_modaem()
            yield "end"
        if len(self.strings) > 0:
            yield "str"
            for str in iter(self.strings):
                yield from str.to_modaem()
            yield "end"
        yield "end"




class Wl0Well:
    """Contains a well, using the WL0 package"""
    ...


class Wl0Collection:
    """Contains the collection of wells using the WL0 package"""
    ...


class Pd0Pond:
    """Contains a pond (circular recharge element), using the PD0 package"""
    ...


class Pd0Collection:
    """Contains the collection of ponds using the PD0 package"""
    ...


class Ls3String:
    """Contains a string of linesinks, using the LS3 package"""
    ...


class LS3Collection:
    """Contains the collection of wells using the LS3 package"""
    ...
