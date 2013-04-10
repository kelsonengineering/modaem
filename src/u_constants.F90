module u_constants

  ! ModAEM 1.8
  ! Copyright(c) 1995-2008 WHPA Inc. and Vic Kelson
  !
  ! This program is free software; you can redistribute it and/or
  ! modify it under the terms of the GNU General Public License
  ! as published by the Free Software Foundation; either version 2
  ! of the License, or(at your option) any later version.
  !
  ! This program is distributed in the hope that it will be useful,
  ! but WITHOUT ANY WARRANTY; without even the implied warranty of
  ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ! GNU General Public License for more details.
  !
  ! You should have received a copy of the GNU General Public License
  ! along with this program; if not, write to the Free Software
  ! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
  !
  ! Contact the author by e-mail at: vic@wittmanhydro.com
  ! Or by regular mail at:
  ! WHPA, Inc
  ! 320 W 8th St
  ! Bloomington, IN 47401

  implicit none

  integer, public, parameter :: AE_INT = selected_int_kind(8)
  integer, public, parameter :: AE_REAL = selected_real_kind(10, 8)
  integer, public, parameter :: AE_DBL = selected_real_kind(10, 8)

  ! Progam-defined types for input file interpreter

  type, public :: DIRECTIVE
    integer(kind=AE_INT) :: iOpCode
    character(len=3) :: sText
  end type DIRECTIVE

  ! OpCode named constants
  ! Error and message codes
  ! Data record found
  integer(kind=AE_INT), public, parameter :: kOpData = 0
  ! Error detected
  integer(kind=AE_INT), public, parameter :: kOpError = -1
  ! End-of-file detected
  integer(kind=AE_INT), public, parameter :: kOpFileEOF = -2
  ! Program control
  ! Begin definition of an AEM model
  integer(kind=AE_INT), public, parameter :: kOpAEM = 1001
  ! End of data for module
  integer(kind=AE_INT), public, parameter :: kOpEND = 1002
  ! End of model input data
  integer(kind=AE_INT), public, parameter :: kOpEOD = 1003
  ! Dimensions an element module
  integer(kind=AE_INT), public, parameter :: kOpDIM = 1004
  ! Dimensions a limear element item
  integer(kind=AE_INT), public, parameter :: kOpSTR = 1005
  ! Tell a module to do its thing
  integer(kind=AE_INT), public, parameter :: kOpRUN = 1006
  ! Generate an overall report
  integer(kind=AE_INT), public, parameter :: kOpRPT = 1009
  ! Create a well
  integer(kind=AE_INT), public, parameter :: kOpWEL = 1010
  ! Close a STR object
  integer(kind=AE_INT), public, parameter :: kOpCLO = 1011
  ! Partial-penetration information follows
  integer(kind=AE_INT), public, parameter :: kOpPPW = 1012
  ! Enter well efficiency information
  integer(kind=AE_INT), public, parameter :: kOpEFF = 1013
  ! Desired head in a well
  integer(kind=AE_INT), public, parameter :: kOpDHD = 1014
  ! Drawdown flags follow
  integer(kind=AE_INT), public, parameter :: kOpDDN = 1015
  ! Standard model modules
  ! Reference information follows
  integer(kind=AE_INT), public, parameter :: kOpREF = 2001
  ! Regional aquifer property data follows
  integer(kind=AE_INT), public, parameter :: kOpBDY = 2003
  ! Solution module
  integer(kind=AE_INT), public, parameter :: kOpSOL = 2004
  ! Grid module
  integer(kind=AE_INT), public, parameter :: kOpGRI = 2005
  ! Data inquiry module
  integer(kind=AE_INT), public, parameter :: kOpOBS = 2006
  ! Here is an island in the aquifer
  integer(kind=AE_INT), public, parameter :: kOpISL = 2007
  ! Data extraction module
  integer(kind=AE_INT), public, parameter :: kOpEXT = 2008
  ! Salt-water interface
  integer(kind=AE_INT), public, parameter :: kOpSWI = 2009

  ! Analytic element modules
  ! Aquifer data module
  integer(kind=AE_INT), public, parameter :: kOpAQU = 3000
  ! Strength-specified wells
  integer(kind=AE_INT), public, parameter :: kOpWL0 = 3001
  ! head-specified wells
  integer(kind=AE_INT), public, parameter :: kOpWL1 = 3002
  ! Strength-specified line-sinks
  integer(kind=AE_INT), public, parameter :: kOpLS0 = 3003
  ! head-specified line-sinks
  integer(kind=AE_INT), public, parameter :: kOpLS1 = 3004
  ! No-flow horizontal barrier
  integer(kind=AE_INT), public, parameter :: kOpHB0 = 3005
  ! Inhomogeneity module
  integer(kind=AE_INT), public, parameter :: kOpIN0 = 3007
  ! Inhomogeneity domain
  integer(kind=AE_INT), public, parameter :: kOpDOM = 3008
  ! Discharge-specifed pond module
  integer(kind=AE_INT), public, parameter :: kOpPD0 = 3009
  ! Linesinks with resistance, routing
  integer(kind=AE_INT), public, parameter :: kOpLS2 = 3010
  ! Linesinks with resistance, routing, next generation
  integer(kind=AE_INT), public, parameter :: kOpLS3 = 3011
  ! Strength-specified area-sinks
  integer(kind=AE_INT), public, parameter :: kOpAS0 = 3012
#ifndef __GPL__
  ! Collector wells
  integer(kind=AE_INT), public, parameter :: kOpCW0 = 30013
#endif
  ! Analysis Modules
  ! Pathline trace module
  integer(kind=AE_INT), public, parameter :: kOpTR0 = 3506

  ! Other commands
  ! Set the window
  integer(kind=AE_INT), public, parameter :: kOpWIN = 5001
  ! Select/inquire heads
  integer(kind=AE_INT), public, parameter :: kOpHEA = 5002
  ! Select/inquire potentials
  integer(kind=AE_INT), public, parameter :: kOpPOT = 5003
  ! Select/inquire streamfunction
  integer(kind=AE_INT), public, parameter :: kOpPSI = 5004
  ! Select/inquire discharge
  integer(kind=AE_INT), public, parameter :: kOpDIS = 5005
  ! Select/inquire velocity
  integer(kind=AE_INT), public, parameter :: kOpVEL = 5006
  ! Select/inquire total flow
  integer(kind=AE_INT), public, parameter :: kOpFLO = 5007
  ! Select/inquire gradient in Phi
  integer(kind=AE_INT), public, parameter :: kOpGRA = 5008
  ! Create a Gaussian RBF in LKG
  integer(kind=AE_INT), public, parameter :: kOpGAU = 5009
  ! Create a Bessel RBF in LKG
  integer(kind=AE_INT), public, parameter :: kOpBES = 5010
  ! Read a SURFER grid of leakances in LKG
  integer(kind=AE_INT), public, parameter :: kOpSUR = 5011
  ! Set a layer leakance in LKG
  integer(kind=AE_INT), public, parameter :: kOpLAY = 5012
  ! Grid of leakage errors
  integer(kind=AE_INT), public, parameter :: kOpLER = 5013
  ! Specify an analysis output file
  integer(kind=AE_INT), public, parameter :: kOpFIL = 5014
  ! Transport parameters
  integer(kind=AE_INT), public, parameter :: kOpTRA = 5015
  ! Tuning parameters
  integer(kind=AE_INT), public, parameter :: kOpTUN = 5016
  ! Timing parameters
  integer(kind=AE_INT), public, parameter :: kOpTIM = 5017
  ! Specify a single point
  integer(kind=AE_INT), public, parameter :: kOpPOI = 5018
  ! Specify a line
  integer(kind=AE_INT), public, parameter :: kOpLIN = 5019
  ! Select/inquire Qx
  integer(kind=AE_INT), public, parameter :: kOpQ_X = 5020
  ! Select/inquire Qy
  integer(kind=AE_INT), public, parameter :: kOpQ_Y = 5021
  ! Set a module option
  integer(kind=AE_INT), public, parameter :: kOpOPT = 5022
  ! Set the input file
  integer(kind=AE_INT), public, parameter :: kOpINP = 5023
  ! Set the output file
  integer(kind=AE_INT), public, parameter :: kOpOUT = 5024
  ! Select/inquire recharge
  integer(kind=AE_INT), public, parameter :: kOpRCH = 5025
  ! Select/inquire saturated thickness
  integer(kind=AE_INT), public, parameter :: kOpSAT = 5026
  ! Select/inquire aquifer properties
  integer(kind=AE_INT), public, parameter :: kOpPRO = 5027
  ! Do routing(LS2/LS3)
  integer(kind=AE_INT), public, parameter :: kOpROU = 5028
  ! Gage value(OBS)
  integer(kind=AE_INT), public, parameter :: kOpGAG = 5029
  ! Select/inquire laplacian in Phi
  integer(kind=AE_INT), public, parameter :: kOpLAP = 5030

  ! Pre-packaged standard directives
  type(DIRECTIVE), public, parameter :: dirAEM = DIRECTIVE(kOpAEM, "AEM")
  type(DIRECTIVE), public, parameter :: dirEND = DIRECTIVE(kOpEND, "END")
  type(DIRECTIVE), public, parameter :: dirEOD = DIRECTIVE(kOpEOD, "EOD")
  type(DIRECTIVE), public, parameter :: dirDIM = DIRECTIVE(kOpDIM, "DIM")
  type(DIRECTIVE), public, parameter :: dirSTR = DIRECTIVE(kOpSTR, "STR")
  type(DIRECTIVE), public, parameter :: dirRUN = DIRECTIVE(kOpRUN, "RUN")
  type(DIRECTIVE), public, parameter :: dirRPT = DIRECTIVE(kOpRPT, "RPT")
  type(DIRECTIVE), public, parameter :: dirWEL = DIRECTIVE(kOpWEL, "WEL")
  type(DIRECTIVE), public, parameter :: dirPPW = DIRECTIVE(kOpPPW, "PPW")
  type(DIRECTIVE), public, parameter :: dirEFF = DIRECTIVE(kOpEFF, "EFF")
  type(DIRECTIVE), public, parameter :: dirDHD = DIRECTIVE(kOpDHD, "DHD")
  type(DIRECTIVE), public, parameter :: dirDDN = DIRECTIVE(kOpDDN, "DDN")
  !
  type(DIRECTIVE), public, parameter :: dirREF = DIRECTIVE(kOpREF, "REF")
  type(DIRECTIVE), public, parameter :: dirBDY = DIRECTIVE(kOpBDY, "BDY")
  type(DIRECTIVE), public, parameter :: dirSOL = DIRECTIVE(kOpSOL, "SOL")
  type(DIRECTIVE), public, parameter :: dirGRI = DIRECTIVE(kOpGRI, "GRI")
  type(DIRECTIVE), public, parameter :: dirOBS = DIRECTIVE(kOpOBS, "OBS")
  type(DIRECTIVE), public, parameter :: dirEXT = DIRECTIVE(kOpEXT, "EXT")
  type(DIRECTIVE), public, parameter :: dirTR0 = DIRECTIVE(kOpTR0, "TR0")
  type(DIRECTIVE), public, parameter :: dirISL = DIRECTIVE(kOpISL, "ISL")
  type(DIRECTIVE), public, parameter :: dirSWI = DIRECTIVE(kOpSWI, "SWI")
  !
  type(DIRECTIVE), public, parameter :: dirAQU = DIRECTIVE(kOpAQU, "AQU")
  type(DIRECTIVE), public, parameter :: dirWL0 = DIRECTIVE(kOpWL0, "WL0")
  type(DIRECTIVE), public, parameter :: dirWL1 = DIRECTIVE(kOpWL1, "WL1")
  type(DIRECTIVE), public, parameter :: dirLS0 = DIRECTIVE(kOpLS0, "LS0")
  type(DIRECTIVE), public, parameter :: dirLS1 = DIRECTIVE(kOpLS1, "LS1")
  type(DIRECTIVE), public, parameter :: dirHB0 = DIRECTIVE(kOpHB0, "HB0")
  type(DIRECTIVE), public, parameter :: dirIN0 = DIRECTIVE(kOpIN0, "IN0")
  type(DIRECTIVE), public, parameter :: dirDOM = DIRECTIVE(kOpDOM, "DOM")
  type(DIRECTIVE), public, parameter :: dirPD0 = DIRECTIVE(kOpPD0, "PD0")
  type(DIRECTIVE), public, parameter :: dirLS2 = DIRECTIVE(kOpLS2, "LS2")
  type(DIRECTIVE), public, parameter :: dirLS3 = DIRECTIVE(kOpLS3, "LS3")
  type(DIRECTIVE), public, parameter :: dirAS0 = DIRECTIVE(kOpAS0, "AS0")
#ifndef __GPL__
  type(DIRECTIVE), public, parameter :: dirCW0 = DIRECTIVE(kOpCW0, "CW0")
#endif

  !
  type(DIRECTIVE), public, parameter :: dirWIN = DIRECTIVE(kOpWIN, "WIN")
  type(DIRECTIVE), public, parameter :: dirHEA = DIRECTIVE(kOpHEA, "HEA")
  type(DIRECTIVE), public, parameter :: dirPOT = DIRECTIVE(kOpPOT, "POT")
  type(DIRECTIVE), public, parameter :: dirPSI = DIRECTIVE(kOpPSI, "PSI")
  type(DIRECTIVE), public, parameter :: dirDIS = DIRECTIVE(kOpDIS, "DIS")
  type(DIRECTIVE), public, parameter :: dirVEL = DIRECTIVE(kOpVEL, "VEL")
  type(DIRECTIVE), public, parameter :: dirFLO = DIRECTIVE(kOpFLO, "FLO")
  type(DIRECTIVE), public, parameter :: dirGRA = DIRECTIVE(kOpGRA, "GRA")
  type(DIRECTIVE), public, parameter :: dirLAP = DIRECTIVE(kOpLAP, "LAP")
  type(DIRECTIVE), public, parameter :: dirGAU = DIRECTIVE(kOpGAU, "GAU")
  type(DIRECTIVE), public, parameter :: dirBES = DIRECTIVE(kOpBES, "BES")
  type(DIRECTIVE), public, parameter :: dirSUR = DIRECTIVE(kOpSUR, "SUR")
  type(DIRECTIVE), public, parameter :: dirLAY = DIRECTIVE(kOpLAY, "LAY")
  type(DIRECTIVE), public, parameter :: dirLER = DIRECTIVE(kOpLER, "LER")
  type(DIRECTIVE), public, parameter :: dirFIL = DIRECTIVE(kOpFIL, "FIL")
  type(DIRECTIVE), public, parameter :: dirTRA = DIRECTIVE(kOpTRA, "TRA")
  type(DIRECTIVE), public, parameter :: dirTUN = DIRECTIVE(kOpTUN, "TUN")
  type(DIRECTIVE), public, parameter :: dirTIM = DIRECTIVE(kOpTIM, "TIM")
  type(DIRECTIVE), public, parameter :: dirPOI = DIRECTIVE(kOpPOI, "POI")
  type(DIRECTIVE), public, parameter :: dirLIN = DIRECTIVE(kOpLIN, "LIN")
  type(DIRECTIVE), public, parameter :: dirQ_X = DIRECTIVE(kOpQ_X, "Q_X")
  type(DIRECTIVE), public, parameter :: dirQ_Y = DIRECTIVE(kOpQ_Y, "Q_Y")
  type(DIRECTIVE), public, parameter :: dirOPT = DIRECTIVE(kOpOPT, "OPT")
  type(DIRECTIVE), public, parameter :: dirINP = DIRECTIVE(kOpINP, "INP")
  type(DIRECTIVE), public, parameter :: dirOUT = DIRECTIVE(kOpOUT, "OUT")
  type(DIRECTIVE), public, parameter :: dirRCH = DIRECTIVE(kOpRCH, "RCH")
  type(DIRECTIVE), public, parameter :: dirSAT = DIRECTIVE(kOpSAT, "SAT")
  type(DIRECTIVE), public, parameter :: dirPRO = DIRECTIVE(kOpPRO, "PRO")
  type(DIRECTIVE), public, parameter :: dirROU = DIRECTIVE(kOpROU, "ROU")
  type(DIRECTIVE), public, parameter :: dirGAG = DIRECTIVE(kOpGAG, "GAG")

  ! Other constants...

  ! Error Constants
  integer(kind=AE_INT), public, parameter :: errOK = 0
  integer(kind=AE_INT), public, parameter :: errRunTime = -1
  integer(kind=AE_INT), public, parameter :: errFatal = -2
  integer(kind=AE_INT), public, parameter :: errNonFatal = -3
  integer(kind=AE_INT), public, parameter :: errInternal = -4
  integer(kind=AE_INT), public, parameter :: errAssertion = -5
  integer(kind=AE_INT), public, parameter :: errInvalidDirective = -101
  integer(kind=AE_INT), public, parameter :: errIllegalValue = -102
  integer(kind=AE_INT), public, parameter :: errAquiferNotSet = -103
  integer(kind=AE_INT), public, parameter :: errNoSolution = -104
  integer(kind=AE_INT), public, parameter :: errNotAllocated = -105
  integer(kind=AE_INT), public, parameter :: errPreviouslyAllocated = -106
  integer(kind=AE_INT), public, parameter :: errSpaceExhausted = -107
  integer(kind=AE_INT), public, parameter :: errUnexpectedEOF = -108
  integer(kind=AE_INT), public, parameter :: errUnexpectedData = -109
  integer(kind=AE_INT), public, parameter :: errMissingArgument = -110
  integer(kind=AE_INT), public, parameter :: errNoData = -111
  integer(kind=AE_INT), public, parameter :: errNotImplemented = -112
  integer(kind=AE_INT), public, parameter :: errNoOutputFile = -113

  ! Error messages

  type, public :: ERRORMSG
    integer(kind=AE_INT) :: Code
    character(len=80) :: Message
  end type ERRORMSG
  !
  type(ERRORMSG), public, parameter, dimension(12) :: msgErrorMessages = &
                    (/ERRORMSG(errInvalidDirective, "Invalid program directive"), &
                    ERRORMSG(errIllegalValue, "Illegal data value"), &
                    ERRORMSG(errAquiferNotSet, "Aquifer parameters are not set"), &
                    ERRORMSG(errNoSolution, "No model solution is present"), &
                    ERRORMSG(errNotAllocated, "Array has not been allocated"), &
                    ERRORMSG(errPreviouslyAllocated, "Array has already been allocated"), &
                    ERRORMSG(errSpaceExhausted, "Buffer space exhausted"), &
                    ERRORMSG(errUnexpectedEOF, "Unexpected end-of-file detected"), &
                    ERRORMSG(errUnexpectedData, "Unexpected data record found"), &
                    ERRORMSG(errMissingArgument, "Missing argument in parameter list"), &
                    ERRORMSG(errNoData, "No data available"), &
                    ERRORMSG(errAssertion, "Assertion error") &
                    /)

  ! Constants for matrix generator
  ! Equation types
  integer(kind=AE_INT), public, parameter :: EQN_HEAD = 0
  integer(kind=AE_INT), public, parameter :: EQN_FLOW = 1
  integer(kind=AE_INT), public, parameter :: EQN_INHO = 2
  integer(kind=AE_INT), public, parameter :: EQN_DISCHARGE = 3
  integer(kind=AE_INT), public, parameter :: EQN_RECHARGE = 4
  integer(kind=AE_INT), public, parameter :: EQN_CONTINUITY = 5
  integer(kind=AE_INT), public, parameter :: EQN_POTENTIALDIFF = 6
  integer(kind=AE_INT), public, parameter :: EQN_TOTALFLOW = 7
  integer(kind=AE_INT), public, parameter :: EQN_BDYGHB = 8
  ! Element types
  integer(kind=AE_INT), public, parameter :: ELEM_AQU = 1
  integer(kind=AE_INT), public, parameter :: ELEM_WL0 = 2
  integer(kind=AE_INT), public, parameter :: ELEM_PD0 = 3
  integer(kind=AE_INT), public, parameter :: ELEM_LS0 = 4
  integer(kind=AE_INT), public, parameter :: ELEM_LS1 = 5
  integer(kind=AE_INT), public, parameter :: ELEM_HB0 = 6
  integer(kind=AE_INT), public, parameter :: ELEM_IN0 = 7
  integer(kind=AE_INT), public, parameter :: ELEM_RCH = 8
  integer(kind=AE_INT), public, parameter :: ELEM_WL1 = 9
  integer(kind=AE_INT), public, parameter :: ELEM_LS2 = 10
  integer(kind=AE_INT), public, parameter :: ELEM_LS3 = 11
  integer(kind=AE_INT), public, parameter :: ELEM_AS0 = 12
  integer(kind=AE_INT), public, parameter :: ELEM_CW0 = 1001

  ! Influence Function types
  integer(kind=AE_INT), public, parameter :: INFLUENCE_P = 1
  integer(kind=AE_INT), public, parameter :: INFLUENCE_W = 2
  integer(kind=AE_INT), public, parameter :: INFLUENCE_F = 3
  integer(kind=AE_INT), public, parameter :: INFLUENCE_G = 4
  integer(kind=AE_INT), public, parameter :: INFLUENCE_Q = 5
  integer(kind=AE_INT), public, parameter :: INFLUENCE_J = 6
  integer(kind=AE_INT), public, parameter :: INFLUENCE_D = 7
  integer(kind=AE_INT), public, parameter :: INFLUENCE_Z = 8
  ! Size and other inquiry options for XXX_GetInfo()
  integer(kind=AE_INT), public, parameter :: SIZE_FWL = 1
  integer(kind=AE_INT), public, parameter :: SIZE_FDP = 2
  integer(kind=AE_INT), public, parameter :: SIZE_FPD = 3
  integer(kind=AE_INT), public, parameter :: SIZE_UNKNOWNS = 4
  integer(kind=AE_INT), public, parameter :: SIZE_EQUATIONS = 5
  integer(kind=AE_INT), public, parameter :: INFO_REGENERATE = 6
  ! Date retrieval selectors
  integer(kind=AE_INT), public, parameter :: VALUE_NEXT = -1
  integer(kind=AE_INT), public, parameter :: VALUE_HEAD = 0
  integer(kind=AE_INT), public, parameter :: VALUE_FLOW = 1
  integer(kind=AE_INT), public, parameter :: VALUE_VELOCITY = 2
  integer(kind=AE_INT), public, parameter :: VALUE_DISCHARGE = 3
  integer(kind=AE_INT), public, parameter :: VALUE_RECHARGE = 4
  integer(kind=AE_INT), public, parameter :: VALUE_POTENTIALDIFF = 5
  integer(kind=AE_INT), public, parameter :: VALUE_TOTALFLOW = 6
  integer(kind=AE_INT), public, parameter :: VALUE_POTENTIAL = 7
  integer(kind=AE_INT), public, parameter :: VALUE_SATDTHICK = 8
  integer(kind=AE_INT), public, parameter :: VALUE_TRANSMISSIVITY = 9
  integer(kind=AE_INT), public, parameter :: VALUE_EXTRACTION = 10
  ! Boundary element flags
  integer(kind=AE_INT), public, parameter :: BDY_HEAD = 0
  integer(kind=AE_INT), public, parameter :: BDY_FLUX = 1
  integer(kind=AE_INT), public, parameter :: BDY_GHB = 2
  integer(kind=AE_INT), public, parameter :: BDY_FREESURF = 3
  ! Flags for the LS2 elements
  integer(kind=AE_INT), public, parameter :: LS2_MODE_GHB = 0
  integer(kind=AE_INT), public, parameter :: LS2_MODE_RIVER = 1
  integer(kind=AE_INT), public, parameter :: LS2_MODE_DRAIN = 2
  ! Flags for the LS3 elements
  integer(kind=AE_INT), public, parameter :: LS3_MODE_GHB = 0
  integer(kind=AE_INT), public, parameter :: LS3_MODE_RIVER = 1
  integer(kind=AE_INT), public, parameter :: LS3_MODE_DRAIN = 2
  ! Flags for the top and bottom AS0 elements
  integer(kind=AE_INT), public, parameter :: AS0_TOP = 0
  integer(kind=AE_INT), public, parameter :: AS0_BOTTOM = 1

  ! Computational constants
  real(kind=AE_REAL), public, parameter :: rPI = 3.14159265359_AE_REAL
  real(kind=AE_REAL), public, parameter :: rPIOVER180 = 3.14159265359_AE_REAL / 180.0_AE_REAL
  real(kind=AE_REAL), public, parameter :: rTINY = 1.0e-20_AE_REAL
  real(kind=AE_REAL), public, parameter :: rHUGE = 1.0e+20_AE_REAL
  real(kind=AE_REAL), public, parameter :: rZERO = 0.0_AE_REAL
  real(kind=AE_REAL), public, parameter :: rHALF = 0.5_AE_REAL
  real(kind=AE_REAL), public, parameter :: rONE_FOURTH = 0.25_AE_REAL
  real(kind=AE_REAL), public, parameter :: rONE_TENTH = 0.10_AE_REAL
  real(kind=AE_REAL), public, parameter :: rTHREE_FOURTHS = 0.75_AE_REAL
  real(kind=AE_REAL), public, parameter :: rONE_EIGHTH = 0.125_AE_REAL
  real(kind=AE_REAL), public, parameter :: rONE = 1.0_AE_REAL
  real(kind=AE_REAL), public, parameter :: rTWO = 2.0_AE_REAL
  real(kind=AE_REAL), public, parameter :: rFOUR = 4.0_AE_REAL
  real(kind=AE_REAL), public, parameter :: rEIGHT = 8.0_AE_REAL
  real(kind=AE_REAL), public, parameter :: rHUNDRED = 100.0_AE_REAL
  complex(kind=AE_REAL), public, parameter :: cHUGE = (1.0e+20_AE_REAL, 1.0e+20_AE_REAL)
  complex(kind=AE_REAL), public, parameter :: cZERO = (0.0_AE_REAL, 0.0_AE_REAL)
  complex(kind=AE_REAL), public, parameter :: cONE = (1.0_AE_REAL, 0.0_AE_REAL)
  complex(kind=AE_REAL), public, parameter :: cI = (0.0_AE_REAL, 1.0_AE_REAL)
  ! Tolerance near a vertex along a linear feature for point moves(in Big-Z space)
  real(kind=AE_REAL), public, parameter :: rVERTEXTOL = 1.0e-9_AE_REAL
  ! Movement factor for vertex movements(to make sure the new position is more than rVERTEXTOL away
  real(kind=AE_REAL), public, parameter :: rMOVE_FACTOR = 1.1_AE_REAL
  ! Toerance near the intersection between a tracing segment and a line-element
  real(kind=AE_REAL), public, parameter :: rCROSSING_TOLERANCE = 0.001_AE_REAL
  ! Regeneration flags
  integer(kind=AE_INT), public, parameter :: REGENERATE_NO = 0
  integer(kind=AE_INT), public, parameter :: REGENERATE_YES = 1

  ! Types and constants for check iterators
  ! The iterator system is implemented in every element module that has unknown
  ! strengths.  The iterator is stored internally inside the module and is
  ! managed using the XXX_ResetIterator(), XXX_GetIterator(), and XXX_SetIterator()
  ! methods.  In module AEM, the ITERATOR_RESULT class provides a mechanism
  ! for inquiries in AEM_Check()

  type, public :: ITERATOR_RESULT
    integer(kind=AE_INT) :: iValueSelector
    integer(kind=AE_INT) :: iElementType
    integer(kind=AE_INT) :: iElementString
    integer(kind=AE_INT) :: iElementVertex
    integer(kind=AE_INT) :: iElementFlag
    complex(kind=AE_REAL), dimension(:), pointer :: cZ
  end type ITERATOR_RESULT


contains


  function Uppercase (s) result(upper_s)
    ! Returns s in uppercase
    character(len=*), intent(in) :: s
    character(len=256) :: upper_s
    integer(kind=AE_INT) :: ic

    upper_s = s
    do ic = 1, len(s)
      if (ichar(s(ic:ic)) > 96 .and. ichar(s(ic:ic)) < 122) then
        upper_s(ic:ic) = char(ichar(s(ic:ic))-32)
      else
        upper_s(ic:ic) = s(ic:ic)
      end if
    end do
  end function Uppercase

end module u_constants
