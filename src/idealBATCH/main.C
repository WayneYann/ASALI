/*##############################################################################################
#                                                                                              #
#     #############       #############       #############       ####                ####     #
#    #             #     #             #     #             #     #    #              #    #    #
#    #    #####    #     #    #########      #    #####    #     #    #              #    #    #
#    #    #   #    #     #    #              #    #   #    #     #    #              #    #    #
#    #    #####    #     #    #              #    #####    #     #    #              #    #    #
#    #             #     #    #########      #             #     #    #              #    #    #
#    #             #     #             #     #             #     #    #              #    #    #
#    #    #####    #      #########    #     #    #####    #     #    #              #    #    #
#    #    #   #    #              #    #     #    #   #    #     #    #              #    #    #
#    #    #   #    #      #########    #     #    #   #    #     #    #########      #    #    #
#    #    #   #    #     #             #     #    #   #    #     #             #     #    #    #
#     ####     ####       #############       ####     ####       #############       ####     #
#                                                                                              #
#   Department of Energy                                                                       #
#   Politecnico di Milano                                                                      #
#   Author: Stefano Rebughini <stefano.rebughini@polimi.it>                                    #
#                                                                                              #
################################################################################################
#                                                                                              #
#   License                                                                                    #
#                                                                                              #
#   This file is part of ASALI.                                                                #
#                                                                                              #
#   ASALI is free software: you can redistribute it and/or modify                              #
#   it under the terms of the GNU General Public License as published by                       #
#   the Free Software Foundation, either version 3 of the License, or                          #
#   (at your option) any later version.                                                        #
#                                                                                              #
#   ASALI is distributed in the hope that it will be useful,                                   #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of                             #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                              #
#   GNU General Public License for more details.                                               #
#                                                                                              #
#   You should have received a copy of the GNU General Public License                          #
#   along with ASALI. If not, see <http://www.gnu.org/licenses/>.                              #
#                                                                                              #
##############################################################################################*/

// C++
#include <string>
#include <iostream>
#include <math.h>
#include <ctime>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <algorithm>

#if ASALI_USE_BZZ == 1
// BzzMath
#define BZZ_COMPILER 101
#include "BzzMath.hpp"
#endif

#if ASALI_USE_SUNDIALS == 1
#define OPENSMOKE_USE_SUNDIALS 1
#endif

#if ASALI_USE_SUNDIALS == 0
#define OPENSMOKE_USE_SUNDIALS 0
#endif

// Eigen
#include <Eigen/Dense>

// Boost
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

// RapidXML
#include "rapidxml.hpp"

// OpenSMOKE++ Definitions
#include "OpenSMOKEpp"

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"

// Reactor utilities
#include "reactors/utilities/Utilities"

// Ode
#include "math/multivalue-ode-solvers/MultiValueSolver"
#include "ode/ODE_Parameters.h"

// Utilities
#include "readInput.h"
#include "functions.h"
#include "vector.h"

// Lumped reactions
#include "lumpedReactions/reactionRates.H"

// Equations
#include "memoryAllocation.H"

#include "ODEsystem.h"
#include "odeInterfaces.h"

int main( int argc, char** argv )
{
    double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

    #include "input.H"
    #include "ODEresolution.H"
    #include "write.H"

    remove("BzzFile.txt");

    double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
    ASALI::CPUtime(tStart,tEnd);

    return 0;
}
