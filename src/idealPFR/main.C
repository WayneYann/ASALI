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
#define BZZ_COMPILER 101
#define OPENSMOKE_USE_BZZMATH 1
#include "BzzMath.hpp"
#endif

#if ASALI_USE_BZZ == 0
#define OPENSMOKE_USE_BZZMATH 0
#endif

// OpenSMOKE++
#include "OpenSMOKEpp"
#include "maps/Maps_CHEMKIN"
#include "reactors/utilities/Utilities"
#include "math/OpenSMOKEVector.h"
#include "math/multivalue-ode-solvers/MultiValueSolver"
#include "math/multivalue-dae-solvers/MultiValueSolver"

// Eigen
#include <Eigen/Dense>

// Boost
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

// RapidXML
#include "rapidxml.hpp"

// Utilities
#include "readInput.h"
#include "functions.h"
#include "vector.h"

// Lumped reactions
#include "lumpedReactions/reactionRates.H"

// Equations
#include "memoryAllocation.H"
#include "equations.h"
#include "odeInterfaces.h"
#include "daeInterfaces.h"

int main( int argc, char** argv )
{
    double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

    #include "input.H"
    #include "resolution.H"
    #include "write.H"

    remove("BzzFile.txt");

    double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
    ASALI::CPUtime(tStart,tEnd);

    return 0;
}
