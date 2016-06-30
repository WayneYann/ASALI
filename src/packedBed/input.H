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

ASALI::check(argc);
std::string*      file = new std::string(argv[1]);
ASALI::readInput  input(*file);

ASALI::equationSystem  *eq;
{
    eq = new ASALI::equationSystem();

    eq->setGasProperties(input.getGasSpecificHeat(),input.getGasConductivity(),
                         input.getGasViscosity(),   input.getGasDiffusivity());
                         
    eq->setCatalystProperties(input.getCatalystSpecificHeat(),
                              input.getCatalystConductivity(),
                              input.getCatalystDensity(),
                              input.getCatalystVoidFraction(),
                              input.getCatalystTortuosity());

    eq->setReactionType(input.getReactionType());
    
    eq->setPressure(input.getPressure());
    eq->setCoolantTemperature(input.getCoolantTemperature());
    eq->setFlowRate(input.getSpecificMassFlowRate());

    eq->setThermodinamicsData(input.getMW(),input.numberOfSpecies());

    eq->setGeometry(input.getTubeDiameter(),
                    input.getParticleDiameter());

    eq->setHeatExchange(input.externalExchange());
    eq->setEnergy(input.energy());
}

{
    ASALI::logo();
    input.recapOnScreen();
    unsigned int WP = 6;
    std::cout.setf(std::ios::scientific);
    std::cout.precision(WP);
}
