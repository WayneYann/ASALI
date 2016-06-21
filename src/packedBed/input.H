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

ASALI::ODESystem  *ode;
{
    ode = new ASALI::ODESystem();

    ode->setGasProperties(input.getGasSpecificHeat(),input.getGasConductivity(),
                          input.getGasViscosity(),   input.getGasDiffusivity());
                         
    ode->setCatalystProperties(input.getCatalystSpecificHeat(),
                               input.getCatalystConductivity(),
                               input.getCatalystDensity(),
                               input.getCatalystVoidFraction(),
                               input.getCatalystTortuosity());

    ode->setReactionType(input.getReactionType());
    
    ode->setPressure(input.getPressure());
    ode->setCoolantTemperature(input.getCoolantTemperature());
    ode->setFlowRate(input.getSpecificMassFlowRate());

    ode->setThermodinamicsData(input.getMW(),input.numberOfSpecies());

    ode->setGeometry(input.getTubeDiameter(),
                     input.getParticleDiameter());

    ode->setHeatExchange(input.externalExchange());
}

ASALI::DAESystem  *dae;
{
    dae = new ASALI::DAESystem();

    dae->setGasProperties(input.getGasSpecificHeat(),input.getGasConductivity(),
                          input.getGasViscosity(),   input.getGasDiffusivity());
                         
    dae->setCatalystProperties(input.getCatalystSpecificHeat(),
                               input.getCatalystConductivity(),
                               input.getCatalystDensity(),
                               input.getCatalystVoidFraction(),
                               input.getCatalystTortuosity());

    dae->setReactionType(input.getReactionType());
    
    dae->setPressure(input.getPressure());
    dae->setCoolantTemperature(input.getCoolantTemperature());
    dae->setFlowRate(input.getSpecificMassFlowRate());

    dae->setThermodinamicsData(input.getMW(),input.numberOfSpecies());

    dae->setGeometry(input.getTubeDiameter(),
                     input.getParticleDiameter());

    dae->setHeatExchange(input.externalExchange());
}

{
    ASALI::logo();
    input.recapOnScreen();
    unsigned int WP = 6;
    std::cout.setf(std::ios::scientific);
    std::cout.precision(WP);
}
