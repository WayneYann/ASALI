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
                              
    eq->setSupportProperties(input.getSupportSpecificHeat(),
                              input.getSupportConductivity(),
                              input.getSupportDensity());
    
    eq->setGrid(input.axialPoints(),input.radialPoints());
    
    eq->setReactionType(input.getReactionType());
    eq->setInertLength(input.getInertLength());

    eq->setPressure(input.getPressure());
    eq->setCoolantTemperature(input.getCoolantTemperature());
    eq->setFeedTemperature(input.getFeedTemperature());
    eq->setFlowRate(input.getSpecificMassFlowRate());
    eq->setMassFraction(input.getMassFraction());

    eq->setThermodinamicsData(input.getMW(),input.numberOfSpecies());

    if ( input.getModels()[0] == true) // Honeycomb
    {
        eq->setHoneyComb(input.getHoneyCombTubeDiameter(),
                         input.getHoneyCombCPSI(),
                         input.getHoneyCombWall(),         // [mills = inch*1000]
                         input.getHoneyCombWashCoat(),
                         input.getHoneyCombLength(),
                         input.getHoneyCombType());
    }

    if ( input.getModels()[1] == true) // Packed bed
    {
        eq->setPackedBed(input.getPackedBedTubeDiameter(),
                         input.getPackedBedParticleDiameter(),
                         input.getPackedBedLength());
    }

    if ( input.getModels()[2] == true) // Micro packed bed
    {
        eq->setMicroBed(input.getMicroBedTubeDiameter(),
                        input.getMicroBedCPSI(),
                        input.getMicroBedParticleDiameter(),
                        input.getMicroBedWall(),         // [mills = inch*1000]
                        input.getMicroBedLength());
    }

}

{
    ASALI::logo();
    input.recapOnScreen();
    unsigned int WP = 6;
    std::cout.setf(std::ios::scientific);
    std::cout.precision(WP);
}
