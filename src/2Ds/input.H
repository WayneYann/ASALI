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
                              
    ode->setSupportProperties(input.getSupportSpecificHeat(),
                              input.getSupportConductivity(),
                              input.getSupportDensity());
    
    ode->setGrid(input.axialPoints(),input.radialPoints());
    
    ode->setReactionType(input.getReactionType());
    ode->setInertLength(input.getInertLength());
    
    ode->setPressure(input.getPressure());
    ode->setCoolantTemperature(input.getCoolantTemperature());
    ode->setFeedTemperature(input.getFeedTemperature());
    ode->setFlowRate(input.getSpecificMassFlowRate());
    ode->setMassFraction(input.getMassFraction());

    ode->setThermodinamicsData(input.getMW(),input.numberOfSpecies());

    if ( input.getModels()[0] == true) // Honeycomb
    {
        ode->setHoneyComb(input.getHoneyCombTubeDiameter(),
                          input.getHoneyCombCPSI(),
                          input.getHoneyCombWall(),         // [mills = inch*1000]
                          input.getHoneyCombWashCoat(),
                          input.getHoneyCombLength(),
                          input.getHoneyCombType());
    }

    if ( input.getModels()[1] == true) // Packed bed
    {
        ode->setPackedBed(input.getPackedBedTubeDiameter(),
                          input.getPackedBedParticleDiameter(),
                          input.getPackedBedLength());
    }

    if ( input.getModels()[2] == true) // Micro packed bed
    {
        ode->setMicroBed(input.getMicroBedTubeDiameter(),
                         input.getMicroBedCPSI(),
                         input.getMicroBedParticleDiameter(),
                         input.getMicroBedWall(),         // [mills = inch*1000]
                         input.getMicroBedLength());
    }

}
ASALI::BVPSystem  *bvp;
{
    bvp = new ASALI::BVPSystem();

    bvp->setGasProperties(input.getGasSpecificHeat(),input.getGasConductivity(),
                          input.getGasViscosity(),   input.getGasDiffusivity());
                         
    bvp->setCatalystProperties(input.getCatalystSpecificHeat(),
                               input.getCatalystConductivity(),
                               input.getCatalystDensity(),
                               input.getCatalystVoidFraction(),
                               input.getCatalystTortuosity());
                              
    bvp->setSupportProperties(input.getSupportSpecificHeat(),
                              input.getSupportConductivity(),
                              input.getSupportDensity());
    
    bvp->setGrid(input.axialPoints(),input.radialPoints());
    
    bvp->setReactionType(input.getReactionType());
    bvp->setInertLength(input.getInertLength());

    bvp->setPressure(input.getPressure());
    bvp->setCoolantTemperature(input.getCoolantTemperature());
    bvp->setFeedTemperature(input.getFeedTemperature());
    bvp->setFlowRate(input.getSpecificMassFlowRate());
    bvp->setMassFraction(input.getMassFraction());

    bvp->setThermodinamicsData(input.getMW(),input.numberOfSpecies());

    if ( input.getModels()[0] == true) // Honeycomb
    {
        bvp->setHoneyComb(input.getHoneyCombTubeDiameter(),
                          input.getHoneyCombCPSI(),
                          input.getHoneyCombWall(),         // [mills = inch*1000]
                          input.getHoneyCombWashCoat(),
                          input.getHoneyCombLength(),
                          input.getHoneyCombType());
    }

    if ( input.getModels()[1] == true) // Packed bed
    {
        bvp->setPackedBed(input.getPackedBedTubeDiameter(),
                          input.getPackedBedParticleDiameter(),
                          input.getPackedBedLength());
    }

    if ( input.getModels()[2] == true) // Micro packed bed
    {
        bvp->setMicroBed(input.getMicroBedTubeDiameter(),
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
