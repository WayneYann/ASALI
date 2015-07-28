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
std::string* file;

file    = new std::string(argv[1]);

ASALI::READinput  input(*file);
ASALI::ODESystem *ode;
{
    std::string folderPath = input.getKineticsPath() + "/";
    boost::filesystem::path path_kinetics(folderPath);
    {
        rapidxml::xml_document<> doc;
        std::vector<char> xml_string;
        OpenSMOKE::OpenInputFileXML(doc, xml_string, path_kinetics / "kinetics.xml");

        thermodynamicsMapXML = new OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>(doc); 
        transportMapXML = new OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>(doc); 
    }

    {
        ASALI::logo();
        input.recapOnScreen();
    }

    {
        unsigned int WP = 6;
        std::cout.setf(std::ios::scientific);
        std::cout.precision(WP);
    }
    
    
    {
        double MW;
        thermodynamicsMapXML->SetTemperature(input.getFeedTemperature());
        thermodynamicsMapXML->SetPressure(input.getFeedPressure());
        MW = thermodynamicsMapXML->MW()[thermodynamicsMapXML->IndexOfSpecies(input.getSpecie())];

        if ( input.getReactorType() == "Monolith" )
        {
            epsi = input.getVoidFraction();
            G = input.getFeedVelocity()*MW*input.getFeedPressure()/(PhysicalConstants::R_J_kmol*input.getFeedTemperature());
        }
        else if ( input.getReactorType() == "PackedBed" )
        {
            epsi = ASALI::voidFraction(input.getParticleDiameter(),input.getShellDiameter());
            G = (input.getFeedVelocity()*MW*input.getFeedPressure()/(PhysicalConstants::R_J_kmol*input.getFeedTemperature()));
        }
    }
    
    
    {
        ode = new ASALI::ODESystem(*thermodynamicsMapXML, *transportMapXML);

        ode->SetReactorType(input.getReactorType());
        ode->SetFluid(input.getSpecie());
        ode->SetTransportLaw(input.getThermophysicalProperties());
        ode->SetReactor(G,input.getSolidTemperature(),input.getParticleDiameter(),input.getShellDiameter(),epsi);
        ode->SetCorrelation(input.getPressureCorrelation(),input.getHeatTransferCorrelation());
    }
}
