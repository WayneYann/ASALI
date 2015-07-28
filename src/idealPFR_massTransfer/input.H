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
ASALI::DAESystem *dae;
ASALI::ODESystem *ode;

{
    std::string folderPath = input.getKineticsPath() + "/";
    std::cout << input.getKineticsPath() << std::endl;
    boost::filesystem::path path_kinetics(folderPath);
    {
        rapidxml::xml_document<> doc;
        std::vector<char> xml_string;
        OpenSMOKE::OpenInputFileXML(doc, xml_string, path_kinetics / "kinetics.xml");

        thermodynamicsMapXML = new OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>(doc); 
        transportMapXML = new OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>(doc); 
        kineticsMapXML = new OpenSMOKE::KineticsMap_CHEMKIN<double>(*thermodynamicsMapXML, doc);
    }
    {
        rapidxml::xml_document<> doc;
        std::vector<char> xml_string;
        OpenSMOKE::OpenInputFileXML(doc, xml_string, path_kinetics / "kinetics.surface.xml");

        thermodynamicsSurfaceMapXML = new OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN<double>(doc); 
        kineticsSurfaceMapXML = new OpenSMOKE::KineticsMap_Surface_CHEMKIN_Lumped<double>(*thermodynamicsSurfaceMapXML, doc);
    }
    std::cout << "Sono qui" << std::endl;
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
        NC = thermodynamicsMapXML->NumberOfSpecies();
        ChangeDimensions(NC,      &x0, true);

        thermodynamicsMapXML->SetTemperature(input.getFeedTemperature());
        thermodynamicsMapXML->SetPressure(input.getFeedPressure());

        if ( input.getInletType() == "mole" )
        {
            OpenSMOKE::OpenSMOKEVectorDouble x(NC);
            for(unsigned int k=0;k<input.getFeedName().size();k++)
                x[thermodynamicsMapXML->IndexOfSpecies(input.getFeedName()[k])] = input.getFeedComposition()[k];

            thermodynamicsMapXML->MassFractions_From_MoleFractions(x0,MW,x);
        }
        else if ( input.getInletType() == "mass" )
        {
            for(unsigned int k=0;k<input.getFeedName().size();k++)
                x0[thermodynamicsMapXML->IndexOfSpecies(input.getFeedName()[k])] = input.getFeedComposition()[k];

            thermodynamicsMapXML->MolecularWeight_From_MassFractions(MW,x0);
        }

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
        dae = new ASALI::DAESystem(*thermodynamicsMapXML,        *kineticsMapXML,
                                   *thermodynamicsSurfaceMapXML, *kineticsSurfaceMapXML,
                                   *transportMapXML);


        dae->SetReactorType(input.getReactorType());
        dae->SetInertSpecie(input.getInertSpecie());
        dae->SetReactor(G,input.getParticleDiameter(),input.getShellDiameter(),epsi);
        dae->SetTemperature(input.getFeedTemperature());
        dae->SetCorrelation(input.getPressureCorrelation(),input.getMassTransferCorrelation());

        ode = new ASALI::ODESystem(*thermodynamicsMapXML,        *kineticsMapXML,
                                   *thermodynamicsSurfaceMapXML, *kineticsSurfaceMapXML,
                                   *transportMapXML);

        ode->SetReactorType(input.getReactorType());
        ode->SetInertSpecie(input.getInertSpecie());
        ode->SetReactor(G,input.getParticleDiameter(),input.getShellDiameter(),epsi);
        ode->SetTemperature(input.getFeedTemperature());
        ode->SetCorrelation(input.getPressureCorrelation(),input.getMassTransferCorrelation());
    }

}
