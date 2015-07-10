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
ASALI::ICSystem  *ic;

{
    std::string folderPath = input.getKineticsPath() + "/";
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

    {
        ASALI::logo();
        input.recapOnScreen();
    }

    {
        unsigned int WP = 6;
        std::cout.setf(std::ios::scientific);
        std::cout.precision(WP);
    }
}

{
    std::string folderName = "results/";
    boost::filesystem::remove_all(folderName);
    boost::filesystem::create_directory(folderName);

    NC      = thermodynamicsMapXML->NumberOfSpecies();
    SURF_NC = thermodynamicsSurfaceMapXML->number_of_site_species();
    SURF_NP = thermodynamicsSurfaceMapXML->number_of_site_phases(0);
    NE      = NC + SURF_NC + SURF_NP + 1 + 1;

    ChangeDimensions(NC,      &x0bulk, true);
    ChangeDimensions(SURF_NC, &x0site, true);

    thermodynamicsMapXML->SetTemperature(input.getTemperature());
    thermodynamicsMapXML->SetPressure(input.getPressure());

    x0site[1] = 1.0;
    for (unsigned int k=1;k<=SURF_NC;k++)
        x0site[k]=std::max(x0site[k],0.);

    std::vector<double> InletCOMPOSITION(input.getFraction().size());
    double sum = 0.;
    for(unsigned int k=0; k<input.getFraction().size(); k++)
        sum = sum + input.getFraction()[k];

    if ( sum != 1. )
    {
        std::cout << "\nASALI::WARNING Sum of inlet composition is not 1!\n" << std::endl;
        std::cout << "\nPress Enter to normalized them..." << std::endl;
        getchar();
        for(unsigned int k=0; k<input.getFraction().size(); k++)
            InletCOMPOSITION[k]=input.getFraction()[k]/sum;
    }
    else
    {
        for(unsigned int k=0; k<input.getFraction().size(); k++)
            InletCOMPOSITION[k]=input.getFraction()[k];
    }

    double MW;
    if ( input.getFeed() == "Mole" )
    {
        OpenSMOKE::OpenSMOKEVectorDouble x(NC);
        for(unsigned int k=0;k<input.getFractionName().size();k++)
            x[thermodynamicsMapXML->IndexOfSpecies(input.getFractionName()[k])] = InletCOMPOSITION[k];

        thermodynamicsMapXML->MassFractions_From_MoleFractions(x0bulk,MW,x);
    }
    else if ( input.getFeed() == "Mass" )
    {
        for(unsigned int k=0;k<input.getFractionName().size();k++)
            x0bulk[thermodynamicsMapXML->IndexOfSpecies(input.getFractionName()[k])] = InletCOMPOSITION[k];

        thermodynamicsMapXML->MolecularWeight_From_MassFractions(MW,x0bulk);
    }


    ChangeDimensions(NE, &x0, true);
    unsigned int counter = 1;
    {
        for (unsigned int j=1;j<=NC;j++)
            x0[counter++] = x0bulk[j];
        x0[counter++] = input.getPressure()*MW*input.getVolume()/(8314.*input.getTemperature());
        for (unsigned int j=1;j<=SURF_NP;j++)
            x0[counter++] = thermodynamicsSurfaceMapXML->matrix_densities_site_phases()[0][j-1];
        for (unsigned int j=1;j<=SURF_NC;j++)
            x0[counter++] = x0site[j];
        x0[counter++] = input.getTemperature();
    }
}

{
    ic = new ASALI::ICSystem(*thermodynamicsMapXML,        *kineticsMapXML,
                             *thermodynamicsSurfaceMapXML, *kineticsSurfaceMapXML,
                             *transportMapXML);

    ic->setTemperature(input.getTemperature());
    ic->setPressure(input.getPressure());
    ic->setVolume(input.getVolume());
    ic->setCatalystLoad(input.getAlfa());
    ic->setHeterogeneusReactions(input.getHeterogenousReactions());
    ic->setResolutionType(input.getResolutionType());

    ode = new ASALI::ODESystem(*thermodynamicsMapXML,        *kineticsMapXML,
                               *thermodynamicsSurfaceMapXML, *kineticsSurfaceMapXML,
                               *transportMapXML);

    ode->setTemperature(input.getTemperature());
    ode->setPressure(input.getPressure());
    ode->setVolume(input.getVolume());
    ode->setCatalystLoad(input.getAlfa());
    ode->setHomogeneusReactions(input.getHomogeneousReactions());
    ode->setHeterogeneusReactions(input.getHeterogenousReactions());
    ode->setEnergy(input.getEnergy());
    ode->setResolutionType(input.getResolutionType());
}
