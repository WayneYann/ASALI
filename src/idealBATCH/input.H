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

ASALI::READinput       input(*file);
ASALI::equationSystem *eq;

{
    std::string folderPath = input.getKineticsPath() + "/";
    boost::filesystem::path path_kinetics(folderPath);
    {
        rapidxml::xml_document<> doc;
        std::vector<char> xml_string;
        OpenSMOKE::OpenInputFileXML(doc, xml_string, path_kinetics / "kinetics.xml");

        thermodynamicsMapXML = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(doc); 
        transportMapXML      = new OpenSMOKE::TransportPropertiesMap_CHEMKIN(doc); 
        kineticsMapXML       = new OpenSMOKE::KineticsMap_CHEMKIN(*thermodynamicsMapXML, doc);
    }
    {
        rapidxml::xml_document<> doc;
        std::vector<char> xml_string;
        OpenSMOKE::OpenInputFileXML(doc, xml_string, path_kinetics / "kinetics.surface.xml");

        thermodynamicsSurfaceMapXML = new OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN(doc); 
        kineticsSurfaceMapXML       = new OpenSMOKE::KineticsMap_Surface_CHEMKIN_Lumped(*thermodynamicsSurfaceMapXML, doc);
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

    {
        std::vector<double> InletCOMPOSITION(input.getCoverage().size());
        double sum = 0.;
        for(unsigned int k=0; k<input.getCoverage().size(); k++)
            sum = sum + input.getCoverage()[k];

        if ( sum != 1. )
        {
            std::cout << "\nASALI::WARNING Sum of coverage composition is not 1!\n" << std::endl;
            std::cout << "\nPress Enter to normalized them..." << std::endl;
            getchar();
            for(unsigned int k=0; k<input.getCoverage().size(); k++)
                InletCOMPOSITION[k]=input.getCoverage()[k]/sum;
        }
        else
        {
            for(unsigned int k=0; k<input.getCoverage().size(); k++)
                InletCOMPOSITION[k]=input.getCoverage()[k];
        }

        double MW;

        OpenSMOKE::OpenSMOKEVectorDouble x(SURF_NC);
        for(unsigned int k=0;k<input.getCoverageName().size();k++)
        {
            x0site[thermodynamicsSurfaceMapXML->IndexOfSpecies(input.getCoverageName()[k])-NC] = InletCOMPOSITION[k];
        }

    }

    double MW;
    {
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
    eq = new ASALI::equationSystem(*thermodynamicsMapXML,        *kineticsMapXML,
                                   *thermodynamicsSurfaceMapXML, *kineticsSurfaceMapXML,
                                   *transportMapXML);

    eq->setTemperature(input.getTemperature());
    eq->setPressure(input.getPressure());
    eq->setVolume(input.getVolume());
    eq->setArea(input.getCatalyticArea());
    eq->setCatalystLoad(input.getAlfa());
    eq->setHomogeneusReactions(input.getHomogeneousReactions());
    eq->setHeterogeneusReactions(input.getHeterogenousReactions());
    eq->setEnergy(input.getEnergy());
    eq->setResolutionType(input.getResolutionType());
}
