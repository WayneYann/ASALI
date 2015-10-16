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
ASALI::BVPSystem *bvp;
dustyGas*         dgm;

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
        dgm = new dustyGas(*thermodynamicsMapXML, *transportMapXML);
    }

    ASALI::logo();
    input.recapOnScreen();

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

    N       = input.getNumberOfPoints();
    NC      = thermodynamicsMapXML->NumberOfSpecies();
    SURF_NC = thermodynamicsSurfaceMapXML->number_of_site_species();
    NE      = (NC + SURF_NC)*N;

    ChangeDimensions(N,       &z,      true);
    ChangeDimensions(NC,      &x0bulk, true);
    ChangeDimensions(SURF_NC, &x0site, true);

    thermodynamicsMapXML->SetTemperature(input.getTemperature());
    thermodynamicsMapXML->SetPressure(input.getPressure());

    x0site[1] = 1.;
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

    double step = input.getReactorLength()/(N-1);
    z[1] = 0.;
    for (unsigned int k=2;k<=z.Size();k++)
        z[k] = z[k-1] + step;

    ChangeDimensions(NE, &x0, true);
    unsigned int counter = 1;
    for (unsigned int k=1;k<=N;k++)
    {
        for (unsigned int j=1;j<=NC;j++)
        {
            if (k==1)
                x0[counter++] = x0bulk[j];
            else
                x0[counter++] = x0bulk[j];
        }
        for (unsigned int j=1;j<=SURF_NC;j++)
        {
            x0[counter++] = x0site[j];
        }
    }
}

{
    ode = new ASALI::ODESystem(*thermodynamicsMapXML,        *kineticsMapXML,
                               *thermodynamicsSurfaceMapXML, *kineticsSurfaceMapXML,
                               *transportMapXML);

    ode->setReactorGeometry(input.getReactorLength(),input.getReactorDiameter());
    ode->setFeedValue(input.getPressure(),input.getTemperature(),x0bulk);
    ode->setInert(input.getInert());

    bvp = new ASALI::BVPSystem(*thermodynamicsMapXML,        *kineticsMapXML,
                               *thermodynamicsSurfaceMapXML, *kineticsSurfaceMapXML,
                               *transportMapXML,             *dgm);

    bvp->setReactions(input.getReactions());
    bvp->setCatalystProperties(input.getVoidFraction(),input.getTourtuosity(),input.getPoreRadius());
    bvp->setReactorGeometry(input.getAlfa(),input.getReactorLength(),input.getReactorDiameter());
    bvp->setFeedValue(input.getPressure(),input.getTemperature(),x0bulk);
    bvp->setInert(input.getInert());
    bvp->setDiffusionModel(input.getDiffusionType());
}

