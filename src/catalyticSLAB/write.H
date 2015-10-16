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

{
    std::string folderName;
    std::string results;

    folderName = "results/" + ASALI::to_string(N);
    boost::filesystem::create_directory(folderName);

    if ( input.getResults() == "mass")
    {
        results = "mass";
    }
    else if ( input.getResults() == "mole")
    {
        results = "mole";
    }
 
    std::ofstream bulk;
    std::string fileNameBulk = folderName + "/bulkSpecies.txt";
    const char *pathBulk = fileNameBulk.c_str();
    remove(pathBulk);
    bulk.open(pathBulk,std::ios::out);

    std::ofstream site;
    std::string fileNameSite = folderName + "/adsorbedSpecies.txt";
    const char *pathSite = fileNameSite.c_str();
    remove(pathSite);
    site.open(pathSite,std::ios::out);

    std::ofstream grid;
    std::string fileNameGrid = folderName + "/grid.txt";
    const char *pathGrid = fileNameGrid.c_str();
    remove(pathGrid);
    grid.open(pathGrid,std::ios::out);

    std::ofstream extra;
    std::string fileNameExtra = folderName + "/information.txt";
    const char *pathExtra = fileNameExtra.c_str();
    remove(pathExtra);
    extra.open(pathExtra,std::ios::out);

    {
        unsigned int WP = 6;
        bulk.setf(std::ios::scientific);
        grid.setf(std::ios::scientific);
        site.setf(std::ios::scientific);
        extra.setf(std::ios::scientific);
        bulk.precision(WP);
        grid.precision(WP);
        site.precision(WP);
        extra.precision(WP);
    }

    OpenSMOKE::OpenSMOKEVectorDouble *omegaBulk;
    OpenSMOKE::OpenSMOKEVectorDouble *teta;

    omegaBulk = new OpenSMOKE::OpenSMOKEVectorDouble[N];
    teta      = new OpenSMOKE::OpenSMOKEVectorDouble[N];

    for (unsigned int i=0;i<N;i++)
    {
        ChangeDimensions(NC,      &omegaBulk[i], true);
        ChangeDimensions(SURF_NC, &teta[i],      true);
    }

    unsigned int counter = 1;
    for (unsigned int i=0;i<N;i++)
    {
        for (unsigned int j=1;j<=NC;j++)
            omegaBulk[i][j] = std::max(0.,yF[counter++]);
        for (unsigned int j=1;j<=SURF_NC;j++)
            teta[i][j] = std::max(0.,yF[counter++]);
    }
    
    for (unsigned int j=1;j<=NC;j++)
        bulk << thermodynamicsMapXML->NamesOfSpecies()[j-1] << "\t";
    bulk << std::endl;

    for (unsigned int j=1;j<=SURF_NC;j++)
        site << thermodynamicsSurfaceMapXML->NamesOfSpecies()[j-1+NC] << "\t";
    site << std::endl;

    if ( results == "mass")
    {
        for (unsigned int k=0;k<N;k++)
        {
            for (unsigned int j=1;j<=NC;j++)
                bulk << omegaBulk[k][j] << "\t";
            bulk << std::endl;

            for (unsigned int j=1;j<=SURF_NC;j++)
                site << teta[k][j] << "\t";
            site << std::endl;
        }
    }
    else if ( results == "mole" )
    {
        for (unsigned int k=0;k<N;k++)
        {
            double dummy;
            OpenSMOKE::OpenSMOKEVectorDouble x(NC);
            thermodynamicsMapXML->MoleFractions_From_MassFractions(x,dummy,omegaBulk[k]);
            for (unsigned int j=1;j<=NC;j++)
                bulk << x[j] << "\t" ;
            bulk << std::endl;

            for (unsigned int j=1;j<=SURF_NC;j++)
                site << teta[k][j] << "\t";
            site << std::endl;
        }
    }

    bulk.close();
    site.close();

    for (unsigned int k=1;k<=N;k++)
    {
        grid << z[k] << std::endl;
    }

    grid.close();


    std::vector<double> IN;
    std::vector<double> OUT;
    {
        thermodynamicsMapXML->SetPressure(input.getPressure());
        thermodynamicsMapXML->SetTemperature(input.getTemperature());
        Eigen::MatrixXd SpecieAtomicComposition;
        Eigen::VectorXd Element;
        std::vector<std::string> ElementsNames;
        ElementsNames = thermodynamicsMapXML->elements();
        SpecieAtomicComposition = thermodynamicsMapXML->atomic_composition();
        Element.resize(NC);
        IN.resize(ElementsNames.size());
        for (unsigned int k=0;k<ElementsNames.size();k++)
        {
            IN[k] = 0.;
            Element = SpecieAtomicComposition.col(k);
            for (unsigned int j=0;j<NC;j++)
            {
                IN[k] += x0bulk[j+1]*Element[j]/thermodynamicsMapXML->MW()[j+1];
            }
        }
    }

    {
        thermodynamicsMapXML->SetPressure(input.getPressure());
        thermodynamicsMapXML->SetTemperature(input.getTemperature());

        Eigen::MatrixXd SpecieAtomicComposition;
        Eigen::VectorXd Element;
        std::vector<std::string> ElementsNames;
        ElementsNames = thermodynamicsMapXML->elements();
        SpecieAtomicComposition = thermodynamicsMapXML->atomic_composition();
        Element.resize(NC);
        OUT.resize(ElementsNames.size());
        for (unsigned int k=0;k<ElementsNames.size();k++)
        {
            OUT[k] = 0.;
            Element = SpecieAtomicComposition.col(k);
            for (unsigned int j=0;j<NC;j++)
            {
                OUT[k] += omegaBulk[N-1][j+1]*Element[j]/thermodynamicsMapXML->MW()[j+1];
            }
        }
    }

    std::cout << "\n############################" << std::endl;
    std::cout << "ATOMIC BALANCE               " << std::endl;
    for (unsigned int k=0;k<IN.size();k++)
    {
        if ( IN[k] != 0. )
            std::cout << thermodynamicsMapXML->elements()[k] << " :\t" << (IN[k] - OUT[k])*100./IN[k] << "\t[%]" << std::endl;
    }
    std::cout << "############################\n" << std::endl;

    extra << "ATOMIC BALANCE               " << std::endl;
    for (unsigned int k=0;k<IN.size();k++)
    {
        if ( IN[k] != 0. )
            extra << thermodynamicsMapXML->elements()[k] << " :\t" << (IN[k] - OUT[k])*100./IN[k] << "\t[%]" << std::endl;
    }

    extra.close();






    delete [] omegaBulk;
    delete [] teta;
}

