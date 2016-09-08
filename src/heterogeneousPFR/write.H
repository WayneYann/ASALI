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

    if ( input.getStart() == "new" || input.getStart() == "latest" )
    {
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
    }
    else if ( input.getStart() == "converter")
    {
        folderName = "results/" + ASALI::to_string(N) + ".converted";
        boost::filesystem::create_directory(folderName);

        if ( input.getRestartResults() == "mass")
        {
            results = "mole";
        }
        else if ( input.getRestartResults() == "mole")
        {
            results = "mass";
        }
    }
    else if ( input.getStart() == "sampling")
    {
        folderName = "results/" + ASALI::to_string(N) + ".sampling";
        boost::filesystem::create_directory(folderName);
        
        if ( input.getRestartResults() == "mass")
        {
            results = "mass";
        }
        else if ( input.getRestartResults() == "mole")
        {
            results = "mole";
        }
    }
    else if ( input.getStart() == "kinetic")
    {
        folderName = "results/" + ASALI::to_string(N) + ".kinetic";
        boost::filesystem::remove_all(folderName);
        boost::filesystem::create_directory(folderName);
        
        if ( input.getRestartResults() == "mass")
        {
            results = "mass";
        }
        else if ( input.getRestartResults() == "mole")
        {
            results = "mole";
        }
    }
    else if ( input.getStart() == "RPA")
    {
        folderName = "results/" + ASALI::to_string(N) + ".RPA";
        boost::filesystem::create_directory(folderName);

        if ( input.getRestartResults() == "mass")
        {
            results = "mass";
        }
        else if ( input.getRestartResults() == "mole")
        {
            results = "mole";
        }
    }
    else if ( input.getStart() == "conversion")
    {
        folderName = "results/" + ASALI::to_string(N) + ".conversion";
        boost::filesystem::create_directory(folderName);

        if ( input.getRestartResults() == "mass")
        {
            results = "mass";
        }
        else if ( input.getRestartResults() == "mole")
        {
            results = "mole";
        }
    }
 
    std::ofstream bulk;
    std::string fileNameBulk = folderName + "/bulkSpecies.txt";
    const char *pathBulk = fileNameBulk.c_str();
    remove(pathBulk);
    bulk.open(pathBulk,std::ios::out);

    std::ofstream wall;
    std::string fileNameWall = folderName + "/wallSpecies.txt";
    const char *pathWall = fileNameWall.c_str();
    remove(pathWall);
    wall.open(pathWall,std::ios::out);

    std::ofstream site;
    std::string fileNameSite = folderName + "/adsorbedSpecies.txt";
    const char *pathSite = fileNameSite.c_str();
    remove(pathSite);
    site.open(pathSite,std::ios::out);

    std::ofstream temp;
    std::string fileNameTemp = folderName + "/temperature.txt";
    const char *pathTemp = fileNameTemp.c_str();
    remove(pathTemp);
    temp.open(pathTemp,std::ios::out);

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

    std::ofstream vel;
    std::string fileNameVel = folderName + "/velocity.txt";
    const char *pathVel = fileNameVel.c_str();
    remove(pathVel);
    vel.open(pathVel,std::ios::out);

    {
        unsigned int WP = 6;
        bulk.setf(std::ios::scientific);
        wall.setf(std::ios::scientific);
        temp.setf(std::ios::scientific);
        grid.setf(std::ios::scientific);
        extra.setf(std::ios::scientific);
        vel.setf(std::ios::scientific);
        site.setf(std::ios::scientific);
        bulk.precision(WP);
        wall.precision(WP);
        temp.precision(WP);
        grid.precision(WP);
        extra.precision(WP);
        vel.precision(WP);
        site.precision(WP);
    }
    
    
    OpenSMOKE::OpenSMOKEVectorDouble *omegaBulk;
    OpenSMOKE::OpenSMOKEVectorDouble *omegaWall;
    OpenSMOKE::OpenSMOKEVectorDouble *teta;

    omegaBulk = new OpenSMOKE::OpenSMOKEVectorDouble[N];
    omegaWall = new OpenSMOKE::OpenSMOKEVectorDouble[N];
    teta      = new OpenSMOKE::OpenSMOKEVectorDouble[N];

    for (unsigned int i=0;i<N;i++)
    {
        ChangeDimensions(NC,      &omegaBulk[i], true);
        ChangeDimensions(NC,      &omegaWall[i], true);
        ChangeDimensions(SURF_NC, &teta[i],      true);
    }

    OpenSMOKE::OpenSMOKEVectorDouble  Tbulk(N);
    OpenSMOKE::OpenSMOKEVectorDouble  Twall(N);

    unsigned int counter = 1;
    for (unsigned int i=0;i<N;i++)
    {
        for (unsigned int j=1;j<=NC;j++)
            omegaBulk[i][j] = std::max(0.,yF[counter++]);
        for (unsigned int j=1;j<=NC;j++)
            omegaWall[i][j] = std::max(0.,yF[counter++]);
        for (unsigned int j=1;j<=SURF_NC;j++)
            teta[i][j] = std::max(0.,yF[counter++]);
        Tbulk[i+1] = std::max(0.,yF[counter++]);
        Twall[i+1] = std::max(0.,yF[counter++]);
    }
    
    for (unsigned int j=1;j<=NC;j++)
        bulk << thermodynamicsMapXML->NamesOfSpecies()[j-1] << "\t";
    bulk << std::endl;

    for (unsigned int j=1;j<=NC;j++)
        wall << thermodynamicsMapXML->NamesOfSpecies()[j-1] << "\t";
    wall << std::endl;

    for (unsigned int j=1;j<=SURF_NC;j++)
        site << thermodynamicsSurfaceMapXML->NamesOfSpecies()[j-1+NC] << "\t";
    site << std::endl;

    temp << "Bulk\tWall" << std::endl;

    if ( results == "mass")
    {
        for (unsigned int k=0;k<N;k++)
        {
            for (unsigned int j=1;j<=NC;j++)
                bulk << omegaBulk[k][j] << "\t";
            bulk << std::endl;

            for (unsigned int j=1;j<=NC;j++)
                wall << omegaWall[k][j] << "\t";
            wall << std::endl;

            for (unsigned int j=1;j<=SURF_NC;j++)
                site << teta[k][j] << "\t";
            site << std::endl;

            temp << Tbulk[k+1] << "\t" << Twall[k+1] << std::endl;
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

            thermodynamicsMapXML->MoleFractions_From_MassFractions(x, dummy,omegaWall[k]);
            for (unsigned int j=1;j<=NC;j++)
                wall << x[j] << "\t" ;
            wall << std::endl;

            for (unsigned int j=1;j<=SURF_NC;j++)
                site << teta[k][j] << "\t";
            site << std::endl;

            temp << Tbulk[k+1] << "\t" << Twall[k+1] << std::endl;
        }
    }

    bulk.close();
    wall.close();
    site.close();
    temp.close();

    for (unsigned int k=1;k<=N;k++)
    {
        grid << z[k] << std::endl;
    }

    grid.close();

    //- Change discretization schemes to CDS
    if ( grow == true) 
    {
        eq->resize(N);
        eq->setGrid(z);
        double errorEnthalpy  = eq->enthalpyBalanceError(yF);
        
        NJ++;
        if ( NJ > 1 )
        {
            if ( fabs(errorEnthalpy) <= input.getEnergyError() )
            {
                discretizationScheme = "CDS";
            }
        }
    }

    if ( input.getStart() == "new" || input.getStart() == "latest" )
    {
        eq->resize(N);
        eq->setGrid(z);
        std::vector<double> errorAtom      = eq->atomicBalanceError(yF);
        double              errorMass      = eq->massBalanceError(yF);
        double              errorEnthalpy  = eq->enthalpyBalanceError(yF);

        extra << "Results written in:        " << results << " fractions\n" << std::endl; 
        extra << "Discretization scheme:     " << discretizationScheme << "\n" << std::endl;  
        extra << "Number of grid point     = " << N << "\n" << std::endl; 
        extra << "ERROR on MASS BALANCE    = " << errorMass << "    [%]\n" << std::endl;
        extra << "ERROR on ENERGY BALANCE  = " << errorEnthalpy << "    [%]\n" << std::endl;
        extra << "\nATOMIC BALANCE               " << std::endl;
        for (unsigned int k=0;k<errorAtom.size();k++)
        {
            extra << thermodynamicsMapXML->elements()[k] << " :\t" << errorAtom[k] << "\t[%]" << std::endl;
        }
        extra << "\n" << std::endl;
        
        extra.close();

        std::cout << "\n#################################################" << std::endl;
        std::cout << "Results written in:        " << results << " fractions\n" << std::endl; 
        std::cout << "Number of grid point     = " << N << "\n" << std::endl; 
        std::cout << "ERROR on MASS BALANCE    = " << errorMass << "    [%]\n" << std::endl;
        std::cout << "ERROR on ENERGY BALANCE  = " << errorEnthalpy << "    [%]\n" << std::endl;
        std::cout << "#################################################" << std::endl;
        std::cout << "\n############################" << std::endl;
        std::cout << "ATOMIC BALANCE               " << std::endl;
        for (unsigned int k=0;k<errorAtom.size();k++)
        {
            std::cout << thermodynamicsMapXML->elements()[k] << " :\t" << errorAtom[k] << "\t[%]" << std::endl;
        }
        std::cout << "############################\n" << std::endl;
        
        
        OpenSMOKE::OpenSMOKEVectorDouble velocity = eq->getVelocity(yF);

        vel << "Velocity [m/s]" << std::endl;
        for (unsigned int i=1;i<=N;i++)
        {
            vel << velocity[i] << std::endl;
        }
        vel.close();
    }
    else if ( input.getStart() == "sampling" )
    {
        vel.close();
        extra.close();
        remove(pathVel);
        remove(pathExtra);
    }
    else if ( input.getStart() == "converter" )
    {
        extra << "Results written in:        " << results << " fractions\n" << std::endl;
        extra << "Number of grid point     = " << N << "\n" << std::endl;

        extra.close();
        vel.close();
        remove(pathVel);
    }
    else if ( input.getStart() == "kinetic" )
    {
        #include "kinetic.H"
    }
    else if ( input.getStart() == "RPA" )
    {
        #include "RPA.H"
    }
    else if ( input.getStart() == "conversion" )
    {
        std::ofstream conv;
        std::string fileNameConv = folderName + "/conversion.txt";
        const char *pathConv = fileNameConv.c_str();
        remove(pathConv);
        conv.open(pathConv,std::ios::out);
        conv.setf(std::ios::scientific);

        for (unsigned int j=1;j<=NC;j++)
        {
            if ( omegaBulk[0][j] > 1e-09 )
            {
                conv << thermodynamicsMapXML->NamesOfSpecies()[j-1] << "\t" << (omegaBulk[0][j] - omegaBulk[N-1][j])/omegaBulk[0][j] << std::endl;
            }
        }

        conv.close();
        remove(pathBulk);
        remove(pathWall);
        remove(pathSite);
        remove(pathTemp);
        remove(pathGrid);
        remove(pathExtra);
        remove(pathVel);
    }

    delete [] omegaBulk;
    delete [] omegaWall;
    delete [] teta;
}

