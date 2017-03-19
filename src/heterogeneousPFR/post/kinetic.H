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
    OpenSMOKE::OpenSMOKEVectorDouble *Rbulk;
    OpenSMOKE::OpenSMOKEVectorDouble *Rwall;
    OpenSMOKE::OpenSMOKEVectorDouble *Rteta;

    Rbulk = new OpenSMOKE::OpenSMOKEVectorDouble[N];
    Rwall = new OpenSMOKE::OpenSMOKEVectorDouble[N];
    Rteta = new OpenSMOKE::OpenSMOKEVectorDouble[N];

    for (unsigned int i=0;i<N;i++)
    {
        ChangeDimensions(NC,      &Rbulk[i], true);
        ChangeDimensions(NC,      &Rwall[i], true);
        ChangeDimensions(SURF_NC, &Rteta[i], true);
    }

    OpenSMOKE::OpenSMOKEVectorDouble xBulk_(NC);
    OpenSMOKE::OpenSMOKEVectorDouble xWall_(NC);
    OpenSMOKE::OpenSMOKEVectorDouble cBulk_(NC);
    OpenSMOKE::OpenSMOKEVectorDouble cWall_(NC);

    unsigned int SURF_NP = thermodynamicsSurfaceMapXML->number_of_site_phases(0);

    for (unsigned int i=0;i<N;i++)
    {
        if ( z[i+1] <= input.getInertLength())
        {
            //-Setting conditions for GAS phase
            thermodynamicsMapXML->SetTemperature(Tbulk[i+1]);
            thermodynamicsMapXML->SetPressure(input.getPressure());
            kineticsMapXML->SetTemperature(Tbulk[i+1]);
            kineticsMapXML->SetPressure(input.getPressure());
            thermodynamicsSurfaceMapXML->SetPressure(input.getPressure());
            thermodynamicsSurfaceMapXML->SetTemperature(Tbulk[i+1]);
            kineticsSurfaceMapXML->SetPressure(input.getPressure());
            kineticsSurfaceMapXML->SetTemperature(Tbulk[i+1]);
            transportMapXML->SetPressure(input.getPressure());
            transportMapXML->SetTemperature(Tbulk[i+1]);

            //-Evaluating GAS phase properties
            double MWbulk_;
            thermodynamicsMapXML->MoleFractions_From_MassFractions(xBulk_,MWbulk_,omegaBulk[i]);

            double cTotBulk_ = input.getPressure()/(PhysicalConstants::R_J_kmol*Tbulk[i+1]);     // [Kmol/m3]
            Product(cTotBulk_, xBulk_, &cBulk_);

            if ( input.getHomogeneousReactions() == true )
            {
                kineticsMapXML->ReactionRates(cBulk_);
                kineticsMapXML->FormationRates(&Rbulk[i]);
                ElementByElementProduct(Rbulk[i], thermodynamicsMapXML->MW(), &Rbulk[i]);
            }
            else
            {
                for (unsigned int j=1;j<=NC;j++)
                    Rbulk[i][j] = 0.;
            }

            for (unsigned int j=1;j<=NC;j++)
            {
                Rwall[i][j] = 0.;
            }
            
            for (unsigned int j=1;j<=SURF_NC;j++)
            {
                Rteta[i][j] = 0.;
            }
        }
        else
        {
            //-Setting conditions for GAS phase
            thermodynamicsMapXML->SetTemperature(Tbulk[i+1]);
            thermodynamicsMapXML->SetPressure(input.getPressure());
            kineticsMapXML->SetTemperature(Tbulk[i+1]);
            kineticsMapXML->SetPressure(input.getPressure());
            thermodynamicsSurfaceMapXML->SetPressure(input.getPressure());
            thermodynamicsSurfaceMapXML->SetTemperature(Tbulk[i+1]);
            kineticsSurfaceMapXML->SetPressure(input.getPressure());
            kineticsSurfaceMapXML->SetTemperature(Tbulk[i+1]);
            transportMapXML->SetPressure(input.getPressure());
            transportMapXML->SetTemperature(Tbulk[i+1]);

            //-Evaluating GAS phase properties
            double MWbulk_;
            thermodynamicsMapXML->MoleFractions_From_MassFractions(xBulk_,MWbulk_,omegaBulk[i]);

            double cTotBulk_ = input.getPressure()/(PhysicalConstants::R_J_kmol*Tbulk[i+1]);     // [Kmol/m3]
            Product(cTotBulk_, xBulk_, &cBulk_);

            if ( input.getHomogeneousReactions() == true )
            {
                kineticsMapXML->ReactionRates(cBulk_);
                kineticsMapXML->FormationRates(&Rbulk[i]);
                ElementByElementProduct(Rbulk[i], thermodynamicsMapXML->MW(), &Rbulk[i]);
            }
            else
            {
                for (unsigned int j=1;j<=NC;j++)
                    Rbulk[i][j] = 0.;
            }
            
            //-Setting conditions for SOLID phase
            thermodynamicsMapXML->SetTemperature(Twall[i+1]);
            thermodynamicsMapXML->SetPressure(input.getPressure());
            kineticsMapXML->SetTemperature(Twall[i+1]);
            kineticsMapXML->SetPressure(input.getPressure());
            thermodynamicsSurfaceMapXML->SetPressure(input.getPressure());
            thermodynamicsSurfaceMapXML->SetTemperature(Twall[i+1]);
            kineticsSurfaceMapXML->SetPressure(input.getPressure());
            kineticsSurfaceMapXML->SetTemperature(Twall[i+1]);
            transportMapXML->SetPressure(input.getPressure());
            transportMapXML->SetTemperature(Twall[i+1]);
            
            //-Evaluating properties for SOLID phase
            double MWwall_;
            thermodynamicsMapXML->MoleFractions_From_MassFractions(xWall_, MWwall_, omegaWall[i]);

            double cTotWall_ = input.getPressure()/(PhysicalConstants::R_J_kmol*Twall[i+1]);
            Product(cTotWall_, xWall_, &cWall_); 
            
            if ( input.getHeterogenousReactions() == true )
            {
                OpenSMOKE::OpenSMOKEVectorDouble fakephase(SURF_NP);
                OpenSMOKE::OpenSMOKEVectorDouble Rphase(SURF_NP);
                OpenSMOKE::OpenSMOKEVectorDouble fakebulk(NC);
                OpenSMOKE::OpenSMOKEVectorDouble dummy;
                for(unsigned int j=1;j<=SURF_NP;j++)
                    fakephase[j] = thermodynamicsSurfaceMapXML->matrix_densities_site_phases()[0][j-1];

                kineticsSurfaceMapXML->ReactionEnthalpiesAndEntropies();
                kineticsSurfaceMapXML->ReactionRates(cWall_, teta[i], fakebulk, fakephase);
                kineticsSurfaceMapXML->FormationRates(&Rwall[i], &Rteta[i], &fakebulk, &Rphase);
                ElementByElementProduct(Rwall[i], thermodynamicsMapXML->MW(), &Rwall[i]);
            }
            else
            {
                for (unsigned int j=1;j<=NC;j++)
                    Rwall[i][j] = 0.;
                for (unsigned int j=1;j<=SURF_NC;j++)
                    Rteta[i][j] = 0.;
            }
        }
    }

    {
        remove(pathBulk);
        remove(pathWall);
        remove(pathSite);
        remove(pathTemp);
        remove(pathVel);
        vel.close();
    }

    {
        std::ofstream bulkReaction;
        std::string fileNameBulkReaction = folderName + "/bulkReactions.txt";
        const char *pathBulkReaction = fileNameBulkReaction.c_str();
        remove(pathBulkReaction);
        bulkReaction.open(pathBulkReaction,std::ios::out);

        std::ofstream wallReaction;
        std::string fileNameWallReaction = folderName + "/wallReactions.txt";
        const char *pathWallReaction = fileNameWallReaction.c_str();
        remove(pathWallReaction);
        wallReaction.open(pathWallReaction,std::ios::out);

        std::ofstream siteReaction;
        std::string fileNameSiteReaction = folderName + "/siteReactions.txt";
        const char *pathSiteReaction = fileNameSiteReaction.c_str();
        remove(pathSiteReaction);
        siteReaction.open(pathSiteReaction,std::ios::out);

        {
            unsigned int WP = 6;
            bulkReaction.setf(std::ios::scientific);
            wallReaction.setf(std::ios::scientific);
            siteReaction.setf(std::ios::scientific);
            bulkReaction.precision(WP);
            wallReaction.precision(WP);
            siteReaction.precision(WP);
        }

        for (unsigned int j=1;j<=NC;j++)
            bulkReaction << thermodynamicsMapXML->NamesOfSpecies()[j-1] << "\t";
        bulkReaction << std::endl;

        for (unsigned int j=1;j<=NC;j++)
            wallReaction << thermodynamicsMapXML->NamesOfSpecies()[j-1] << "\t";
        wallReaction << std::endl;

        for (unsigned int j=1;j<=SURF_NC;j++)
            siteReaction << thermodynamicsSurfaceMapXML->NamesOfSpecies()[j-1+NC] << "\t";
        siteReaction << std::endl;

        for (unsigned int k=0;k<N;k++)
        {
            for (unsigned int j=1;j<=NC;j++)
                bulkReaction << Rbulk[k][j] << "\t";
            bulkReaction << std::endl;

            for (unsigned int j=1;j<=NC;j++)
                wallReaction << Rwall[k][j] << "\t";
            wallReaction << std::endl;

            for (unsigned int j=1;j<=SURF_NC;j++)
                siteReaction << Rteta[k][j] << "\t";
            siteReaction << std::endl;
        }

        bulkReaction.close();
        wallReaction.close();
        siteReaction.close();
        vel.close();
        remove(pathVel);

        extra << "Unit dimensions:\n" << std::endl;
        extra << "BULK reaction rates:   Kg/m3/s" << std::endl;
        extra << "WALL reaction rates:   Kg/m2/s" << std::endl;
        extra << "SITE reaction rates: Kmol/m2/s" << std::endl;
        extra.close();
    }

    delete [] Rbulk;
    delete [] Rwall;
    delete [] Rteta;
}
