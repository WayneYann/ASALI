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
    vel.close();
    extra.close();
    remove(pathVel);
    remove(pathExtra);
    remove(pathBulk);
    remove(pathWall);
    remove(pathSite);
    remove(pathTemp);
}

{
    OpenSMOKE::OpenSMOKEVectorDouble Qbulk(N);
    OpenSMOKE::OpenSMOKEVectorDouble Qwall(N);

    OpenSMOKE::OpenSMOKEVectorDouble xBulk_(NC);
    OpenSMOKE::OpenSMOKEVectorDouble xWall_(NC);
    OpenSMOKE::OpenSMOKEVectorDouble cBulk_(NC);
    OpenSMOKE::OpenSMOKEVectorDouble cWall_(NC);
    OpenSMOKE::OpenSMOKEVectorDouble Rbulk_(NC);
    OpenSMOKE::OpenSMOKEVectorDouble Rwall_(NC);
    OpenSMOKE::OpenSMOKEVectorDouble Rteta_(SURF_NC);

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
                kineticsMapXML->FormationRates(&Rbulk_);
                Qbulk[i+1] = kineticsMapXML->HeatRelease(Rbulk_);
            }
            else
            {
                Qbulk[i+1] = 0.;
            }
            
            Qwall[i+1] = 0.;
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
                kineticsMapXML->FormationRates(&Rbulk_);
                Qbulk[i+1] = kineticsMapXML->HeatRelease(Rbulk_);
            }
            else
            {
                Qbulk[i+1] = 0.;
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
                kineticsSurfaceMapXML->ArrheniusKineticConstants();
                kineticsSurfaceMapXML->ReactionRates(cWall_, teta[i], fakebulk, fakephase);
                kineticsSurfaceMapXML->FormationRates(&Rwall_, &Rteta_, &fakebulk, &Rphase);
                Qwall[i+1] = kineticsSurfaceMapXML->HeatRelease(Rwall_, Rteta_, dummy);
            }
            else
            {
                Qwall[i+1] = 0.;
            }
        }
    }

    {
        std::ofstream heat;
        std::string fileNameHeat = folderName + "/heat.txt";
        const char *pathHeat = fileNameHeat.c_str();
        remove(pathHeat);
        heat.open(pathHeat,std::ios::out);

        {
            unsigned int WP = 6;
            heat.setf(std::ios::scientific);
            heat.precision(WP);
        }

        heat << "Bulk" << "\t" << "Wall" << std::endl;

        for (unsigned int k=0;k<N;k++)
        {
            heat << Qbulk[k+1] << "\t" << Qwall[k+1] << std::endl;
        }

        heat.close();
    }
}
