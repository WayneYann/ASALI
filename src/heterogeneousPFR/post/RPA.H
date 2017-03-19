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
    unsigned int index     = thermodynamicsSurfaceMapXML->IndexOfSpecies(*options);
    unsigned int Nreaction = kineticsSurfaceMapXML->NumberOfReactions();
    unsigned int SURF_NP   = thermodynamicsSurfaceMapXML->number_of_site_phases(0);

    OpenSMOKE::OpenSMOKEVectorDouble *stoichiometricCoeff;
    OpenSMOKE::OpenSMOKEVectorDouble *rNet;
    OpenSMOKE::OpenSMOKEVectorDouble *rProd;
    OpenSMOKE::OpenSMOKEVectorDouble *rCons;
    OpenSMOKE::OpenSMOKEVectorDouble *ratioProd;
    OpenSMOKE::OpenSMOKEVectorDouble *ratioCons;

    OpenSMOKE::OpenSMOKEVectorDouble rProdSum(N);
    OpenSMOKE::OpenSMOKEVectorDouble rConsSum(N);
    OpenSMOKE::OpenSMOKEVectorDouble xBulk_(NC);
    OpenSMOKE::OpenSMOKEVectorDouble xWall_(NC);
    OpenSMOKE::OpenSMOKEVectorDouble cBulk_(NC);
    OpenSMOKE::OpenSMOKEVectorDouble cWall_(NC);

    stoichiometricCoeff = new OpenSMOKE::OpenSMOKEVectorDouble[Nreaction];
    rNet                = new OpenSMOKE::OpenSMOKEVectorDouble[N];
    rProd               = new OpenSMOKE::OpenSMOKEVectorDouble[N];
    rCons               = new OpenSMOKE::OpenSMOKEVectorDouble[N];
    ratioProd           = new OpenSMOKE::OpenSMOKEVectorDouble[N];
    ratioCons           = new OpenSMOKE::OpenSMOKEVectorDouble[N];

    for (unsigned int i=0;i<Nreaction;i++)
    {
        ChangeDimensions(NC + SURF_NC, &stoichiometricCoeff[i], true);
    }

    for (unsigned int i=0;i<N;i++)
    {
        ChangeDimensions(int(Nreaction/2), &rNet[i],      true);
        ChangeDimensions(int(Nreaction/2), &rProd[i],     true);
        ChangeDimensions(int(Nreaction/2), &rCons[i],     true);
        ChangeDimensions(int(Nreaction/2), &ratioProd[i], true);
        ChangeDimensions(int(Nreaction/2), &ratioCons[i], true);
    }

    //- Stoichiometry
    {
        kineticsSurfaceMapXML->stoichiometry().BuildStoichiometricMatrix();
        for (unsigned int i=0; i<kineticsSurfaceMapXML->stoichiometry().stoichiometric_matrix_products().outerSize();++i)
        {
            for (Eigen::SparseMatrix<double>::InnerIterator it(kineticsSurfaceMapXML->stoichiometry().stoichiometric_matrix_products(),i); it;++it)
            {
                stoichiometricCoeff[it.index()][i+1]=it.value();
            }
        }
        for (unsigned int i=0; i<kineticsSurfaceMapXML->stoichiometry().stoichiometric_matrix_reactants().outerSize();++i)
        {
            for (Eigen::SparseMatrix<double>::InnerIterator it(kineticsSurfaceMapXML->stoichiometry().stoichiometric_matrix_reactants(),i); it;++it)
            {
                stoichiometricCoeff[it.index()][i+1]=-it.value();
            }
        }
    }

    //- Reaction rates
    {
        for (unsigned int i=0;i<N;i++)
        {
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
            
            OpenSMOKE::OpenSMOKEVectorDouble R(Nreaction);
            OpenSMOKE::OpenSMOKEVectorDouble Rwall(NC);
            OpenSMOKE::OpenSMOKEVectorDouble Rteta(SURF_NC);
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
                kineticsSurfaceMapXML->GetReactionRates(&R);
            }
            
            unsigned int counter = 0.;
            for (unsigned int k=1;k<=int(Nreaction/2);k++)
            {
                rNet[i][k] = stoichiometricCoeff[counter][index]*R[counter+1] - stoichiometricCoeff[counter+1][index]*R[counter+2];
                counter = counter + 2;

                if ( rNet[i][k] > 0 )
                {
                    rProd[i][k] = std::fabs(rNet[i][k]);
                    rCons[i][k] = 0.;
                }
                else
                {
                    rProd[i][k] = 0.;
                    rCons[i][k] = std::fabs(rNet[i][k]);
                }
                
                rProdSum[i] += rProd[i][k];
                rConsSum[i] += rCons[i][k];
            }

            for (unsigned int k=1;k<=int(Nreaction/2);k++)
            {
                if ( rProdSum[i] == 0. )
                {
                    ratioProd[i][k] = 0.;
                }
                else if ( rProd[i][k] == 0. )
                {
                    ratioProd[i][k] = 0.;
                }
                else
                {
                    ratioProd[i][k] = rProd[i][k]/rProdSum[i];
                }
                
                if ( rConsSum[i] == 0. )
                {
                    ratioCons[i][k] = 0.;
                }
                else if ( rCons[i][k] == 0. )
                {
                    ratioCons[i][k] = 0.;
                }
                else
                {
                    ratioCons[i][k] = rCons[i][k]/rConsSum[i];
                }
            }
        }
    }

    {
        remove(pathBulk);
        remove(pathWall);
        remove(pathSite);
        remove(pathTemp);
        remove(pathVel);
        remove(pathExtra);
        vel.close();
    }

    {
        std::string folderNameProd = "results/" + ASALI::to_string(N) + ".RPA/production";
        std::string folderNameCons = "results/" + ASALI::to_string(N) + ".RPA/consumption";
        boost::filesystem::create_directory(folderNameProd);
        boost::filesystem::create_directory(folderNameCons);

        std::ofstream prod;
        std::string fileNameProd = folderNameProd + "/" + *options + ".txt";
        const char *pathProd = fileNameProd.c_str();
        remove(pathProd);
        prod.open(pathProd,std::ios::out);

        std::ofstream cons;
        std::string fileNameCons = folderNameCons + "/" + *options + ".txt";
        const char *pathCons = fileNameCons.c_str();
        remove(pathCons);
        cons.open(pathCons,std::ios::out);

        {
            unsigned int WP = 3;
            prod.setf(std::ios::fixed);
            cons.setf(std::ios::fixed);
            prod.precision(WP);
            cons.precision(WP);
        }

        for (unsigned int k=1;k<=int(Nreaction/2);k++)
        {
            prod << k << "-" << k+1 << "\t";
            cons << k << "-" << k+1 << "\t";
        }
        prod << std::endl;
        cons << std::endl;

        for (unsigned int i=0;i<N;i++)
        {
            for (unsigned int k=1;k<=int(Nreaction/2);k++)
            {
                prod << ratioProd[i][k] << "\t";
                cons << ratioCons[i][k] << "\t";
            }
            prod << std::endl;
            cons << std::endl;
        }

        prod.close();
        cons.close();
        vel.close();
        remove(pathVel);
    }

}
