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
    boost::filesystem::create_directory("results");
    std::string folderName = "results/" + input.getModelName()[j];
    boost::filesystem::create_directory(folderName);

    std::ofstream mass;
    std::string fileNameMass = folderName + "/mass-fraction.txt";
    const char *pathMass = fileNameMass.c_str();
    remove(pathMass);
    mass.open(pathMass,std::ios::out);

    std::ofstream mole;
    std::string fileNameMole = folderName + "/mole-fraction.txt";
    const char *pathMole = fileNameMole.c_str();
    remove(pathMole);
    mole.open(pathMole,std::ios::out);

    std::ofstream temp;
    std::string fileNameWall = folderName + "/temperture.txt";
    const char *pathWall = fileNameWall.c_str();
    remove(pathWall);
    temp.open(pathWall,std::ios::out);

    std::ofstream grid;
    std::string fileNameGrid = folderName + "/grid.txt";
    const char *pathGrid = fileNameGrid.c_str();
    remove(pathGrid);
    grid.open(pathGrid,std::ios::out);

    {
        unsigned int WP = 6;
        mass.setf(std::ios::scientific);
        mole.setf(std::ios::scientific);
        temp.setf(std::ios::scientific);
        grid.setf(std::ios::scientific);
        mass.precision(WP);
        mole.precision(WP);
        temp.precision(WP);
        grid.precision(WP);
    }

    if (input.getModelName()[j] == "heterogeneous")
    {
        unsigned int NC = input.numberOfSpecies();
        std::vector<double> Tb;
        std::vector<double> omegab;
        std::vector<double> Tw;
        std::vector<double> omegaw;
        std::vector<double> z;
        eq->GetDaeProfile(z,Tb,Tw,omegab,omegaw);

        OpenSMOKE::OpenSMOKEVectorDouble *omegaBulk;
        OpenSMOKE::OpenSMOKEVectorDouble *omegaWall;

        omegaBulk = new OpenSMOKE::OpenSMOKEVectorDouble[z.size()];
        omegaWall = new OpenSMOKE::OpenSMOKEVectorDouble[z.size()];

        for (unsigned int i=0;i<z.size();i++)
        {
            ChangeDimensions(NC, &omegaBulk[i], true);
            ChangeDimensions(NC, &omegaWall[i], true);
        }

        {
            unsigned int conter = 0;
            for (unsigned int k=0;k<z.size();k++)
            {
                for (unsigned int i=1;i<=NC;i++)
                {
                    omegaBulk[k][i] = omegab[conter++];
                    omegaWall[k][i] = omegaw[conter++];
                }
            }
        }

        for (unsigned int s=1;s<=NC;s++)
        {
            mass << input.getSpecieName()[s-1] << "-gas\t";
            mass << input.getSpecieName()[s-1] << "-cat\t";
            mole << input.getSpecieName()[s-1] << "-gas\t";
            mole << input.getSpecieName()[s-1] << "-cat\t";
        }
        mass << std::endl;
        mole << std::endl;
        
        temp << "gas\tcat" << std::endl;

        for (unsigned int a=0;a<z.size();a++)
        {
            OpenSMOKE::OpenSMOKEVectorDouble xb(NC);
            OpenSMOKE::OpenSMOKEVectorDouble xw(NC);
            xb = ASALI::FromMassToMoleFraction(omegaBulk[a],input.getMW());
            xw = ASALI::FromMassToMoleFraction(omegaWall[a],input.getMW());
            for (unsigned int s=1;s<=NC;s++)
            {
                mass << omegaBulk[a][s] << "\t";
                mass << omegaWall[a][s] << "\t";
                mole << xb[s] << "\t";
                mole << xw[s] << "\t";
            }
            mass << std::endl;
            mole << std::endl;

            temp << Tb[a] << "\t" << Tw[a] << std::endl;
            
            grid << z[a] << std::endl;
        }
    }
    else if (input.getModelName()[j] == "pseudoHomogeneous")
    {
        unsigned int NC = input.numberOfSpecies();
        std::vector<double> T;
        std::vector<double> comp;
        std::vector<double> z;
        eq->GetOdeProfile(z,T,comp);

        OpenSMOKE::OpenSMOKEVectorDouble *omega;

        omega = new OpenSMOKE::OpenSMOKEVectorDouble[z.size()];

        for (unsigned int i=0;i<z.size();i++)
        {
            ChangeDimensions(NC, &omega[i], true);
        }

        {
            unsigned int conter = 0;
            for (unsigned int k=0;k<z.size();k++)
            {
                for (unsigned int i=1;i<=NC;i++)
                {
                    omega[k][i] = comp[conter++];
                }
            }
        }

        for (unsigned int s=1;s<=NC;s++)
        {
            mass << input.getSpecieName()[s-1] << "\t";
            mole << input.getSpecieName()[s-1] << "\t";
        }
        mass << std::endl;
        mole << std::endl;

        for (unsigned int a=0;a<z.size();a++)
        {
            OpenSMOKE::OpenSMOKEVectorDouble x(NC);
            x = ASALI::FromMassToMoleFraction(omega[a],input.getMW());
            for (unsigned int s=1;s<=NC;s++)
            {
                mass << omega[a][s] << "\t";
                mole << x[s] << "\t";
            }
            mass << std::endl;
            mole << std::endl;

            temp << T[a] << std::endl;
            
            grid << z[a] << std::endl;
        }
    }

    mass.close();
    mole.close();
    temp.close();
    grid.close();
}

