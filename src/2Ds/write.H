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
    std::string fileNameWall = folderName + "/temperature.txt";
    const char *pathWall = fileNameWall.c_str();
    remove(pathWall);
    temp.open(pathWall,std::ios::out);

    std::ofstream grid;
    std::string fileNameGrid = folderName + "/axial-coordinate.txt";
    const char *pathGrid = fileNameGrid.c_str();
    remove(pathGrid);
    grid.open(pathGrid,std::ios::out);

    std::ofstream rad;
    std::string fileNameRad = folderName + "/radial-coordinate.txt";
    const char *pathRad = fileNameRad.c_str();
    remove(pathRad);
    rad.open(pathRad,std::ios::out);

    std::ofstream twoD;
    std::string fileName2 = folderName + "/2D-temp.txt";
    const char *path2 = fileName2.c_str();
    remove(path2);
    twoD.open(path2,std::ios::out);


    {
        unsigned int WP = 6;
        mass.setf(std::ios::scientific);
        mole.setf(std::ios::scientific);
        temp.setf(std::ios::scientific);
        twoD.setf(std::ios::scientific);
        grid.setf(std::ios::scientific);
        rad.setf(std::ios::scientific);
        mass.precision(WP);
        mole.precision(WP);
        temp.precision(WP);
        twoD.precision(WP);
        grid.precision(WP);
        rad.precision(WP);
    }

    if (input.getModelName()[j] == "honeyComb")
    {
        unsigned int Na = input.axialPoints();
        unsigned int Nr = input.radialPoints();
        unsigned int NC = input.numberOfSpecies();
        OpenSMOKE::OpenSMOKEVectorDouble **omegaBulk;
        OpenSMOKE::OpenSMOKEVectorDouble  *TWall;
        OpenSMOKE::OpenSMOKEVectorDouble  *omega;
        OpenSMOKE::OpenSMOKEVectorDouble   Tg(Na);
        OpenSMOKE::OpenSMOKEVectorDouble   z(Na);
        OpenSMOKE::OpenSMOKEVectorDouble   d(Nr);

        omegaBulk   = new OpenSMOKE::OpenSMOKEVectorDouble*[Na];
        TWall       = new OpenSMOKE::OpenSMOKEVectorDouble[Na];
        omega       = new OpenSMOKE::OpenSMOKEVectorDouble[Na];

        for (unsigned int a=0;a<Na;a++)
        {
            omegaBulk[a] = new OpenSMOKE::OpenSMOKEVectorDouble[Nr];

            for (unsigned int r=0;r<Nr;r++)
            {
                ChangeDimensions(NC, &omegaBulk[a][r],     true);
            }
            
            ChangeDimensions(Nr, &TWall[a],           true);
            ChangeDimensions(NC, &omega[a],           true);
        }


        {
            unsigned int counter = 1;
            for (unsigned int a=0;a<Na;a++)
            {
                for (unsigned int r=0;r<Nr;r++)
                {
                    for(unsigned int s=1;s<=NC;s++)
                        omegaBulk[a][r][s] = yF[counter++];
                    for(unsigned int s=1;s<=NC;s++)
                        counter++;
                    counter++;
                    TWall[a][r+1] = yF[counter++];
                }
            }
        }

        {
            for (unsigned int a=0;a<Na;a++)
            {
                OpenSMOKE::OpenSMOKEVectorDouble v(Nr);
                for(unsigned int s=1;s<=NC;s++)
                {
                    for (unsigned int r=0;r<Nr;r++)
                    {
                        v[r+1] = omegaBulk[a][r][s];
                    }

                    omega[a][s] = v.SumElements()/v.Size();
                }
                
                for (unsigned int r=0;r<Nr;r++)
                {
                    v[r+1] = TWall[a][r+1];
                }

                Tg[a+1] = v.SumElements()/v.Size();
            }
        }
        
        {
            z[1] = 0.;
            for (unsigned int a=2;a<=Na;a++)
            {
                z[a] = z[a-1] + (input.getHoneyCombLength() + input.getInertLength())/double(Na-1);
            }
        }

        {
            d[1] = 0.;
            for (unsigned int r=2;r<=Nr;r++)
            {
                d[r] = d[r-1] + (input.getHoneyCombTubeDiameter()*0.5)/double(Nr-1);
            }
        }

        for (unsigned int s=1;s<=NC;s++)
            mass << input.getSpecieName()[s-1] << "\t";
        mass << std::endl;

        for (unsigned int s=1;s<=NC;s++)
            mole << input.getSpecieName()[s-1] << "\t";
        mole << std::endl;

        for (unsigned int a=0;a<Na;a++)
        {
            for (unsigned int s=1;s<=NC;s++)
                mass << omega[a][s] << "\t";
            mass << std::endl;

            OpenSMOKE::OpenSMOKEVectorDouble x(NC);
            x = ASALI::fromMassFractionToMoleFraction(omega[a],input.getMW());
            for (unsigned int s=1;s<=NC;s++)
                mole << x[s] << "\t";
            mole << std::endl;

            temp << Tg[a+1] << std::endl;
            
            grid << z[a+1] << std::endl;
            
            for (unsigned int r=0;r<Nr;r++)
                twoD << TWall[a][r+1] << "\t";
            twoD << std::endl;
        }
        
        for (unsigned int r=0;r<Nr;r++)
        {
            rad << d[r+1] << std::endl;
        }
    }
    else if (input.getModelName()[j] == "packedBed")
    {
        unsigned int Na = input.axialPoints();
        unsigned int Nr = input.radialPoints();
        unsigned int NC = input.numberOfSpecies();
        OpenSMOKE::OpenSMOKEVectorDouble **omegaBulk;
        OpenSMOKE::OpenSMOKEVectorDouble  *TBulk;
        OpenSMOKE::OpenSMOKEVectorDouble  *omega;
        OpenSMOKE::OpenSMOKEVectorDouble   Tg(Na);
        OpenSMOKE::OpenSMOKEVectorDouble   z(Na);
        OpenSMOKE::OpenSMOKEVectorDouble   d(Nr);

        omegaBulk   = new OpenSMOKE::OpenSMOKEVectorDouble*[Na];
        TBulk       = new OpenSMOKE::OpenSMOKEVectorDouble[Na];
        omega       = new OpenSMOKE::OpenSMOKEVectorDouble[Na];

        for (unsigned int a=0;a<Na;a++)
        {
            omegaBulk[a] = new OpenSMOKE::OpenSMOKEVectorDouble[Nr];

            for (unsigned int r=0;r<Nr;r++)
            {
                ChangeDimensions(NC, &omegaBulk[a][r],     true);
            }
            
            ChangeDimensions(Nr, &TBulk[a],           true);
            ChangeDimensions(NC, &omega[a],           true);
        }


        {
            unsigned int counter = 1;
            for (unsigned int a=0;a<Na;a++)
            {
                for (unsigned int r=0;r<Nr;r++)
                {
                    for(unsigned int s=1;s<=NC;s++)
                        omegaBulk[a][r][s] = yF[counter++];
                    TBulk[a][r+1] = yF[counter++];
                }
            }
        }

        {
            for (unsigned int a=0;a<Na;a++)
            {
                OpenSMOKE::OpenSMOKEVectorDouble v(Nr);
                for(unsigned int s=1;s<=NC;s++)
                {
                    for (unsigned int r=0;r<Nr;r++)
                    {
                        v[r+1] = omegaBulk[a][r][s];
                    }

                    omega[a][s] = v.SumElements()/v.Size();
                }
                
                for (unsigned int r=0;r<Nr;r++)
                {
                    v[r+1] = TBulk[a][r+1];
                }

                Tg[a+1] = v.SumElements()/v.Size();
            }
        }
        
        {
            z[1] = 0.;
            for (unsigned int a=2;a<=Na;a++)
            {
                z[a] = z[a-1] + (input.getPackedBedLength() + input.getInertLength())/double(Na-1);
            }
        }

        {
            d[1] = 0.;
            for (unsigned int r=2;r<=Nr;r++)
            {
                d[r] = d[r-1] + (input.getPackedBedTubeDiameter()*0.5)/double(Nr-1);
            }
        }

        for (unsigned int s=1;s<=NC;s++)
            mass << input.getSpecieName()[s-1] << "\t";
        mass << std::endl;

        for (unsigned int s=1;s<=NC;s++)
            mole << input.getSpecieName()[s-1] << "\t";
        mole << std::endl;

        for (unsigned int a=0;a<Na;a++)
        {
            for (unsigned int s=1;s<=NC;s++)
                mass << omega[a][s] << "\t";
            mass << std::endl;

            OpenSMOKE::OpenSMOKEVectorDouble x(NC);
            x = ASALI::fromMassFractionToMoleFraction(omega[a],input.getMW());
            for (unsigned int s=1;s<=NC;s++)
                mole << x[s] << "\t";
            mole << std::endl;

            temp << Tg[a+1] << std::endl;
            
            grid << z[a+1] << std::endl;

            for (unsigned int r=0;r<Nr;r++)
                twoD << TBulk[a][r+1] << "\t";
            twoD << std::endl;
        }
        
        for (unsigned int r=0;r<Nr;r++)
        {
            rad << d[r+1] << std::endl;
        }
    }
    else if (input.getModelName()[j] == "microBed")
    {
        unsigned int Na = input.axialPoints();
        unsigned int Nr = input.radialPoints();
        unsigned int NC = input.numberOfSpecies();
        OpenSMOKE::OpenSMOKEVectorDouble **omegaBulk;
        OpenSMOKE::OpenSMOKEVectorDouble  *TBulk;
        OpenSMOKE::OpenSMOKEVectorDouble  *omega;
        OpenSMOKE::OpenSMOKEVectorDouble   Tg(Na);
        OpenSMOKE::OpenSMOKEVectorDouble   z(Na);
        OpenSMOKE::OpenSMOKEVectorDouble   d(Nr);

        omegaBulk   = new OpenSMOKE::OpenSMOKEVectorDouble*[Na];
        TBulk       = new OpenSMOKE::OpenSMOKEVectorDouble[Na];
        omega       = new OpenSMOKE::OpenSMOKEVectorDouble[Na];

        for (unsigned int a=0;a<Na;a++)
        {
            omegaBulk[a] = new OpenSMOKE::OpenSMOKEVectorDouble[Nr];

            for (unsigned int r=0;r<Nr;r++)
            {
                ChangeDimensions(NC, &omegaBulk[a][r],     true);
            }
            
            ChangeDimensions(Nr, &TBulk[a],           true);
            ChangeDimensions(NC, &omega[a],           true);
        }


        {
            unsigned int counter = 1;
            for (unsigned int a=0;a<Na;a++)
            {
                for (unsigned int r=0;r<Nr;r++)
                {
                    for(unsigned int s=1;s<=NC;s++)
                        omegaBulk[a][r][s] = yF[counter++];
                    TBulk[a][r+1] = yF[counter++];
                    counter++;
                }
            }
        }

        {
            for (unsigned int a=0;a<Na;a++)
            {
                OpenSMOKE::OpenSMOKEVectorDouble v(Nr);
                for(unsigned int s=1;s<=NC;s++)
                {
                    for (unsigned int r=0;r<Nr;r++)
                    {
                        v[r+1] = omegaBulk[a][r][s];
                    }

                    omega[a][s] = v.SumElements()/v.Size();
                }
                
                for (unsigned int r=0;r<Nr;r++)
                {
                    v[r+1] = TBulk[a][r+1];
                }

                Tg[a+1] = v.SumElements()/v.Size();
            }
        }
        
        {
            z[1] = 0.;
            for (unsigned int a=2;a<=Na;a++)
            {
                z[a] = z[a-1] + (input.getMicroBedLength() + input.getInertLength())/double(Na-1);
            }
        }

        {
            d[1] = 0.;
            for (unsigned int r=2;r<=Nr;r++)
            {
                d[r] = d[r-1] + (input.getMicroBedTubeDiameter()*0.5)/double(Nr-1);
            }
        }

        for (unsigned int s=1;s<=NC;s++)
            mass << input.getSpecieName()[s-1] << "\t";
        mass << std::endl;

        for (unsigned int s=1;s<=NC;s++)
            mole << input.getSpecieName()[s-1] << "\t";
        mole << std::endl;

        for (unsigned int a=0;a<Na;a++)
        {
            for (unsigned int s=1;s<=NC;s++)
                mass << omega[a][s] << "\t";
            mass << std::endl;

            OpenSMOKE::OpenSMOKEVectorDouble x(NC);
            x = ASALI::fromMassFractionToMoleFraction(omega[a],input.getMW());
            for (unsigned int s=1;s<=NC;s++)
                mole << x[s] << "\t";
            mole << std::endl;

            temp << Tg[a+1] << std::endl;
            
            grid << z[a+1] << std::endl;
            for (unsigned int r=0;r<Nr;r++)
                twoD << TBulk[a][r+1] << "\t";
            twoD << std::endl;
        }
        
        for (unsigned int r=0;r<Nr;r++)
        {
            rad << d[r+1] << std::endl;
        }
    }
    mass.close();
    mole.close();
    temp.close();
    grid.close();
    twoD.close();
    rad.close();
}

{
    std::string folderName = "results/" + input.getModelName()[j];

    std::ofstream info;
    std::string fileNameInfo = folderName + "/geometry.txt";
    const char *pathInfo = fileNameInfo.c_str();
    remove(pathInfo);
    info.open(pathInfo,std::ios::out);
    
    if (input.getModelName()[j] == "honeyComb")
    {
        info << "Tube diameter:      " << input.getHoneyCombTubeDiameter() << " m" << std::endl;
        info << "CPSI:               " << input.getHoneyCombCPSI() << std::endl;
        info << "Wall thickness:     " << input.getHoneyCombWall() << " mills" << std::endl;
        info << "Washcoat thickness: " << input.getHoneyCombWashCoat() << " m" << std::endl;
        info << "Catalytic lenght:   " << input.getHoneyCombLength() << " m" << std::endl;
        info << "Inert lenght:       " << input.getInertLength() << " m" << std::endl;
        info << "Type:               " << input.getHoneyCombType() << std::endl;
    }
    else if (input.getModelName()[j] == "packedBed")
    {
        info << "Tube diameter:      " << input.getPackedBedTubeDiameter() << " m" << std::endl;
        info << "Particle diameter:  " << input.getPackedBedParticleDiameter() << " m" << std::endl;
        info << "Catalytic lenght:   " << input.getPackedBedLength() << " m" << std::endl;
        info << "Inert lenght:       " << input.getInertLength() << " m" << std::endl;
    }
    else if (input.getModelName()[j] == "microBed")
    {
        info << "Tube diameter:      " << input.getMicroBedTubeDiameter() << " m" << std::endl;
        info << "CPSI:               " << input.getMicroBedCPSI() << std::endl;
        info << "Particle diameter:  " << input.getMicroBedParticleDiameter() << " m" << std::endl;
        info << "Washcoat thickness: " << input.getMicroBedWall() << " m" << std::endl;
        info << "Catalytic lenght:   " << input.getMicroBedLength() << " m" << std::endl;
        info << "Inert lenght:       " << input.getInertLength() << " m" << std::endl;
    }

    info.close();
}

