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

OpenSMOKE::OpenSMOKEVectorDouble ReactionRate(const OpenSMOKE::OpenSMOKEVectorDouble omega,const double T)
{
    OpenSMOKE::OpenSMOKEVectorDouble R(NC_);
    
    if ( reactionType_ == "O-xylene-to-phthalic" )
    {
        OpenSMOKE::OpenSMOKEVectorDouble s(NC_);
        s[1] = -3.;
        s[2] = -1.;
        s[3] =  1.;
        s[4] =  3.;
        s[5] =  0.;

        OpenSMOKE::OpenSMOKEVectorDouble p(NC_);

        double MWmix = mixtureMolecularWeight(omega);
        for (unsigned int i=1;i<=NC_;i++)
        {
            p[i] = omega[i]*MWmix*p_/(1.e05*MW_[i]);
        }

        double k   = std::exp(19.837 - 13636./T);
        double kp  = k*p[1]*8.314*T/1e05;
        double De  = diff_*epsiC_/tauC_;
        double eff = std::min(1.,1./(Lgeo_*std::sqrt(kp*rhoC_/De)));

        for (unsigned int i=1;i<=NC_;i++)
        {
            R[i] = eff*s[i]*k*p[1]*p[2]*rhoC_*MW_[i]/3600.; //[Kg/m3cat/s]
        }
    }
    else if ( reactionType_ == "O-xylene-to-phthalic-complex" )
    {
        unsigned int NR = 3;
        OpenSMOKE::OpenSMOKEVectorDouble r(NR);
        OpenSMOKE::OpenSMOKEVectorDouble eff(NR);
        OpenSMOKE::OpenSMOKEVectorDouble p(NC_);
        {
            double k1    = std::exp(-13588./T + 19.837);
            double k2    = std::exp(-15803./T + 20.860);
            double k3    = std::exp(-14394./T + 18.970);
            double MWmix = mixtureMolecularWeight(omega);

            for (unsigned int i=1;i<=NC_;i++)
            {
                p[i] = omega[i]*MWmix*p_/(1.e05*MW_[i]);
            }

            unsigned int iO2 = 1;
            unsigned int iXY = 2;
            unsigned int iPA = 3;
            r[1] = k1*p[iO2]*p[iXY];                //[Kmol/Kgcat/s]
            r[2] = k2*p[iPA]*p[iO2];                //[Kmol/Kgcat/s]
            r[3] = k3*p[iO2]*p[iXY];                //[Kmol/Kgcat/s]

            double kp1  = k1*p[iO2]*8.314*T/1e05;
            double kp2  = k2*p[iO2]*8.314*T/1e05;
            double kp3  = k3*p[iO2]*8.314*T/1e05;
            double De   = diff_*epsiC_/tauC_;
            eff[1]      = std::min(1.,1./(Lgeo_*std::sqrt(kp1*rhoC_/De)));
            eff[2]      = std::min(1.,1./(Lgeo_*std::sqrt(kp2*rhoC_/De)));
            eff[3]      = std::min(1.,1./(Lgeo_*std::sqrt(kp3*rhoC_/De)));
        }

        OpenSMOKE::OpenSMOKEVectorDouble *s;
        s = new OpenSMOKE::OpenSMOKEVectorDouble[NC_];
        {
            for (unsigned int i=0;i<NC_;i++)
            {
                ChangeDimensions(NR, &s[i], true);
                
                for (unsigned int j=1;j<=NR;j++)
                {
                    s[i][j] = 0.;
                }
            }

            // XYLENE + 3O2 -> PHTHALIC + 3H2O [1 + 0 -> 2 + 5]
            {
                s[1][1] = -1.;
                s[0][1] = -3.;
                s[2][1] =  1.;
                s[5][1] =  3.;
            }

            // PHTHALIC + (15/2)O2 -> 8CO2 + 2H2O [2 + 0 -> 3 + 5]
            {
                s[2][2] = -1.;
                s[0][2] = -15./2.;
                s[3][2] =  8.;
                s[5][2] =  2.;
            }

            // XYLENE + (21/2)O2 -> 8CO2 + 5H2O [1 + 0 -> 3 + 5]
            {
                s[1][3] = -1.;
                s[0][3] = -21./2.;
                s[3][3] =  9.;
                s[5][3] =  5.;
            }

            for (unsigned int i=1;i<=NC_;i++)
            {
                R[i] = 0.;
                for (unsigned int j=1;j<=NR;j++)
                {
                    R[i] = R[i] + s[i-1][j]*r[j]*eff[j]*rhoC_*MW_[i]/3600.;    //[Kmol/Kgcat/h] --> [Kg/m3cat/s]
                }
            }
        }
        delete [] s;
    }
    else if ( "ethylene-partial-oxidation" )
    {
        unsigned int NR = 2;
        OpenSMOKE::OpenSMOKEVectorDouble r(NR);
        OpenSMOKE::OpenSMOKEVectorDouble k(NR);
        OpenSMOKE::OpenSMOKEVectorDouble K(NR);
        OpenSMOKE::OpenSMOKEVectorDouble p(NC_);
        OpenSMOKE::OpenSMOKEVectorDouble *c;
        c = new OpenSMOKE::OpenSMOKEVectorDouble[NC_];
        {
            double MWmix = mixtureMolecularWeight(omega);

            for (unsigned int i=1;i<=NC_;i++)
            {
                p[i] = omega[i]*MWmix*p_/MW_[i];
                
                ChangeDimensions(NR, &c[i-1], true);
                
                for (unsigned int j=1;j<=NR;j++)
                {
                    c[i-1][j] = 0.;
                }
            }

            unsigned int iO2  = 1;
            unsigned int iE   = 0;
            unsigned int iCO2 = 3;

            //- Reaction 1
            {
                c[iO2][1] = 0.5;
                c[iE][1]  = 0.6;

                k[1]      = 6.275*1e06*std::exp(-74.9*1e03/(8.314*T));
                
                K[1]      = 1.985*1e02*std::exp(2400./T);
            }

            //- Reaction 2
            {
                c[iO2][2] = 0.5;
                c[iE][2]  = 0.5;

                k[2]      = 1.206*1e07*std::exp(-89.8*1e03/(8.314*T));
                
                K[2]      = 1.080*1e02*std::exp(1530./T);
            }

            r[1] = k[1]*std::pow(p[iO2+1],c[iO2][1])*std::pow(p[iE+1],c[iE][1])/(1. + K[1]*p[iCO2+1]);                //[mol/Kgcat/s]
            r[2] = k[2]*std::pow(p[iO2+1],c[iO2][2])*std::pow(p[iE+1],c[iE][2])/(1. + K[1]*p[iCO2+1]);                //[mol/Kgcat/s]
            
            r[1] = r[1]*1e-03;                //[mol/Kgcat/s] --> [Kmol/Kgcat/s]
            r[2] = r[2]*1e-03;                //[mol/Kgcat/s] --> [Kmol/Kgcat/s]
        }
        delete [] c;

        OpenSMOKE::OpenSMOKEVectorDouble *s;
        s = new OpenSMOKE::OpenSMOKEVectorDouble[NC_];
        {
            for (unsigned int i=0;i<NC_;i++)
            {
                ChangeDimensions(NR, &s[i], true);
                
                for (unsigned int j=1;j<=NR;j++)
                {
                    s[i][j] = 0.;
                }
            }

            // C2H4 + 0.5O2 -> C2H4O [0 + 1 -> 2]
            {
                s[0][1] = -1.0;
                s[1][1] = -0.5;
                s[2][1] =  1.0;
            }

            // C2H4 + 3O2 -> 2CO2 + 2H2O [0 + 1 -> 3 + 5]
            {
                s[0][2] = -1.;
                s[1][2] = -3.;
                s[3][2] =  2.;
                s[5][2] =  2.;
            }

            for (unsigned int i=1;i<=NC_;i++)
            {
                R[i] = 0.;
                for (unsigned int j=1;j<=NR;j++)
                {
                    R[i] = R[i] + s[i-1][j]*r[j]*rhoC_*MW_[i];    //[Kmol/Kgcat/s] --> [Kg/m3cat/s]
                }
            }
        }
        delete [] s;
    }

    return R;
}

double ReactionHeat(const OpenSMOKE::OpenSMOKEVectorDouble omega,const double T)
{
    double Q = 0.;
    if ( reactionType_ == "O-xylene-to-phthalic" )
    {
        OpenSMOKE::OpenSMOKEVectorDouble p(NC_);
        double MWmix = mixtureMolecularWeight(omega);

        for (unsigned int i=1;i<=NC_;i++)
        {
            p[i] = omega[i]*MWmix*p_/(1.e05*MW_[i]);
        }

        double k   = std::exp(19.837 - 13636./T);
        double kp  = k*p[1]*8.314*T/1e05;
        double De  = diff_*epsiC_/tauC_;
        double eff = std::min(1.,1./(Lgeo_*std::sqrt(kp*rhoC_/De)));

        double dH = -1285409.*1.e03;      //[J/Kmol]

        Q = eff*dH*k*p[1]*p[2]*rhoC_/3600.;   //[J/m3cat/s]
    }
    else if ( reactionType_ == "O-xylene-to-phthalic-complex" )
    {
        unsigned int NR = 3;
        OpenSMOKE::OpenSMOKEVectorDouble r(NR);
        OpenSMOKE::OpenSMOKEVectorDouble eff(NR);
        OpenSMOKE::OpenSMOKEVectorDouble p(NC_);
        {
            double k1    = std::exp(-13588./T + 19.837);
            double k2    = std::exp(-15803./T + 20.860);
            double k3    = std::exp(-14394./T + 18.970);
            double MWmix = mixtureMolecularWeight(omega);

            for (unsigned int i=1;i<=NC_;i++)
            {
                p[i] = omega[i]*MWmix*p_/(1e05*MW_[i]);
            }

            unsigned int iO2 = 1;
            unsigned int iXY = 2;
            unsigned int iPA = 3;
            r[1] = k1*p[iO2]*p[iXY];                //[Kmol/Kgcat/s]
            r[2] = k2*p[iPA]*p[iO2];                //[Kmol/Kgcat/s]
            r[3] = k3*p[iO2]*p[iXY];                //[Kmol/Kgcat/s]

            double kp1  = k1*p[iO2]*8.314*T/1e05;
            double kp2  = k2*p[iO2]*8.314*T/1e05;
            double kp3  = k3*p[iO2]*8.314*T/1e05;
            double De   = diff_*epsiC_/tauC_;
            eff[1]      = std::min(1.,1./(Lgeo_*std::sqrt(kp1*rhoC_/De)));
            eff[2]      = std::min(1.,1./(Lgeo_*std::sqrt(kp2*rhoC_/De)));
            eff[3]      = std::min(1.,1./(Lgeo_*std::sqrt(kp3*rhoC_/De)));
        }
        
        OpenSMOKE::OpenSMOKEVectorDouble DHr(NR);
        {
            DHr[1] = -1.285*1e06*1e03;                        //[J/Kmol]
            DHr[3] = -4.564*1e06*1e03;                        //[J/Kmol]
            DHr[2] = (DHr[3]/r[3] - DHr[1]/r[1])*r[2];        //[J/Kmol]
        }

        for (unsigned int i=1;i<=NR;i++)
        {
            Q = Q + DHr[i]*r[i]*eff[i]*rhoC_/3600.;                //[J/Kgcat/h] --> [W/m3cat]
        }
    }
    else if ( "ethylene-partial-oxidation" )
    {
        unsigned int NR = 2;
        OpenSMOKE::OpenSMOKEVectorDouble r(NR);
        OpenSMOKE::OpenSMOKEVectorDouble k(NR);
        OpenSMOKE::OpenSMOKEVectorDouble K(NR);
        OpenSMOKE::OpenSMOKEVectorDouble p(NC_);
        OpenSMOKE::OpenSMOKEVectorDouble *c;
        c = new OpenSMOKE::OpenSMOKEVectorDouble[NC_];
        {
            double MWmix = mixtureMolecularWeight(omega);

            for (unsigned int i=1;i<=NC_;i++)
            {
                p[i] = omega[i]*MWmix*p_/MW_[i];

                ChangeDimensions(NR, &c[i-1], true);

                for (unsigned int j=1;j<=NR;j++)
                {
                    c[i-1][j] = 0.;
                }
            }

            unsigned int iO2  = 1;
            unsigned int iE   = 0;
            unsigned int iCO2 = 3;

            //- Reaction 1
            {
                c[iO2][1] = 0.5;
                c[iE][1]  = 0.6;

                k[1]      = 6.275*1e06*std::exp(-74.9*1e03/(8.314*T));
                
                K[1]      = 1.985*1e02*std::exp(2400./T);
            }

            //- Reaction 2
            {
                c[iO2][2] = 0.5;
                c[iE][2]  = 0.5;

                k[2]      = 1.206*1e07*std::exp(-89.8*1e03/(8.314*T));
                
                K[2]      = 1.080*1e02*std::exp(1530./T);
            }

            r[1] = k[1]*std::pow(p[iO2+1],c[iO2][1])*std::pow(p[iE+1],c[iE][1])/(1. + K[1]*p[iCO2+1]);                //[mol/Kgcat/s]
            r[2] = k[2]*std::pow(p[iO2+1],c[iO2][2])*std::pow(p[iE+1],c[iE][2])/(1. + K[1]*p[iCO2+1]);                //[mol/Kgcat/s]
            
            r[1] = r[1]*1e-03;                //[mol/Kgcat/s] --> [Kmol/Kgcat/s]
            r[2] = r[2]*1e-03;                //[mol/Kgcat/s] --> [Kmol/Kgcat/s]
        }
        delete [] c;
        
        OpenSMOKE::OpenSMOKEVectorDouble DHr(NR);
        {
            DHr[1] = -106.7*1e03/1e-03;               //[J/Kmol]
            DHr[2] = -1323.*1e03/1e-03;               //[J/Kmol]
        }

        for (unsigned int i=1;i<=NR;i++)
        {
            Q = Q + DHr[i]*r[i]*rhoC_;                //[J/Kgcat/s] --> [W/m3cat]
        }
    }

    return Q;
}

void runAway(const double T)
{
    if ( reactionType_ == "O-xylene-to-phthalic" )
    {
        if ( T > (750. + 273.15) )
        {
            std::cout << "\n######################################" << std::endl;
            std::cout << "#                                    #" << std::endl;
            std::cout << "#           ##############           #" << std::endl;
            std::cout << "#          #              #          #" << std::endl;
            std::cout << "#         #    ##    ##    #         #" << std::endl;
            std::cout << "#         #    ##    ##    #         #" << std::endl;
            std::cout << "#         #                #         #" << std::endl;
            std::cout << "#          #              #          #" << std::endl;
            std::cout << "#            #   #   #   #           #" << std::endl;
            std::cout << "#            #   #   #   #           #" << std::endl;
            std::cout << "#            #############           #" << std::endl;
            std::cout << "#                                    #" << std::endl;
            std::cout << "#              RUN-AWAY              #" << std::endl;
            std::cout << "#                                    #" << std::endl;
            std::cout << "######################################\n" << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    else if ( reactionType_ == "O-xylene-to-phthalic-complex" )
    {
        if ( T > (750. + 273.15) )
        {
            std::cout << "\n######################################" << std::endl;
            std::cout << "#                                    #" << std::endl;
            std::cout << "#           ##############           #" << std::endl;
            std::cout << "#          #              #          #" << std::endl;
            std::cout << "#         #    ##    ##    #         #" << std::endl;
            std::cout << "#         #    ##    ##    #         #" << std::endl;
            std::cout << "#         #                #         #" << std::endl;
            std::cout << "#          #              #          #" << std::endl;
            std::cout << "#            #   #   #   #           #" << std::endl;
            std::cout << "#            #   #   #   #           #" << std::endl;
            std::cout << "#            #############           #" << std::endl;
            std::cout << "#                                    #" << std::endl;
            std::cout << "#              RUN-AWAY              #" << std::endl;
            std::cout << "#                                    #" << std::endl;
            std::cout << "######################################\n" << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    else if ( "ethylene-partial-oxidation" )
    {
        if ( T > (750. + 273.15) )
        {
            std::cout << "\n######################################" << std::endl;
            std::cout << "#                                    #" << std::endl;
            std::cout << "#           ##############           #" << std::endl;
            std::cout << "#          #              #          #" << std::endl;
            std::cout << "#         #    ##    ##    #         #" << std::endl;
            std::cout << "#         #    ##    ##    #         #" << std::endl;
            std::cout << "#         #                #         #" << std::endl;
            std::cout << "#          #              #          #" << std::endl;
            std::cout << "#            #   #   #   #           #" << std::endl;
            std::cout << "#            #   #   #   #           #" << std::endl;
            std::cout << "#            #############           #" << std::endl;
            std::cout << "#                                    #" << std::endl;
            std::cout << "#              RUN-AWAY              #" << std::endl;
            std::cout << "#                                    #" << std::endl;
            std::cout << "######################################\n" << std::endl;
            exit(EXIT_FAILURE);
        }
    }
}
