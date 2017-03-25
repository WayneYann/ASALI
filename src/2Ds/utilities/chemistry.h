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

OpenSMOKE::OpenSMOKEVectorDouble ReactionRate(const OpenSMOKE::OpenSMOKEVectorDouble omega,const double T, const double epsi)
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
    else if ( reactionType_ == "ethylene-partial-oxidation" )
    {
        unsigned int NR = 2;
        OpenSMOKE::OpenSMOKEVectorDouble r(NR);
        OpenSMOKE::OpenSMOKEVectorDouble k(NR);
        OpenSMOKE::OpenSMOKEVectorDouble K(NC_);
        OpenSMOKE::OpenSMOKEVectorDouble p(NC_);
        {
            double MWmix = mixtureMolecularWeight(omega);

            for (unsigned int i=1;i<=NC_;i++)
            {
                p[i] = omega[i]*MWmix*p_/(1.e05*MW_[i]);
            }

            unsigned int iO2  = 1 + 1;
            unsigned int iE   = 0 + 1;
            unsigned int iCO2 = 3 + 1;
            unsigned int iH2O = 5 + 1;

            //- Reaction 1
            {
                for (unsigned int i=1;i<=NC_;i++)
                {
                    K[i] = 0;
                }
                
                K[iO2]  = 2.2;
                K[iCO2] = 24.;
                K[iH2O] = 50.;
                K[iE]   = 0.01*std::exp(3000./T);
                
                k[1]    = 0.76*1e06*std::exp(-8800./T);

                double D= (1. + K[iE]*p[iE]+std::sqrt(K[iO2]*p[iO2])+K[iCO2]*p[iCO2]+K[iH2O]*p[iH2O]);

                r[1]    = k[1]*K[iE]*std::sqrt(K[iO2])*p[iE]*std::sqrt(p[iO2])/std::pow(D,2.);                            //[mol/Kgcat/s]
            }

            //- Reaction 2
            {
                for (unsigned int i=1;i<=NC_;i++)
                {
                    K[i] = 0;
                }
                
                K[iO2]  = 10.4;
                K[iCO2] = 89.0;
                K[iH2O] = 40.0;
                K[iE]   = 0.0021*std::exp(3500./T);
                
                k[2]    = 10.5*1e09*std::exp(-12800./T);

                double D= (1. + K[iE]*p[iE]+std::sqrt(K[iO2]*p[iO2])+K[iCO2]*p[iCO2]+K[iH2O]*p[iH2O]);

                r[2]    = k[2]*K[iE]*std::sqrt(K[iO2])*p[iE]*std::sqrt(p[iO2])/std::pow(D,2.);                            //[mol/Kgcat/s]
            }

            r[1] = r[1]*1e-03;                //[mol/Kgcat/s] --> [Kmol/Kgcat/s]
            r[2] = r[2]*1e-03;                //[mol/Kgcat/s] --> [Kmol/Kgcat/s]
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
                    R[i] = R[i] + s[i-1][j]*r[j]*MW_[i]*rhoC_;    //[Kmol/Kgcat/s] --> [Kg/m3cat/s]
                }
            }
        }
        delete [] s;
    }
    else if (  reactionType_ == "heat-generation" )
    {
        for (unsigned int i=1;i<=NC_;i++)
        {
            R[i] = 0.;
        }
    }
    else if ( reactionType_ == "O-xylene-to-phthalic-with-tolualdehyde" )
    {
        unsigned int NR = 5;
        OpenSMOKE::OpenSMOKEVectorDouble r(NR);
        OpenSMOKE::OpenSMOKEVectorDouble p(NC_);
        {
            double k1 = 6.917*std::exp(-61.379e03/(8.314*(T + 273.15))); // [kmol/(kg s bar)]
            double k2 = 2.497*std::exp(-46.422e03/(8.314*(T + 273.15)));
            double k3 = 0.419*std::exp(-51.170e03/(8.314*(T + 273.15)));
            double k4 = 0.000*std::exp(-54.476e03/(8.314*(T + 273.15)));
            double k5 = 57.00*std::exp(-57.907e03/(8.314*(T + 273.15)));
            double MWmix = mixtureMolecularWeight(omega);

            for (unsigned int i=1;i<=NC_;i++)
            {
                p[i] = omega[i]*MWmix*p_/(1.e05*MW_[i]);  // [bar]
            }

            unsigned int iO2 = 1;  // Oxygen
            unsigned int iXY = 2;  // o-Xylene
            unsigned int iPA = 3;  // Phthalic anhydride
            unsigned int iTA = 4;  // o-Tolualdehyde
            unsigned int iPD = 5;  // Phthalide

			double kcpo2 = 0.722e-05;

			double a = kcpo2/(kcpo2 + p[iXY]*(k1 + 6.5*k3 + 3*k4) + k2*p[iTA] + k5*p[iPD]);

            r[1] = k1*a*p[iXY];                //[Kmol/Kgcat/s]
            r[2] = k2*a*p[iTA];                //[Kmol/Kgcat/s]
            r[3] = k3*a*p[iXY];                //[Kmol/Kgcat/s]
            r[4] = k4*a*p[iXY];                //[Kmol/Kgcat/s]
            r[5] = k5*a*p[iPD];                //[Kmol/Kgcat/s]
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

            // 0/ Oxygen
            // 1/ o-Xylene              C8H10
            // 2/ Phthalic anhydride    C8H4O3
            // 3/ o-Tolualdehyde        C8H8O
            // 4/ Phthalide             C8H6O2
            // 5/ Carbon dioxide
            // 6/ Nitrogen
            // 7/ Water

            // C8H10 + O2 -> C8H8O + H2O [1 + 0 -> 3 + 7]
            {
                s[1][1] = -1.;
                s[0][1] = -1.;
                s[3][1] =  1.;
                s[7][1] =  1.;
            }

            // C8H8O + O2 -> C8H6O2 + H2O [3 + 0 -> 4 + 7]
            {
                s[3][2] = -1.;
                s[0][2] = -1.;
                s[4][2] =  1.;
                s[7][2] =  1.;
            }

            // C8H10 + (21/2)O2 -> 8CO2 + 5H2O [1 + 0 -> 5 + 7]
            {
                s[1][3] = -1.;
                s[0][3] = -21./2.;
                s[5][3] =  9.;
                s[7][3] =  5.;
            }
            
            // C8H10 + 3O2 -> C8H4O3 + 3H2O [1 + 0 -> 2 + 7]
            {
                s[1][4] = -1.;
                s[0][4] = -3.;
                s[2][4] =  1.;
                s[7][4] =  3.;
            }

            // C8H6O2 + O2 -> C8H4O3 + H2O [4 + 0 -> 2 + 7]
            {
                s[4][5] = -1.;
                s[0][5] = -1.;
                s[2][5] =  1.;
                s[7][5] =  1.;
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

double ReactionHeat(const OpenSMOKE::OpenSMOKEVectorDouble omega,const double T,const double epsi)
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
    else if ( reactionType_ == "ethylene-partial-oxidation" )
    {
        unsigned int NR = 2;
        OpenSMOKE::OpenSMOKEVectorDouble r(NR);
        OpenSMOKE::OpenSMOKEVectorDouble k(NR);
        OpenSMOKE::OpenSMOKEVectorDouble K(NC_);
        OpenSMOKE::OpenSMOKEVectorDouble p(NC_);
        {
            double MWmix = mixtureMolecularWeight(omega);

            for (unsigned int i=1;i<=NC_;i++)
            {
                p[i] = omega[i]*MWmix*p_/(1.e05*MW_[i]);
            }

            unsigned int iO2  = 1 + 1;
            unsigned int iE   = 0 + 1;
            unsigned int iCO2 = 3 + 1;
            unsigned int iH2O = 5 + 1;

            //- Reaction 1
            {
                for (unsigned int i=1;i<=NC_;i++)
                {
                    K[i] = 0;
                }
                
                K[iO2]  = 2.2;
                K[iCO2] = 24.;
                K[iH2O] = 50.;
                K[iE]   = 0.01*std::exp(3000./T);
                
                k[1]    = 0.76*1e06*std::exp(-8800./T);

                double D= (1. + K[iE]*p[iE]+std::sqrt(K[iO2]*p[iO2])+K[iCO2]*p[iCO2]+K[iH2O]*p[iH2O]);

                r[1]    = k[1]*K[iE]*std::sqrt(K[iO2])*p[iE]*std::sqrt(p[iO2])/std::pow(D,2.);                            //[mol/Kgcat/s]
            }

            //- Reaction 2
            {
                for (unsigned int i=1;i<=NC_;i++)
                {
                    K[i] = 0;
                }
                
                K[iO2]  = 10.4;
                K[iCO2] = 89.0;
                K[iH2O] = 40.0;
                K[iE]   = 0.0021*std::exp(3500./T);
                
                k[2]    = 10.5*1e09*std::exp(-12800./T);

                double D= (1. + K[iE]*p[iE]+std::sqrt(K[iO2]*p[iO2])+K[iCO2]*p[iCO2]+K[iH2O]*p[iH2O]);

                r[2]    = k[2]*K[iE]*std::sqrt(K[iO2])*p[iE]*std::sqrt(p[iO2])/std::pow(D,2.);                            //[mol/Kgcat/s]
            }

            r[1] = r[1]*1e-03;                //[mol/Kgcat/s] --> [Kmol/Kgcat/s]
            r[2] = r[2]*1e-03;                //[mol/Kgcat/s] --> [Kmol/Kgcat/s]
        }
        
        OpenSMOKE::OpenSMOKEVectorDouble DHr(NR);
        {
            DHr[1] = -25.0*4.186*1e03/1e-03;               //[J/Kmol]
            DHr[2] = -317.*4.186*1e03/1e-03;               //[J/Kmol]
        }

        for (unsigned int i=1;i<=NR;i++)
        {
            Q = Q + DHr[i]*r[i]*rhoC_;                //[J/Kgcat/s] --> [W/m3cat]
        }
    }
    else if (  reactionType_ == "heat-generation" )
    {
        Q = Qext_;                                    //[W/m3cat]
    }
    else if ( reactionType_ == "O-xylene-to-phthalic-with-tolualdehyde" )
    {
        unsigned int NR = 5;
        OpenSMOKE::OpenSMOKEVectorDouble r(NR);
        OpenSMOKE::OpenSMOKEVectorDouble p(NC_);
        {
            double k1 = 6.917*std::exp(-61.379e03/(8.314*(T + 273.15))); // [kmol/(kg s bar)]
            double k2 = 2.497*std::exp(-46.422e03/(8.314*(T + 273.15)));
            double k3 = 0.419*std::exp(-51.170e03/(8.314*(T + 273.15)));
            double k4 = 0.000*std::exp(-54.476e03/(8.314*(T + 273.15)));
            double k5 = 57.00*std::exp(-57.907e03/(8.314*(T + 273.15)));
            double MWmix = mixtureMolecularWeight(omega);

            for (unsigned int i=1;i<=NC_;i++)
            {
                p[i] = omega[i]*MWmix*p_/(1.e05*MW_[i]);  // [bar]
            }

            unsigned int iO2 = 1;  // Oxygen
            unsigned int iXY = 2;  // o-Xylene
            unsigned int iPA = 3;  // Phthalic anhydride
            unsigned int iTA = 4;  // o-Tolualdehyde
            unsigned int iPD = 5;  // Phthalide

			double kcpo2 = 0.722e-05;

			double a = kcpo2/(kcpo2 + p[iXY]*(k1 + 6.5*k3 + 3*k4) + k2*p[iTA] + k5*p[iPD]);

            r[1] = k1*a*p[iXY];                //[Kmol/Kgcat/s]
            r[2] = k2*a*p[iTA];                //[Kmol/Kgcat/s]
            r[3] = k3*a*p[iXY];                //[Kmol/Kgcat/s]
            r[4] = k4*a*p[iXY];                //[Kmol/Kgcat/s]
            r[5] = k5*a*p[iPD];                //[Kmol/Kgcat/s]
        }
        
        OpenSMOKE::OpenSMOKEVectorDouble DHr(NR);
        {
            DHr[1] = -109*4.186*1e03/1e-03;                  //[J/Kmol]
            DHr[2] = -99*4.186*1e03/1e-03;                   //[J/Kmol]
            DHr[3] = -550*4.186*1e03/1e-03;                  //[J/Kmol]
            DHr[4] = -296*4.186*1e03/1e-03;                  //[J/Kmol]
            DHr[5] = -88*4.186*1e03/1e-03;                   //[J/Kmol]
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
    else if ( reactionType_ == "ethylene-partial-oxidation" )
    {
        if ( T > (900. + 273.15) )
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
