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
    //-Assigned variables
    {
        unsigned int counter = 1;
        for (unsigned int a=0;a<Na_;a++)
        {
            for (unsigned int r=0;r<Nr_;r++)
            {
                for(unsigned int s=1;s<=NC_;s++)
                    omegaBulk_[a][r][s] = yOS_[counter++];
                for(unsigned int s=1;s<=NC_;s++)
                    omegaWall_[a][r][s] = yOS_[counter++];
                TBulk_[a][r+1] = yOS_[counter++];
                TWall_[a][r+1] = yOS_[counter++];
            }
        }
        t_ = t;
    }

    //- Properties estimation
    {
        for (unsigned int a=0;a<Na_;a++)
        {
            for (unsigned int r=0;r<Nr_;r++)
            {
                rho_[a][r+1]  = p_*mixtureMolecularWeight(omegaBulk_[a][r])/(8314.*TBulk_[a][r+1]);
                v_[a][r+1]    = G_/rho_[a][r+1];
                h_[a][r+1]    = heatTransfer(zH_[a+1]);
                kMat_[a][r+1] = massTransfer(zH_[a+1],rho_[a][r+1]);
                if ( zH_[a + 1] < Linert_ )
                {
                    R_[a][r]      = 0.;
                    Q_[a][r+1]    = 0.;
                }
                else
                {
                    R_[a][r]      = ReactionRate(omegaWall_[a][r],TWall_[a][r+1],epsiHC_);
                    Q_[a][r+1]    = ReactionHeat(omegaWall_[a][r],TWall_[a][r+1],epsiHC_);
                }
            }
        }

        kr_ = RadialConduction(reactorType_);
        ka_ = AxialConduction(reactorType_);
        hw_ = ExternalHeatTransfer(reactorType_);
    }
    
    //- Derivative
    {
        for (unsigned int r=0;r<Nr_;r++)
        {
            OpenSMOKE::OpenSMOKEVectorDouble variable(Na_);
            OpenSMOKE::OpenSMOKEVectorDouble dvariable(Na_);
            
            //- 1st - temperature
            {
                for (unsigned int a=0;a<Na_;a++)
                {
                    variable[a+1] = TBulk_[a][r+1];
                }
                
                dvariable = FirstOrderDerivite(variable,zH_);

                for (unsigned int a=0;a<Na_;a++)
                {
                    d1stTBulk_[a][r+1] = dvariable[a+1];
                }
            }
            
            //- 1st - specie
            for (unsigned int s=1;s<=NC_;s++)
            {
                for (unsigned int a=0;a<Na_;a++)
                {
                    variable[a+1] = omegaBulk_[a][r][s];
                }
                
                dvariable = FirstOrderDerivite(variable,zH_);

                for (unsigned int a=0;a<Na_;a++)
                {
                    d1stomegaBulk_[a][r][s] = dvariable[a+1];
                }
            }
            
            //- 2nd - axial
            {
                for (unsigned int a=0;a<Na_;a++)
                {
                    variable[a+1] = TWall_[a][r+1];
                }
                
                dvariable = SecondOrderDerivite(variable,zH_,ka_,"axial");

                for (unsigned int a=0;a<Na_;a++)
                {
                    d2ndTWallaxial_[a][r+1] = dvariable[a+1];
                }
            }
        }
        
        for (unsigned int a=0;a<Na_;a++)
        {
            OpenSMOKE::OpenSMOKEVectorDouble variable(Nr_);
            OpenSMOKE::OpenSMOKEVectorDouble dvariable(Nr_);

            //- 2nd - radial
            {
                for (unsigned int r=0;r<Nr_;r++)
                {
                    variable[r+1] = TWall_[a][r+1];
                }
                
                dvariable = SecondOrderDerivite(variable,rH_,kr_,"radial");

                for (unsigned int r=0;r<Nr_;r++)
                {
                    d2ndTWallradial_[a][r+1] = dvariable[r+1];
                }
            }
        }
    }


    //- Equations
    {
        for (unsigned int a=0;a<Na_;a++)
        {
            for (unsigned int r=0;r<Nr_;r++)
            {
                runAway(TWall_[a][r+1]);
                if ( a == 0 )
                {
                    for (unsigned int s=1;s<=NC_;s++)
                    {
                        domegaBulk_[a][r][s] = omegaIN_[s] - omegaBulk_[a][r][s];
                    }

                    for (unsigned int s=1;s<=NC_;s++)
                    {
                        domegaWall_[a][r][s] = rho_[a][r+1]*avH_*kMat_[a][r+1]*(omegaBulk_[a][r][s]- omegaWall_[a][r][s]) + R_[a][r][s]*epsiHC_;
                    }

                    dTBulk_[a][r+1] = Tin_ - TBulk_[a][r+1];
                    dTWall_[a][r+1] = TWall_[a+1][r+1] - TWall_[a][r+1];
                }
                else if ( a == (Na_ - 1))
                {
                    for (unsigned int s=1;s<=NC_;s++)
                    {
                        domegaBulk_[a][r][s] = -v_[a][r+1]*d1stomegaBulk_[a][r][s]
                                               -avH_*kMat_[a][r+1]*(omegaBulk_[a][r][s]- omegaWall_[a][r][s]);
                    }

                    for (unsigned int s=1;s<=NC_;s++)
                    {
                        domegaWall_[a][r][s] = rho_[a][r+1]*avH_*kMat_[a][r+1]*(omegaBulk_[a][r][s]- omegaWall_[a][r][s]) + R_[a][r][s]*epsiHC_;
                    }

                    if ( energy_ == true )
                    {
                        dTBulk_[a][r+1] = -v_[a][r+1]*d1stTBulk_[a][r+1]
                                          -avH_*h_[a][r+1]*(TBulk_[a][r+1] - TWall_[a][r+1])/(rho_[a][r+1]*cpG_*epsiH_);
                    }
                    else
                    {
                        dTBulk_[a][r+1] = 0.;
                    }
                                      
                    dTWall_[a][r+1] = TWall_[a][r+1] - TWall_[a-1][r+1];
                }
                else
                {
                    for (unsigned int s=1;s<=NC_;s++)
                    {
                        domegaBulk_[a][r][s] = -v_[a][r+1]*d1stomegaBulk_[a][r][s]
                                               -avH_*kMat_[a][r+1]*(omegaBulk_[a][r][s]- omegaWall_[a][r][s]);
                    }

                    for (unsigned int s=1;s<=NC_;s++)
                    {
                        domegaWall_[a][r][s] = rho_[a][r+1]*avH_*kMat_[a][r+1]*(omegaBulk_[a][r][s]- omegaWall_[a][r][s]) + R_[a][r][s]*epsiHC_;
                    }
                    
                    if ( r == 0 )
                    {
                        if ( energy_ == true )
                        {
                            dTBulk_[a][r+1] = -v_[a][r+1]*d1stTBulk_[a][r+1]
                                              -avH_*h_[a][r+1]*(TBulk_[a][r+1] - TWall_[a][r+1])/(rho_[a][r+1]*cpG_*epsiH_);
                        }
                        else
                        {
                            dTBulk_[a][r+1] = 0.;
                        }

                        dTWall_[a][r+1] = TWall_[a][r+1+1] - TWall_[a][r+1];
                    }
                    else if ( r == (Nr_ - 1))
                    {
                        if ( energy_ == true )
                        {
                            dTBulk_[a][r+1] = -v_[a][r+1]*d1stTBulk_[a][r+1]
                                              -avH_*h_[a][r+1]*(TBulk_[a][r+1] - TWall_[a][r+1])/(rho_[a][r+1]*cpG_*epsiH_);
                                              
                            dTWall_[a][r+1] = (kr_/hw_)*(TWall_[a][r+1] - TWall_[a][r+1-1])/(rH_[r+1] - rH_[r+1-1]) + (TWall_[a][r+1] - Tw_);
                        }
                        else
                        {
                            dTBulk_[a][r+1] = 0.;
                            dTWall_[a][r+1] = (TWall_[a][r+1] - TWall_[a][r+1-1]);
                        }
                    }
                    else
                    {
                        if ( energy_ == true )
                        {
                            dTBulk_[a][r+1] = -v_[a][r+1]*d1stTBulk_[a][r+1]
                                              -avH_*h_[a][r+1]*(TBulk_[a][r+1] - TWall_[a][r+1])/(rho_[a][r+1]*cpG_*epsiH_);
                                              
                            dTWall_[a][r+1] =  d2ndTWallaxial_[a][r+1]/(rhoS_*cpS_*(1. - epsiH_))
                                              +d2ndTWallradial_[a][r+1]/(rhoS_*cpS_*(1. - epsiH_))
                                              +avH_*h_[a][r+1]*(TBulk_[a][r+1] - TWall_[a][r+1])/(rhoS_*cpS_*(1. - epsiH_))
                                              -Q_[a][r+1]*epsiHC_/(rhoS_*cpS_*(1. - epsiH_));
                        }
                        else
                        {
                            dTBulk_[a][r+1] = 0.;
                            dTWall_[a][r+1] = 0.;
                        }
                    }
                }
            }
        }
    }

    //-Assigned residuals
    {
        unsigned int counter = 1;
        for (unsigned int a=0;a<Na_;a++)
        {
            for (unsigned int r=0;r<Nr_;r++)
            {
                for(unsigned int s=1;s<=NC_;s++)
                    dyOS_[counter++] = domegaBulk_[a][r][s];
                for(unsigned int s=1;s<=NC_;s++)
                    dyOS_[counter++] = domegaWall_[a][r][s];
                dyOS_[counter++] = dTBulk_[a][r+1];
                dyOS_[counter++] = dTWall_[a][r+1];
            }
        }
    }
}
