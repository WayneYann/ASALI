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
        s[4] =  0.;
        s[5] =  0.;
        s[6] =  3.;
        s[7] =  0.;
        s[8] =  0.;

        OpenSMOKE::OpenSMOKEVectorDouble p(NC_);

        double MWmix = mixtureMolecularWeight(omega);
        for (unsigned int i=1;i<=NC_;i++)
        {
            p[i] = omega[i]*MWmix*p_/(1.e05*MW_[i]);
        }

        double k = std::exp(19.837 - 13636./T);
        
        for (unsigned int i=1;i<=NC_;i++)
        {
            R[i] = s[i]*k*p[1]*p[2]*rhoC_*MW_[i]/3600.; //[Kg/m3cat/s]
        }
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

        double k  = std::exp(19.837 - 13636./T);

        double dH = -1285409.*1.e03;      //[J/Kmol]

        Q = dH*k*p[1]*p[2]*rhoC_/3600.;   //[J/m3cat/s]
    }
    return Q;
}

void runAway(const double T)
{
    if ( reactionType_ == "O-xylene-to-phthalic" )
    {
        if ( T > (420. + 273.14) )
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
