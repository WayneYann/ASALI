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
        for(unsigned int s=1;s<=NC_;s++)
        {
            omegaBulk_[s] = yOS_[counter++];
        }
        
        for(unsigned int s=1;s<=NC_;s++)
        {
            omegaWall_[s] = yOS_[counter++];
        }
        
        TBulk_ = yOS_[counter++];
        
        TWall_ = yOS_[counter++];
        
        t_ = t;
    }

    //- Properties estimation
    {
        rho_  = p_*mixtureMolecularWeight(omegaBulk_)/(8314.*TBulk_);
        v_    = G_/rho_;
        R_    = ReactionRate(omegaWall_,TWall_);
        Q_    = ReactionHeat(omegaWall_,TWall_);
        U_    = ExternalHeatTransfer(extHeat_);
        h_    = heatTransfer();
        kMat_ = massTransfer(rho_);
    }

    //- Equations
    {
        for(unsigned int s=1;s<=NC_;s++)
        {
            domegaBulk_[s] = -kMat_*av_*(omegaBulk_[s] - omegaWall_[s])*rho_ /(rho_*v_);
        }

        for(unsigned int s=1;s<=NC_;s++)
        {
            domegaWall_[s] = kMat_*av_*(omegaBulk_[s] - omegaWall_[s])*rho_  + R_[s]*(1. - epsi_);
        }

        dTBulk_ = (-h_*av_*(TBulk_ - TWall_)- U_*aex_*(TBulk_ - Tw_))/(rho_*cpG_*v_);
        
        dTWall_ = h_*av_*(TBulk_ - TWall_) - Q_*(1. - epsi_);
    }

    //-Assigned residuals
    {
        unsigned int counter = 1;
        for(unsigned int s=1;s<=NC_;s++)
        {
            dyOS_[counter++] = domegaBulk_[s];
        }
        for(unsigned int s=1;s<=NC_;s++)
        {
            dyOS_[counter++] = domegaWall_[s];
        }
        dyOS_[counter++] = dTBulk_;
        dyOS_[counter++] = dTWall_;
    }
}
