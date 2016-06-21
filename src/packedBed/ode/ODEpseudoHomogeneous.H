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
        TBulk_ = yOS_[counter++];
        t_ = t;
    }

    //- Properties estimation
    {
        rho_  = p_*mixtureMolecularWeight(omegaBulk_)/(8314.*TBulk_);
        v_    = G_/rho_;
        R_    = ReactionRate(omegaBulk_,TBulk_);
        Q_    = ReactionHeat(omegaBulk_,TBulk_);
        U_    = ExternalHeatTransfer(extHeat_);
    }

    //- Equations
    {
        for(unsigned int s=1;s<=NC_;s++)
        {
            domegaBulk_[s] = R_[s]*(1. - epsi_)/(rho_*v_);
        }
        
        dTBulk_ = (-Q_*(1. - epsi_) - U_*aex_*(TBulk_ - Tw_))/(rho_*cpG_*v_);
    }

    //-Assigned residuals
    {
        unsigned int counter = 1;
        for(unsigned int s=1;s<=NC_;s++)
        {
            dyOS_[counter++] = domegaBulk_[s];
        }
        dyOS_[counter++] = dTBulk_;
    }
}
