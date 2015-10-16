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
        for (unsigned int i=0;i<NP_;i++)
        {
            for (unsigned int j=1;j<=NC_;j++)
                omegaBulk_[i][j] = yOS_[counter++];
            for (unsigned int j=1;j<=SURF_NC_;j++)
                teta_[i][j] = yOS_[counter++];
        }
        t_ = t;
    }

    unsigned int counter = 1;
    for (unsigned int i=0;i<NP_;i++)
    {
        //- Equations
        if(i==0)
        {
            //-INLET conditions
            for (unsigned int j=1;j<=NC_;j++)
            {
                if ( thermodynamicsMapXML->NamesOfSpecies()[j-1] != inert_ )
                {
                    dyOS_[counter++] = x0bulk_[j] - omegaBulk_[i][j];// ALGEBRAIC
                }
                else
                {
                    dyOS_[counter++] = 1e06*(1. - omegaBulk_[i].SumElements());
                }
            }
            for (unsigned int j=1;j<=SURF_NC_;j++)
            {
                if ( thermodynamicsSurfaceMap_.NamesOfSpecies()[j-1+thermodynamicsSurfaceMap_.number_of_gas_species()] != "Rh(s)" )
                {
                    dyOS_[counter++] = 0.; //DIFFERENTIAL
                }
                else
                {
                    dyOS_[counter++] = 1e06*(1. - teta_[i].SumElements());
                }
            }
        }
        else if (i==(NP_-1))
        {
            //-OUTLET conditions
            for (unsigned int j=1;j<=NC_;j++)
            {
                if ( thermodynamicsMapXML->NamesOfSpecies()[j-1] != inert_ )
                {
                    dyOS_[counter++] = 0.; //DIFFERENTIAL
                }
                else
                {
                    dyOS_[counter++] = 1e06*(1. - omegaBulk_[i].SumElements());
                }
            }
            for (unsigned int j=1;j<=SURF_NC_;j++)
            {
                if ( thermodynamicsSurfaceMap_.NamesOfSpecies()[j-1+thermodynamicsSurfaceMap_.number_of_gas_species()] != "Rh(s)" )
                {
                    dyOS_[counter++] = 0.; //DIFFERENTIAL
                }
                else
                {
                    dyOS_[counter++] = 1e06*(1. - teta_[i].SumElements());
                }
            }
        }
        else
        {
            for (unsigned int j=1;j<=NC_;j++)
            {
                if ( thermodynamicsMapXML->NamesOfSpecies()[j-1] != inert_ )
                {
                    dyOS_[counter++] = 0.; //DIFFERENTIAL
                }
                else
                {
                    dyOS_[counter++] = 1e06*(1. - omegaBulk_[i].SumElements());// ALGEBRAIC
                }
            }
            for (unsigned int j=1;j<=SURF_NC_;j++)
            {
                if ( thermodynamicsSurfaceMap_.NamesOfSpecies()[j-1+thermodynamicsSurfaceMap_.number_of_gas_species()] != "Rh(s)" )
                {
                    dyOS_[counter++] = 0.; //DIFFERENTIAL
                }
                else
                {
                    dyOS_[counter++] = 1e06*(1. - teta_[i].SumElements());
                }
            }
        }
    }
}
