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
    //- Assigned variable
    {
        p_ = yOS_[1];
        T_ = yOS_[2];
        t_ = t;
    }

    //-Setting conditions for GAS phase
    {
        thermodynamicsMap_.SetTemperature(T_);
        thermodynamicsMap_.SetPressure(p_);
        transportMap_.SetPressure(p_);
        transportMap_.SetTemperature(T_);
    }

    //- Evaluting properties
    {
        for (unsigned int k=0;k<=NC_;k++)
            x_[k] = 0;

        x_[thermodynamicsMap_.IndexOfSpecies(name_)] = 1.;

        thermodynamicsMap_.MolecularWeight_From_MoleFractions(MW_,x_);

        //- Density [Kg/m3]
        rho_ = p_*MW_/(PhysicalConstants::R_J_kmol*T_);

        //- Velocity [m/s]
        v_   = G_/rho_;


        if ( law_ == "OpenSMOKE")
        {
            //- Specific heat [J/Kg/K]
            thermodynamicsMap_.cpMolar_Mixture_From_MoleFractions(cp_, x_);
            cp_ = cp_/MW_;

            //- Thermal conductivity [W/m/K]
            transportMap_.ThermalConductivity(k_,x_);

            //- Viscosity [kg/m/s]
            transportMap_.DynamicViscosity(mu_, x_);
        }
        else if ( law_ == "sutherland")
        {
            //- Specific heat [J/Kg/K]
            cp_ = specificHeatConstantPressure(T_);
            cv_ = specificHeatConstantVolume(T_);
            cp_ = cp_/MW_;
            cv_ = cv_/MW_;

            //- Viscosity [kg/m/s]
            mu_ = viscosity(T_);

            //- Thermal conductivity [W/m/K]
            k_ = termalConducitivity(mu_,cv_,MW_);
        }

        //- Specific area [1/m]
        av();

        //- Prandt number
        Pr_ = mu_*cp_/k_;

        //- Fanning friction factor
        FrictionFactor(ReynoldsForFrictionFactor());

        //- Nusselt number
        Nusselt(ReynoldsForHeatTransfer());

        //- Heat transfer coefficient [W/m2/K]
        h();
    }

    if ( type_ == "Monolith" )
    {
        dyOS_[1] = -rho_*v_*v_*fm_/Dt_;
    }
    else if ( type_ == "PackedBed" )
    {
        dyOS_[1] = -rho_*v_*v_*fm_/Dp_;
    }

    dyOS_[2] = -h_*av_*(T_ - Tw_)/(rho_*v_*cp_/epsi_);
}
