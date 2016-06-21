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

double mixtureMolecularWeight(const OpenSMOKE::OpenSMOKEVectorDouble omega)
{
    double MWmix_ = 0.;
    for (unsigned int i=1;i<=NC_;i++)
    {
        MWmix_ = MWmix_ + omega[i]/MW_[i];
    }
    return 1./MWmix_;
}

double massTransfer(const double rho)
{
    double Re = G_*Dp_/(mu_*epsi_);
    double Sc = mu_*rho/diff_;
    double jM = 0.357/(epsi_*std::pow(Re,0.359));
    double Sh = jM*std::pow(Sc,(1./3.))*Re;
    double kMat = Sh*diff_/Dp_;
    return kMat;
}

double heatTransfer()
{
    double Re = G_*Dp_/(mu_*epsi_);
    double jH = 0.357/(epsi_*std::pow(Re,0.359));
    double Nu = jH*std::pow(Pr_,(1./3.))*Re;
    double h  = Nu*kG_/Dp_;
    return h;
}

double ExternalHeatTransfer(const bool check)
{
    double U = 0.;
    if ( check == true )
    {
        double Re  = G_*Dp_/mu_;
        double krf = G_*cpG_*Dp_/(8.6*(1.+19.4*std::pow((Dp_/Dt_),2.)));
        double A   = kC_/kG_;
        double p1  = ((1./3.)*std::pow((1.-1./A),2.))/(std::log(A-0.577*(A-1))-0.423*(1.-1./A)) - 2./(3.*A);
        double p2  = ((0.072)*std::pow((1.-1./A),2.))/(std::log(A-0.925*(A-1))-0.075*(1.-1./A)) - 2./(3.*A);
        double p   = 0.;
        if ( epsi_ <= 0.26 )
        {
            p = p2;
        }
        else if ( epsi_ >= 0.476 )
        {
            p = p1;
        }
        else
        {
            p = p2 + (p1-p2)*(epsi_-0.26)/0.216;
        }

        double krs   = kG_*(1.-epsi_)*(1.+2.66*std::sqrt(Dp_/Dt_))/(2./(3.*A) + p);
        double Bi    = 5.73*std::sqrt(Dt_/Dp_)*std::pow(Re,-0.262);
        double h     = 0.393*G_*cpG_*std::pow(Re,-0.384)/(epsi_*std::pow(Pr_,(2./3.)));        
        double Ns    = 1.5*(1-epsi_)*h*Dt_*Dt_/(krs*Dp_);
        double kr    = krf + krs*(1+4/Bi)/(1+8/Ns);
        double hw    = Bi*2*kr/Dt_;
               U     = 1./(1./hw + (Dt_/(6.*kr))*(Bi+3.)/(Bi+4.));
    }
    else
    {
        U = 0.;
    }

    return U;
}
