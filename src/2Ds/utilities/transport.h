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

double AxialConduction(const std::string type)
{
    double condAx = 0.;
    if ( type == "honeyComb" )
    {
        if ( typeH_ == "washcoated" )
        {
            condAx = kS_*(epsiHS_ + (kC_/kS_)*epsiHC_);
        }
        else if  ( typeH_ == "extruded" )
        {
            condAx = kC_*epsiHS_;
        }
    }
    else if ( type == "packedBed" )
    {
        double Re     = G_*DpP_/mu_;
        double B      = 1.25*std::pow(((1. - epsiP_)/epsiP_),(10./9.));
        double A      = kG_/kC_;
        double teta   = (1. - A)*B*std::log(1./(B*A))/std::pow((1. - B*A),2.) - (B - 1.)/(1. - B*A) - (B + 1.)/2.;
        double kr0    = kG_*((1. - std::sqrt(1. - epsiP_)) + 2.*std::sqrt(1.-epsiP_)*teta/(1. - B*A));
        double Peinf  = 8.*(2. - std::pow((1. - 2.*DpP_/DtP_),2.));
        double Pe     = 1./(1./Peinf + (2./3.)*epsiP_/(Re*Pr_));
        double kr     = (1./Pe)*Re*Pr_*kG_;
               condAx = kr0 + kr;
    }
    else if ( type == "microBed" )
    {
        condAx = kC_*(1. - epsiMh_);
    }

    return condAx;
}

double RadialConduction(const std::string type)
{
    double condR = 0.;
    if ( type == "honeyComb" )
    {
        if ( typeH_ == "washcoated" )
        {
            condR = kS_/(1.-std::sqrt(epsiH_ + epsiHC_) + (std::sqrt(epsiH_+epsiHC_)-std::sqrt(epsiH_))/(1.-std::sqrt(epsiH_+epsiHC_)+(kC_/kS_)*std::sqrt(epsiH_+epsiHC_))+std::sqrt(epsiH_)/(1.-std::sqrt(epsiH_+epsiHC_)+(kC_/kS_)*(std::sqrt(epsiH_+epsiHC_)-std::sqrt(epsiH_)) + (kG_/kS_)*std::sqrt(epsiH_))); 
        }
        else if  ( typeH_ == "extruded" )
        {
            condR = kC_/(1.-std::sqrt(epsiH_) + (std::sqrt(epsiH_)-std::sqrt(epsiH_))/(1.-std::sqrt(epsiH_))+std::sqrt(epsiH_)/(1.-std::sqrt(epsiH_)+(kG_/kC_)*std::sqrt(epsiH_)));
        }
    }
    else if ( type == "packedBed" )
    {
        double Re    = G_*DpP_/mu_;
        double B     = 1.25*std::pow(((1. - epsiP_)/epsiP_),(10./9.));
        double A     = kG_/kC_;
        double teta  = (1. - A)*B*std::log(1./(B*A))/std::pow((1. - B*A),2.) - (B - 1.)/(1. - B*A) - (B + 1.)/2.;
        double kr0   = kG_*((1. - std::sqrt(1. - epsiP_)) + 2.*std::sqrt(1.-epsiP_)*teta/(1. - B*A));
        double Peinf = 8.*(2. - std::pow((1. - 2.*DpP_/DtP_),2.));
        double Pe    = 1./(1./Peinf + (2./3.)*epsiP_/(Re*Pr_));
        double kr    = (1./Pe)*Re*Pr_*kG_;
               condR = kr0 + kr;
    }
    else if ( type == "microBed" )
    {
        condR = kS_/(1.-std::sqrt(epsiMh_) + std::sqrt(epsiMh_)/(1.-std::sqrt(epsiMh_)+(kG_/kS_)*std::sqrt(epsiMh_)));
    }

    return condR;
}

double massTransfer(const double z,const double rho)
{
    double Re   = G_*DcH_/mu_;
    double Sc   = mu_*rho/diff_;
    double zS   = std::max(1e-06,z/(DcH_*Re*Sc));
    double Sh   = 3.692 + 6.877*(std::pow((1000*zS),-0.488)*std::exp(-57.2*zS));
    double kMat = Sh*diff_/DcH_;
    return kMat;
}

double heatTransfer(const double z)
{
    double Re   = G_*DcH_/mu_;
    double zS   = std::max(1e-06,z/(DcH_*Re*Pr_));
    double Nu   = 2.997 + 8.827*(std::pow((1000*zS),-0.545)*std::exp(-48.2*zS));
    double h    = Nu*kG_/DcH_;
    return h;
}

double ExternalHeatTransfer(const std::string type)
{
    double U = 0.;
    if ( type == "honeyComb" )
    {
        U = 500; //[W/m/K]
    }
    else if ( type == "packedBed" )
    {
        double Re    = G_*DpP_/mu_;
        double B     = 1.25*std::pow(((1. - epsiP_)/epsiP_),(10./9.));
        double A     = kG_/kC_;
        double teta  = (1. - A)*B*std::log(1./(B*A))/std::pow((1. - B*A),2.) - (B - 1.)/(1. - B*A) - (B + 1.)/2.;
        double kr0   = kG_*((1. - std::sqrt(1. - epsiP_)) + 2.*std::sqrt(1.-epsiP_)*teta/(1. - B*A));
        double Nw0   = (1.3 + 5.*DpP_/DtP_)*(kr0/kG_);
        double Nw    = Nw0 + 0.19*std::pow(Re,0.75)*std::pow(Pr_,(1./3.));
               U     = Nw*kG_/DpP_;
    }
    else if ( type == "microBed" )
    {
        double Re  = G_*DpM_/mu_;
        double krf = G_*cpG_*DpM_/(8.6*(1.+19.4*std::pow((DpM_/DcM_),2.)));
        double A   = kC_/kG_;
        double p1  = ((1./3.)*std::pow((1.-1./A),2.))/(std::log(A-0.577*(A-1))-0.423*(1.-1./A)) - 2./(3.*A);
        double p2  = ((0.072)*std::pow((1.-1./A),2.))/(std::log(A-0.925*(A-1))-0.075*(1.-1./A)) - 2./(3.*A);
        double p   = 0.;
        if ( epsiP_ <= 0.26 )
        {
            p = p2;
        }
        else if ( epsiP_ >= 0.476 )
        {
            p = p1;
        }
        else
        {
            p = p2 + (p1-p2)*(epsiP_-0.26)/0.216;
        }

        double krs   = kG_*(1.-epsiP_)*(1.+2.66*std::sqrt(DpM_/DcM_))/(2./(3.*A) + p);
        double Bi    = 5.73*std::sqrt(DcM_/DpM_)*std::pow(Re,-0.262);
        double h     = 0.393*G_*cpG_*std::pow(Re,-0.384)/(epsiP_*std::pow(Pr_,(2./3.)));        
        double Ns    = 1.5*(1-epsiP_)*h*DcM_*DcM_/(krs*DpM_);
        double kr    = krf + krs*(1+4/Bi)/(1+8/Ns);
        double hw    = Bi*2*kr/DcM_;
               U     = 1./(1./hw + (DcM_/(6.*kr))*(Bi+3.)/(Bi+4.));
    }

    return U;
}
