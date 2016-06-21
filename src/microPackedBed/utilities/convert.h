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

void ConvertsToMeter(double &L, std::string &m)
{
    if (m == "m")
        L=L;
    else if(m == "micron")
        L*=0.000001;
    else if(m == "mm")
        L*=0.001;
    else if(m == "cm")
        L*=0.01;
    else if(m == "dm")
        L*=0.1;
    else if(m == "Km")
        L*=1000.;
    else
    {
        std::cout << m << " is not an accepted unit!" << std::endl;
        exit(-1);
    }
}

void ConvertsToOneOnMeter(double &L, std::string &m)
{
    if (m == "1/m")
        L=L;
    else if(m == "1/micron")
        L/=0.000001;
    else if(m == "1/mm")
        L/=0.001;
    else if(m == "1/cm")
        L/=0.01;
    else if(m == "1/dm")
        L/=0.1;
    else if(m == "1/Km")
        L/=1000.;
    else
    {
        std::cout << m << " is not an accepted unit!" << std::endl;
        exit(-1);
    }
}

void ConvertsToKg(double &mass, std::string &m)
{
    if (m == "Kg")
        mass=mass;
    else if(m == "g")
        mass*=1.e-03;
    else if(m == "mg")
        mass*=1.e-06;
    else
    {
        std::cout << m << " is not an accepted unit!" << std::endl;
        exit(-1);
    }
}

void ConvertsToSecond(double &t, std::string &m)
{
    if (m == "s")
        t=t;
    else if(m == "h")
        t*=3600.;
    else if(m == "min")
        t*=60.;
    else if(m == "d")
        t=t*24.*3600.;
    else if(m == "y")
        t=t*24.*365.*3600.;
    else if(m == "ms")
        t*=1.e-03;
    else
    {
        std::cout << m << " is not an accepted unit!" << std::endl;
        exit(-1);
    }
}


void ConvertsToPascal(double &P, std::string &m)
{
    if (m == "Pa")
        P=P;
    else if(m == "bar")
        P*=1.e05;
    else if(m == "atm")
        P*=101325.;
    else if(m == "mmHg")
        P*=101325./760.;
    else
    {
        std::cout << m << " is not an accepted unit!" << std::endl;
        exit(-1);
    }
}

void ConvertsToNm3perSecond(double &Q, std::string &m)
{
    if (m == "Nm3/s")
        Q=Q;
    else if(m == "Nm3/min")
        Q = Q/60.;
    else if(m == "Nm3/h")
        Q = Q/3600.;
    else if(m == "Nl/s")
        Q = Q*1.e-03;
    else if(m == "Nl/min")
        Q = Q*1.e-03/60.;
    else if(m == "Nl/h")
        Q = Q*1.e-03/3600.;
    else
    {
        std::cout << m << " is not an accepted unit!" << std::endl;
        exit(-1);
    }
}

void ConvertsToMeterPerSecond(double &v, std::string &m)
{
    if (m == "m/s")
        v=v;
    else if(m == "m/min")
        v = v/60.;
    else if(m == "m/h")
        v = v/3600.;
    else
    {
        std::cout << m << " is not an accepted unit!" << std::endl;
        exit(-1);
    }
}

void FromCelsiusToKelvin(double &T, std::string &m)
{
    if (m == "°C")
        T = T + 273.15;
    else
    {
        std::cout << m << " is not an accepted unit!" << std::endl;
        exit(-1);
    }
}

void FromKelvinToCelsius(double &T, std::string &m)
{
    if (m == "K")
        T = T - 273.15;
    else
    {
        std::cout << m << " is not an accepted unit!" << std::endl;
        exit(-1);
    }
}
