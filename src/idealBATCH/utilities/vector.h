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

void FromCppToOS(const std::vector<double> &v, OpenSMOKE::OpenSMOKEVectorDouble &b)
{
    ChangeDimensions(v.size(), &b,true);
    double* ptb = b.GetHandle();
    for(int i=0;i<v.size();i++)
        *ptb++ = v[i];
}

void FromOSToCpp(const OpenSMOKE::OpenSMOKEVectorDouble &b, std::vector<double> &v)
{
    v.resize(b.Size());
    const double* ptb = b.GetHandle();
    for(int i=0;i<v.size();i++)
        v[i] = *ptb++;
}

#if ASALI_USE_BZZ == 1
void FromBzzToOS(const BzzVector &bzz, OpenSMOKE::OpenSMOKEVectorDouble &os)
{
    memcpy(os.GetHandle(), bzz.GetHandle(), bzz.Size()*sizeof(double));
}

void FromOSToBzz(const OpenSMOKE::OpenSMOKEVectorDouble&os, BzzVector &bzz)
{
    memcpy(bzz.GetHandle(), os.GetHandle(), bzz.Size()*sizeof(double));
}
#endif
