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

#ifndef OpenSMOKE_DaeInterfaces_H
#define OpenSMOKE_DaeInterfaces_H

#include "math/OpenSMOKEVector.h"
#include "dae/OpenSMOKE_DaeSystemObject.h"


#if ASALI_USE_SUNDIALS == 1
#include "dae/OpenSMOKE_IDA_Sundials_Interface.h"
#include "dae/OpenSMOKE_IDA_Sundials.h"
#endif

namespace OpenSMOKE
{
    #if ASALI_USE_SUNDIALS == 1

        class DAESystem_IDA_Template : public OpenSMOKE::OpenSMOKE_DaeSystemObject
        {
            DEFINE_DAESOLVERINTERFACE_IDA_Sundials(DAESystem_IDA_Template)

            ASALI::BVPSystem* bvp_;
            OpenSMOKE::OpenSMOKEVectorDouble y_;
            OpenSMOKE::OpenSMOKEVectorDouble dy_;

        public:

            void SetDaeSystem(ASALI::BVPSystem* bvp)
            {
                bvp_ = bvp;
                ChangeDimensions(bvp_->NumberOfEquations(), &y_, true);
                ChangeDimensions(bvp_->NumberOfEquations(), &dy_, false);
            }

            int GetSystemFunctions(const double t, double* y, double* dy)
            {
                y_.CopyFrom(y);
                int flag = bvp_->Equations(t, y_, dy_);
                dy_.CopyTo(dy);
                return(flag);
            }

            int GetAnalyticalJacobian(const double t, double* y, double* J)
            {
                return(0);
            }

            int GetWriteFunction(const double t, double *y)
            {
                y_.CopyFrom(y);
                int flag = bvp_->Print(t, y_);
                return 0;
            }
        };
        COMPLETE_DAESOLVERINTERFACE_IDA_Sundials(DAESystem_IDA_Template)

    #endif
}

#endif    // OpenSMOKE_DaeInterfaces_H
