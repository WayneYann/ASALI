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

#ifndef OpenSMOKE_OdeInterfaces_H
#define OpenSMOKE_OdeInterfaces_H

#include "math/OpenSMOKEVector.h"
#include "ode/OpenSMOKE_OdeSystemObject.h"


#if ASALI_USE_SUNDIALS == 1
#include "ode/OpenSMOKE_CVODE_Sundials_Interface.h"
#include "ode/OpenSMOKE_CVODE_Sundials.h"
#endif

namespace OpenSMOKE
{
    #if ASALI_USE_SUNDIALS == 1

        class ODESystem_CVODE_Template : public OpenSMOKE::OpenSMOKE_OdeSystemObject
        {
            DEFINE_ODESOLVERINTERFACE_CVODE_Sundials(ODESystem_CVODE_Template)

            ASALI::ODESystem* ode_;
            OpenSMOKE::OpenSMOKEVectorDouble y_;
            OpenSMOKE::OpenSMOKEVectorDouble dy_;

        public:

            void SetOdeSystem(ASALI::ODESystem* ode)
            {
                ode_ = ode;
                ChangeDimensions(ode_->NumberOfEquations(), &y_, true);
                ChangeDimensions(ode_->NumberOfEquations(), &dy_, false);
            }

            int GetSystemFunctions(const double t, double* y, double* dy)
            {
                y_.CopyFrom(y);
                int flag = ode_->Equations(t, y_, dy_);
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
                int flag = ode_->Print(t, y_);
                return 0;
            }
        };
        COMPLETE_ODESOLVERINTERFACE_CVODE_Sundials(ODESystem_CVODE_Template)

    #endif
}
#endif    // OpenSMOKE_OdeInterfaces_H
