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
#include "math/external-ode-solvers/OpenSMOKE_OdeSystemObject.h"

#if OPENSMOKE_USE_BZZMATH == 1
#include "BzzMath.hpp"
#endif

#include "math/native-ode-solvers/MultiValueSolver"

namespace OpenSMOKE
{
    class ODESystem_OpenSMOKE
    {
    public:

            ODESystem_OpenSMOKE() {};

            void SetReactor(ASALI::equationSystem* reactor)
            {
                reactor_ = reactor;
            }

        protected:

            unsigned int ne_;

            void MemoryAllocation()
            {
                OpenSMOKE::ChangeDimensions(ne_, &y_, true);
                OpenSMOKE::ChangeDimensions(ne_, &dy_, false);
            }

            virtual void Equations(const Eigen::VectorXd &Y, const double t, Eigen::VectorXd &DY)
            {
                y_.CopyFrom(Y.data());
                reactor_->OdeEquations(t, y_, dy_);
                dy_.CopyTo(DY.data());
            }

            virtual void Jacobian(const Eigen::VectorXd &Y, const double t, Eigen::MatrixXd &J) {};

            void Print(const double t, const Eigen::VectorXd &Y)
            {
                y_.CopyFrom(Y.data());
                reactor_->OdePrint(t, y_);
            }

        private:

            ASALI::equationSystem* reactor_;
            OpenSMOKE::OpenSMOKEVectorDouble  y_;
            OpenSMOKE::OpenSMOKEVectorDouble dy_;
    };

    #if OPENSMOKE_USE_BZZMATH == 1
    class ODESystem_BzzOde : public BzzOdeSystemObject 
    {
        public:
            ODESystem_BzzOde(ASALI::equationSystem& reactor) :
            reactor_(reactor)
            {
                ChangeDimensions(reactor_.NumberOfEquations(), &y_, true);
                ChangeDimensions(reactor_.NumberOfEquations(), &dy_, false);
            }
        
            virtual void GetSystemFunctions(BzzVector &Y, double t, BzzVector &DY)
            {
                y_.CopyFrom(Y.GetHandle());
                reactor_.OdeEquations(t, y_, dy_);
                dy_.CopyTo(DY.GetHandle());
                MyPrint(Y,t);
            }

            void MyPrint(BzzVector &Y, double t)
            {
                y_.CopyFrom(Y.GetHandle());
                reactor_.OdePrint(t, y_);
            }

            virtual void ObjectBzzPrint(void) {};

        private:

            ASALI::equationSystem& reactor_;
            OpenSMOKE::OpenSMOKEVectorDouble  y_;
            OpenSMOKE::OpenSMOKEVectorDouble dy_;
    };
    #endif
}

#endif    // OpenSMOKE_OdeInterfaces_H
