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
if ( input.getSolver() == "OpenSMOKE")
{
    unsigned int NE = eq->NumberOfEquations();

    Eigen::VectorXd y0(NE);
    Eigen::VectorXd yf(NE);

    y0(0) = input.getFeedPressure();
    y0(1) = input.getFeedTemperature();

    {
        // Create the solver
        typedef OdeSMOKE::KernelDense<OpenSMOKE::ODESystem_OpenSMOKE> denseOde;
        typedef OdeSMOKE::MethodGear<denseOde> methodGear;
        OdeSMOKE::MultiValueSolver<methodGear> ode_solver;
        ode_solver.SetReactor(eq);

        // Set initial conditions
        ode_solver.SetInitialConditions(0., y0);

        // Set linear algebra options
        ode_solver.SetFullPivoting(true);

        // Set relative and absolute tolerances
        ode_solver.SetAbsoluteTolerances(1e-12);
        ode_solver.SetRelativeTolerances(1e-07);


        // Solve the system
        OdeSMOKE::OdeStatus status = ode_solver.Solve(input.getReactorLength());
        ode_solver.Solution(yf);
        ChangeDimensions(NE, &yF, true);
        yF.CopyFrom(yf.data());
    }
}
else if ( input.getSolver() == "BzzMath")
{
#if ASALI_USE_BZZ == 1
    unsigned int NE = eq->NumberOfEquations();

    BzzVector y0(NE);
    BzzVector yf(NE);

    y0[1] = input.getFeedPressure();
    y0[2] = input.getFeedTemperature();

    {
        // Create the solver
        OpenSMOKE::ODESystem_BzzOde ode(*eq);
        BzzOdeStiffObject ode_solver(y0, 0., &ode);

        // Set relative and absolute tolerances
        ode_solver.SetTolAbs(1e-12);
        ode_solver.SetTolRel(1e-07);

        // Solve the system
        yf = ode_solver(input.getReactorLength());
        ChangeDimensions(NE, &yF, true);
        yF.CopyFrom(yf.GetHandle());
    }
#endif
}
