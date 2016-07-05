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
    eq->resize(N);
    eq->setGrid(z);
    eq->setDiscretizationScheme(discretizationScheme);
}

if ( input.getSolver() == "OpenSMOKE")
{
    unsigned int NE = eq->NumberOfEquations();

    Eigen::VectorXd                      yMin(NE);
    Eigen::VectorXd                      yMax(NE);
    Eigen::VectorXd                      y0(NE);
    Eigen::VectorXd                      yf(NE);
    std::vector<OpenSMOKE::EquationType> algebraic(NE);

    {
        unsigned int counter = 0;
        for (unsigned int k=1;k<=N;k++)
        {
            for (unsigned int j=1;j<=NC;j++)
            {
                yMax(counter) = 1.;
                yMin(counter) = 0.;
                counter++;
            }
            for (unsigned int j=1;j<=NC;j++)
            {
                yMax(counter) = 1.;
                yMin(counter) = 0.;
                counter++;
            }
            for (unsigned int j=1;j<=SURF_NC;j++)
            {
                yMax(counter) = 1.;
                yMin(counter) = 0.;
                counter++;
            }
            
            yMax(counter) = 5000.;
            yMin(counter) = 0.;
            counter++;
            
            yMax(counter) = 5000.;
            yMin(counter) = 0.;
            counter++;
        }
    }
    
    {
        ASALI::ODEstart();
        x0.CopyTo(y0.data());
        // Create the solver
        typedef OdeSMOKE::KernelBand<OpenSMOKE::ODESystem_OpenSMOKE> bandOde;
        typedef OdeSMOKE::MethodGear<bandOde> methodGear;
        OdeSMOKE::MultiValueSolver<methodGear> ode_solver;
        ode_solver.SetReactor(eq);

        // Set initial conditions
        ode_solver.SetInitialConditions(0., y0);
        ode_solver.SetBandSizes(eq->BandDimension(),eq->BandDimension());

        // Set relative and absolute tolerances
        ode_solver.SetAbsoluteTolerances(input.getAbsTol());
        ode_solver.SetRelativeTolerances(input.getRelTol());

        // Set minimum and maximum values
        if ( input.getConstraints() == true)
        {
            ode_solver.SetMinimumValues(yMin);
            ode_solver.SetMaximumValues(yMax);
        }

        // Solve the system
        OdeSMOKE::OdeStatus status = ode_solver.Solve(1e10);
        ode_solver.Solution(y0);
        ASALI::ODEend();
    }

    //- DAE
    {
        ASALI::DAEstart();
        typedef DaeSMOKE::KernelBand<OpenSMOKE::BVPSystem_OpenSMOKE> bandBvp;
        typedef DaeSMOKE::MethodGear<bandBvp> methodGear;
        DaeSMOKE::MultiValueSolver<methodGear> bvp_solver;
        bvp_solver.SetReactor(eq);

        // Set initial conditions
        eq->AlgebraicEquations(algebraic);
        bvp_solver.SetInitialConditions(0., y0, algebraic);
        bvp_solver.SetBandSizes(eq->BandDimension(),eq->BandDimension());

        // Set relative and absolute tolerances
        bvp_solver.SetAbsoluteTolerances(input.getAbsTol());
        bvp_solver.SetRelativeTolerances(input.getRelTol());

        // Set minimum and maximum values
        if ( input.getConstraints() == true)
        {
            bvp_solver.SetMinimumValues(yMin);
            bvp_solver.SetMaximumValues(yMax);
        }

        DaeSMOKE::DaeStatus status = bvp_solver.Solve(input.getIntegrationTime());

        bvp_solver.Solution(yf);
        ChangeDimensions(NE, &yF, true);
        yF.CopyFrom(yf.data());
        ASALI::DAEend();
    }
}
else if ( input.getSolver() == "BzzMath")
{
#if ASALI_USE_BZZ == 1
    unsigned int NE = eq->NumberOfEquations();

    BzzVector                            yMin(NE);
    BzzVector                            yMax(NE);
    BzzVector                            y0(NE);
    BzzVector                            yf(NE);
    BzzVectorInt                         integer(NE);
    std::vector<OpenSMOKE::EquationType> algebraic(NE);

    {
        unsigned counter = 1;
        for (unsigned int k=1;k<=N;k++)
        {
            for (unsigned int j=1;j<=NC;j++)
            {
                yMax[counter] = 1.;
                yMin[counter] = 0.;
                counter++;
            }
            for (unsigned int j=1;j<=NC;j++)
            {
                yMax[counter] = 1.;
                yMin[counter] = 0.;
                counter++;
            }
            for (unsigned int j=1;j<=SURF_NC;j++)
            {
                yMax[counter] = 1.;
                yMin[counter] = 0.;
                counter++;
            }

            yMax[counter] = 5000.;
            yMin[counter] = 0.;
            counter++;
            
            yMax[counter] = 5000.;
            yMin[counter] = 0.;
            counter++;
        }
    }
    
    {
        ASALI::ODEstart();
        x0.CopyTo(y0.GetHandle());

        // Create the solver
        OpenSMOKE::ODESystem_BzzOde ode(*eq);
        BzzOdeSparseStiffObject ode_solver(y0, 0., &ode,eq->BlockDimension());

        // Set relative and absolute tolerances
        ode_solver.SetTolAbs(input.getAbsTol());
        ode_solver.SetTolRel(input.getRelTol());

        // Set minimum and maximum values
        if ( input.getConstraints() == true)
        {
            ode_solver.SetMinimumConstraints(yMin);
            ode_solver.SetMaximumConstraints(yMax);
        }

        // Solve the system
        y0 = ode_solver(1e10);
        ASALI::ODEend();
    }

    //- DAE
    {
        ASALI::DAEstart();
        
        eq->AlgebraicEquations(algebraic);
        for (unsigned int i = 1; i <= NE; i++)
        {
            integer[i] = (algebraic[i-1] == OpenSMOKE::EQUATION_TYPE_DIFFERENTIAL) ? 1 : 0;
        }

        // Create the solver
        OpenSMOKE::BVPSystem_BzzBvp bvp(*eq);
        BzzDaeSparseObject bvp_solver(y0, 0.,integer, &bvp,eq->BlockDimension());
        
        // Set relative and absolute tolerances
        if ( input.getConstraints() == true)
        {
            bvp_solver.SetMinimumConstraints(yMin);
            bvp_solver.SetMaximumConstraints(yMax);
        }

        // Set minimum and maximum values
        bvp_solver.SetTolAbs(input.getAbsTol());
        bvp_solver.SetTolRel(input.getRelTol());

        yf = bvp_solver(input.getIntegrationTime());
        ChangeDimensions(NE, &yF, true);
        yF.CopyFrom(yf.GetHandle());
        ASALI::DAEend();
    }
#endif
}
