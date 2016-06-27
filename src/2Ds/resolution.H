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
    if (input.getModelName()[j] == "honeyComb")
    {
        eq->resize(input.getModelName()[j]);
        ChangeDimensions(eq->NumberOfEquations(), &x0, true);
        unsigned int counter = 1;
        for (unsigned int a=0;a<input.axialPoints();a++)
        {
            for (unsigned int r=0;r<input.radialPoints();r++)
            {
                for(unsigned int s=1;s<=input.numberOfSpecies();s++)
                    x0[counter++] = input.getMassFraction()[s-1];
                for(unsigned int s=1;s<=input.numberOfSpecies();s++)
                    x0[counter++] = input.getMassFraction()[s-1];
                x0[counter++] = input.getCoolantTemperature();
                x0[counter++] = input.getCoolantTemperature();
            }
        }
    }
    else if (input.getModelName()[j] == "packedBed")
    {
        eq->resize(input.getModelName()[j]);
        ChangeDimensions(eq->NumberOfEquations(), &x0, true);
        unsigned int counter = 1;
        for (unsigned int a=0;a<input.axialPoints();a++)
        {
            for (unsigned int r=0;r<input.radialPoints();r++)
            {
                for(unsigned int s=1;s<=input.numberOfSpecies();s++)
                    x0[counter++] = input.getMassFraction()[s-1];
                x0[counter++] = input.getCoolantTemperature();
            }
        }
    }
    else if (input.getModelName()[j] == "microBed")
    {
        eq->resize(input.getModelName()[j]);
        ChangeDimensions(eq->NumberOfEquations(), &x0, true);
        unsigned int counter = 1;
        for (unsigned int a=0;a<input.axialPoints();a++)
        {
            for (unsigned int r=0;r<input.radialPoints();r++)
            {
                for(unsigned int s=1;s<=input.numberOfSpecies();s++)
                    x0[counter++] = input.getMassFraction()[s-1];
                x0[counter++] = input.getCoolantTemperature();
                x0[counter++] = input.getCoolantTemperature();
            }
        }
    }
}

if ( input.getSolver() == "OpenSMOKE")
{
    unsigned int NE = eq->NumberOfEquations();

    Eigen::VectorXd                      yMin(NE);
    Eigen::VectorXd                      yMax(NE);
    Eigen::VectorXd                      y0(NE);
    Eigen::VectorXd                      yf(NE);
    std::vector<OpenSMOKE::EquationType> algebraic(NE);

    if (input.getModelName()[j] == "honeyComb")
    {
        unsigned int counter = 0;
        for (unsigned int a=0;a<input.axialPoints();a++)
        {
            for (unsigned int r=0;r<input.radialPoints();r++)
            {
                for(unsigned int s=1;s<=input.numberOfSpecies();s++)
                {
                    yMax(counter) = 1.;
                    yMin(counter) = 0.;
                    counter++;
                }
                for(unsigned int s=1;s<=input.numberOfSpecies();s++)
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
    }
    else if (input.getModelName()[j] == "packedBed")
    {
        unsigned int counter = 0;
        for (unsigned int a=0;a<input.axialPoints();a++)
        {
            for (unsigned int r=0;r<input.radialPoints();r++)
            {
                for(unsigned int s=1;s<=input.numberOfSpecies();s++)
                {
                    yMax(counter) = 1.;
                    yMin(counter) = 0.;
                    counter++;
                }

                yMax(counter) = 5000.;
                yMin(counter) = 0.;
                counter++;
            }
        }
    }
    else if (input.getModelName()[j] == "microBed")
    {
        unsigned int counter = 0;
        for (unsigned int a=0;a<input.axialPoints();a++)
        {
            for (unsigned int r=0;r<input.radialPoints();r++)
            {
                for(unsigned int s=1;s<=input.numberOfSpecies();s++)
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
    }

    {
        eq->start(input.getModelName()[j],"ode");
        x0.CopyTo(y0.data());
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

        // Set minimum and maximum values
        ode_solver.SetMinimumValues(yMin);
        ode_solver.SetMaximumValues(yMax);

        // Solve the system
        OdeSMOKE::OdeStatus status = ode_solver.Solve(1e10);
        ode_solver.Solution(y0);
        eq->end(input.getModelName()[j],"ode");
    }

    //- DAE
    {
        eq->start(input.getModelName()[j],"bvp");
        typedef DaeSMOKE::KernelDense<OpenSMOKE::BVPSystem_OpenSMOKE> denseBvp;
        typedef DaeSMOKE::MethodGear<denseBvp> methodGear;
        DaeSMOKE::MultiValueSolver<methodGear> bvp_solver;
        bvp_solver.SetReactor(eq);

        eq->AlgebraicEquations(algebraic);
        bvp_solver.SetInitialConditions(0., y0, algebraic);

        bvp_solver.SetFullPivoting(true);


        bvp_solver.SetAbsoluteTolerances(1.e-12);
        bvp_solver.SetRelativeTolerances(1.e-07);


        bvp_solver.SetMinimumValues(yMin);
        bvp_solver.SetMaximumValues(yMax);


        DaeSMOKE::DaeStatus status = bvp_solver.Solve(1e06);

        bvp_solver.Solution(yf);
        ChangeDimensions(NE, &yF, true);
        yF.CopyFrom(yf.data());
        eq->end(input.getModelName()[j],"bvp");
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

    if (input.getModelName()[j] == "honeyComb")
    {
        unsigned int counter = 1;
        for (unsigned int a=0;a<input.axialPoints();a++)
        {
            for (unsigned int r=0;r<input.radialPoints();r++)
            {
                for(unsigned int s=1;s<=input.numberOfSpecies();s++)
                {
                    yMax[counter] = 1.;
                    yMin[counter] = 0.;
                    counter++;
                }
                for(unsigned int s=1;s<=input.numberOfSpecies();s++)
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
    }
    else if (input.getModelName()[j] == "packedBed")
    {
        unsigned int counter = 1;
        for (unsigned int a=0;a<input.axialPoints();a++)
        {
            for (unsigned int r=0;r<input.radialPoints();r++)
            {
                for(unsigned int s=1;s<=input.numberOfSpecies();s++)
                {
                    yMax[counter] = 1.;
                    yMin[counter] = 0.;
                    counter++;
                }

                yMax[counter] = 5000.;
                yMin[counter] = 0.;
                counter++;
            }
        }
    }
    else if (input.getModelName()[j] == "microBed")
    {
        unsigned int counter = 1;
        for (unsigned int a=0;a<input.axialPoints();a++)
        {
            for (unsigned int r=0;r<input.radialPoints();r++)
            {
                for(unsigned int s=1;s<=input.numberOfSpecies();s++)
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
    }

    {
        eq->start(input.getModelName()[j],"ode");
        x0.CopyTo(y0.GetHandle());

        // Create the solver
        OpenSMOKE::ODESystem_BzzOde ode(*eq);
        BzzOdeStiffObject ode_solver(y0, 0., &ode);

        // Set relative and absolute tolerances
        ode_solver.SetTolAbs(1e-12);
        ode_solver.SetTolRel(1e-07);

        // Set minimum and maximum values
        ode_solver.SetMinimumConstraints(yMin);
        ode_solver.SetMaximumConstraints(yMax);

        y0 = ode_solver(1e10);
        eq->end(input.getModelName()[j],"ode");
    }

    //- DAE
    {
        eq->start(input.getModelName()[j],"bvp");
        
        eq->AlgebraicEquations(algebraic);
        for (unsigned int i = 1; i <= NE; i++)
        {
            integer[i] = (algebraic[i-1] == OpenSMOKE::EQUATION_TYPE_DIFFERENTIAL) ? 1 : 0;
        }

        // Create the solver
        OpenSMOKE::BVPSystem_BzzBvp bvp(*eq);
        BzzDaeObject bvp_solver(y0, 0.,integer, &bvp);
        
        // Set relative and absolute tolerances
        bvp_solver.SetTolAbs(1e-12);
        bvp_solver.SetTolRel(1e-07);

        // Set minimum and maximum values
        bvp_solver.SetMinimumConstraints(yMin);
        bvp_solver.SetMaximumConstraints(yMax);

        yf = bvp_solver(1e06);
        ChangeDimensions(NE, &yF, true);
        yF.CopyFrom(yf.GetHandle());
        eq->end(input.getModelName()[j],"bvp");
    }
#endif
}

{
    eq->cleanMemory();
}
