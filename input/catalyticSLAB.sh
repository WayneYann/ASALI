touch catalyticSLAB.txt

echo '/*##############################################################################################' >> catalyticSLAB.txt
echo '#                                                                                              #' >> catalyticSLAB.txt
echo '#     #############       #############       #############       ####                ####     #' >> catalyticSLAB.txt
echo '#    #             #     #             #     #             #     #    #              #    #    #' >> catalyticSLAB.txt
echo '#    #    #####    #     #    #########      #    #####    #     #    #              #    #    #' >> catalyticSLAB.txt
echo '#    #    #   #    #     #    #              #    #   #    #     #    #              #    #    #' >> catalyticSLAB.txt
echo '#    #    #####    #     #    #              #    #####    #     #    #              #    #    #' >> catalyticSLAB.txt
echo '#    #             #     #    #########      #             #     #    #              #    #    #' >> catalyticSLAB.txt
echo '#    #             #     #             #     #             #     #    #              #    #    #' >> catalyticSLAB.txt
echo '#    #    #####    #      #########    #     #    #####    #     #    #              #    #    #' >> catalyticSLAB.txt
echo '#    #    #   #    #              #    #     #    #   #    #     #    #              #    #    #' >> catalyticSLAB.txt
echo '#    #    #   #    #      #########    #     #    #   #    #     #    #########      #    #    #' >> catalyticSLAB.txt
echo '#    #    #   #    #     #             #     #    #   #    #     #             #     #    #    #' >> catalyticSLAB.txt
echo '#     ####     ####       #############       ####     ####       #############       ####     #' >> catalyticSLAB.txt
echo '#                                                                                              #' >> catalyticSLAB.txt
echo '#   Department of Energy                                                                       #' >> catalyticSLAB.txt
echo '#   Politecnico di Milano                                                                      #' >> catalyticSLAB.txt
echo '#   Author: Stefano Rebughini <stefano.rebughini@polimi.it>                                    #' >> catalyticSLAB.txt
echo '#                                                                                              #' >> catalyticSLAB.txt
echo '################################################################################################' >> catalyticSLAB.txt
echo ' ' >> catalyticSLAB.txt
echo 'Kinetics path           kinetics/kinetics' >> catalyticSLAB.txt
echo ' ' >> catalyticSLAB.txt
echo 'Pressure                1 bar' >> catalyticSLAB.txt
echo ' ' >> catalyticSLAB.txt
echo 'Temperature             500 °C' >> catalyticSLAB.txt
echo ' ' >> catalyticSLAB.txt
echo 'Reactor' >> catalyticSLAB.txt
echo '{' >> catalyticSLAB.txt
echo '  length                0.5  cm' >> catalyticSLAB.txt
echo '  diameter              24   mm' >> catalyticSLAB.txt
echo '}' >> catalyticSLAB.txt
echo ' ' >> catalyticSLAB.txt
echo 'Catalyst' >> catalyticSLAB.txt
echo '{' >> catalyticSLAB.txt
echo '  alfa                  100  1/m' >> catalyticSLAB.txt
echo '}' >> catalyticSLAB.txt
echo ' ' >> catalyticSLAB.txt
echo 'Solid' >> catalyticSLAB.txt
echo '{' >> catalyticSLAB.txt
echo '  void fraction         0.86' >> catalyticSLAB.txt
echo '  tortuosity            8.00' >> catalyticSLAB.txt
echo '  pore radius           1e-06 m' >> catalyticSLAB.txt
echo '} ' >> catalyticSLAB.txt
echo ' ' >> catalyticSLAB.txt
echo 'Mole fractions' >> catalyticSLAB.txt
echo '{' >> catalyticSLAB.txt
echo '  CO                    0.1' >> catalyticSLAB.txt
echo '  O2                    0.4' >> catalyticSLAB.txt
echo '  N2                    0.5' >> catalyticSLAB.txt
echo '  inert                 N2' >> catalyticSLAB.txt
echo '}' >> catalyticSLAB.txt
echo ' ' >> catalyticSLAB.txt
echo 'Solver options' >> catalyticSLAB.txt
echo '{' >> catalyticSLAB.txt
echo '  Constraints          true' >> catalyticSLAB.txt
echo ' ' >> catalyticSLAB.txt
echo '  Absolute tollerance  1e-12' >> catalyticSLAB.txt
echo ' ' >> catalyticSLAB.txt
echo '  Relative tollerance  1e-09' >> catalyticSLAB.txt
echo '' >> catalyticSLAB.txt
echo '  Reactions          on' >> catalyticSLAB.txt
echo '' >> catalyticSLAB.txt
echo '  Number of points   800' >> catalyticSLAB.txt
echo '' >> catalyticSLAB.txt
echo '  Diffusion model    DustyGas' >> catalyticSLAB.txt
echo '' >> catalyticSLAB.txt
echo '  Integration time            1000 s' >> catalyticSLAB.txt
echo ' ' >> catalyticSLAB.txt
echo '  Results                     mole' >> catalyticSLAB.txt
echo '}' >> catalyticSLAB.txt
echo ' ' >> catalyticSLAB.txt
echo 'Numerical solvers' >> catalyticSLAB.txt
echo '{' >> catalyticSLAB.txt
echo '  ODE       BzzMath' >> catalyticSLAB.txt
echo '  DAE       BzzMath' >> catalyticSLAB.txt
echo '}' >> catalyticSLAB.txt
