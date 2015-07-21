touch idealPFR.txt

echo '/*##############################################################################################' >> idealPFR.txt
echo '#                                                                                              #' >> idealPFR.txt
echo '#     #############       #############       #############       ####                ####     #' >> idealPFR.txt
echo '#    #             #     #             #     #             #     #    #              #    #    #' >> idealPFR.txt
echo '#    #    #####    #     #    #########      #    #####    #     #    #              #    #    #' >> idealPFR.txt
echo '#    #    #   #    #     #    #              #    #   #    #     #    #              #    #    #' >> idealPFR.txt
echo '#    #    #####    #     #    #              #    #####    #     #    #              #    #    #' >> idealPFR.txt
echo '#    #             #     #    #########      #             #     #    #              #    #    #' >> idealPFR.txt
echo '#    #             #     #             #     #             #     #    #              #    #    #' >> idealPFR.txt
echo '#    #    #####    #      #########    #     #    #####    #     #    #              #    #    #' >> idealPFR.txt
echo '#    #    #   #    #              #    #     #    #   #    #     #    #              #    #    #' >> idealPFR.txt
echo '#    #    #   #    #      #########    #     #    #   #    #     #    #########      #    #    #' >> idealPFR.txt
echo '#    #    #   #    #     #             #     #    #   #    #     #             #     #    #    #' >> idealPFR.txt
echo '#     ####     ####       #############       ####     ####       #############       ####     #' >> idealPFR.txt
echo '#                                                                                              #' >> idealPFR.txt
echo '#   Department of Energy                                                                       #' >> idealPFR.txt
echo '#   Politecnico di Milano                                                                      #' >> idealPFR.txt
echo '#   Author: Stefano Rebughini <stefano.rebughini@polimi.it>                                    #' >> idealPFR.txt
echo '#                                                                                              #' >> idealPFR.txt
echo '################################################################################################' >> idealPFR.txt
echo ' ' >> idealPFR.txt
echo 'Kinetics path CH4_UBI/kinetics' >> idealPFR.txt
echo ' ' >> idealPFR.txt
echo 'Pressure    1 atm' >> idealPFR.txt
echo ' ' >> idealPFR.txt
echo 'Temperature 1097.15 K' >> idealPFR.txt
echo ' ' >> idealPFR.txt
echo 'Volume      5.4e-08 m3' >> idealPFR.txt
echo ' ' >> idealPFR.txt
echo 'Reactor' >> idealPFR.txt
echo '{' >> idealPFR.txt
echo ' type	        annular' >> idealPFR.txt
echo ' length	        1.0 cm' >> idealPFR.txt
echo ' inner diameter	0.1 cm' >> idealPFR.txt
echo ' outer diameter  0.15 cm' >> idealPFR.txt
echo '}' >> idealPFR.txt
echo 'Catalyst' >> idealPFR.txt
echo '{' >> idealPFR.txt
echo '  mass                  500 mg' >> idealPFR.txt
echo '  dispersion            0.2' >> idealPFR.txt
echo '  Rh fraction           0.02' >> idealPFR.txt
echo '  deactivation factor   1.0' >> idealPFR.txt
echo '}' >> idealPFR.txt
echo ' ' >> idealPFR.txt
echo 'Mole fractions' >> idealPFR.txt
echo '{' >> idealPFR.txt
echo '  CH4                   0.27' >> idealPFR.txt
echo '  O2                    0.1593' >> idealPFR.txt
echo '  N2                    0.5707' >> idealPFR.txt
echo '  inert                 N2' >> idealPFR.txt
echo '}' >> idealPFR.txt
echo ' ' >> idealPFR.txt

echo 'Solver options' >> idealPFR.txt
echo '{' >> idealPFR.txt
echo '  Constraints          false' >> idealPFR.txt
echo ' ' >> idealPFR.txt
echo '  Absolute tollerance' >> idealPFR.txt
echo '  (' >> idealPFR.txt
echo '    specie           1e-15' >> idealPFR.txt
echo '  )' >> idealPFR.txt
echo ' ' >> idealPFR.txt
echo '  Relative tollerance' >> idealPFR.txt
echo '  (' >> idealPFR.txt
echo '    specie           1e-08' >> idealPFR.txt
echo '  )' >> idealPFR.txt
echo ' ' >> idealPFR.txt
echo '  Energy equation    on' >> idealPFR.txt
echo ' ' >> idealPFR.txt
echo '  Reactions' >> idealPFR.txt
echo '  (' >> idealPFR.txt
echo '    Homogeneous      off' >> idealPFR.txt
echo '    Heterogeneous    on' >> idealPFR.txt
echo '  )' >> idealPFR.txt
echo ' ' >> idealPFR.txt
echo '}' >> idealPFR.txt
echo ' ' >> idealPFR.txt
echo 'Numerical solvers' >> idealPFR.txt
echo '{' >> idealPFR.txt
echo '     ODE BzzMath' >> idealPFR.txt
echo '     DAE BzzMath' >> idealPFR.txt
echo '}' >> idealPFR.txt
