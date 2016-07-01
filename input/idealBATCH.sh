touch idealBATCH.txt

echo '/*##############################################################################################' >> idealBATCH.txt
echo '#                                                                                              #' >> idealBATCH.txt
echo '#     #############       #############       #############       ####                ####     #' >> idealBATCH.txt
echo '#    #             #     #             #     #             #     #    #              #    #    #' >> idealBATCH.txt
echo '#    #    #####    #     #    #########      #    #####    #     #    #              #    #    #' >> idealBATCH.txt
echo '#    #    #   #    #     #    #              #    #   #    #     #    #              #    #    #' >> idealBATCH.txt
echo '#    #    #####    #     #    #              #    #####    #     #    #              #    #    #' >> idealBATCH.txt
echo '#    #             #     #    #########      #             #     #    #              #    #    #' >> idealBATCH.txt
echo '#    #             #     #             #     #             #     #    #              #    #    #' >> idealBATCH.txt
echo '#    #    #####    #      #########    #     #    #####    #     #    #              #    #    #' >> idealBATCH.txt
echo '#    #    #   #    #              #    #     #    #   #    #     #    #              #    #    #' >> idealBATCH.txt
echo '#    #    #   #    #      #########    #     #    #   #    #     #    #########      #    #    #' >> idealBATCH.txt
echo '#    #    #   #    #     #             #     #    #   #    #     #             #     #    #    #' >> idealBATCH.txt
echo '#     ####     ####       #############       ####     ####       #############       ####     #' >> idealBATCH.txt
echo '#                                                                                              #' >> idealBATCH.txt
echo '#   Department of Energy                                                                       #' >> idealBATCH.txt
echo '#   Politecnico di Milano                                                                      #' >> idealBATCH.txt
echo '#   Author: Stefano Rebughini <stefano.rebughini@polimi.it>                                    #' >> idealBATCH.txt
echo '#                                                                                              #' >> idealBATCH.txt
echo '################################################################################################' >> idealBATCH.txt
echo ' ' >> idealBATCH.txt
echo 'Kinetics path CH4_UBI/kinetics' >> idealBATCH.txt
echo ' ' >> idealBATCH.txt
echo 'Pressure    1 atm' >> idealBATCH.txt
echo ' ' >> idealBATCH.txt
echo 'Temperature 1097.15 K' >> idealBATCH.txt
echo ' ' >> idealBATCH.txt
echo 'Volume      5.4e-08 m3' >> idealBATCH.txt
echo ' ' >> idealBATCH.txt
echo 'Catalyst' >> idealBATCH.txt
echo '{' >> idealBATCH.txt
echo '	alfa	14630	1/m' >> idealBATCH.txt
echo '	area	6.007495e-05 m2' >> idealBATCH.txt
echo '}' >> idealBATCH.txt
echo ' ' >> idealBATCH.txt
echo 'Mole fractions' >> idealBATCH.txt
echo '{' >> idealBATCH.txt
echo '  CH4                   0.27' >> idealBATCH.txt
echo '  O2                    0.1593' >> idealBATCH.txt
echo '  N2                    0.5707' >> idealBATCH.txt
echo '  inert                 N2' >> idealBATCH.txt
echo '}' >> idealBATCH.txt
echo ' ' >> idealBATCH.txt
echo 'Coverage' >> idealBATCH.txt
echo '{' >> idealBATCH.txt
echo '	Rh(s)  1.0' >> idealBATCH.txt
echo '}' >> idealBATCH.txt
echo ' ' >> idealBATCH.txt
echo 'Solver options' >> idealBATCH.txt
echo '{' >> idealBATCH.txt
echo '  Constraints          false' >> idealBATCH.txt
echo ' ' >> idealBATCH.txt
echo '  Constant             pressure' >> idealBATCH.txt
echo ' ' >> idealBATCH.txt
echo '  Absolute tollerance  1e-12' >> idealBATCH.txt
echo '  Relative tollerance  1e-07' >> idealBATCH.txt
echo ' ' >> idealBATCH.txt
echo '  Energy equation    on' >> idealBATCH.txt
echo ' ' >> idealBATCH.txt
echo '  Reactions' >> idealBATCH.txt
echo '  (' >> idealBATCH.txt
echo '    Homogeneous      off' >> idealBATCH.txt
echo '    Heterogeneous    on' >> idealBATCH.txt
echo '  )' >> idealBATCH.txt
echo ' ' >> idealBATCH.txt
echo '  Integration time            100 s' >> idealBATCH.txt
echo ' ' >> idealBATCH.txt
echo '  Numerical solvers           OpenSMOKE' >> idealBATCH.txt
echo '}' >> idealBATCH.txt
