touch microPACKEDbed.txt

echo '/*##############################################################################################' >> microPACKEDbed.txt
echo '#                                                                                              #' >> microPACKEDbed.txt
echo '#     #############       #############       #############       ####                ####     #' >> microPACKEDbed.txt
echo '#    #             #     #             #     #             #     #    #              #    #    #' >> microPACKEDbed.txt
echo '#    #    #####    #     #    #########      #    #####    #     #    #              #    #    #' >> microPACKEDbed.txt
echo '#    #    #   #    #     #    #              #    #   #    #     #    #              #    #    #' >> microPACKEDbed.txt
echo '#    #    #####    #     #    #              #    #####    #     #    #              #    #    #' >> microPACKEDbed.txt
echo '#    #             #     #    #########      #             #     #    #              #    #    #' >> microPACKEDbed.txt
echo '#    #             #     #             #     #             #     #    #              #    #    #' >> microPACKEDbed.txt
echo '#    #    #####    #      #########    #     #    #####    #     #    #              #    #    #' >> microPACKEDbed.txt
echo '#    #    #   #    #              #    #     #    #   #    #     #    #              #    #    #' >> microPACKEDbed.txt
echo '#    #    #   #    #      #########    #     #    #   #    #     #    #########      #    #    #' >> microPACKEDbed.txt
echo '#    #    #   #    #     #             #     #    #   #    #     #             #     #    #    #' >> microPACKEDbed.txt
echo '#     ####     ####       #############       ####     ####       #############       ####     #' >> microPACKEDbed.txt
echo '#                                                                                              #' >> microPACKEDbed.txt
echo '#   Department of Energy                                                                       #' >> microPACKEDbed.txt
echo '#   Politecnico di Milano                                                                      #' >> microPACKEDbed.txt
echo '#   Author: Stefano Rebughini <stefano.rebughini@polimi.it>                                    #' >> microPACKEDbed.txt
echo '#                                                                                              #' >> microPACKEDbed.txt
echo '################################################################################################' >> microPACKEDbed.txt
echo ' ' >> microPACKEDbed.txt
echo 'Kinetics path CH4_UBI/kinetics' >> microPACKEDbed.txt
echo ' ' >> microPACKEDbed.txt
echo 'Volumetric flow rate    10 Nl/min' >> microPACKEDbed.txt
echo ' ' >> microPACKEDbed.txt
echo 'Pressure                1 bar' >> microPACKEDbed.txt
echo ' ' >> microPACKEDbed.txt
echo 'Temperature' >> microPACKEDbed.txt
echo '{' >> microPACKEDbed.txt
echo '  solid                 500 °C' >> microPACKEDbed.txt
echo '  gas                    40 °C' >> microPACKEDbed.txt
echo '  wall                   10 °C' >> microPACKEDbed.txt
echo '}' >> microPACKEDbed.txt
echo ' ' >> microPACKEDbed.txt
echo 'Reactor' >> microPACKEDbed.txt
echo '{' >> microPACKEDbed.txt
echo '  inert length                 0.5  cm' >> microPACKEDbed.txt
echo '  catalytic length             26   mm' >> microPACKEDbed.txt
echo '  honeycomb diameter           24   mm' >> microPACKEDbed.txt
echo '  channel diameter             1.0  mm' >> microPACKEDbed.txt
echo '  particle diameter            0.3  mm' >> microPACKEDbed.txt
echo '  void fraction                0.75' >> microPACKEDbed.txt
echo '  gas-to-particle correlation  Petrovic' >> microPACKEDbed.txt
echo '}' >> microPACKEDbed.txt
echo ' ' >> microPACKEDbed.txt
echo 'Catalyst' >> microPACKEDbed.txt
echo '{' >> microPACKEDbed.txt
echo '  alfa                    40000 1/m' >> microPACKEDbed.txt
echo '}' >> microPACKEDbed.txt
echo ' ' >> microPACKEDbed.txt
echo 'Honeycomb' >> microPACKEDbed.txt
echo '{' >> microPACKEDbed.txt
echo '  density               2300 Kg/m3' >> microPACKEDbed.txt
echo '  conductivity          2.5  W/m/K' >> microPACKEDbed.txt
echo '  specific heat         925 J/Kg/K' >> microPACKEDbed.txt
echo '} ' >> microPACKEDbed.txt
echo ' ' >> microPACKEDbed.txt
echo 'PackedBed' >> microPACKEDbed.txt
echo '{' >> microPACKEDbed.txt
echo '  density               2300 Kg/m3' >> microPACKEDbed.txt
echo '  conductivity          2.5  W/m/K' >> microPACKEDbed.txt
echo '  specific heat         925 J/Kg/K' >> microPACKEDbed.txt
echo '} ' >> microPACKEDbed.txt
echo ' ' >> microPACKEDbed.txt
echo 'Mole fractions' >> microPACKEDbed.txt
echo '{' >> microPACKEDbed.txt
echo '  CH4                   0.27' >> microPACKEDbed.txt
echo '  O2                    0.1593' >> microPACKEDbed.txt
echo '  N2                    0.5707' >> microPACKEDbed.txt
echo '  inert                 N2' >> microPACKEDbed.txt
echo '}' >> microPACKEDbed.txt
echo ' ' >> microPACKEDbed.txt
echo 'Solver options' >> microPACKEDbed.txt
echo '{' >> microPACKEDbed.txt
echo '  Constraints        false' >> microPACKEDbed.txt
echo ' ' >> microPACKEDbed.txt
echo '  Absolute tollerance' >> microPACKEDbed.txt
echo '  (' >> microPACKEDbed.txt
echo '    specie           1e-15' >> microPACKEDbed.txt
echo '    temperature      1e-08' >> microPACKEDbed.txt
echo '  )' >> microPACKEDbed.txt
echo ' ' >> microPACKEDbed.txt
echo '  Relative tollerance' >> microPACKEDbed.txt
echo '  (' >> microPACKEDbed.txt
echo '    specie           1e-08' >> microPACKEDbed.txt
echo '    temperature      1e-08' >> microPACKEDbed.txt
echo '  )' >> microPACKEDbed.txt
echo ' ' >> microPACKEDbed.txt
echo '  Energy equation    on' >> microPACKEDbed.txt
echo ' ' >> microPACKEDbed.txt
echo '  Reactions' >> microPACKEDbed.txt
echo '  (' >> microPACKEDbed.txt
echo '    Homogeneous      off' >> microPACKEDbed.txt
echo '    Heterogeneous    on' >> microPACKEDbed.txt
echo '  )' >> microPACKEDbed.txt
echo ' ' >> microPACKEDbed.txt
echo '  Grid' >> microPACKEDbed.txt
echo '  (' >> microPACKEDbed.txt
echo '    Growing                   on' >> microPACKEDbed.txt
echo '    Minimum number of points  10' >> microPACKEDbed.txt
echo '    Maximum number of points  800' >> microPACKEDbed.txt
echo '    Add number of points      300' >> microPACKEDbed.txt
echo '    Resolution type           new' >> microPACKEDbed.txt
echo '    Discretization scheme     BDS' >> microPACKEDbed.txt
echo '  )' >> microPACKEDbed.txt
echo ' ' >> microPACKEDbed.txt
echo '  Diffusion in gas on' >> microPACKEDbed.txt
echo ' ' >> microPACKEDbed.txt
echo '  External heat exchange on' >> microPACKEDbed.txt
echo ' ' >> microPACKEDbed.txt
echo '  Accepted errors' >> microPACKEDbed.txt
echo '  (' >> microPACKEDbed.txt
echo '    grid               1e-08' >> microPACKEDbed.txt
echo '    energy             1e-05' >> microPACKEDbed.txt
echo '  )' >> microPACKEDbed.txt
echo ' ' >> microPACKEDbed.txt
echo '  Integration time            100 s' >> microPACKEDbed.txt
echo ' ' >> microPACKEDbed.txt
echo '  Results                     mole' >> microPACKEDbed.txt
echo '}' >> microPACKEDbed.txt
echo ' ' >> microPACKEDbed.txt
echo 'Numerical solvers' >> microPACKEDbed.txt
echo '{' >> microPACKEDbed.txt
echo '  ODE       BzzMath' >> microPACKEDbed.txt
echo '  DAE       BzzMath' >> microPACKEDbed.txt
echo '}' >> microPACKEDbed.txt
