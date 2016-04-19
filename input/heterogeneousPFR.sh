touch heterogeneousPFR.txt

echo '/*##############################################################################################' >> heterogeneousPFR.txt
echo '#                                                                                              #' >> heterogeneousPFR.txt
echo '#     #############       #############       #############       ####                ####     #' >> heterogeneousPFR.txt
echo '#    #             #     #             #     #             #     #    #              #    #    #' >> heterogeneousPFR.txt
echo '#    #    #####    #     #    #########      #    #####    #     #    #              #    #    #' >> heterogeneousPFR.txt
echo '#    #    #   #    #     #    #              #    #   #    #     #    #              #    #    #' >> heterogeneousPFR.txt
echo '#    #    #####    #     #    #              #    #####    #     #    #              #    #    #' >> heterogeneousPFR.txt
echo '#    #             #     #    #########      #             #     #    #              #    #    #' >> heterogeneousPFR.txt
echo '#    #             #     #             #     #             #     #    #              #    #    #' >> heterogeneousPFR.txt
echo '#    #    #####    #      #########    #     #    #####    #     #    #              #    #    #' >> heterogeneousPFR.txt
echo '#    #    #   #    #              #    #     #    #   #    #     #    #              #    #    #' >> heterogeneousPFR.txt
echo '#    #    #   #    #      #########    #     #    #   #    #     #    #########      #    #    #' >> heterogeneousPFR.txt
echo '#    #    #   #    #     #             #     #    #   #    #     #             #     #    #    #' >> heterogeneousPFR.txt
echo '#     ####     ####       #############       ####     ####       #############       ####     #' >> heterogeneousPFR.txt
echo '#                                                                                              #' >> heterogeneousPFR.txt
echo '#   Department of Energy                                                                       #' >> heterogeneousPFR.txt
echo '#   Politecnico di Milano                                                                      #' >> heterogeneousPFR.txt
echo '#   Author: Stefano Rebughini <stefano.rebughini@polimi.it>                                    #' >> heterogeneousPFR.txt
echo '#                                                                                              #' >> heterogeneousPFR.txt
echo '################################################################################################' >> heterogeneousPFR.txt
echo ' ' >> heterogeneousPFR.txt
echo 'Kinetics path           CH4_UBI/kinetics' >> heterogeneousPFR.txt
echo ' ' >> heterogeneousPFR.txt
echo 'Volumetric flow rate    10 Nl/min' >> heterogeneousPFR.txt
echo ' ' >> heterogeneousPFR.txt
echo 'Pressure                1 bar' >> heterogeneousPFR.txt
echo ' ' >> heterogeneousPFR.txt
echo 'Temperature' >> heterogeneousPFR.txt
echo '{' >> heterogeneousPFR.txt
echo '  solid                 500 °C' >> heterogeneousPFR.txt
echo '  gas                    40 °C' >> heterogeneousPFR.txt
echo '  wall                  100 °C' >> heterogeneousPFR.txt
echo '}' >> heterogeneousPFR.txt
echo ' ' >> heterogeneousPFR.txt
echo 'Reactor' >> heterogeneousPFR.txt
echo '{' >> heterogeneousPFR.txt
echo '  type                      honeyComb' >> heterogeneousPFR.txt
echo '  inert length              0.5  cm' >> heterogeneousPFR.txt
echo '  catalytic length          26   mm' >> heterogeneousPFR.txt
echo '  hydraulic diameter        24   mm' >> heterogeneousPFR.txt
echo '  channel diameter          1.1  mm' >> heterogeneousPFR.txt
echo '  channel shape             square' >> heterogeneousPFR.txt
echo '  void fraction             0.75' >> heterogeneousPFR.txt
echo '  gas-to-solid correlation  massTransfer' >> heterogeneousPFR.txt
echo '}' >> heterogeneousPFR.txt
echo ' ' >> heterogeneousPFR.txt
echo 'Catalyst' >> heterogeneousPFR.txt
echo '{' >> heterogeneousPFR.txt
echo '  mass                  500 mg' >> heterogeneousPFR.txt
echo '  dispersion            0.2' >> heterogeneousPFR.txt
echo '  Rh fraction           0.02' >> heterogeneousPFR.txt
echo '  deactivation factor   1.0' >> heterogeneousPFR.txt
echo '}' >> heterogeneousPFR.txt
echo ' ' >> heterogeneousPFR.txt
echo 'Solid' >> heterogeneousPFR.txt
echo '{' >> heterogeneousPFR.txt
echo '  density               2300 Kg/m3' >> heterogeneousPFR.txt
echo '  conductivity          2.5  W/m/K' >> heterogeneousPFR.txt
echo '  specific heat         925 J/Kg/K' >> heterogeneousPFR.txt
echo '} ' >> heterogeneousPFR.txt
echo ' ' >> heterogeneousPFR.txt
echo 'Mole fractions' >> heterogeneousPFR.txt
echo '{' >> heterogeneousPFR.txt
echo '  CH4                   0.27' >> heterogeneousPFR.txt
echo '  O2                    0.1593' >> heterogeneousPFR.txt
echo '  N2                    0.5707' >> heterogeneousPFR.txt
echo '  inert                 N2' >> heterogeneousPFR.txt
echo '}' >> heterogeneousPFR.txt
echo ' ' >> heterogeneousPFR.txt
echo 'Solver options' >> heterogeneousPFR.txt
echo '{' >> heterogeneousPFR.txt
echo '  Constraints          false' >> heterogeneousPFR.txt
echo ' ' >> heterogeneousPFR.txt
echo '  Absolute tollerance' >> heterogeneousPFR.txt
echo '  (' >> heterogeneousPFR.txt
echo '    specie           1e-15' >> heterogeneousPFR.txt
echo '    temperature      1e-08' >> heterogeneousPFR.txt
echo '  )' >> heterogeneousPFR.txt
echo ' ' >> heterogeneousPFR.txt
echo '  Relative tollerance' >> heterogeneousPFR.txt
echo '  (' >> heterogeneousPFR.txt
echo '    specie           1e-08' >> heterogeneousPFR.txt
echo '    temperature      1e-08' >> heterogeneousPFR.txt
echo '  )' >> heterogeneousPFR.txt
echo ' ' >> heterogeneousPFR.txt
echo '  Energy equation    on' >> heterogeneousPFR.txt
echo ' ' >> heterogeneousPFR.txt
echo '  Reactions' >> heterogeneousPFR.txt
echo '  (' >> heterogeneousPFR.txt
echo '    Homogeneous      off' >> heterogeneousPFR.txt
echo '    Heterogeneous    on' >> heterogeneousPFR.txt
echo '  )' >> heterogeneousPFR.txt
echo ' ' >> heterogeneousPFR.txt
echo '  Grid' >> heterogeneousPFR.txt
echo '  (' >> heterogeneousPFR.txt
echo '    Growing                   on' >> heterogeneousPFR.txt
echo '    Minimum number of points  10' >> heterogeneousPFR.txt
echo '    Maximum number of points  800' >> heterogeneousPFR.txt
echo '    Add number of points      300' >> heterogeneousPFR.txt
echo '    Resolution type           new' >> heterogeneousPFR.txt
echo '    Discretization scheme     BDS' >> heterogeneousPFR.txt
echo '  )' >> heterogeneousPFR.txt
echo ' ' >> heterogeneousPFR.txt
echo '  Diffusion in gas off' >> heterogeneousPFR.txt
echo ' ' >> heterogeneousPFR.txt
echo '  External heat exchange off' >> heterogeneousPFR.txt
echo ' ' >> heterogeneousPFR.txt
echo '  Accepted errors' >> heterogeneousPFR.txt
echo '  (' >> heterogeneousPFR.txt
echo '    grid               1e-08' >> heterogeneousPFR.txt
echo '    energy             1e-05' >> heterogeneousPFR.txt
echo '  )' >> heterogeneousPFR.txt
echo ' ' >> heterogeneousPFR.txt
echo '  Integration time            100 s' >> heterogeneousPFR.txt
echo ' ' >> heterogeneousPFR.txt
echo '  Results                     mole' >> heterogeneousPFR.txt
echo '}' >> heterogeneousPFR.txt
echo ' ' >> heterogeneousPFR.txt
echo 'Numerical solvers' >> heterogeneousPFR.txt
echo '{' >> heterogeneousPFR.txt
echo '  ODE       BzzMath' >> heterogeneousPFR.txt
echo '  DAE       BzzMath' >> heterogeneousPFR.txt
echo '}' >> heterogeneousPFR.txt
