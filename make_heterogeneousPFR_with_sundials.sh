source mybashrc

mkdir exe
rm -f exe/heterogeneousPFR_with_sundials.sh

cd src
g++ -O3 -m64 -Wno-write-strings -fpermissive heterogeneousPFR/main.C $SUPPORT -I$MKL_INCLUDE -IheterogeneousPFR -IheterogeneousPFR/utilities -IheterogeneousPFR/ode -IheterogeneousPFR/bvp -I. -I$EIGEN -I$BOOST/include/ -I$OPENSMOKE/ -I$RAPIDXML -I$SUNDIALS_INCLUDE -L$BOOST/lib/ -lgfortran -L$SUNDIALS_LIBS_PATH $SUNDIALS_LIBS -L$MKL_LIBS_PATH $MKL_LIBS -Wl,--no-as-needed -ldl -lpthread -lm -lboost_program_options -lboost_system -lboost_filesystem -lboost_regex -lboost_date_time -o ../exe/heterogeneousPFR_with_sundialss.sh
cd ..
