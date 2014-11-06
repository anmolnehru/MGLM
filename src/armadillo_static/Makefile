all:
	g++ -static legr_core.cpp -o legr_dti -O3 liblapack.a  librefblas.a -Iarmadillo-3.930.2/include/ -lgfortran -lpthread -lboost_filesystem -lboost_system
armabin:
	g++ -static test_armabin.cpp -o test_armabin -O3 liblapack.a  librefblas.a -Iarmadillo-3.930.2/include/ -lgfortran -lpthread -lboost_filesystem -lboost_system
