all:
	g++ -g -static legr_core.cpp -o legr_dti -O3 liblapack.a  librefblas.a -Iarmadillo-3.930.2/include/ -lgfortran -lpthread -lboost_filesystem -lboost_system
multithread:
	g++ -g -static legr_core_multi_thread.cpp -o legr_dti_parallel -fopenmp -O3 liblapack.a  librefblas.a -Iarmadillo-3.930.2/include/ -lgfortran -lpthread -lboost_filesystem -lboost_system
mglm:
	g++ -static legr_core.cpp -o mglm_dti -O3 liblapack.a  librefblas.a -Iarmadillo-3.930.2/include/ -lgfortran -lpthread -lboost_filesystem -lboost_system
armabin:
	g++ -static test_armabin.cpp -o test_armabin -O3 liblapack.a  librefblas.a -Iarmadillo-3.930.2/include/ -lgfortran -lpthread -lboost_filesystem -lboost_system
show_arma_mx:
	g++ -static show_arma_mx.cpp -o show_arma_mx -O2 liblapack.a librefblas.a -Iarmadillo-3.930.2/include/ -lgfortran -lpthread -lboost_system -lboost_filesystem
armabin2asciifull:
	g++ -static armabin2asciifull.cpp -o armabin2asciifull -O2 liblapack.a librefblas.a -Iarmadillo-3.930.2/include/ -lgfortran -lpthread -lboost_system -lboost_filesystem
test_armabin:
	g++ -static test_armabin.cpp -o test_armabin -O2 liblapack.a librefblas.a -Iarmadillo-3.930.2/include/ -lgfortran -lpthread -lboost_system -lboost_filesystem
armabinequaltotxt:
	g++ -static armabinequaltotxt.cpp -o armabinequaltotxt -O2 liblapack.a librefblas.a -Iarmadillo-3.930.2/include/ -lgfortran -lpthread -lboost_system -lboost_filesystem

clean:

	rm *.o *.BIN 
