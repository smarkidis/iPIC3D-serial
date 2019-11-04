CPP = g++
OPTFLAGS=  -O2 




ipic3D: iPIC3D.cpp Particles3Dcomm.o Particles3D.o ConfigFile.o
	${CPP} ${OPTFLAGS} -o iPIC3D  ${INC_MPI} \
	iPIC3D.cpp Particles3Dcomm.o Particles3D.o ConfigFile.o ${LIB_MPI}  

iPIC3D.o: iPIC3D.cpp
	${CPP} ${OPTFLAGS} -c iPIC3D.cpp 

ConfigFile.o: ./ConfigFile/src/ConfigFile.cpp
	${CPP} ${OPTFLAGS} -c ./ConfigFile/src/ConfigFile.cpp

Particles3Dcomm.o: ./particles/Particles3Dcomm.cpp
	${CPP} ${OPTFLAGS} -c ./particles/Particles3Dcomm.cpp


Particles3D.o: ./particles/Particles3D.cpp 
	${CPP} ${OPTFLAGS} -c ./particles/Particles3D.cpp

clean:
	rm -rf *.o iPIC3D
