DYLD_LIBRARY_PATH = 

CPP = mpicxx
OPTFLAGS=  -O2 -DMPICH_IGNORE_CXX_SEEK
# include files
INC_HDF5 = 
# libs
LIB_HDF5 = 

HDF5LIBS = -lhdf5 -lhdf5_hl 
MPELIB = -lmpe


ipic3D: iPIC3D.cpp Particles3Dcomm.o Particles3D.o TrackedParticles.o ConfigFile.o
	${CPP} ${OPTFLAGS} -o iPIC3D ${INC_HDF5} ${INC_MPI} \
	iPIC3D.cpp Particles3Dcomm.o Particles3D.o TrackedParticles.o ConfigFile.o ${LIB_HDF5} ${LIB_MPI} ${HDF5LIBS} ${MPELIB}

iPIC3D.o: iPIC3D.cpp
	${CPP} ${OPTFLAGS} ${INC_HDF5} ${INC_MPI} -c iPIC3D.cpp 

ConfigFile.o: ./ConfigFile/src/ConfigFile.cpp
	${CPP} ${OPTFLAGS} -c ./ConfigFile/src/ConfigFile.cpp

Particles3Dcomm.o: ./particles/Particles3Dcomm.cpp
	${CPP} ${OPTFLAGS} ${INC_HDF5} -c ./particles/Particles3Dcomm.cpp

TrackedParticles.o: ./particles/TrackedParticles.cpp
	${CPP} ${OPTFLAGS} ${INC_HDF5} -c ./particles/TrackedParticles.cpp

Particles3D.o: ./particles/Particles3D.cpp 
	${CPP} ${OPTFLAGS} ${INC_HDF5} -c ./particles/Particles3D.cpp

clean:
	rm -rf *.o iPIC3D
