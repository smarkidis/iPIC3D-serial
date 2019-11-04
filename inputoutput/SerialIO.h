#ifndef SerialIO_H
#define SerialIO_H



#include <iostream>
#include <sstream>
#include <string>
#include <fstream>


using std::string;
using std::stringstream;
using std::ofstream;
using std::cout;
using std::endl;



/** write to binary file Tracked particles */
inline void writeTrackedParticles(int myrank, int cycle,Particles3D *part){
    // stream file to be opened and managed
    string temp;
    stringstream ss;
    stringstream cc;
    ss << myrank;
    cc << cycle;
    temp = "./data/TrackedParticles_proc" + ss.str() + "_cycle"+ cc.str() ;
    temp += ".bin";
    ofstream my_file(temp.c_str(),ios::binary);
    // first data is NOP = number of particles
    int nop = ((int) part->getNOP());
    my_file.write((char*)&nop,sizeof(nop));
    // then ID, X, Y, Z, U, V, W, Q
    double dvalue;
    int ivalue;
    for (int i=0; i < ((int)part->getNOP()); i++){
        ivalue = ((int) part->getParticleID(i));
        my_file.write((char*)&ivalue,sizeof(ivalue));
        dvalue = part->getX(i);
        my_file.write((char*)&dvalue,sizeof(dvalue));
        dvalue = part->getY(i);
        my_file.write((char*)&dvalue,sizeof(dvalue));
        dvalue = part->getZ(i);
        my_file.write((char*)&dvalue,sizeof(dvalue));
        dvalue = part->getU(i);
        my_file.write((char*)&dvalue,sizeof(dvalue));
        dvalue = part->getV(i);
        my_file.write((char*)&dvalue,sizeof(dvalue));
        dvalue = part->getW(i);
        my_file.write((char*)&dvalue,sizeof(dvalue));
        dvalue = part->getQ(i);
        my_file.write((char*)&dvalue,sizeof(dvalue));
    }
    my_file.close();
}



/** write to an ASCII file Electrostatic Potential */
inline void writePHIascii1D(string filename, int myrank, int cycle, Grid *grid,Field* field){
      // stream file to be opened and managed
      string temp;
      stringstream ss;
      stringstream cc;
      ss << myrank;
      cc << cycle;
      temp = "./data/" + filename +"_proc" + ss.str()  + "_cycle"+ cc.str() ;
      temp += ".txt";
      cout << "Opening file: " << temp << endl;
      ofstream my_file(temp.c_str());
	  int nxc = grid->getNXC();
	 
      for (int i=1; i < grid->getNXC()-1; i++){
        my_file <<  grid->getXC(i,0,0) << "\t " << (field->getPHI())[i][0][0]<< endl;
	  }
      my_file.close();
      
     
}
inline void writeRHOascii1D(string filename, int myrank, int cycle, Grid *grid,Field* field){
      // stream file to be opened and managed
      string temp;
      stringstream ss;
      stringstream cc;
      ss << myrank;
      cc << cycle;
      temp = "./data/" + filename +"_proc" + ss.str() + "_cycle"+ cc.str() ;
      temp += ".txt";
      cout << "Opening file: " << temp << endl;
      ofstream my_file(temp.c_str());
	  int nxc = grid->getNXC();
	 
      for (int i=0; i < grid->getNXN(); i++){
        my_file <<  grid->getXN(i,0,0) << "\t " << field->getRHOn(i,0,0)<< endl;
	  }
      my_file.close();
}
inline void writeRHOascii2D(string filename, int myrank, int cycle, Grid *grid,Field* field){
      // stream file to be opened and managed
      string temp;
      stringstream ss;
      stringstream cc;
      ss << myrank;
      cc << cycle;
      temp = "./data/" + filename +"_proc" + ss.str() + "_cycle"+ cc.str() ;
      temp += ".txt";
      cout << "Opening file: " << temp << endl;
      ofstream my_file(temp.c_str());
	  int nxc = grid->getNXC();
	 
       for (int i=1; i < grid->getNXN()-1; i++){
	    for (int j=1; j < grid->getNYN()-1; j++){
          my_file << field->getRHOn(i,j,0)<< "\t";
	    }
		my_file << endl;
	  }
      my_file.close();
}

inline void writeEyascii2D(string filename, int myrank, int cycle, Grid *grid,Field* field){
    // stream file to be opened and managed
    string temp;
    stringstream ss;
    stringstream cc;
    ss << myrank;
    cc << cycle;
    temp = "./data/" + filename +"_proc" + ss.str() + "_cycle"+ cc.str() ;
    temp += ".txt";
    cout << "Opening file: " << temp << endl;
    ofstream my_file(temp.c_str());
    int nxc = grid->getNXC();
    
    for (int i=1; i < grid->getNXN()-1; i++){
	    for (int j=1; j < grid->getNYN()-1; j++){
            my_file << field->getEy(i,j,0)<< "\t";
	    }
		my_file << endl;
    }
    my_file.close();
}

inline void VTK_Write_Scalars(int cycle, Grid *grid,Field* field){
    // stream file to be opened and managed
    string filename = "rhoe";
    string temp;
    stringstream cc;
    cc << cycle;
    temp = "./data/" + filename + "_"+ cc.str() ;
    temp += ".vtk";
    cout << "Opening file: " << temp << endl;
    
    int nxn = grid->getNXN();
    int nyn = grid->getNYN();
    
    double dx = grid->getDX();
    double dy = grid->getDY();
    
    ofstream my_file(temp.c_str());
    my_file << "# vtk DataFile Version 1.0" << endl;
    my_file << "Electron Density from iPIC" << endl;
    my_file << "ASCII" << endl;
    my_file << "DATASET STRUCTURED_POINTS"  << endl;
    my_file << "DIMENSIONS " << (nxn-2) << " " << (nyn-2) << " " << 1 << endl;
    my_file << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
    my_file << "SPACING " << dx << " " << dy << " " << 1.0 << endl;
    my_file << "POINT_DATA " << (nxn-2)*(nyn-2) << endl;
    my_file << "SCALARS rhoe float" << endl;
    my_file << "LOOKUP_TABLE default" << endl;
    
    
	    for (int j=1; j < nyn-1; j++){
            for (int i=1; i < nxn-1; i++){
            my_file << field->getRHOns(i,j,1,0) << endl;
	    }
    }
    my_file.close();
    
    filename = "rhoi";
    temp = "./data/" + filename + "_"+ cc.str() ;
    temp += ".vtk";
    cout << "Opening file: " << temp << endl;
    ofstream my_file2(temp.c_str());
    my_file2 << "# vtk DataFile Version 1.0" << endl;
    my_file2 << "Ion Density from iPIC" << endl;
    my_file2 << "ASCII" << endl;
    my_file2 << "DATASET STRUCTURED_POINTS"  << endl;
    my_file2 << "DIMENSIONS " << (nxn-2) << " " << (nyn-2) << " " << 1 << endl;
    my_file2 << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
    my_file2 << "SPACING " << dx << " " << dy << " " << 1.0 << endl;
    my_file2 << "POINT_DATA " << (nxn-2)*(nyn-2) << endl;
    my_file2 << "SCALARS rhoi float" << endl;
    my_file2 << "LOOKUP_TABLE default" << endl;
    
    
	    for (int j=1; j < grid->getNYN()-1; j++){
            for (int i=1; i < grid->getNXN()-1; i++){
            my_file2 << field->getRHOns(i,j,1,1) << endl;
	    }
    }
    my_file2.close();
    
    filename = "rho_net";
    temp = "./data/" + filename + "_"+ cc.str() ;
    temp += ".vtk";
    cout << "Opening file: " << temp << endl;
    ofstream my_file1(temp.c_str());
    my_file1 << "# vtk DataFile Version 1.0" << endl;
    my_file1 << "Net Charge Density from iPIC" << endl;
    my_file1 << "ASCII" << endl;
    my_file1 << "DATASET STRUCTURED_POINTS"  << endl;
    my_file1 << "DIMENSIONS " << (nxn-2) << " " << (nyn-2) << " " << 1 << endl;
    my_file1 << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
    my_file1 << "SPACING " << dx << " " << dy << " " << 1.0 << endl;
    my_file1 << "POINT_DATA " << (nxn-2)*(nyn-2) << endl;
    my_file1 << "SCALARS rhonet float" << endl;
    my_file1 << "LOOKUP_TABLE default" << endl;
    
    
	    for (int j=1; j < grid->getNYN()-1; j++){
            for (int i=1; i < grid->getNXN()-1; i++){
            my_file1 << field->getRHOn(i,j,0) << endl;
	    }
    }
    my_file1.close();
    
    
}


inline void VTK_Write_Vectors(int cycle, Grid *grid,Field* field){
    // stream file to be opened and managed
    string filename = "E";
    string temp;
    stringstream cc;
    cc << cycle;
    temp = "./data/" + filename + "_"+ cc.str() ;
    temp += ".vtk";
    cout << "Opening file: " << temp << endl;
    
    int nxn = grid->getNXN();
    int nyn = grid->getNYN();
    
    double dx = grid->getDX();
    double dy = grid->getDY();
    
    
    
    ofstream my_fileE(temp.c_str());
    my_fileE << "# vtk DataFile Version 1.0" << endl;
    my_fileE << "E iPIC" << endl;
    my_fileE << "ASCII" << endl;
    my_fileE << "DATASET STRUCTURED_POINTS"  << endl;
    my_fileE << "DIMENSIONS " << (nxn-2) << " " << (nyn-2) << " " << 1 << endl;
    my_fileE << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
    my_fileE << "SPACING " << dx << " " << dy << " " << 1.0 << endl;
    my_fileE << "POINT_DATA " << (nxn-2)*(nyn-2) << endl;
    my_fileE << "VECTORS E float" << endl;
    
    double Ex = 0, Ey = 0, Ez = 0;
    for (int j=1; j < nyn-1; j++){
        for (int i=1; i < nxn-1; i++){
            Ex = field->getEx(i,j,1);
            if (fabs(Ex) < 1E-8)
                Ex = 0.0;
            Ey = field->getEy(i,j,1);
            if (fabs(Ey) < 1E-8)
                Ey = 0.0;
            Ez = field->getEz(i,j,1);
            if (fabs(Ez) < 1E-8)
                Ez = 0.0;
            my_fileE << Ex << " " << Ey <<  " " << Ez <<  endl;
	    }
    }
    my_fileE.close();
    
    filename = "B";
    temp = "./data/" + filename + "_"+ cc.str() ;
    temp += ".vtk";
    cout << "Opening file: " << temp << endl;
    ofstream my_file2(temp.c_str());
    my_file2 << "# vtk DataFile Version 1.0" << endl;
    my_file2 << "B iPIC" << endl;
    my_file2 << "ASCII" << endl;
    my_file2 << "DATASET STRUCTURED_POINTS"  << endl;
    my_file2 << "DIMENSIONS " << (nxn-2) << " " << (nyn-2) << " " << 1 << endl;
    my_file2 << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
    my_file2 << "SPACING " << dx << " " << dy << " " << 1.0 << endl;
    my_file2 << "POINT_DATA " << (nxn-2)*(nyn-2) << endl;
    my_file2 << "VECTORS B float" << endl;
    
    
    for (int j=1; j < nyn-1; j++){
        for (int i=1; i < nxn-1; i++){
            my_file2 << field->getBx(i,j,1) << " " << field->getBy(i,j,1) <<  " " << field->getBz(i,j,1) <<  endl;
	    }
    }
    my_file2.close();
    
    
    filename = "Je";
    temp = "./data/" + filename + "_"+ cc.str() ;
    temp += ".vtk";
    cout << "Opening file: " << temp << endl;
    ofstream my_file1(temp.c_str());
    my_file1 << "# vtk DataFile Version 1.0" << endl;
    my_file1 << "Je from iPIC" << endl;
    my_file1 << "ASCII" << endl;
    my_file1 << "DATASET STRUCTURED_POINTS"  << endl;
    my_file1 << "DIMENSIONS " << (nxn-2) << " " << (nyn-2) << " " << 1 << endl;
    my_file1 << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
    my_file1 << "SPACING " << dx << " " << dy << " " << 1.0 << endl;
    my_file1 << "POINT_DATA " << (nxn-2)*(nyn-2) << endl;
    my_file1 << "VECTORS Je float" << endl;
    
    double Vx, Vy, Vz;
    for (int j=1; j < nyn-1; j++){
        for (int i=1; i < nxn-1; i++){
            Vx=(field->getJxs(i,j,1,0));
            Vy=(field->getJys(i,j,1,0));
            Vz=(field->getJzs(i,j,1,0));
            my_file1 << Vx << " " << Vy <<  " " << Vz <<  endl;
	    }
    }
    my_file1.close();
    
    filename = "Ji";
    temp = "./data/" + filename + "_"+ cc.str() ;
    temp += ".vtk";
    cout << "Opening file: " << temp << endl;
    ofstream my_file3(temp.c_str());
    my_file3 << "# vtk DataFile Version 1.0" << endl;
    my_file3 << "Ion Vel from iPIC" << endl;
    my_file3 << "ASCII" << endl;
    my_file3 << "DATASET STRUCTURED_POINTS"  << endl;
    my_file3 << "DIMENSIONS " << (nxn-2) << " " << (nyn-2) << " " << 1 << endl;
    my_file3 << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
    my_file3 << "SPACING " << dx << " " << dy << " " << 1.0 << endl;
    my_file3 << "POINT_DATA " << (nxn-2)*(nyn-2) << endl;
    my_file3 << "VECTORS Ji float" << endl;
    
    
    for (int j=1; j < nyn-1; j++){
        for (int i=1; i < nxn-1; i++){
            Vx=(field->getJxs(i,j,1,1));
            Vy=(field->getJys(i,j,1,1));
            Vz=(field->getJzs(i,j,1,1));
            my_file3 << Vx << " " << Vy <<  " " << Vz <<  endl;
	    }
    }
    my_file3.close();
    
    
}





inline void writeRHOascii2Del(string filename, int myrank, int cycle, Grid *grid,Field* field){
      // stream file to be opened and managed
      string temp;
      stringstream ss;
      stringstream cc;
      ss << myrank;
      cc << cycle;
      temp = "./data/" + filename +"_proc" + ss.str() + "_cycle"+ cc.str() ;
      temp += ".txt";
      cout << "Opening file: " << temp << endl;
      ofstream my_file(temp.c_str());
	  int nxc = grid->getNXC();
	 
      for (int i=1; i < grid->getNXN()-1; i++){
	    for (int j=1; j < grid->getNYN()-1; j++){
          my_file << field->getRHOns(i,j,0,0)<< "\t";
	    }
		my_file << endl;
	  }
      my_file.close();
}
inline void writeVyascii2Del(string filename, int myrank, int cycle, Grid *grid,Field* field){
    // stream file to be opened and managed
    string temp;
    stringstream ss;
    stringstream cc;
    ss << myrank;
    cc << cycle;
    temp = "./data/" + filename +"_proc" + ss.str() + "_cycle"+ cc.str() ;
    temp += ".txt";
    cout << "Opening file: " << temp << endl;
    ofstream my_file(temp.c_str());
    int nxc = grid->getNXC();
    double v = 0.0;
    for (int i=1; i < grid->getNXN()-1; i++){
	    for (int j=1; j < grid->getNYN()-1; j++){
            if (field->getRHOns(i,j,0,0)==0)
                v = 0.0;
            else
                v = ((field->getJys(i,j,0,0))/(field->getRHOns(i,j,0,0)));
            my_file << v << "\t";
	    }
		my_file << endl;
    }
    my_file.close();
}
inline void writeVxascii2Del(string filename, int myrank, int cycle, Grid *grid,Field* field){
    // stream file to be opened and managed
    string temp;
    stringstream ss;
    stringstream cc;
    ss << myrank;
    cc << cycle;
    temp = "./data/" + filename +"_proc" + ss.str() + "_cycle"+ cc.str() ;
    temp += ".txt";
    cout << "Opening file: " << temp << endl;
    ofstream my_file(temp.c_str());
    int nxc = grid->getNXC();
    
    double v;
    
    for (int i=1; i < grid->getNXN()-1; i++){
	    for (int j=1; j < grid->getNYN()-1; j++){
            if (field->getRHOns(i,j,0,0)==0)
                v = 0.0;
            else
                v = ((field->getJxs(i,j,0,0))/(field->getRHOns(i,j,0,0)));
            my_file << v << "\t";
	    }
		my_file << endl;
    }
    my_file.close();
}
inline void writeVxascii2Dion(string filename, int myrank, int cycle, Grid *grid,Field* field){
    // stream file to be opened and managed
    string temp;
    stringstream ss;
    stringstream cc;
    ss << myrank;
    cc << cycle;
    temp = "./data/" + filename +"_proc" + ss.str() + "_cycle"+ cc.str() ;
    temp += ".txt";
    cout << "Opening file: " << temp << endl;
    ofstream my_file(temp.c_str());
    int nxc = grid->getNXC();
    double v;
    
    for (int i=1; i < grid->getNXN()-1; i++){
	    for (int j=1; j < grid->getNYN()-1; j++){
            if (field->getRHOns(i,j,0,1)==0)
                v = 0.0;
            else
                v = ((field->getJxs(i,j,0,1))/(field->getRHOns(i,j,0,1)));
            my_file << v << "\t";
	    }
		my_file << endl;
    }
    my_file.close();
}
inline void writeVyascii2Dion(string filename, int myrank, int cycle, Grid *grid,Field* field){
    // stream file to be opened and managed
    string temp;
    stringstream ss;
    stringstream cc;
    ss << myrank;
    cc << cycle;
    temp = "./data/" + filename +"_proc" + ss.str() + "_cycle"+ cc.str() ;
    temp += ".txt";
    cout << "Opening file: " << temp << endl;
    ofstream my_file(temp.c_str());
    int nxc = grid->getNXC();
    double v = 0;
    for (int i=1; i < grid->getNXN()-1; i++){
	    for (int j=1; j < grid->getNYN()-1; j++){
            if (field->getRHOns(i,j,0,1)==0)
                v = 0.0;
            else
                v = ((field->getJys(i,j,0,1))/(field->getRHOns(i,j,0,1)));
            my_file << v << "\t";
	    }
		my_file << endl;
    }
    my_file.close();
}
inline void writeRHOascii2Dion(string filename, int myrank, int cycle, Grid *grid,Field* field){
      // stream file to be opened and managed
      string temp;
      stringstream ss;
      stringstream cc;
      ss << myrank;
      cc << cycle;
      temp = "./data/" + filename +"_proc" + ss.str() + "_cycle"+ cc.str() ;
      temp += ".txt";
      cout << "Opening file: " << temp << endl;
      ofstream my_file(temp.c_str());
	  int nxc = grid->getNXC();
	 
       for (int i=1; i < grid->getNXN()-1; i++){
	    for (int j=1; j < grid->getNYN()-1; j++){
          my_file << field->getRHOns(i,j,0,1)<< "\t";
	    }
		my_file << endl;
	  }
      my_file.close();
}
/** write an ASCII with Ex */
inline void writeExascii1D(string filename, int myrank, int cycle, Grid *grid,Field* field){
      // stream file to be opened and managed
      string temp;
      stringstream ss;
      stringstream cc;
      ss << myrank;
      cc << cycle;
      temp = "./data/" + filename + "_cycle"+ cc.str() +"_proc" + ss.str();
      temp += ".txt";
      cout << "Opening file: " << temp << endl;
      ofstream my_file(temp.c_str());
      for (int i=0; i < (grid->getNXN()); i++){
         my_file <<  field->getRHOns(i,0,0,0)  << endl;
      }
      my_file.close();
}
/** write an ASCII with the heat flux */
inline void writeJx1D(string filename, int myrank, int cycle, Grid *grid,Field* field){
      // stream file to be opened and managed
      string temp;
      stringstream ss;
      stringstream cc;
      ss << myrank;
      cc << cycle;
      temp = "./data/" + filename +"_proc" + ss.str() + "_cycle"+ cc.str() ;
      temp += ".txt";
      cout << "Opening file: " << temp << endl;
      ofstream my_file(temp.c_str());
      for (int i=0; i < (grid->getNXN()); i++){
         my_file <<  field->getJx(i,0,0)  << endl;
      }
      my_file.close();
}
/** write an ASCII with the heat flux */
inline void writeHeatFluxascii1D(string filename, int myrank, int cycle, Grid *grid,Field* field){
      // stream file to be opened and managed
      string temp;
      stringstream ss;
      stringstream cc;
      ss << myrank;
      cc << cycle;
      temp = "./data/" + filename + "_proc" + ss.str() +  "_cycle"+ cc.str() ;
      temp += ".txt";
      cout << "Opening file: " << temp << endl;
      ofstream my_file(temp.c_str());
      for (int i=0; i < (grid->getNXN()); i++){
         my_file <<  field->getJy(i,0,0)  << endl;
      }
      my_file.close();
}
/** write to an ASCII file  particles position */
inline void writeParticles1D(string filename, int myrank, int cycle,Particles *part){
      // stream file to be opened and managed
      string temp;
      stringstream ss;
      stringstream cc;
      ss << myrank;
      cc << cycle;
      temp = "./data/" + filename +"_proc" + ss.str() + "_cycle"+ cc.str() ;
      temp += ".txt";
      cout << "Opening file: " << temp << endl;
      ofstream my_file(temp.c_str());
      for (int i=0; i < part->getNOP(); i++){
        my_file << part->getY(i) << "\t" << part->getV(i) << endl;
      }
      my_file.close();
}
/** write to an ASCII file  particles position */
inline void writeParticles1D(string filename, int myrank, int cycle,Particles *part,int p){
      // stream file to be opened and managed
      string temp;
      stringstream ss;
      stringstream cc;
      ss << myrank;
      cc << cycle;
      temp = "./data/" + filename + "_proc" + ss.str() + "_cycle"+ cc.str() ;
      temp += ".txt";
      cout << "Opening file: " << temp << endl;
      ofstream my_file(temp.c_str());
      for (int i=0; i < part->getNOP(); i=i+p){
           my_file << part->getY(i) << "\t" << part->getV(i) << endl;
      }
      my_file.close();
}



#endif
