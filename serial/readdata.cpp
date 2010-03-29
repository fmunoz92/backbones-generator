
#include "readdata.h"

// Las pongo aca porque creo que no tiene utilidad en ningun otro lado. Esta bien?
const float deg2rad_ratio = M_PI/180;
inline float deg2rad ( float deg) {
	return deg * deg2rad_ratio;
}

void readdata(std::ifstream &filer, std::vector<float> &cosfi, std::vector<float> &sinfi, std::vector<float> &cossi, std::vector<float> &sinsi, AnglesMapping *angles_mapping ) 
{

        float fi = 0.0f;
        float si = 0.0f;
	unsigned int i = 0;
	while(filer.good()) { 	
		
		if(filer >> fi && filer >> si) {
		  cosfi.push_back(cos(deg2rad(fi)));
		  sinfi.push_back(sin(deg2rad(fi)));
		  cossi.push_back(cos(deg2rad(si)));
		  sinsi.push_back(sin(deg2rad(si)));
		  angles_mapping->set_mapping(fi, si);
		}
		++i;
	}
	
	// Nota a futuro: se deberia lanzar una Excepcion si el formato del archivo fuera equivocado. 
}



