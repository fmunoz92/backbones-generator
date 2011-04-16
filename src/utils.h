#include "petu.h"

//This function puts 0 in all the atoms coordinates
void clearatm(ATOM* patm, int nres);

void setr(float rn, float rca, float rc, float scal_1_4, float scal_1_5);

//This function calculates the gyration radius of the CA atoms
//It was meant to filter long chains. Now it is replaced by the volume filter
FilterResultType calcRdG(ATOM* patm, int nres, float rgmax);

//This function filter long chains.
//Calculates the CA-CA distances of the recently added CA with the rest.
//Then compares this distances with dmax2, wich cames from the analisys of protein
//structures experimetally determined
FilterResultType islong(const ATOM* patm, int at, float dmax2);

FilterResultType isclash(const ATOM* patm, int at);

FilterResultType volumen_en_rango(int nres, Volume vol_parcial);

