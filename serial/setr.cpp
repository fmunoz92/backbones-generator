#include "setr.h"


void setr( float rn, float rca, float rc, float scal_1_4, float scal_1_5)
{
  /* La matriz de distancias minimas al cuadrado*/

	int i,j;
	float radio[3];

	radio[0] = rn;
	radio[1] = rca;
	radio[2] = rc;

	for(i = 0;i <= 2;i++) 
	{
		for(j = 0;j <= 2;j++) 
		{
			r[i][j][0] = (radio[i] + radio[j]) * (radio[i] + radio[j]) * scal_1_4 * scal_1_4;
			r[i][j][1] = (radio[i] + radio[j]) * (radio[i] + radio[j]) * scal_1_5 * scal_1_5;
			r[i][j][2] = (radio[i] + radio[j]) * (radio[i] + radio[j]);      
		}
	}  

}
