/*
    Spherical Grid Volume Approximation: An algorithm for approximating the
    volume of a compound object by keeping track of its members' positions in space.
    Copyright (C) 2009, Rodrigo Castano, FuDePAN

    This file is part of Spherical Grid Volume Approximation.

    Spherical Grid Volume Approximation is free software: you can redistribute
    it and/or modify it under the terms of the GNU General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Spherical Grid Volume Approximation is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "grillado.h"

// Se asigna memoria a la matriz y se inicializa con 0s.
// Para que no hayan intersecciones diagonales
// la suma de dos radios tiene que ser menor a la
// diagonal. Por Pitagoras, si x es la diagonal,
// x^2 = D^2 + D^2, luego buscamos que 2*R < sqrt(2)*D.
// Que es lo mismo que sqrt(2)*R < D.
Grillado::Grillado(size_t M, size_t N, size_t Z, Length R, Length D) throw(radiusOrDistanceParamException, sizeParamException, bad_alloc)
    : v(M), w(N), z(Z), r(R), d(D), vol_esfera(calc_vol(R)), vol_int(calc_int(R, D)), intersecciones(0), esferas(0)
{
    const Length paramDiagonal = sqrt(2.0f) * R;
    const Length cota = (2.0f) * R;
    if (!((paramDiagonal) < D && (D < cota)))
    {
        throw radiusOrDistanceParamException();
    }
    if (!(v > 2 || w > 2 || z > 2))
    {
        throw sizeParamException();
    }
    matriz = new unsigned int** [M];
    for (size_t m = 0; m < M; m++)
    {
        matriz[m] = new unsigned int* [N];
        for (size_t n = 0; n < N; n++)
        {
            matriz[m][n] = new unsigned int[Z];
            for (size_t x = 0; x < Z; x++)
            {
                matriz[m][n][x] = 0;
            }
        }
    }

}


/*  Se muestra el volumen de las esferas, el volumen de las intersecciones
    y la cantidad de elementos en cada posicion del grillado.
*/
/*  Shows the spheres' volume, the intersections' volume and the number of elements in each
    position of the grid.
*/
void Grillado::info_grillado() const
{
    for (size_t m = 0; m < v; m++)
    {
        for (size_t n = 0; n < w; n++)
        {
            for (size_t Z = 0; Z < z; Z++)
            {
                cerr << "[" << matriz[m][n][Z] << "]";
            }
            cerr << endl;
        }
        cerr << endl << endl;
    }
    cerr << "Sphere volume: " << vol_esfera << endl;
    cerr << "Intersection volume: " << vol_int << endl;
}

// Se libera la memoria correspondiente a la matriz
// All memory resources are released.
Grillado::~Grillado()
{
    for (size_t m = 0; m < v; m++)
    {
        for (size_t n = 0; n < w; n++)
        {
            delete [] matriz[m][n];
        }
        delete [] matriz[m];
    }
    delete [] matriz;
}


// Se agrega un elemento a la esfera correspondiente a las coordenadas ingresadas
// como argumento.
// Adds an element to the sphere determined by coordinates x, y, z.
void Grillado::agregar_esfera(const esferaId& id)
{
    const GridCoord coord_x = id.get_m();
    const GridCoord coord_y = id.get_n();
    const GridCoord coord_z = id.get_z();
    if (matriz[coord_x][coord_y][coord_z] == 0)
    {
        aumentar_vol_parcial(coord_x, coord_y, coord_z);
    }
    matriz[coord_x][coord_y][coord_z]++;
}

// Se quita un elemento de la esfera correspondiente a las coordenadas ingresadas
// como argumento.
// Removes an element from the sphere determined by coordinates x, y, z.
void Grillado::sacar_esfera(const esferaId& id) throw(gridPositionException)
{
    const GridCoord coord_x = id.get_m();
    const GridCoord coord_y = id.get_n();
    const GridCoord coord_z = id.get_z();
    if (matriz[coord_x][coord_y][coord_z] == 0)
    {
        throw gridPositionException();
    }
    matriz[coord_x][coord_y][coord_z]--;
    if (matriz[coord_x][coord_y][coord_z] == 0)
    {
        reducir_vol_parcial(coord_x, coord_y, coord_z);
    }
}

// Reduce el volumen parcial del grillado, se encapsula parte del comportamiento
// de sacar_esfera.
// Reduces the accumulated volume of the grid. It encapsulates part of the behaviour in
// sacar_esfera (removes a sphere)
void Grillado::reducir_vol_parcial(int coord_x, int coord_y, int coord_z)
{
    --esferas;
    intersecciones -= calcular_intersecciones(coord_x, coord_y, coord_z);
}

// Aumenta el volumen parcial del grillado, se encapsula parte del comportamiento
// de agregar_esfera.
// Increases the accumulated volume of the grid. It encapsulates part of the behaviour in
// sacar_esfera (removes a sphere)
void Grillado::aumentar_vol_parcial(int coord_x, int coord_y, int coord_z)
{
    ++esferas;
    intersecciones += calcular_intersecciones(coord_x, coord_y, coord_z);
}

// Asigna a coord_x, coord_y, coord_z los indices del grillado correspondientes a las
// coordenadas x, y, z.
// Sets coord_x, coord_y, coord_z with the grid coordinates corresponding to x, y, z.
void Grillado::calcular_coord(Coord x, Coord y, Coord Z, GridCoord& coord_x, GridCoord& coord_y, GridCoord& coord_z) const
{
    coord_x = modulo(static_cast<int>(round(x / d)), static_cast<int>(v));
    coord_y = modulo(static_cast<int>(round(y / d)), static_cast<int>(w));
    coord_z = modulo(static_cast<int>(round(Z / d)), static_cast<int>(z));
}

/*  Se calcula la cantidad de intersecciones de la esfera ubicada en el grillado en
    coord_x, coord_y, coord_z con esferas vecinas que contengan algun elemento.
    Solo se tienen en cuenta las esferas cuyo centro pertenece a una recta paralela a
    alguno de los ejes que ademas contenga al centro de la esfera cuyas intersecciones
    se estan calculando.
    Es decir no se toman en cuenta intersecciones "diagonales".
    Nota: estas intersecciones diagonales no pueden existir en la actual implementacion ya
    que las combinaciones de dist y radio que las generarian son consideradas
    invalidas por el constructor de la clase.
*/
/*
    Returns the number of intersections with direct neighbours of the sphere determined by coord_x,
    coord_y, coord_z.
*/
unsigned int Grillado::calcular_intersecciones(GridCoord coord_x, GridCoord coord_y, GridCoord coord_z) const
{
    const GridCoord x_plus  = modulo((static_cast<int>(coord_x) + 1) , v);
    const GridCoord x_minus = modulo((static_cast<int>(coord_x) - 1) , v);
    const GridCoord y_minus = modulo((static_cast<int>(coord_y) - 1) , w);
    const GridCoord y_plus  = modulo((static_cast<int>(coord_y) + 1) , w);
    const GridCoord z_plus  = modulo((static_cast<int>(coord_z) + 1) , z);
    const GridCoord z_minus = modulo((static_cast<int>(coord_z) - 1) , z);

    unsigned int inters = matriz[coord_x][coord_y][z_plus] != 0;
    inters += matriz[coord_x][coord_y][z_minus] != 0;
    inters += matriz[x_plus][coord_y][coord_z] != 0;
    inters += matriz[x_minus][coord_y][coord_z] != 0;
    inters += matriz[coord_x][y_plus][coord_z] != 0;
    inters += matriz[coord_x][y_minus][coord_z] != 0;
    return inters;
}

// Sets every grid position to zero. Sets the partial volume to zero.
// Asigna cero a cada posicion del grillado y al volumen parcial.
void Grillado::reset()
{
    for (size_t m = 0; m < v; m++)
    {
        for (size_t n = 0; n < w; n++)
        {
            for (size_t Z = 0; Z < z; Z++)
            {
                matriz[m][n][Z] = 0;
            }
        }
    }
    esferas = 0;
    intersecciones = 0;
}
