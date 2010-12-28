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



#ifndef GRILLADO_H
#define GRILLADO_H

#include <cmath>
#include <iostream>
#include <stdexcept>

typedef float Coord;
typedef float Volume;
typedef float Length;
typedef unsigned int GridCoord;
static const Length PI = atan(1.0) * 4.0f;

using namespace std;

class esferaId
{
public:
    esferaId() {};
    esferaId(GridCoord a, GridCoord b, GridCoord c);
    esferaId(const esferaId& s);
    inline GridCoord get_m() const;
    inline GridCoord get_n() const;
    inline GridCoord get_z() const;
private:
    GridCoord v;
    GridCoord w;
    GridCoord z;
};

class Grillado
{
private:
    const size_t v;
    const size_t w;
    const size_t z;
    const Length r;
    const Length d;
    const Volume vol_esfera;
    const Volume vol_int;

    unsigned int*** matriz;
    unsigned int intersecciones;
    unsigned int esferas;

    inline static Volume calc_vol(Length radio);

    inline static Volume calc_int(Length radio, Length dist);

    void reducir_vol_parcial(int coord_x, int coord_y, int coord_z) ;
    void aumentar_vol_parcial(int coord_x, int coord_y, int coord_z) ;
    void calcular_coord(Coord coord_x, Coord coord_y, Coord coord_z,
                        GridCoord& x, GridCoord& y, GridCoord& z) const;
    unsigned int calcular_intersecciones(GridCoord coord_x, GridCoord coord_y, GridCoord coord_z) const;
    inline esferaId obtener_id(Coord x, Coord y, Coord z) const;

public:
    struct gridPositionException : public domain_error
    {
        gridPositionException()
            : domain_error("The desired position is already empty") {}

    };
    struct sizeParamException : public invalid_argument
    {
        sizeParamException()
            : invalid_argument("M<3 || N<3 || Z<3") {}

    };

    struct radiusOrDistanceParamException : public invalid_argument
    {
        radiusOrDistanceParamException()
            : invalid_argument("Either 1) sqrt(2)*radio >= dist: there may be diagonal intersections between spheres (spheres which differ in more than one coordinate) or 2) dist>=2*R: there are no intersections between spheres.") {}

    };

    Grillado(size_t a, size_t b, size_t c, Length R, Length D) throw
    (radiusOrDistanceParamException, sizeParamException, bad_alloc);

    ~Grillado() ;

    void info_grillado() const;

    inline esferaId agregar_esfera(Coord x, Coord y, Coord z) ;

    void agregar_esfera(const esferaId& id);

    inline void sacar_esfera(Coord x, Coord y, Coord z) throw(gridPositionException) ;

    void sacar_esfera(const esferaId& id) throw(gridPositionException);

    // Devuelve el volumen parcial.
    inline Volume obtener_vol_parcial() const;
    void reset();
};

// Devuelve el volumen parcial.
// Returns the accumulated volume
inline Volume Grillado::obtener_vol_parcial() const
{
    return Volume(esferas) * vol_esfera - Volume(intersecciones) * vol_int;
}

// Se calcula el volumen de las esferas.
// Calculates the volume of a sphere.
inline Volume Grillado::calc_vol(Length radio)
{
    return ((4.0f) * PI * radio * radio * radio / (3.0f));
}

// Devuelve "a (mod m)". El operador % devuelve numeros negativos en algunos casos.
// Returns "a (mod m)". The % operator returns negative values in some cases.
inline GridCoord modulo(int a, int m)
{
    a = (a % m);
    if (a < 0)
    {
        a += m;
    }
    return static_cast<unsigned int>(a);
}

// Se calcula el volumen de las intersecciones entre esferas.
// Calculates the volume of the intersection between neighbour spheres.
inline Volume Grillado::calc_int(Length radio, Length dist)
{
    return (2.0f) * PI * (((2.0f) * radio * radio * radio / (3.0f)) - (dist * radio * radio / (2.0f)) + (dist * dist * dist / (24.0f)));
}





inline GridCoord esferaId::get_m() const
{
    return v;
}
inline GridCoord esferaId::get_n() const
{
    return w;
}
inline GridCoord esferaId::get_z() const
{
    return z;
}

// obtener_id debe estar antes que agregar_esfera y sacar_esfera, para que se pueda satisfacer que sea inline.
inline esferaId Grillado::obtener_id(Coord x, Coord y, Coord z) const
{
    GridCoord temp[3];
    calcular_coord(x, y, z, temp[0], temp[1], temp[2]);
    return esferaId(temp[0], temp[1], temp[2]);
}

inline esferaId Grillado::agregar_esfera(Coord x, Coord y, Coord Z)
{
    esferaId result(obtener_id(x, y, Z));
    agregar_esfera(result);
    return result;
}

inline void Grillado::sacar_esfera(Coord x, Coord y, Coord Z) throw(gridPositionException)
{
    sacar_esfera(obtener_id(x, y, Z));
}

inline esferaId::esferaId(GridCoord input_v, GridCoord input_w, GridCoord input_z) : v(input_v), w(input_w), z(input_z) {}

inline esferaId::esferaId(const esferaId& s) : v(s.get_m()), w(s.get_n()), z(s.get_z()) {}

#endif // GRILLADO_H
