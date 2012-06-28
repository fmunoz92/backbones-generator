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

class esferaId
{
public:
    esferaId();
    esferaId(GridCoord input_v, GridCoord input_w, GridCoord input_z);
    esferaId(const esferaId& s);

    GridCoord get_m() const;
    GridCoord get_n() const;
    GridCoord get_z() const;
private:
    GridCoord v;
    GridCoord w;
    GridCoord z;
};

class Grillado
{
public:
    struct gridPositionException : public std::domain_error
    {
        gridPositionException()
            : domain_error("The desired position is already empty") {}
    };

    struct sizeParamException : public std::invalid_argument
    {
        sizeParamException()
            : invalid_argument("M<3 || N<3 || Z<3") {}
    };

    struct radiusOrDistanceParamException : public std::invalid_argument
    {
        radiusOrDistanceParamException()
            : std::invalid_argument("Either 1) sqrt(2)*radio >= dist: there may be diagonal intersections between spheres (spheres which differ in more than one coordinate) or 2) dist>=2*R: there are no intersections between spheres.") {}

    };

    Grillado(size_t a, size_t b, size_t c, Length R = 4.0f, Length D = 5.7f) throw
    (radiusOrDistanceParamException, sizeParamException, std::bad_alloc);

    ~Grillado();

    void info_grillado() const;

    inline esferaId agregar_esfera(Coord x, Coord y, Coord z) ;

    void agregar_esfera(const esferaId& id);

    inline void sacar_esfera(Coord x, Coord y, Coord z) throw(gridPositionException) ;

    void sacar_esfera(const esferaId& id) throw(gridPositionException);

    // Devuelve el volumen parcial.
    inline Volume obtener_vol_parcial() const;
    void reset();
private:
    inline static Volume calc_vol(Length radio);

    inline static Volume calc_int(Length radio, Length dist);

    void reducir_vol_parcial(int coord_x, int coord_y, int coord_z);
    void aumentar_vol_parcial(int coord_x, int coord_y, int coord_z);
    void calcular_coord(Coord coord_x, Coord coord_y, Coord coord_z,
                        GridCoord& x, GridCoord& y, GridCoord& z) const;
    unsigned int calcular_intersecciones(GridCoord coord_x, GridCoord coord_y, GridCoord coord_z) const;
    inline esferaId obtener_id(Coord x, Coord y, Coord z) const;

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
};

#define GRILLADO_INLINE_H
#include "grillado_inline.h"
#undef GRILLADO_INLINE_H

#endif // GRILLADO_H
