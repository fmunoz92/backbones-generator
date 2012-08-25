#ifndef GRILLADO_INLINE_H
#error Internal header file, DO NOT include this.
#endif

// Se calcula el volumen de las intersecciones entre esferas.
// Calculates the volume of the intersection between neighbour spheres.
inline Volume Grillado::calc_int(Length radio, Length dist)
{
    return (2.0f) * PI * (((2.0f) * radio * radio * radio / (3.0f)) - (dist * radio * radio / (2.0f)) + (dist * dist * dist / (24.0f)));
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
        a += m;

    return static_cast<unsigned int>(a);
}
