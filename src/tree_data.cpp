#include "tree_data.h"

TreeData::TreeData(int nRes, size_t cols, size_t rows, size_t depth, std::istream& input_file) :
    nres(nRes),
    // Maximun gyration radius and maximun CA-CA distance.
    // Both equations constructed from database analisys.
    rgmax(2.72 * mili::cubic_root(nres) + 5.0),
    dmax2(mili::square(8.0 * mili::cubic_root(nres) + 25.0)),
    atm(Atoms(nres * 3)),
    cont(0),
    hubo_algun_exito(false),
    grilla(new Grillado(cols, rows, depth)),
    angles_mapping(new prot_filer::AnglesMapping(nres)),
    angles_data(new prot_filer::AnglesData(nres, angles_mapping.get())),
    fragment_ids(1000),
    output_file("traj.xtc")
{
    readdata(input_file);
}

void TreeData::readdata(std::istream& filer)
{
    float fi = 0.0f;
    float si = 0.0f;
    unsigned int i = 0;

    while (filer.good())
    {
        if (filer >> fi && filer >> si)
        {
            cosfi.push_back(cos(mili::deg2rad(fi)));
            sinfi.push_back(sin(mili::deg2rad(fi)));
            cossi.push_back(cos(mili::deg2rad(si)));
            sinsi.push_back(sin(mili::deg2rad(si)));
            angles_mapping->set_mapping(fi, si);
        }
        ++i;
    }
    // Nota a futuro: se deberia lanzar una Excepcion si el formato del archivo fuera equivocado.
}
