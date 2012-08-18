#include "backbones-generator/tree_data.h"

TreeData::TreeData(int nRes, size_t cols, size_t rows, size_t depth, std::istream& input_file, std::string& output_file)
    : nres(nRes),
      rgmax(2.72 * mili::cubic_root(nres) + 5.0),              // Maximun gyration radius and maximun CA-CA distance.
      dmax2(mili::square(8.0 * mili::cubic_root(nres) + 25.0)),// Both equations constructed from database analisys.
      atm(Atoms(nres * 3)),
      cont(0),
      hubo_algun_exito(false),
      grilla(cols, rows, depth),
      angles_mapping(nres),
      angles_data(nres, &angles_mapping),
      output_file(output_file)
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

            angles_mapping.set_mapping(fi, si);
        }
        ++i;
    }
    // Nota a futuro: se deberia lanzar una Excepcion si el formato del archivo fuera equivocado.
}
