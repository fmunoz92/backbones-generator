// g++ -DMILI_NAMESPACE -lprot-filer test.cpp && ./a.out

#include <iostream>
#include "prot-filer/format_filer.h"
#include "prot-filer/cached_reader.h"

using namespace std;
using namespace prot_filer;

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        cerr << "Wrong number of arguments.\nUse read_xtc <xtc_file>.\n";
    }
    else
    {
        string file = string(argv[1]);
        Coord3DReader* reader = Coord3DReaderFactory::get_instance()->create("xtc");
        reader->open(file);
        CachedReader<FullCache, Coord3DReader, Protein> cached_reader(reader);

        Protein* protein = NULL;

        unsigned int i = 0;
        while ((protein = cached_reader.read(i)) != NULL)
        {
            for (unsigned int c = 0; c < protein->items(); ++c)
            {
                Coord3d cr = (*protein)[c];
                cout << cr.x << " " << cr.y << " " << cr.z << endl;
            }
            cout << "--------------------------------------" << endl;
            ++i;
        }
        cout << "residue size: " << cached_reader.get_reader().get_atom_number() / 3 << endl;
        cout << "size: " << i << endl;
        cached_reader.close();
        Coord3DReaderFactory::destroy_instance();
    }
    return 0;
}

