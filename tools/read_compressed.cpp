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
        AnglesReader* reader = AnglesReaderFactory::get_instance()->create("compressed");
        reader->open(file);
        CachedReader<FullCache, AnglesReader, IncompleteAnglesData> cached_reader(reader);

        IncompleteAnglesData* angles_data = NULL;

        unsigned int i = 0;
        while ((angles_data = cached_reader.read(i)) != NULL)
        {
            for (unsigned int c = 0; c < angles_data->nres - 1; ++c)
            {
                AngleIdPair p = angles_data->angles[c];
                cout << p.fi << " " << p.si  << endl;
            }
            cout << "--------------------------------------" << endl;
            ++i;
        }
        cout << "residue size: " << cached_reader.get_reader().get_atom_number() / 3 << endl;
        cout << "size: " << i << endl;
        cached_reader.close();
        AnglesReaderFactory::destroy_instance();
    }
    return 0;
}
