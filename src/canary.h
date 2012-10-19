#ifndef CANARY_H
#define CANARY_H

#include <iostream>
#include <sstream>
#include <stdlib.h>

static const int UNASSIGNED_CANARY_VALUE = -1;

template <int line>
class _Canary
{

  char buffer[sizeof(float)*16];
  static size_t depth;
  static int canary_value;
public:
  _Canary()
  {
    if (canary_value == UNASSIGNED_CANARY_VALUE)
    {
      const char* envvar = getenv("CANARY_VALUE");
      if (envvar == NULL)
        canary_value = 0;
      else
      {
        std::stringstream ss(envvar);
        ss >> canary_value;
      }
    }

    for (unsigned i = 0; i < sizeof(buffer); i++)
      buffer[i] = char(canary_value);

    depth++;
  }

  ~_Canary()
  {
    for (size_t i = 0; i < sizeof(buffer); i++)
      if (buffer[i] != char(canary_value))
      {
        std::cerr << "Canary at line " << line << " and depth " << depth << " violated: value is not " << canary_value << "\n";
        exit(EXIT_FAILURE);
      }

    depth--;
  }
};
template <int line> size_t _Canary<line>::depth;
template <int line> int _Canary<line>::canary_value = UNASSIGNED_CANARY_VALUE;

#define Canary  _Canary<__LINE__>

#endif

