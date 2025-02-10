#ifndef LATTICE_CONSTANTS_H
#define LATTICE_CONSTANTS_H

#include <vector>
#include <array>
#include <map>
#include <cmath>

extern std::vector<int> ex;
extern std::vector<int> nex;
extern std::vector<int> ey;
extern std::vector<int> ney;
extern std::vector<int> ez;
extern std::vector<int> nez;
extern std::vector<std::vector<int>> drns;
extern std::array <std::array <int,3>, 26> adrns;

extern std::array<int, 3> ax;
extern std::array<int, 3> ay;
extern std::array<int, 3> az;
extern std::array<int, 3> nx;
extern std::array<int, 3> ny;
extern std::array<int, 3> nz;
extern std::array<int, 3> axay;
extern std::array<int, 3> axaz;
extern std::array<int, 3> axny;
extern std::array<int, 3> axnz;
extern std::array<int, 3> nxay;
extern std::array<int, 3> nxaz;
extern std::array<int, 3> nxny;
extern std::array<int, 3> nxnz;
extern std::array<int, 3> ayaz;
extern std::array<int, 3> aynz;
extern std::array<int, 3> nyaz;
extern std::array<int, 3> nynz;
extern std::array<int, 3> axayaz;
extern std::array<int, 3> axaynz;
extern std::array<int, 3> axnyaz;
extern std::array<int, 3> axnynz;
extern std::array<int, 3> nxayaz;
extern std::array<int, 3> nxaynz;
extern std::array<int, 3> nxnyaz;
extern std::array<int, 3> nxnynz;

extern std::map<int, std::array<double, 3>> Or2Dir;
extern std::map<std::array<double, 3>, int> Dir2Or;

#endif // LATTICE_CONSTANTS_H
