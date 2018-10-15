#include <vector>
#include <random>

#include "Parameters.h"
#include "agents.h"

using namespace std;

int Cell::totalAliveCount;

template<> const int _Cell<EC>::cellTypeID(0);
template<> const int _Cell<pmn>::cellTypeID(1);
template<> const int _Cell<mono>::cellTypeID(2);
template<> const int _Cell<TH0>::cellTypeID(3);
template<> const int _Cell<TH1>::cellTypeID(4);
template<> const int _Cell<TH2>::cellTypeID(5);
template<> const int _Cell<pmn_marrow>::cellTypeID(6);
template<> const int _Cell<mono_marrow>::cellTypeID(7);
template<> const int _Cell<TH0_germ>::cellTypeID(8);
template<> const int _Cell<TH1_germ>::cellTypeID(9);
template<> const int _Cell<TH2_germ>::cellTypeID(10);

Cell::Cell(int x, int y) {
    xLoc = x;
    yLoc = y;
    gp = Grid::getInstance()(xLoc,yLoc);
    gp->numCells++;
    totalAliveCount++;
}

EC::EC(int x, int y, int iid) : Cell(x,y) {
    gp->numCells--; // Undo the numCells increment because ECs don't count toward the total
    id = iid;
    ec_activation = 0;
    ec_roll = 0;
    ec_stick = 0;
    ec_migrate = 0;
}

pmn::pmn(int x, int y, int age) : Cell(x,y), MortalCell(age) {
    wbc_roll = 1;
    wbc_stick = 0;
    wbc_migrate = 0;
    pmn_pcd = 10;
}

mono::mono(int x, int y, int age, int iTNFr, int iIL1r) : Cell(x,y), MortalCell(age) {
    TNFr = iTNFr;
    IL_1r = iIL1r;
    specific_c["TNFr"] = 0;
    specific_c["IL1r"] = 0;
    activation = 0;
    wbc_roll = 1;
    wbc_stick = 0;
    wbc_migrate = 0;
}

TH0::TH0(int x, int y, int age) : Cell(x,y), MortalCell(age) {
    activation = 0;
    proTH1 = 0;
    proTH2 = 0;
    rTH1 = 0;
    rTH2 = 0;
}

TH1::TH1(int x, int y, int age) : Cell(x,y), MortalCell(age) {}

TH2::TH2(int x, int y, int age) : Cell(x,y), MortalCell(age) {}

pmn_marrow::pmn_marrow(int x, int y) : Cell(x,y) {}

mono_marrow::mono_marrow(int x, int y) : Cell(x,y) {}

TH0_germ::TH0_germ(int x, int y) : Cell(x,y) {}

TH1_germ::TH1_germ(int x, int y) : Cell(x,y) {}

TH2_germ::TH2_germ(int x, int y) : Cell(x,y) {}
