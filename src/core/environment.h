#ifndef environment_h
#define environment_h

#include <vector>

#include "Parameters.h" // For enums

using namespace std;

class Cytokine {
public:
    Cytokine() {};
    void init();
    float val;
    float tempVal;
    double diffuseFactor; // Should be static
    double evapFactor; // Should be static
    //double scale; // Should be static
};

class GridPoint {
public:
    GridPoint() {};
    void init();
    int x;
    int y;
    
    float oxy;
    float infection;    //infectious vector

    int numCells;
    int numCellsOfEachType[11]; // Order: EC, pmn, mono, TH0, TH1, TH2, pmn_marrow, mono_marrow, TH0_germ, TH1_germ, TH2_germ
    
    Cytokine c[NUM_DIFFUSIBLE_CYTOKINES];
    
    vector<GridPoint*> neighbors; // Pointers to neighboring GridPoints (can't have vector of references)
    void setNeighbors();
};

class Grid {
public:
    static Grid& getInstance() {
        static Grid instance;
        return instance;
    }
    
    void init(); // Calls GridPoint.init() for all GridPoints. Important so that old values don't persist between episodes.
    vector< vector<GridPoint> > field;
    GridPoint* operator ()(int x, int y) { return &field[x][y]; }
    
    Grid(Grid const&) = delete;
    void operator=(Grid const&) = delete;
    
private:
    Grid();
};

class IntroCompartment {
public:
    DiffusibleCytokine id;
    float rate; // Rate at which material is passed onto Grid
    float reserve;
    
    IntroCompartment(DiffusibleCytokine id_arg, float rate_arg);
//--- added by santiago10
//--- LPS injection 
    IntroCompartment():reserve(0){}
    void set(DiffusibleCytokine id_arg, float rate_arg) { id = id_arg; rate = rate_arg; }
    
    void meter(int time);
    void add(float amount);
};

#endif /* environment_h */
