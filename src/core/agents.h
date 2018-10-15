#ifndef AGENTS_H
#define AGENTS_H

#include <vector>
#include <array>
#include <map>
#include <iostream>
#include <random> // mt19937

#include "environment.h" // For Grid, GridPoint

using namespace std;

extern mt19937 generator;
extern uniform_int_distribution<int> distribution8, distribution50;
extern map<string, array<array<float, NUM_COEFFICIENTS>, NUM_COEFFICIENTS>> coefficientMatrix;

typedef pair<int, array<float, NUM_COEFFICIENTS>> Rule;

class Cell {
public:
    GridPoint* gp; // Pointer to the GridPoint on which this Cell sits (which may change)
    void setGP(int x, int y) { gp = Grid::getInstance()(x,y); }
    
    Cell() { cout << "This shouldn't be called!" << endl; }
    Cell(int x, int y);
    int xLoc;
    int yLoc;
    map<string, float> specific_c;
    static int totalAliveCount;
    virtual void step() = 0;
    inline bool operator==(const Cell& other) { return (this == &other); }
};

template<class Derived>
class _Cell : virtual public Cell {
public:
    static map<string, vector<Rule>> rules;
    static int aliveCount;  // Cells of this type currently alive
    static const int cellTypeID; // Arbitrary ID for Cells of this type
    
    _Cell() {
        Derived::aliveCount++;
        gp->numCellsOfEachType[Derived::cellTypeID]++;
    }
    ~_Cell(){}
    void update(string name, bool sequential);
    static void addRule(string ruleName);
};

template<class Derived>
map<string, vector<Rule>> _Cell<Derived>::rules;

template<class Derived>
int _Cell<Derived>::aliveCount(0);

// Populates _Cell<Derived>.rules based on coefficientMatrix.
// This should be called any time the (global) coefficientMatrix is changed.
// Note it doesn't create a Rule if all elements are zero.
template<class Derived>
void _Cell<Derived>::addRule(string ruleName) {
    vector<Rule> updates;
    // For each row (update)
    for(int i = 0; i < NUM_COEFFICIENTS; i++) {
        // Check if the row has any non-zero elements
        bool foundNonZeroEntry = false;
        for(int j = 0; j < NUM_COEFFICIENTS; j++) {
            if(coefficientMatrix[ruleName][i][j] != 0) {
                foundNonZeroEntry = true;
                break;
            }
        }
        // If a non-zero element exists in this row, create a Rule.
        if(foundNonZeroEntry)
            updates.push_back(Rule(i, coefficientMatrix[ruleName][i]));
    }
    if(!updates.empty())
        rules[ruleName] = updates;
}

template<class Derived>
class MortalCell : virtual public _Cell<Derived> {
public:
    int age;
    bool dead;
    
    MortalCell(int age_arg) {
        dead = false;
        if(age_arg == -1) {
            age = distribution50(generator);
        }
        else {
            age = age_arg;
        }
    }
    
    void die() {
        dead = true;
        this->gp->numCells--;
        this->gp->numCellsOfEachType[Derived::cellTypeID]--;
        Derived::aliveCount--;
        Cell::totalAliveCount--;
    }
};

template<class Derived>
class MobileCell : virtual public _Cell<Derived> {
public:
    int orientation;
    
    MobileCell():Cell() { orientation = distribution8(generator); }
    
    void wiggle();
    void move();
    void sniff(DiffusibleCytokine gradient); // TBD
};

struct EC : public _Cell<EC> {
public:
   	int id;
    int ec_activation;
    int ec_roll;         //rolling
    int ec_stick;        //sticking
    int ec_migrate;      //migration
    EC(int x, int y, int id);
    void step();
    void ECfunction(float oxyHeal);
    void inj_function(int infectSpread, int numInfectRepeat);
    void activate();
    void patch_inj_spread(float oxyHeal);
};

struct pmn : MortalCell<pmn>, MobileCell<pmn> {
public:
    int wbc_roll;        //selectins
    int wbc_stick;       //integrens
    int wbc_migrate;     //diapedesis
    float pmn_pcd;
    pmn(int x, int y, int age = -1);
    void step();
    void pmn_burst();
    void pmn_sniff();
};

struct mono : MortalCell<mono>, MobileCell<mono> {
public:
    int wbc_roll;        //selectins
    int wbc_stick;       //integrens
    int wbc_migrate;     //diapedesis
    float TNFr; // petersen33 mono-specific variable
    float IL_1r; // petersen33 mono-specific variable
    float activation;
    mono(int x, int y, int age = -1, int iTNFr = 0, int iIL1r = 0);
    void step();
    void mono_sniff();
};

struct TH0 : MortalCell<TH0>, MobileCell<TH0> {
public:
    float activation;
    float proTH1;
    float proTH2;
    float rTH1;           //random holder for pro-TH1
    float rTH2;           //random holder for pro-TH2
    TH0(int x, int y, int age = -1);
    void step();
};

struct TH1 : MortalCell<TH1>, MobileCell<TH1> {
public:
    TH1(int x, int y, int age = -1);
    void step();
};

struct TH2 : MortalCell<TH2>, MobileCell<TH2> {
public:
    TH2(int x, int y, int age = -1);
    void step();
};

struct pmn_marrow : public _Cell<pmn_marrow> {
public:
    pmn_marrow(int x, int y);
    void step();
};

struct mono_marrow : public _Cell<mono_marrow> {
public:
    mono_marrow(int x, int y);
    void step();
};

struct TH0_germ : public _Cell<TH0_germ> {
public:
   	TH0_germ(int x, int y);
    void step();
};

struct TH1_germ : public _Cell<TH1_germ> {
public:
   	TH1_germ(int x, int y);
    void step();
};

struct TH2_germ : public _Cell<TH2_germ> {
public:
   	TH2_germ(int x, int y);
    void step();
};

#endif
