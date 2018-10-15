#ifndef simulationFunctions_h
#define simulationFunctions_h
#include <algorithm>

extern normal_distribution<float> infusionDistribution, bolusDistribution;

// Cell data
extern vector<EC> ecArray;
extern list<pmn> pmnArray;
extern list<mono> monoArray;
extern list<TH0> TH0array;
extern list<TH1> TH1array;
extern list<TH2> TH2array;
extern list<pmn_marrow> pmn_marrowArray;
extern list<mono_marrow> mono_marrowArray;
extern list<TH0_germ> TH0_germArray;
extern list<TH1_germ> TH1_germArray;
extern list<TH2_germ> TH2_germArray;

// Aggregate data
extern float timeperiod,system_oxy,oxyDeficit,totalInfection,total_TNF,total_sTNFr,total_IL10,
total_IL6,total_sIL6r,total_GCSF,total_proTH1,total_proTH2,oxyheal,total_IFNg,total_PAF,
total_IL1,total_IL4,total_IL8,total_IL12,total_sIL1r,total_IL1ra,total_cytotox,total_endotoxin;

// Control data
extern float multipliers[NUM_DIFFUSIBLE_CYTOKINES];

// Random number generator data
extern mt19937 generator;
extern uniform_int_distribution<int> distribution_xDim;
extern uniform_int_distribution<int> distribution_yDim;
extern uniform_int_distribution<int> distribution10k;
extern uniform_int_distribution<int> distribution1000;
extern uniform_int_distribution<int> distribution100;
extern uniform_int_distribution<int> distribution50;
extern uniform_int_distribution<int> distribution10;
extern uniform_int_distribution<int> distribution9;
extern uniform_int_distribution<int> distribution8;
extern uniform_int_distribution<int> distribution3;
extern uniform_int_distribution<int> distribution2;

void getAhead(int orient, int x, int y, int *xl, int *xm, int *xr, int *yl, int *ym, int *yr);
void adjustOrientation(int* orientation, int leftOrRight);
void updateSystemOxy();
void diffuse();
void evaporate();
void applyAntibiotics();
void recur_injury();
void addAllRules();
void initialize();
void injure_infectionFRD(int inj_number);
void stepAllCells();
void infuseEndotoxin(bool bolus);

/*template<typename T>
void erase_element(vector<T> &v, const T &item) {
    //v.erase(v.begin() + (&item - &v.front())); // Also works
    v.erase(std::find(v.begin(), v.end(), item));
}*/

template<typename T>
void erase_dead(list<T> &v) {
    v.erase(std::remove_if(v.begin(), v.end(), [](T cell) { return cell.dead; } ), v.end());
}

template<typename T>
void stepCells(list<T>& v) {
    vector<T*> cells;
    for(typename list<T>::iterator it = v.begin(); it != v.end(); ++it)
    {
        cells.push_back(&(*it));
    }
    
    shuffle(cells.begin(), cells.end(), generator);

    for(typename vector<T*>::iterator it = cells.begin(); it != cells.end(); ++it) {
        (*it)->step();
    }
    //v.erase(std::remove_if(v.begin(), v.end(), [](T cell) { return cell.dead; } ), v.end()); // Can't incldue this now because stepCells isn't just called on MortalCells, which are the only ones with the variable "dead."
}

template<class Derived>
void MobileCell<Derived>::wiggle() { //Should always be followed by move() to match NetLogo
    int tempOrient = orientation;
    int dir = distribution3(generator);
    if(dir == 0) { tempOrient--; }
    if(dir == 2) { tempOrient++; }
    if(tempOrient > 7) { tempOrient=0; }
    if(tempOrient < 0) { tempOrient=7; }
    orientation = tempOrient;
}

template<class Derived>
void MobileCell<Derived>::move() {
    int oldx = Cell::xLoc;
    int oldy = Cell::yLoc;
    int newx = oldx;
    int newy = oldy;

    bool wrapped = false; // Whether the move results in wrapping around the Grid.
    
    if(orientation == 0){ //Move North
        newy = oldy + 1;
        if(newy >= yDim) { wrapped = true; newy=0; }
    }
    
    if(orientation == 1){ //Move Northeast
        newy = oldy + 1;
        newx = oldx + 1;
        if(newy >= yDim) { wrapped = true; newy=0; }
        if(newx >= yDim) { wrapped = true; newx=0; }
    }
    
    if(orientation == 2){ //Move East
        newx=oldx+1;
        if(newx>=xDim){wrapped = true; newx=0;}
    }
    
    if(orientation == 3){ //Move Southeast
        newy=oldy-1;
        newx=oldx+1;
        if(newx>=xDim){wrapped = true; newx=0;}
        if(newy<0){wrapped = true; newy=yDim-1;}
    }
    
    if(orientation == 4){ //Move South
        newy=oldy-1;
        if(newy<0){wrapped = true; newy=yDim-1;}
    }
    
    if(orientation == 5){ //Move Southwest
        newy=oldy-1;
        newx=oldx-1;
        if(newy<0){wrapped = true; newy=yDim-1;}
        if(newx<0){wrapped = true; newx=xDim-1;}
    }
    
    if(orientation == 6){ //Move West
        newx=oldx-1;
        if(newx<0){wrapped = true; newx=xDim-1;}
    }
    
    if(orientation == 7){ //Move Northwest
        newx=oldx-1;
        newy=oldy+1;
        if(newx<0){wrapped = true; newx=xDim-1;}
        if(newy>=yDim){wrapped = true; newy=0;}
    }

    // Bail if boundary is hard walls and the move results in wrapping around the Grid.
    if(boundary_mode == HARD_WALLS && wrapped)
        return;
    
    // Actually do the move
    if(Grid::getInstance()(newx,newy)->numCells < cellCapacity){
        Cell::xLoc = newx;
        Cell::yLoc = newy;
        Cell::setGP(newx, newy);
        Grid::getInstance()(oldx,oldy)->numCells--;
        Grid::getInstance()(newx,newy)->numCells++;
        Grid::getInstance()(oldx,oldy)->numCellsOfEachType[this->cellTypeID]--;
        Grid::getInstance()(newx,newy)->numCellsOfEachType[this->cellTypeID]++;
    }
}

#endif /* simulationFunctions_h */
