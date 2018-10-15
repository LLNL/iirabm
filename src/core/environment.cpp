#include "environment.h"

void Cytokine::init() {
    val = 0;
    tempVal = 0;
    diffuseFactor = 0;
    evapFactor = 0;
}

// Initializes everything except neighbors, which must be done after all GridPoints have otherwise been initialized.
void GridPoint::init() {
    oxy = 100;
    infection = 0;
    numCells = 0;
    for(int i = 0; i < 5; i++) {
        numCellsOfEachType[i] = 0;
    }
    
    for(int i = 0; i < NUM_DIFFUSIBLE_CYTOKINES; i++) {
        c[i].init();
    }
    
    c[cytotox].diffuseFactor = 0.4;
    c[endotoxin].diffuseFactor = 1.0;
    c[PAF].diffuseFactor = 0.6;
    c[TNF].diffuseFactor = 0.6;
    c[sTNFr].diffuseFactor = 0.8;
    c[IL1].diffuseFactor = 0.6;
    c[sIL1r].diffuseFactor = 0.8;
    c[IL1ra].diffuseFactor = 0.8;
    c[IFNg].diffuseFactor = 0.8;
    c[IL4].diffuseFactor = 0.8;
    c[IL8].diffuseFactor = 0.6;
    c[IL10].diffuseFactor = 0.8;
    c[IL12].diffuseFactor = 0.8;
    c[GCSF].diffuseFactor = 1.0;
    c[IL6].diffuseFactor = 0.8;
    c[sIL6r].diffuseFactor = 0.8;
    
    c[endotoxin].evapFactor = 0.7;
    c[PAF].evapFactor = 0.7;
    c[cytotox].evapFactor = 0.7;
    c[TNF].evapFactor = 0.8;
    c[IL1].evapFactor = 0.8;
    c[sTNFr].evapFactor = 0.9;
    c[IL1ra].evapFactor = 0.9;
    c[sIL1r].evapFactor = 0.9;
    c[IFNg].evapFactor = 0.8;
    c[IL8].evapFactor = 0.7;
    c[IL10].evapFactor = 0.95;
    c[IL12].evapFactor = 0.8;
    c[IL4].evapFactor = 0.95;
    c[GCSF].evapFactor = 0.95;
    c[IL6].evapFactor = 0.95;
    c[sIL6r].evapFactor = 0.95;
}

void GridPoint::setNeighbors() {

    int numNeighbors;
    int* directions_x;
    int* directions_y;

    // Directions read first in W->E direciton then N->S direction, starting with NW GridPoint.
    if(neighborhood_mode == MOORE) {
        if(interaction_length_mode == ONE_STEP) {
            numNeighbors = 8;

            // Directions: NW, N, NE, W, E, SW, S, SE
            directions_x = new int[8]   { -1,  0,  1,
                                          -1,      1,
                                          -1,  0,  1};

            directions_y = new int[8]   {  1,  1,  1,
                                           0,      0,
                                          -1, -1, -1 };    
        }
        else if(interaction_length_mode == TWO_STEPS) {
            numNeighbors = 24;

            directions_x = new int[24]   { -2, -1,  0,  1,  2,
                                           -2, -1,  0,  1,  2, 
                                           -2, -1,      1,  2,
                                           -2, -1,  0,  1,  2, 
                                           -2, -1,  0,  1,  2 };

            directions_y = new int[24]   {  2,  2,  2,  2,  2,
                                            1,  1,  1,  1,  1, 
                                            0,  0,      0,  0,
                                           -1, -1, -1, -1, -1, 
                                           -2, -2, -2, -2, -2 };
        }
        else
            cout << "Interaction length mode must be ONE_STEP (0) or TWO_STEPS (1)." << endl;
    }
    else if(neighborhood_mode == VON_NEUMANN) {
        if(interaction_length_mode == ONE_STEP) {
            numNeighbors = 4;

            directions_x = new int[4]   {       0,    
                                           -1,      1,
                                                0,    };

            directions_y = new int[4]   {       1,    
                                            0,      0,
                                               -1,    };

        }
        else if(interaction_length_mode == TWO_STEPS) {
            numNeighbors = 12;

            directions_x = new int[12]   {          0,        
                                               -1,  0,  1,
                                           -2, -1,      1,  2,
                                               -1,  0,  1,    
                                                    0         };

            directions_y = new int[12]   {          2,        
                                                1,  1,  1,
                                            0,  0,      0,  0,
                                               -1, -1, -1,    
                                                   -2         };
        }
        else
            cout << "Interaction length mode must be ONE_STEP (0) or TWO_STEPS (1)." << endl;
    }
    else
        cout << "Neighborhood mode must be MOORE (0) or VON_NEUMANN (1)." << endl;
    
    Grid& grid = Grid::getInstance();
    
    neighbors.clear(); // Important so that we don't add more neighbors each time initialize() is called.
    for(int i = 0; i < numNeighbors; i++) {
        int new_x = x + directions_x[i];
        int new_y = y + directions_y[i];
        
        // Skip adding a neighbor if looking beyond the boundary and boundaries are hard walls.
        if(boundary_mode == HARD_WALLS && (new_x < 0 || new_x >= xDim || new_y < 0 || new_y >= yDim))
            continue;

        // Othewrise, the Grid is toroidal.
        if(new_x < 0)       new_x += xDim;
        if(new_x >= xDim)   new_x -= xDim;
        if(new_y < 0)       new_y += yDim;
        if(new_y >= yDim)   new_y -= yDim;
        
        neighbors.push_back(grid(new_x, new_y));
    }
}

// Constructor does nothing. Grid::init() resizes and initializes the field, so that the Grid can be reset after changing xDim/yDim.
Grid::Grid() {}

// Initialize/reset the grid. After this is called, xDim/yDim should not change.
void Grid::init() {
    // Resize field to xDim, yDim
    field.clear();
    field.resize(yDim, vector<GridPoint>(xDim, GridPoint()));

    // Initialize the GridPoints in the field.
    for(int i = 0; i < xDim; i++) {
        for(int j = 0; j < yDim; j++) {
            field[i][j].x = i;
            field[i][j].y = j;
            field[i][j].init();
        }
    }

    // Set each GridPoint's neighbors. This must be done after all have been initialized.
    for(int i = 0; i < xDim; i++) {
        for(int j = 0; j < yDim; j++) {
            field[i][j].setNeighbors();
        }
    }
}

IntroCompartment::IntroCompartment(DiffusibleCytokine id_arg, float rate_arg) : id(id_arg), rate(rate_arg), reserve(0) {}

void IntroCompartment::add(float amount) { reserve += amount; }

void IntroCompartment::meter(int time) {
    float totalAmountToTransfer = reserve*(1 - exp(-rate*(time+1)));
    float amountToTransferPerCell = totalAmountToTransfer/(xDim*yDim);
    
    if(reserve < 0.01) {
        reserve = 0;
        return;
    }
    
    for(int y = 0; y < yDim; y++) {
        for(int x = 0; x < xDim; x++) {
            Grid::getInstance()(x,y)->c[id].val += amountToTransferPerCell;
        }
    }
    
    reserve -= totalAmountToTransfer;
}
