#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>
#include <map>
#include <vector>
#include <array>
#include <list>
#include <random> // mt19937

#include "../SepsisIOmacros.h"

using namespace std;

// Cytokine enums
enum DiffusibleCytokine { // Cytokines (and receptors)
    cytotox = 0, // O2rads and enzymes
    endotoxin,
    PAF,
    TNF,
    sTNFr,
    IL1,
    sIL1r,
    IL1ra,
    IFNg,
    IL4,
    IL8,
    IL10,
    IL12,
    GCSF,
    IL6,
    sIL6r,
    NUM_DIFFUSIBLE_CYTOKINES
};
enum NonDiffusibleCytokine { // Cytokine (receptors) specific to mono
    TNFr = NUM_DIFFUSIBLE_CYTOKINES,
    IL1r,
    NUM_CYTOKINES
};
enum NonCytokine { // Non-cytokines
    oxy = NUM_CYTOKINES,
    infection,
    constant,
    NUM_COEFFICIENTS
};

// Behavior mode enums

// Boundary_Modes control how movement at the edges is handled.
// 		Affects: Cell movement, "sniffing," GridPoint neighbors, diffusion, infection/oxy spreading.
enum Boundary_Modes {
	TOROIDAL = 0, 	// Edges wrap around. Default behavior.
	HARD_WALLS		// Edges are hard walls. Cytokines diffuse around and Cells cannot move/see through.
};

// Neighborhood_Modes control how neighborhoods are defined.
// 		Affects: Cell movement, "sniffing," GridPoint neighbors, diffusion, infection/oxy spreading.
enum Neighborhood_Modes {
	MOORE = 0, 		// Up, down, left, right, diagonals. Neighborhood size = 4*neighborhood_length*(neighborhood_length + 1).
	VON_NEUMANN 	// Up, down, left right. Neighborhood size = 4*neighborhood_length.
};

// Interaction_Length_Modes control how far objects can interact.
// 		Affects: GridPoint neighbors, diffusion, infection/oxy spreading. Does NOT affect Cell movement, "sniffing."
enum Interaction_Length_Modes {
	ONE_STEP = 0, 	// Single-step in the Moore neighborhood. Neighborhood size = 8 (Moore) or 4 (von Neumann).
	TWO_STEPS 		// Two steps in the neighborhood. Neighborhood size = 24 (Moore) or 12 (von Neumann).
};

// Mechanism mode parameters
extern bool new_pmn;
extern bool break_ties;
extern bool big_shuffle;

// Mechanism parameters
extern float proTH2_coefficient;
extern float proTH2_sIL6r_coefficient;

// Behavior mode parameters
extern int boundary_mode;
extern int neighborhood_mode;
extern int interaction_length_mode;

// Physiological parameters
extern int inj_number;
extern int numInfectRepeat;
extern float oxyHeal;
extern int infectSpread;
extern int numRecurInj;
extern int seed;

// Experimental parameters
extern int numTimeSteps; //5760=28 days
extern int validation;
extern int injuryStep;
extern float antibioticMultiplier;
extern float endotoxinInfusionStop;

// Structural parameters
extern int xDim;
extern int yDim;
extern int cellCapacity;

// Configuration parameters
extern int verbosity;
extern vector<string> coefficientNames; // Names of coefficients.
extern map<string, int> coefficientMap; // Map from coefficient name to its index/position along a row or column.

// LPS injection parameters
extern float bolusSize; // Amount to inject at time 0
extern float infusionSize; // Total amount to inject from time 1 to time endotoxinInfusionStop
extern float infusionSizeMultiplier; // Fold increase in infusionSize compared to bolusSize. Overwrites infusionSize (to value bolusSize*infusionSizeMultiplier) if infusionSizeMultiplier is positive and infusionSize == 0.
extern float meteringRate; // Exponential rate constant used in metering function

// Debug parameters
extern float mGCSF,mPAF,mTNF,mSTNFR,mIL1,mSIL1R,mIL1RA,mIFNg,mIL4,mIL8,mIL10,mIL12,mCytotox,mEndotoxin,mIL6,mSIL6R;
extern float TNF_scale;

// Cytokine matrix parameters
extern map<string, array<array<float, NUM_COEFFICIENTS>, NUM_COEFFICIENTS>> coefficientMatrix; // Map from rule name to square matrix

// Functions
inline string retname(string name) { return name; }
void setDefaultMatrixParameters();
void setMatrixToZeros();
void initializeCoefficientMap();

//Parameters:
//xDim,yDim => The dimensions of the grid in units of endothelial cells;
//	the model has been calibrated for 101x101
//antibioticMultiplier=>infection variable in each cell is multiplied by this factor each
//	time antibiotics are applied; set to <0 for no antibiotic application.
//numTimeSteps => the number of time steps per simulation
//inj_number => size of initial infectious injury; inj)numberMax is the maximum size of the
//			injury when performing a parameter sweep
//oxyheal => the amount of oxygen healing when 30<oxy<60
//seed => the random number seed
//injuryStep => the amount of steps between recurring small injuries.  Each time step is 7 mins
//infectSpread => measure of how much infection an infected cell spreads to neighboring cells
//numInfectRepeat => how many neighboring cells an infected cell will infect
//cellCapacity => the maximum number of immune cells occupying the same space as an endothelial cell
#endif


