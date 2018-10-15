#include "Parameters.h"

using namespace std;

// Mechanism mode parameters
bool new_pmn = false; // If true, use new pmn aging mechanism
bool break_ties = true; // If true, ties during chemotactic movement are broken randomly
bool big_shuffle = true; // If true, all Cell types are shuffled together before executing randomly

// Mechanism parameters
float proTH2_coefficient = 1.0;
float proTH2_sIL6r_coefficient = 1.0;

// Behavior mode parameters (cast enums as ints). Values of zero are default.
int boundary_mode           = 0; // TOROIDAL
int neighborhood_mode       = 0; // MOORE
int interaction_length_mode = 0; // ONE_STEP

// Physiological parameters
int inj_number = 37;
int numInfectRepeat = 1;
float oxyHeal = 0.2;
int infectSpread = 10;
int numRecurInj = 2;

// Experimental parameters
int numTimeSteps = 515; //5760 = 28 days
int validation = 0; // If 1 (true), numTimeSteps will be overwritten.
int injuryStep = 205;
float antibioticMultiplier = 0.2;
float endotoxinInfusionStop = 3.0F*60.0F/7.0F; // ~3 hr

// Structural parameters
int xDim = 101;
int yDim = 101;
int cellCapacity = 28;

// Configuration parameters
int verbosity = 0;
vector<string> coefficientNames = { "cytotox", "endotoxin", "PAF", "TNF", "sTNFr", "IL1", "sIL1r", "IL1ra", "IFNg", "IL4", "IL8", "IL10", "IL12", "GCSF", "IL6", "sIL6r", "TNFr", "IL1r", "oxy", "infection", "constant" };
map<string, int> coefficientMap;

// LPS injection parameters
float bolusSize = 0;
float infusionSize = 0;
float infusionSizeMultiplier = 0;

float meteringRate = 0.5; // arbitrary default value

// Debug parameters
float mGCSF,mPAF,mTNF,mSTNFR,mIL1,mSIL1R,mIL1RA,mIFNg,mIL4,mIL8,mIL10,mIL12,mCytotox,mEndotoxin,mIL6,mSIL6R;
float TNF_scale = 1;

// Random number generator data/parameters
int seed = 32;

// Cytokine matrix parameters
// Each ruleset has a square matrix. But _Cell has a vector of Updates instead of a square matrix, so that all-zero rows in the matrix are not included.
map<string, array<array<float, NUM_COEFFICIENTS>, NUM_COEFFICIENTS>> coefficientMatrix;

// Sets all elements of coefficientMatrix to zero
void setMatrixToZeros() {
    const int numRules = 10;
    string ruleNames[numRules] = {"EC_activation", "EC_midhealthy","EC_unhealthy", "pmn_primed", "pmn_burst", "mono_function", "mono_activation", "mono_unactivated", "TH1_threshold", "TH2_threshold"};
    
    // Start by filling in each matrix with zeroes.
    for(int r = 0; r < numRules; r++) {
        for(int i = 0; i < NUM_COEFFICIENTS; i++) {
            for(int j = 0; j < NUM_COEFFICIENTS; j++) {
                coefficientMatrix[ruleNames[r]][i][j] = 0;
            }
        }
    }
}

// TBD: Refactor to setDefaultMatrixCoefficients()
// Set all matrix coefficients to hard-coded default values
void setDefaultMatrixParameters() {
    setMatrixToZeros();
    
    // Fill in the hard-coded values
    // Note: This depends on column ordering. Index by row/column enum (setting a single matrix value at a time) values to prevent this.
    coefficientMatrix["EC_activation"][PAF]     = {{     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1   }};
    coefficientMatrix["EC_activation"][IL8]     = {{     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1   }};
    
    coefficientMatrix["EC_midhealthy"][PAF]     = {{     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1   }};
    
    coefficientMatrix["EC_unhealthy"][PAF]      = {{     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1   }};
    
    coefficientMatrix["pmn_primed"][IL1ra]      = {{     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1   }};
    
    
    coefficientMatrix["pmn_burst"][TNF]         = {{     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1   }};
    coefficientMatrix["pmn_burst"][IL1]         = {{     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1   }};
    
    coefficientMatrix["mono_function"][IL1ra]   = {{     0,0,0,0,0,0.5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0   }};
    coefficientMatrix["mono_function"][sTNFr]   = {{     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.5,0,0,0,0   }};
    coefficientMatrix["mono_function"][sIL1r]   = {{     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.5,0,0,0   }};
    
    coefficientMatrix["mono_activation"][GCSF]  = {{     0,1,1,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0   }};
    coefficientMatrix["mono_activation"][IL8]   = {{     0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0   }};
    coefficientMatrix["mono_activation"][IL12]  = {{     0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0   }};
    coefficientMatrix["mono_activation"][IL10]  = {{     0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0   }};
    coefficientMatrix["mono_activation"][IL1]   = {{     0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0   }};
    coefficientMatrix["mono_activation"][TNF]   = {{     0,1,1,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0   }};
    coefficientMatrix["mono_activation"][IL6]   = {{     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1   }};
    
    coefficientMatrix["mono_unactivated"][IL10] = {{     0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0   }};

    coefficientMatrix["TH1_threshold"][IFNg]    = {{     0,0,0,1,0,1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0   }};
    
    coefficientMatrix["TH2_threshold"][IL4]     = {{     0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0   }};
    coefficientMatrix["TH2_threshold"][IL10]    = {{     0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0   }};
}

void initializeCoefficientMap()
{
    for(int i = 0; i < coefficientNames.size(); i++)
    {
        coefficientMap[coefficientNames[i]] = i;
    }
}
