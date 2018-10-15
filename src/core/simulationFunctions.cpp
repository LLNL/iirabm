#include <vector>
#include <iostream>
#include <random>

#include "Parameters.h"
#include "agents.h"
#include "simulationFunctions.h"

using namespace std;

// Cell data
vector<EC> ecArray;
list<pmn> pmnArray;
list<mono> monoArray;
list<TH0> TH0array;
list<TH1> TH1array;
list<TH2> TH2array;
list<pmn_marrow> pmn_marrowArray;
list<mono_marrow> mono_marrowArray;
list<TH0_germ> TH0_germArray;
list<TH1_germ> TH1_germArray;
list<TH2_germ> TH2_germArray;

// Aggregate data
float timeperiod, system_oxy,oxyDeficit,totalInfection,total_TNF,total_sTNFr,total_IL10,
      total_IL6,total_sIL6r,total_GCSF,total_proTH1,total_proTH2,total_IFNg,total_PAF,
      total_IL1,total_IL4,total_IL8,total_IL12,total_sIL1r,total_IL1ra,total_cytotox,total_endotoxin;

// Control data
float multipliers[NUM_DIFFUSIBLE_CYTOKINES];

// Random number generator data
mt19937 generator;
uniform_int_distribution<int> distribution_xDim(0,xDim - 1);
uniform_int_distribution<int> distribution_yDim(0,yDim - 1);
uniform_int_distribution<int> distribution10k(0,9999);
uniform_int_distribution<int> distribution1000(0,999);
uniform_int_distribution<int> distribution100(0,99);
uniform_int_distribution<int> distribution50(0,49);
uniform_int_distribution<int> distribution10(0,9);
uniform_int_distribution<int> distribution9(0,8);
uniform_int_distribution<int> distribution8(0,7);
uniform_int_distribution<int> distribution3(0,2);
uniform_int_distribution<int> distribution2(0,1);

void stepAllCells() {
    // Create a concatenated array
    vector<Cell*> cells;
    for(list<TH0>::iterator it = TH0array.begin(); it != TH0array.end(); ++it) {
        cells.push_back(&(*it));
    }
    for(vector<EC>::iterator it = ecArray.begin(); it != ecArray.end(); ++it) {
        cells.push_back(&(*it));
    }
    for(list<pmn>::iterator it = pmnArray.begin(); it != pmnArray.end(); ++it) {
        cells.push_back(&(*it));
    }
    for(list<mono>::iterator it = monoArray.begin(); it != monoArray.end(); ++it) {
        cells.push_back(&(*it));
    }
    for(list<TH1>::iterator it = TH1array.begin(); it != TH1array.end(); ++it) {
        cells.push_back(&(*it));
    }
    for(list<TH2>::iterator it = TH2array.begin(); it != TH2array.end(); ++it) {
        cells.push_back(&(*it));
    }
    for(list<pmn_marrow>::iterator it = pmn_marrowArray.begin(); it != pmn_marrowArray.end(); ++it) {
        cells.push_back(&(*it));
    }
    for(list<mono_marrow>::iterator it = mono_marrowArray.begin(); it != mono_marrowArray.end(); ++it) {
        cells.push_back(&(*it));
    }
    for(list<TH1_germ>::iterator it = TH1_germArray.begin(); it != TH1_germArray.end(); ++it) {
        cells.push_back(&(*it));
    }
    for(list<TH2_germ>::iterator it = TH2_germArray.begin(); it != TH2_germArray.end(); ++it) {
        cells.push_back(&(*it));
    }
    for(list<TH0_germ>::iterator it = TH0_germArray.begin(); it != TH0_germArray.end(); ++it) {
        cells.push_back(&(*it));
    }
    
    // Shuffle the array
    shuffle(cells.begin(), cells.end(), generator);
    
    // Step all Cells
    for(vector<Cell*>::iterator it = cells.begin(); it != cells.end(); ++it) {
        (*it)->step();
    }
    
    // Erase dead MortalCells
    erase_dead<TH0>(TH0array);
    erase_dead<pmn>(pmnArray);
    erase_dead<mono>(monoArray);
    erase_dead<TH1>(TH1array);
    erase_dead<TH2>(TH2array);
}

void diffuse() {	
	//float nFactor = 1.0/8.0; // Replaced with 1/gp->neighbors.size(), which is variable with hard boundaries.
    
    // Reset all tempVals to 0
    for(int y = 0; y < yDim; y++) {
        for(int x = 0; x < xDim; x++) {
            Cytokine (&c)[NUM_DIFFUSIBLE_CYTOKINES] = Grid::getInstance()(x,y)->c;
            for(int i = 0; i < NUM_DIFFUSIBLE_CYTOKINES; i++) {
                c[i].tempVal = 0;
            }
        }
    }
    
    // Accumulate tempVals for all Cells
    for(int y = 0; y < yDim; y++) {
        for(int x = 0; x < xDim; x++) {
            GridPoint* gp = Grid::getInstance()(x,y);
            for(int i = 0; i < NUM_DIFFUSIBLE_CYTOKINES; i++) {
                float amountToGive = gp->c[i].val * gp->c[i].diffuseFactor / gp->neighbors.size();
                for(int n = 0; n < gp->neighbors.size(); n++) {
                    gp->neighbors[n]->c[i].tempVal += amountToGive;
                }
            }
        }
    }
    
    // Adjust cytokine levels for all Cells
    for(int y = 0; y < yDim; y++) {
        for(int x = 0; x < xDim; x++) {
            Cytokine (&c)[NUM_DIFFUSIBLE_CYTOKINES] = Grid::getInstance()(x,y)->c;
            for(int i = 0; i < NUM_DIFFUSIBLE_CYTOKINES; i++) {
                c[i].val -= c[i].diffuseFactor * c[i].val;
                c[i].val += c[i].tempVal;
            }
        }
    }
}

void evaporate() {
    //bool print = false;
    for(int y = 0; y < yDim; y++) {
        for(int x = 0; x < xDim; x++) {
            Cytokine (&c)[NUM_DIFFUSIBLE_CYTOKINES] = Grid::getInstance()(x,y)->c;
            for(int i = 0; i < NUM_DIFFUSIBLE_CYTOKINES; i++) {
                c[i].val *= c[i].evapFactor;
                
                //if(i==endotoxin && x==30 && y==30 && c[endotoxin].val>0 && c[endotoxin].val < 0.01) cout << "Flooring to 0: " << c[endotoxin].val << "\n";
                //if(i==endotoxin && !print && x==30 && y==30 && c[endotoxin].val==0) print = true;
                //if(i==endotoxin && print && c[endotoxin].val>0.01) cout << "Endotoxin("<<x<<","<<y<<"): " << c[endotoxin].val << "\n";
                
                if(i == TNF && c[i].val/TNF_scale < 0.01)
                    c[i].val = 0;
                
                if(c[i].val < 0.01)
                    c[i].val = 0;
            }
        }
    }
}

void applyAntibiotics() {
    for(int y = 0; y < yDim; y++) {
        for(int x = 0; x < xDim; x++) {
            Grid::getInstance()(x,y)->infection *= antibioticMultiplier;
        }
    }
}

void updateSystemOxy() {
	system_oxy=0;
	oxyDeficit=0;
	totalInfection=0;
	total_TNF=0;
	total_sTNFr=0;
	total_IL10=0;
	total_GCSF=0;
	total_proTH1=0;
	total_proTH2=0;
	total_IFNg=0;
	total_PAF=0;
	total_IL1=0;
	total_IL4=0;
	total_IL8=0;
	total_IL12=0;
	total_sIL1r=0;
	total_IL1ra=0;
    total_cytotox=0;
    total_endotoxin=0;
    total_IL6=0;
    total_sIL6r=0;
    for(int y = 0; y < yDim; y++) {
        for(int x = 0; x < xDim; x++) {
            GridPoint* gp = Grid::getInstance()(x,y);
        
            system_oxy+=(gp->oxy/100);
            totalInfection+=(gp->infection/100);
            total_TNF+=(gp->c[TNF].val/100);
            if(gp->c[TNF].val>mTNF){mTNF=gp->c[TNF].val;}
            total_sTNFr+=(gp->c[sTNFr].val/100);
            if(gp->c[sTNFr].val>mSTNFR){mSTNFR=gp->c[sTNFr].val;}
            total_IL10+=(gp->c[IL10].val/100);
            if(gp->c[IL10].val>mIL10){mIL10=gp->c[IL10].val;}
            total_GCSF+=(gp->c[GCSF].val/100);
            if(gp->c[GCSF].val>mGCSF){mGCSF=gp->c[GCSF].val;}
            total_IFNg+=(gp->c[IFNg].val/100);
            if(gp->c[IFNg].val>mIFNg){mIFNg=gp->c[IFNg].val;}
            total_PAF+=(gp->c[PAF].val/100);
            if(gp->c[PAF].val>mPAF){mPAF=gp->c[PAF].val;}
            total_IL1+=(gp->c[IL1].val/100);
            if(gp->c[IL1].val>mIL1){mIL1=gp->c[IL1].val;}
            total_IL4+=(gp->c[IL4].val/100);
            if(gp->c[IL4].val>mIL4){mIL4=gp->c[IL4].val;}
            total_IL8+=(gp->c[IL8].val/100);
            if(gp->c[IL8].val>mIL8){mIL8=gp->c[IL8].val;}
            total_IL12+=(gp->c[IL12].val/100);
            if(gp->c[IL12].val>mIL12){mIL12=gp->c[IL12].val;}
            total_sIL1r+=(gp->c[sIL1r].val/100);
            if(gp->c[sIL1r].val>mSIL1R){mSIL1R=gp->c[sIL1r].val;}
            total_IL1ra+=(gp->c[IL1ra].val/100);
            if(gp->c[IL1ra].val>mIL1RA){mIL1RA=gp->c[IL1ra].val;}
            total_cytotox+=(gp->c[cytotox].val/100);
            if(gp->c[cytotox].val>mCytotox){mCytotox=gp->c[cytotox].val;}
            total_endotoxin+=(gp->c[endotoxin].val/100);
            if(gp->c[endotoxin].val>mEndotoxin){mEndotoxin=gp->c[endotoxin].val;}
            total_IL6+=(gp->c[IL6].val/100);
            if(gp->c[IL6].val>mIL6){mIL6=gp->c[IL6].val;}
            total_sIL6r+=(gp->c[sIL6r].val/100);
            if(gp->c[sIL6r].val>mSIL6R){mSIL6R=gp->c[sIL6r].val;}
        }
	}
    for(list<TH0>::iterator it = TH0array.begin(); it != TH0array.end(); ++it) {
		total_proTH1+=((*it).proTH1/100);
		total_proTH2+=((*it).proTH2/100);
	}
	oxyDeficit=(xDim*yDim)-system_oxy;
}

void adjustOrientation(int* orientation, int leftOrRight) {
//if leftOrRight=-1, adjust orientation left, if =1, adjust orientation right
	int tempOrient;
	tempOrient=*orientation;
	tempOrient+=leftOrRight;
	if(tempOrient>7){tempOrient=0;}
	if(tempOrient<0){tempOrient=7;}
	*orientation=tempOrient;
}

void getAhead(int orientation, int x, int y, int *xl, int *xm, int *xr, int *yl, int *ym, int *yr) {
	int txl,txm,txr,tyl,tym,tyr;

	if(orientation == 0){
		txl=x-1;
		txm=x;
		txr=x+1;
		tyl=y+1;
		tym=y+1;
		tyr=y+1;
	}
	
	if(orientation == 1){
		txl=x;
		txm=x+1;
		txr=x+1;
		tyl=y+1;
		tym=y+1;
		tyr=y;
	}
	
	if(orientation == 2){
		txl=x+1;
		txm=x+1;
		txr=x+1;
		tyl=y+1;
		tym=y;
		tyr=y-1;
	}
	
	if(orientation == 3){
		txl=x+1;
		txm=x+1;
		txr=x;
		tyl=y;
		tym=y-1;
		tyr=y-1;
	}
	
	if(orientation == 4){
		txl=x+1;
		txm=x;
		txr=x-1;	
		tyl=y-1;
		tym=y-1;
		tyr=y-1;
	}
	
	if(orientation == 5){
		txl=x;
		txm=x-1;
		txr=x-1;
		tyl=y-1;
		tym=y-1;
		tyr=y;
	}
	
	if(orientation == 6){
		txl=x-1;
		txm=x-1;
		txr=x-1;
		tyl=y-1;
		tym=y;
		tyr=y+1;
	}
	
	if(orientation == 7){
		txl=x-1;
		txm=x-1;
		txr=x;
		tyl=y;
		tym=y+1;
		tyr=y+1;
	}
	
	if(boundary_mode == HARD_WALLS) {
		// Value of -1 means that GridPoint will not be "seen."
		if(txl < 0 || txl >= xDim) {txl=-1;}
		if(txm < 0 || txm >= xDim) {txm=-1;}
		if(txr < 0 || txr >= xDim) {txr=-1;}
		if(tyl < 0 || tyl >= yDim) {tyl=-1;}
		if(tym < 0 || tym >= yDim) {tym=-1;}
		if(tyr < 0 || tyr >= yDim) {tyr=-1;}
	}
	else { // Any toroidal mode
		if(txl < 0)		{txl += xDim;}
		if(txl >= xDim)	{txl -= xDim;}
		if(txm < 0)		{txm += xDim;}
		if(txm >= xDim)	{txm -= xDim;}
		if(txr < 0)		{txr += xDim;}
		if(txr >= xDim)	{txr -= xDim;}
		if(tyl < 0)		{tyl += yDim;}
		if(tyl >= yDim)	{tyl -= yDim;}
		if(tym < 0)		{tym += yDim;}
		if(tym >= yDim)	{tym -= yDim;}
		if(tyr < 0)		{tyr += yDim;}
		if(tyr >= yDim)	{tyr -= yDim;}
	}

	*xl=txl;
	*xm=txm;
	*xr=txr;
	*yl=tyl;
	*ym=tym;
	*yr=tyr;
}

void injure_infectionFRD(int inj_number) { //Fixed Radius Disk
    int radius = inj_number;
    if(radius == 0) return; // Special case for radius = 0. Otherwise, 1 Grid Point will be infected due to the "<=" below.
    for(int y = 0; y < yDim; y++) {
        for(int x = 0; x < xDim; x++) {
            if(pow(x-xDim/2,2) + pow(y-yDim/2,2) <= pow(radius,2)) {
                Grid::getInstance()(x,y)->infection = 100;
            }
        }
    }
}

void recur_injury() {
	int x = distribution_xDim(generator);
	int y = distribution_yDim(generator);
    Grid::getInstance()(x,y)->infection=100;
}

// Replace all rules
void addAllRules() {
    EC::rules.clear();
    EC::addRule("EC_activation");
    EC::addRule("EC_midhealthy");
    EC::addRule("EC_unhealthy");
    
    pmn::rules.clear();
    pmn::addRule("pmn_primed");
    pmn::addRule("pmn_burst");
    
    mono::rules.clear();
    mono::addRule("mono_function");
    mono::addRule("mono_activation");
    mono::addRule("mono_unactivated");
    
    TH1::rules.clear();
    TH1::addRule("TH1_threshold");
    
    TH2::rules.clear();
    TH2::addRule("TH2_threshold");
}

// Note: all parameters (namely xDim and yDim) EXCEPT the cytokine update matrix are expected to be set BEFORE reaching this function.
// 		The cytokine update matrix (if not default) is expected to be set/edited AFTER reaching this function.
void initialize() {
	int i,j,k,xTemp,yTemp;
    
    // Reset all Cell counts
    Cell::totalAliveCount = 0;
    EC::aliveCount = 0;
    pmn::aliveCount = 0;
    mono::aliveCount = 0;
    TH0::aliveCount = 0;
    TH1::aliveCount = 0;
    TH2::aliveCount = 0;
    pmn_marrow::aliveCount = 0;
    mono_marrow::aliveCount = 0;
    TH0_germ::aliveCount = 0;
    TH1_germ::aliveCount = 0;
    TH2_germ::aliveCount = 0;

    // Reset the xDim/yDim random number generators in case xDim/yDim have changed.
    distribution_xDim = uniform_int_distribution<int>(0,xDim - 1);
	distribution_yDim = uniform_int_distribution<int>(0,yDim - 1);
    
    // Initializes Grid to size (xDim, yDim) and initializes GridPoint values and neighbors.
    Grid::getInstance().init();
    
    // petersen33: This is now called only at module import
    //initializeCoefficientMap(); // Initialize map from coefficient names to position along rows/columns.
    
    // petersen33: This is now called only at module import, or when resetMatrix() is called.
    // setDefaultMatrixParameters(); // Start with default matrix parameters.

    // petersen33: This is now called only at module import, or when setMatrix() or editMatrix() are called.
    // addAllRules(); // Actually add the rules to each Cell type.
    
    // Reset all Cell arrays
	ecArray.clear();
	pmnArray.clear();
	monoArray.clear();
	TH0array.clear();
	TH1array.clear();
	TH2array.clear();
	pmn_marrowArray.clear();
	mono_marrowArray.clear();
	TH0_germArray.clear();
	TH1_germArray.clear();
	TH2_germArray.clear();

    //petersen33: Are these needed?
	system_oxy=xDim*yDim;
	oxyDeficit=0;
	totalInfection=0;
	total_TNF=0;
	total_sTNFr=0;
	total_IL10=0;
	total_GCSF=0;
	total_proTH1=0;
	total_proTH2=0;
    total_IL6=0;
    total_sIL6r=0;
    
    // Reset all multipliers
    for(int i = 0; i < NUM_DIFFUSIBLE_CYTOKINES; i++) {
        if(i == endotoxin || i == cytotox) {
            multipliers[i] = -1; // endotoxin and cytotox don't have multipliers
        }
        else {
            multipliers[i] = 1;
        }
    }
    
    mGCSF=0;
    mPAF=0;
    mTNF=0;
    mSTNFR=0;
    mIL1=0;
    mSIL1R=0;
    mIL1RA=0;
    mIFNg=0;
    mIL4=0;
    mIL8=0;
    mIL10=0;
    mIL12=0;
	
    // petersen33: Switched order of i and j indices to match the formula "id=y*xDim+xLoc" and the convention N=up (y+1), W=left (x-1), etc.
	k=0; //initialization
	for(j=0;j<yDim;j++){
        for(i=0;i<xDim;i++){       //Initialize EC grid
			ecArray.push_back(EC(i,j,k));
			k++;
		}
	}

	for(i=0;i<500;i++){
		xTemp=distribution_xDim(generator);
		yTemp=distribution_yDim(generator);
		int age = new_pmn ? 2*distribution100(generator) : distribution50(generator);
		pmnArray.push_back(pmn(xTemp,yTemp,age));
	}
	
	for(i=0;i<50;i++){
		xTemp=distribution_xDim(generator);
		yTemp=distribution_yDim(generator);
		monoArray.push_back(mono(xTemp,yTemp));
	}
	
	for(i=0;i<50;i++){
		xTemp=distribution_xDim(generator);
		yTemp=distribution_yDim(generator);
		TH1array.push_back(TH1(xTemp,yTemp));
	}
	
	for(i=0;i<50;i++){
		xTemp=distribution_xDim(generator);
		yTemp=distribution_yDim(generator);
		TH2array.push_back(TH2(xTemp,yTemp));
	}
	
	for(i=0;i<100;i++){
		xTemp=distribution_xDim(generator);
		yTemp=distribution_yDim(generator);
		pmn_marrowArray.push_back(pmn_marrow(xTemp,yTemp));
	}
	
	for(i=0;i<100;i++){
		xTemp=distribution_xDim(generator);
		yTemp=distribution_yDim(generator);
//		cout<<xTemp<<" "<<yTemp<<"\n";
		mono_marrowArray.push_back(mono_marrow(xTemp,yTemp));
	}
		
	for(i=0;i<100;i++){
		xTemp=distribution_xDim(generator);
		yTemp=distribution_yDim(generator);
		TH0_germArray.push_back(TH0_germ(xTemp,yTemp));
	}
	
	for(i=0;i<100;i++){
		xTemp=distribution_xDim(generator);
		yTemp=distribution_yDim(generator);
		TH1_germArray.push_back(TH1_germ(xTemp,yTemp));
	}
	
	for(i=0;i<100;i++){
		xTemp=distribution_xDim(generator);
		yTemp=distribution_yDim(generator);
		TH2_germArray.push_back(TH2_germ(xTemp,yTemp));
	}
}

