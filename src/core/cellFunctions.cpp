#include <vector>
#include <iostream>
#include <random>

#include "Parameters.h"
#include "agents.h"
#include "simulationFunctions.h"

using namespace std;

void heal(int index);
float updateWrapper(float val, float multiplier);

float operator * (Cytokine& c, float x) {
    return c.val*x;
} // Linker errors when this was in agents.h

template<class Derived>
void _Cell<Derived>::update(string name, bool sequential) {
    vector<Rule>& updates = rules[name];
    
    int numToUpdate = updates.size();
    
    float preMultiplierValue[numToUpdate];
    
    for(int i = 0; i < numToUpdate; i++) {
        int target = updates[i].first;
        preMultiplierValue[i] = this->gp->c[target].val + std::inner_product(begin(this->gp->c), end(this->gp->c), begin(updates[i].second), updates[i].second.back());
        
        if(this->specific_c.find("TNFr") != this->specific_c.end()) {
            preMultiplierValue[i] += updates[i].second[TNFr]*this->specific_c["TNFr"];
        }
        if(this->specific_c.find("IL1r") != this->specific_c.end()) {
            preMultiplierValue[i] += updates[i].second[IL1r]*this->specific_c["IL1r"];
        }
        
        if(sequential) {
            this->gp->c[target].val = max(updateWrapper(preMultiplierValue[i], multipliers[target]),float(0));
        }
    }
    
    if(!sequential) {
        for(int i = 0; i < numToUpdate; i++) {
            int target = updates[i].first;
            this->gp->c[target].val = max(updateWrapper(preMultiplierValue[i], multipliers[target]),float(0));
        }
    }
}

float updateWrapper(float valueBeforeMultiplier, float multiplier) {
    if(multiplier <= 1.0) {
        return valueBeforeMultiplier*multiplier;
    }
    else {
        return valueBeforeMultiplier + (multiplier - 1.0);
    }
}

void EC::step() {
    inj_function(infectSpread, numInfectRepeat);
    ECfunction(oxyHeal);
}

void EC::ECfunction(float oxyHeal){
    if((gp->c[endotoxin].val >= 1) || (gp->oxy < 60)) {
        ec_activation = 1;
    }
    if(ec_activation == 1) {
        activate();
    }
    patch_inj_spread(oxyHeal);
}

void EC::inj_function(int infectSpread, int numInfectRepeat) {
    gp->oxy = max(float(0), gp->oxy - gp->infection);
    gp->c[endotoxin].val += (gp->infection/10);
    
    uniform_int_distribution<int> distribution_neighbors(0, gp->neighbors.size() - 1);
    for(int i = 1; i <= numInfectRepeat; i++) {
        if(gp->infection >= 100){
            int temp = distribution_neighbors(generator);
            gp->neighbors[temp]->infection += infectSpread;
            gp->infection = 100;
        }
    }
    if(gp->infection>0) {
        gp->infection=max(float(0),gp->infection - gp->c[cytotox].val + float(0.1)); // petersen33: What is the +0.1 for?
    }
}

void EC::activate() {
    ec_roll++;
    ec_stick++;
    
    //gp->c[PAF].val = updateWrapper(gp->c[PAF].val+1, multipliers[PAF]);
    //gp->c[IL8].val = updateWrapper(gp->c[IL8].val+1, multipliers[IL8]);
    update("EC_activation", false);
}

void EC::patch_inj_spread(float oxyHeal) {
    gp->oxy = gp->oxy - gp->c[cytotox].val;
    
    if(gp->oxy >= 60){
        gp->oxy = min(float(100),gp->oxy + oxyHeal);
    }
    
    if((gp->oxy<60)&&(gp->oxy>30)) { //ischemia
        ec_roll++;
        gp->oxy-=0.5;
        //		gp->c[PAF].val = updateWrapper(gp->c[PAF].val + 1, multipliers[PAF]);
        
        update("EC_midhealthy", false);
        
        float amount_to_reduce = 0.05*8/gp->neighbors.size(); // Assuming 0.05 is calibrated to MOORE + ONE_STEP (8 neighbors) and that ischemia spreads uniformly across all neighbors.
        for(int i = 0; i < gp->neighbors.size(); i++) {
            gp->neighbors[i]->oxy -= amount_to_reduce;
            if(gp->neighbors[i]->oxy < 0) { gp->neighbors[i]->oxy = 0; } // added by petersen33
        }
    }
    if(gp->oxy<=30) { //infarction
        ec_stick++;
        gp->oxy -= 2;
        //gp->c[PAF].val = updateWrapper(gp->c[PAF].val+1, multipliers[PAF]);
        update("EC_unhealthy", false);
        float amount_to_reduce = 0.25*8/gp->neighbors.size(); // Assuming 0.05 is calibrated to MOORE + ONE_STEP (8 neighbors) and that ischemia spreads uniformly across all neighbors.
        for(int i = 0; i < gp->neighbors.size(); i++) {
            gp->neighbors[i]->oxy -= amount_to_reduce;
            if(gp->neighbors[i]->oxy<0){gp->neighbors[i]->oxy=0;}
        }
    }
    if(gp->oxy<0) {
        gp->oxy=0;
    }
}

void pmn::step() {
    int id = yLoc*xDim+xLoc;
    int roll = ecArray[id].ec_roll;

    // Convert all IL6 to sIL6r
    gp->c[sIL6r].val += gp->c[IL6].val;
    gp->c[IL6].val = 0.0;

    if(wbc_migrate>0) {
        pmn_burst();
    }
    else{
        if((roll > 3) && (wbc_roll == 1)) {
            pmn_sniff();
        }
        else {
            pmn_sniff();
            pmn_sniff();
        }
        // petersen33: Note there is a difference here. Before, when using ecArray[id].IL1ra, "id" was referring to the old site. Now, since "gp" changes in pmn_sniff(), the following lines adjust a *different* Cell's IL1ra than before.
        // petersen33: This bug can be reintroduced by placing the if(gp->c[TNF].val + gp->c[PAF].val > 1) update function before the calls to pmn_sniff().
        id=yLoc*xDim+xLoc; // petersen33: This is needed here because of the usage of id in ecArray[id].ec_stick >= 100.
        if(gp->c[TNF].val/TNF_scale + gp->c[PAF].val > 1){
            wbc_stick = gp->c[IL1].val;
            //gp->c[IL1ra].val = updateWrapper(gp->c[IL1ra].val + 1, multipliers[IL1ra]);
            update("pmn_primed", false);
        }
        if((wbc_stick >= 1)&&(ecArray[id].ec_stick >= 100)) {
            wbc_migrate = max(float(0),(gp->c[TNF].val/TNF_scale + gp->c[IL1].val + gp->c[GCSF].val - gp->c[IL10].val));
        }
        if(new_pmn) age -= 4;
        else age--;
        if(age < 0) {
            die();
            if(gp->numCells < 0) cout << "Cell Error PMNFUNC " << gp->numCells << endl;
            return;
        }
    }
}

void pmn::pmn_burst() {
    int id = yLoc*xDim+xLoc;
    float tnf, gcsf, ifng;
    gp->c[cytotox].val = max(float(10), gp->c[TNF].val/TNF_scale);
    gp->oxy = 100;
    ecArray[id].ec_roll = 0;
    ecArray[id].ec_stick = 0;
    ecArray[id].ec_migrate = 0;
    //gp->c[TNF].val = updateWrapper(gp->c[TNF].val+1, multipliers[TNF]);
    //gp->c[IL1].val = updateWrapper(gp->c[IL1].val+1, multipliers[IL1]);
    update("pmn_burst", false);
    if(new_pmn) {
        tnf = gp->c[TNF].val/TNF_scale;
        gcsf = gp->c[GCSF].val;
        ifng = gp->c[IFNg].val;
        if(tnf > 3) {
            age -= 5;
        }
        if(tnf < 3 && tnf > 1 && (gcsf > 0 || ifng > 0)) {
            age -= 2;
        }
        if(tnf < 3 && tnf > 1 && gcsf == 0 && ifng == 0) {
            age -= 3;
        }
        if(tnf < 1) {
            age-= 1;
        }   
    }
    else { // Old pmn aging mechanism
        age = pmn_pcd;
        pmn_pcd = pmn_pcd - 1 + max(float(0),(gp->c[TNF].val + gp->c[IFNg].val + gp->c[GCSF].val - gp->c[IL10].val)/100);
    }
    if(age < 0) {
        die();
        return;
    }
}

void pmn::pmn_sniff(){
    int x,y,xl,xm,xr,yl,ym,yr,flag;
    float pmnahead,pmnright,pmnleft;
    x = xLoc;
    y = yLoc;
    
    flag=0;  //Flag=-1 for left, 0 for middle, and 1 for right
    
    getAhead(orientation, x, y, &xl, &xm, &xr, &yl, &ym, &yr);
    
    // Value of -1 means the GridPoint can't be seen, and the sniffed value is set to 0.
    if(xr == -1 || yr == -1) pmnright = 0;
    else pmnright = Grid::getInstance()(xr,yr)->c[IL8].val;
    if(xm == -1 || ym == -1) pmnahead = 0;
    else pmnahead = Grid::getInstance()(xm,ym)->c[IL8].val;
    if(xl == -1 || yl == -1) pmnleft = 0;
    else pmnleft = Grid::getInstance()(xl,yl)->c[IL8].val;
    
    if(break_ties) {
        // Ties are broken randomly.
        
        // Three-way tie
        if(pmnright == pmnahead && pmnright == pmnleft)
            flag = distribution3(generator) - 1; // -1, 0, or 1
        
        // Two-way ties
        else if (pmnleft == pmnahead && pmnleft > pmnright)
            flag = distribution2(generator) - 1; // -1 or 0
        else if (pmnahead == pmnright && pmnahead > pmnleft)
            flag = distribution2(generator); // 0 or 1
        else if (pmnleft == pmnright && pmnleft > pmnahead)
            flag = distribution2(generator) == 0 ? -1 : 1; // -1 or 1
        
        // Now we're sure there are no ties.
        else if(pmnleft > pmnahead && pmnleft > pmnright)
            flag = -1;
        else if(pmnahead > pmnleft && pmnahead > pmnright)
            flag = 0;
        else if(pmnright > pmnleft && pmnright > pmnahead)
            flag = 1;
        else
            cout << "Should not get here!" << endl;
    }
    else {
        // Ties favor right, then left, then middle
        if((pmnright >= pmnahead) && (pmnright >= pmnleft)) flag = 1;
        else if(pmnleft >= pmnahead) flag = -1;
    }
    
    adjustOrientation(&orientation,flag);

    // Bail early if the chosen flag is an unreachable GridPoint. (This can happen due to how ties are handled.)
    if(flag == -1 && (xl == -1 || yl == -1)) return;
    if(flag ==  0 && (xm == -1 || ym == -1)) return;
    if(flag ==  1 && (xr == -1 || yr == -1)) return;
    
    if(flag==-1){
        if(Grid::getInstance()(xl,yl)->numCells < cellCapacity)
        {
            xLoc=xl;
            yLoc=yl;
            Grid::getInstance()(x,y)->numCells--;
            Grid::getInstance()(xl,yl)->numCells++;
            Grid::getInstance()(x,y)->numCellsOfEachType[this->cellTypeID]--;
            Grid::getInstance()(xl,yl)->numCellsOfEachType[this->cellTypeID]++;
        }
    }
    else if(flag==0){
        if(Grid::getInstance()(xm,ym)->numCells < cellCapacity)
        {
            xLoc=xm;
            yLoc=ym;
            Grid::getInstance()(x,y)->numCells--;
            Grid::getInstance()(xm,ym)->numCells++;
            Grid::getInstance()(x,y)->numCellsOfEachType[this->cellTypeID]--;
            Grid::getInstance()(xm,ym)->numCellsOfEachType[this->cellTypeID]++;
        }
    }
    else if(flag==1){
        if(Grid::getInstance()(xr,yr)->numCells < cellCapacity)
        {
            xLoc=xr;
            yLoc=yr;
            Grid::getInstance()(x,y)->numCells--;
            Grid::getInstance()(xr,yr)->numCells++;
            Grid::getInstance()(x,y)->numCellsOfEachType[this->cellTypeID]--;
            Grid::getInstance()(xr,yr)->numCellsOfEachType[this->cellTypeID]++;
        }
    }
    else
        cout << "PMN can't move!" << endl;
    setGP(xLoc,yLoc);
}

void mono::step() {
    int id = yLoc*xDim+xLoc;
    
    if(gp->c[sTNFr].val <= 100) {
        specific_c["TNFr"] = min(float(100),(gp->c[TNF].val/TNF_scale + gp->c[sTNFr].val));
    }
    else {
        specific_c["TNFr"] = min(float(100),max(float(0), gp->c[TNF].val/TNF_scale - gp->c[sTNFr].val));
    }
    specific_c["IL1r"] = min(float(100) , max(float(0), gp->c[IL1].val - gp->c[IL1ra].val - gp->c[sIL1r].val));
    //gp->c[IL1ra].val = updateWrapper(gp->c[IL1ra].val + gp->c[IL1].val/2, multipliers[IL1ra]);
    //gp->c[sTNFr].val = updateWrapper(gp->c[sTNFr].val + specific_c["TNFr"]/2, multipliers[sTNFr]);
    //gp->c[sIL1r].val = updateWrapper(gp->c[sIL1r].val + specific_c["IL1r"]/2, multipliers[sIL1r]);
    update("mono_function", false);
    
    activation = gp->c[endotoxin].val + gp->c[PAF].val + gp->c[IFNg].val - gp->c[IL10].val;
    if(activation > 0) {
        /*gp->c[GCSF].val = updateWrapper(gp->c[GCSF].val + gp->c[endotoxin].val + gp->c[PAF].val + gp->c[TNF].val + gp->c[IFNg].val, multipliers[GCSF]);
         gp->c[IL8].val = updateWrapper(gp->c[IL8].val + gp->c[TNF].val + gp->c[IL1].val, multipliers[IL8]);
         gp->c[IL12].val = updateWrapper(gp->c[IL12].val + gp->c[TNF].val + gp->c[IL1].val, multipliers[IL12]);
         gp->c[IL10].val = updateWrapper(gp->c[IL10].val + gp->c[TNF].val + gp->c[IL1].val, multipliers[IL10]);
         gp->c[IL1].val = updateWrapper(gp->c[IL1].val + gp->c[endotoxin].val + gp->c[PAF].val + specific_c["IL1r"] + gp->c[TNF].val, multipliers[IL1]);
         gp->c[TNF].val = updateWrapper(gp->c[TNF].val + gp->c[endotoxin].val + gp->c[PAF].val + specific_c["TNFr"] + gp->c[IFNg].val, multipliers[TNF]);*/
        update("mono_activation", false);
        if((wbc_stick == 1) && (ecArray[id].ec_stick >= 100)) {
            wbc_migrate = 1;
        }
        if(wbc_roll == 1) {
            wbc_stick = 1;
        }
        wbc_roll=1;
    }
    if(activation<0) {
        //gp->c[IL10].val = updateWrapper(gp->c[TNF].val + gp->c[IL1].val + gp->c[IL10].val, multipliers[IL10]);
        update("mono_unactivated", false);
    }
    if(wbc_migrate == 1) {
        heal(id);
    }
    if(wbc_roll == 1) {
        mono_sniff();
    }
    else{
        mono_sniff();
        mono_sniff();
    }
    age--;
    if(age < 0) {
        die();
        return;
    }
    if(activation > 20) {
        activation = 20;
    }
}

void mono::mono_sniff() {
    int x,y,xl,xm,xr,yl,ym,yr,flag;
    float pafahead,pafright,pafleft;
    x=xLoc;
    y=yLoc;
    
    flag=0;

    getAhead(orientation, x, y, &xl, &xm, &xr, &yl, &ym, &yr);

    // Value of -1 means the GridPoint can't be seen, and the sniffed value is set to 0.    
    if(xr == -1 || yr == -1) pafright = 0;
    else pafright = Grid::getInstance()(xr,yr)->c[PAF].val;
    if(xm == -1 || ym == -1) pafahead = 0;
    else pafahead = Grid::getInstance()(xm,ym)->c[PAF].val;
    if(xl == -1 || yl == -1) pafleft = 0;
    else pafleft = Grid::getInstance()(xl,yl)->c[PAF].val;
    
    if(break_ties) {
        // Ties are broken randomly.
        
        // Three-way tie
        if(pafright == pafahead && pafright == pafleft)
            flag = distribution3(generator) - 1; // -1, 0, or 1
        
        // Two-way ties
        else if (pafleft == pafahead && pafleft > pafright)
            flag = distribution2(generator) - 1; // -1 or 0
        else if (pafahead == pafright && pafahead > pafleft)
            flag = distribution2(generator); // 0 or 1
        else if (pafleft == pafright && pafleft > pafahead)
            flag = distribution2(generator) == 0 ? -1 : 1; // -1 or 1
        
        // Now we're sure there are no ties.
        else if(pafleft > pafahead && pafleft > pafright)
            flag = -1;
        else if(pafahead > pafleft && pafahead > pafright)
            flag = 0;
        else if(pafright > pafleft && pafright > pafahead)
            flag = 1;
        else
            cout << "Should not get here!" << endl;
    }
    else {
        // Ties favor right, then left, then middle
        if((pafright >= pafahead) && (pafright >= pafleft)) flag = 1;
        else if(pafleft >= pafahead) flag = -1;
    }

    adjustOrientation(&orientation,flag);
    
    // Bail early if the chosen flag is an unreachable GridPoint. (This can happen due to how ties are handled.)
    if(flag == -1 && (xl == -1 || yl == -1)) return;
    if(flag ==  0 && (xm == -1 || ym == -1)) return;
    if(flag ==  1 && (xr == -1 || yr == -1)) return;

    if(flag==-1){
        if(Grid::getInstance()(xl,yl)->numCells < cellCapacity) {
            xLoc=xl;
            yLoc=yl;
            Grid::getInstance()(x,y)->numCells--;
            Grid::getInstance()(xl,yl)->numCells++;
            Grid::getInstance()(x,y)->numCellsOfEachType[this->cellTypeID]--;
            Grid::getInstance()(xl,yl)->numCellsOfEachType[this->cellTypeID]++;
        }
    }
    
    else if(flag==0){
        if(Grid::getInstance()(xm,ym)->numCells < cellCapacity) {
            xLoc=xm;
            yLoc=ym;
            Grid::getInstance()(x,y)->numCells--;
            Grid::getInstance()(xm,ym)->numCells++;
            Grid::getInstance()(x,y)->numCellsOfEachType[this->cellTypeID]--;
            Grid::getInstance()(xm,ym)->numCellsOfEachType[this->cellTypeID]++;
        }
    }
    
    else if(flag==1){
        if(Grid::getInstance()(xr,yr)->numCells < cellCapacity) {
            xLoc=xr;
            yLoc=yr;
            Grid::getInstance()(x,y)->numCells--;
            Grid::getInstance()(xr,yr)->numCells++;
            Grid::getInstance()(x,y)->numCellsOfEachType[this->cellTypeID]--;
            Grid::getInstance()(xr,yr)->numCellsOfEachType[this->cellTypeID]++;
        }
    }
    else
        cout << "Mono can't move!" << endl;
    setGP(xLoc,yLoc);
}

void TH0::step() {
    if(gp->c[IL12].val + gp->c[IL4].val + gp->c[sIL6r].val > 0) {
        proTH1 = (gp->c[IL12].val + gp->c[IFNg].val)*100;
        proTH2 = proTH2_coefficient*(gp->c[IL10].val + gp->c[IL4].val + proTH2_sIL6r_coefficient*gp->c[sIL6r].val)*100;
        if((proTH1 > 0) && (proTH2 > 0)) {
            rTH1 = distribution10k(generator) % int(ceil(proTH1));
            rTH2 = distribution10k(generator) % int(ceil(proTH2));
            if(rTH1 > rTH2) {
                activation++;
            }
            if(rTH1 < rTH2) {
                activation--;
            }
        }
        if(proTH1 == 0) {
            activation--;
        }
        if(proTH2 == 0) {
            activation++;
        }
    }
    wiggle();
    move();
    age--;
    if(age < 0) {
        die();
        return;
    }
    if(activation >= 10) {
        TH1array.push_back(TH1(xLoc,yLoc,age));
        die();
        return;
    }
    if(activation <= -10) {
        TH2array.push_back(TH2(xLoc,yLoc,age));
        die();
        return;
    }
}

void TH1::step() {
    if(gp->c[IL12].val > 0) {
        //gp->c[IFNg].val = updateWrapper(2*gp->c[IFNg].val + gp->c[TNF].val + gp->c[IL1].val + gp->c[IL12].val, multipliers[IFNg]);
        update("TH1_threshold", false);
    }
    wiggle();
    move();
    age--;
    if(age < 0) {
        die();
        return;
    }
}

void TH2::step() {
    if(gp->c[IL10].val > 0) {
        //gp->c[IL4].val = updateWrapper(gp->c[IL4].val + gp->c[IL10].val, multipliers[IL4]);
        //gp->c[IL10].val = updateWrapper(gp->c[IL10].val*2, multipliers[IL10]);
        update("TH2_threshold", false);
    }
    wiggle();
    move();
    age--;
    if(age < 0) {
        die();
        return;
    }
}

void pmn_marrow::step() {
    int x,y,temp,n,i;
    int divisor = new_pmn ? 50 : 100;
    n = int(1 + total_GCSF/divisor);
    
    for(i = 0; i < n; i++) {
        temp = distribution10(generator);
        if(temp < 1) {
            x = distribution_xDim(generator);
            y = distribution_yDim(generator);
            int age = new_pmn ? 200 : 50;
            pmnArray.push_back(pmn(x,y,age));
        }
    }
}

void mono_marrow::step() {
    int x,y,temp;
    temp = distribution100(generator);
    if(temp < 1) {
        x = distribution_xDim(generator);
        y = distribution_yDim(generator);
        monoArray.push_back(mono(x,y,50,1000,1000));
    }
}

void TH0_germ::step() {
    int x,y,temp;
    x = xLoc;
    y = yLoc;
    temp = distribution100(generator);
    if(temp < 1) {
        TH0array.push_back(TH0(x,y,100));
    }
}

void TH1_germ::step() {
    int x,y,temp;
    x = xLoc;
    y = yLoc;
    temp = distribution100(generator);
    if(temp < 1) {
        TH1array.push_back(TH1(x,y,100));
    }
}

void TH2_germ::step() {
    int x,y,temp;
    x = xLoc;
    y = yLoc;
    temp = distribution100(generator);
    if(temp < 1) {
        TH2array.push_back(TH2(x,y,100));
    }
}

void heal(int index){
    ecArray[index].gp->oxy = 100;
    ecArray[index].ec_roll = 0;
    ecArray[index].ec_stick = 0;
    ecArray[index].ec_migrate = 0;
    ecArray[index].gp->infection = 0;
    ecArray[index].ec_activation = 0;
}
