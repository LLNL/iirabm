#define BOOST_PYTHON_MAX_ARITY 21

#include <iostream>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include "boost/python/extract.hpp"

#include "core/Parameters.h"
#include "core/agents.h"
#include "core/environment.h"
#include "core/simulationFunctions.h"
#include "simulationControl.h"

using namespace std;
using namespace boost::python;

Simulation::Simulation(int inj_number_arg, float oxyHeal_arg, int infectSpread_arg, int numInfectRepeat_arg, int numRecurInj_arg, int seed_arg=seed, bool new_pmn_arg=false, bool break_ties_arg=true, bool big_shuffle_arg=true)
{
    new_pmn = new_pmn_arg;
    break_ties = break_ties_arg;
    big_shuffle = big_shuffle_arg;
    generator.seed(seed_arg);
    initialize();
    i = 0; // Time
    actionIsApplied = false;
    
    inj_number = inj_number_arg;
    oxyHeal = oxyHeal_arg;
    infectSpread = infectSpread_arg;
    numInfectRepeat = numInfectRepeat_arg;
    numRecurInj = numRecurInj_arg;    
}

// This method implements the controls. It sets the multipliers for this time step only, then reverts them back to 1.0 the following time step. (If they are given another action the following time step, the 1s will be overwritten.)
void Simulation::receiveAction(boost::python::numpy::ndarray action)
{
    for(int i = 0; i < NUM_DIFFUSIBLE_CYTOKINES; i++) {
        if(i == cytotox || i == endotoxin) {
            continue; // cytotox and endotoxin don't have multipliers and are not part of the intervention
        }
        int j = i-2; // HACK: Multiplier indices from Python start at PAF=0, not PAF=2
        multipliers[i] = extract<float>(action[j]);
    }
    
    actionIsApplied = true;
}

void Simulation::clearIntervention()
{
    for(int i = 0; i < NUM_DIFFUSIBLE_CYTOKINES; i++) {
        if(i == cytotox || i == endotoxin) {
            continue; // cytotox and endotoxin don't have multipliers and are not part of the intervention
        }
        multipliers[i] = 1;
    }
    
    actionIsApplied = false;
}

void Simulation::doStep()
{
    if(verbosity){cout << i << " " << total_IL1 << endl;}
    //cout<<" i == "<<i<<" numTimeSteps: "<<numTimeSteps<<endl;

    if(i == 0) {
        injure_infectionFRD(inj_number);
    }

//--- added by santiago10
//--- LPS injection 
    if(i == 0) { ic.add(bolusSize*xDim*yDim); }
    else if(i <= endotoxinInfusionStop) { ic.add(infusionSize*xDim*yDim/endotoxinInfusionStop); }
    if(ic.reserve > 0)
       ic.meter(i);
    
    updateSystemOxy();
    
    // Execute all Cells
    if(big_shuffle) stepAllCells();
    else {
        stepCells<TH0>(TH0array);
        erase_dead<TH0>(TH0array);
        for(int i = 0; i < ecArray.size(); i++) ecArray[i].step(); // Can't use stepCells<EC> because ecArray is a vector
        stepCells<pmn>(pmnArray);
        erase_dead<pmn>(pmnArray);
        stepCells<mono>(monoArray);
        erase_dead<mono>(monoArray);
        stepCells<TH1>(TH1array);
        erase_dead<TH1>(TH1array);
        stepCells<TH2>(TH2array);
        erase_dead<TH2>(TH2array);
        stepCells<pmn_marrow>(pmn_marrowArray);
        stepCells<mono_marrow>(mono_marrowArray);
        stepCells<TH1_germ>(TH1_germArray);
        stepCells<TH2_germ>(TH2_germArray);
        stepCells<TH0_germ>(TH0_germArray);
    }
    
    // Degradation
    evaporate();		
    
    if((i+1) % (injuryStep-1) == 0) {
        for(int ii=1; ii <= numRecurInj; ii++) {
            recur_injury();
        }
    }
    
    // Diffusion
    diffuse();
    
    if((i+1-102) % injuryStep == 0 && antibioticMultiplier > 0) {
        applyAntibiotics();
    }
    
    i++;
    
    // At the end of the time step, clear the intervention (revert multipliers to 1.0). If an action is given every time step, this function will not affect simulation output.
    if(actionIsApplied) {
        clearIntervention();
    }
}

boost::python::numpy::ndarray Simulation::getState(int channel)
{
    boost::python::list variables; // If channel == -1 (default), list of tuples. Otherwise, a list of scalars.
    
    for(int y = 0; y < yDim; y++) {
        for(int x = 0; x < xDim; x++) {
            GridPoint* gp = Grid::getInstance()(x,y);
            // channel = -1: All values (oxy, infection, cytokines)
            // channel = 0, 1: oxy, infection
            // channel = 2..15: cytokine corresponding to channel-2
            if(channel >= 2 + NUM_DIFFUSIBLE_CYTOKINES)
                cout << "Error: Bad channel input." << endl;
            switch(channel) {
                case -1:
                    variables.append(boost::python::make_tuple(
                                                               // Previous order:
                                                               // oxy, infection, TNF, endotoxin, sTNFr, IL10, GCSF, IFNg, PAF, IL1, IL4, IL8, IL12, sIL1r, IL1ra
                                                               gp->oxy,
                                                               gp->infection,
                                                               gp->c[0].val,
                                                               gp->c[1].val,
                                                               gp->c[2].val,
                                                               gp->c[3].val,
                                                               gp->c[4].val,
                                                               gp->c[5].val,
                                                               gp->c[6].val,
                                                               gp->c[7].val,
                                                               gp->c[8].val,
                                                               gp->c[9].val,
                                                               gp->c[10].val,
                                                               gp->c[11].val,
                                                               gp->c[12].val,
                                                               gp->c[13].val,
                                                               gp->c[14].val,
                                                               gp->c[15].val
                                                               ));
                    break;
                case 0:
                    variables.append(gp->oxy);
                    break;
                case 1:
                    variables.append(gp->infection);
                    break;
                default: // Cases 2 to NUM_DIFFUSIBLE_CYTOKINES+1
                    variables.append(gp->c[channel - 2].val);
                    break;
            }
        }
    }
    return boost::python::numpy::array(variables); // Convert list to ndarray
}

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(getState_overloads, Simulation::getState, 0, 1)

boost::python::numpy::ndarray Simulation::getTotals()
{
    // Order: oxyDeficit, totalInfection, total_CYT
    boost::python::list totals;
    
    // TBD: Convert totals to an array
    totals.append(oxyDeficit);
    totals.append(totalInfection);
    totals.append(total_cytotox);
    totals.append(total_endotoxin);
    totals.append(total_PAF);
    totals.append(total_TNF);
    totals.append(total_sTNFr);
    totals.append(total_IL1);
    totals.append(total_sIL1r);
    totals.append(total_IL1ra);
    totals.append(total_IFNg);
    totals.append(total_IL4);
    totals.append(total_IL8);
    totals.append(total_IL10);
    totals.append(total_IL12);
    totals.append(total_GCSF);
    totals.append(total_IL6);
    totals.append(total_sIL6r);

    return boost::python::numpy::array(totals); // Convert list to ndarray
}

boost::python::numpy::ndarray Simulation::getCellCounts()
{
    // Order: pmn, mono, TH0, TH1, TH2
    boost::python::list counts;
    
    counts.append(pmn::aliveCount);
    counts.append(mono::aliveCount);
    counts.append(TH0::aliveCount);
    counts.append(TH1::aliveCount);
    counts.append(TH2::aliveCount);
    
    return boost::python::numpy::array(counts); // Convert list to ndarray
}

// Returns a xDim*yDim by 5 ndarray of the number of Cells at each location
boost::python::numpy::ndarray Simulation::getCellLocations()
{
    int numCellTypes = 5; // Order: pmn, mono, TH0, TH1, TH2
    
    int cellLocations[xDim*yDim][numCellTypes];
    for(int x = 0; x < xDim; x++) {
        for(int y = 0; y < yDim; y++) {
            int id = y*xDim + x;
            for(int i = 0; i < numCellTypes; i++) {
                cellLocations[id][i] = Grid::getInstance()(x,y)->numCellsOfEachType[i+1]; // The +1 is to skip ECs. I.e. the cellTypeIDs of MobileCells are 1..5.7
            }
        }
    }
    
    // Need to start with list of tuples since this version of Boost.Python does not support multidimensional ndarrays
    boost::python::list cellLocations_PyList;
    for(int i = 0; i < xDim*yDim; i++) {
        cellLocations_PyList.append(boost::python::make_tuple(
                                                              cellLocations[i][0],
                                                              cellLocations[i][1],
                                                              cellLocations[i][2],
                                                              cellLocations[i][3],
                                                              cellLocations[i][4]
                                                              ));
    }
    
    return boost::python::numpy::array(cellLocations_PyList);
}

// Returns an (xDim*yDim,) ndarray of which Cell type is at each location (zero means no Cells), using an arbitrary priority.
boost::python::numpy::ndarray Simulation::getPriorityCell()
{
    // Priority: pmn, mono, TH0, TH1, TH2
    int numMobileCells = 5;
    boost::python::list priorityCell;
    for(int y = 0; y < yDim; y++) {
        for(int x = 0; x < xDim; x++) {
            int priority;
            int *count = Grid::getInstance()(x,y)->numCellsOfEachType;
            int max = *std::max_element(count + 1, count + numMobileCells + 1); // +1 to skip ECs.
            if(max == 0) priority = 0; // Zero means no Cells
            else priority = std::distance(count + 1, std::max_element(count + 1, count + numMobileCells + 1)) + 1; // First three +1s to skip ECs. Last +1 because zero means no Cells.
            priorityCell.append(priority);
        }
    }

    return boost::python::numpy::array(priorityCell);
}

// Edits individual entries of they cytokine matrix. matrix_dict keys are 3-tuples of strings (ruleName, target, effector); values are floats.
void editMatrix(boost::python::dict matrix_dict)
{
    namespace bp = boost::python;
    
    bp::list keys = matrix_dict.keys();
    for(int i = 0; i < bp::len(keys); i++) {
        //string ruleName = bp::extract<string>(bp::extract<bp::tuple>(keys[i])()[0]);
        string ruleName = bp::extract<string>(keys[i][0]);
        string target   = bp::extract<string>(keys[i][1]);
        string effector = bp::extract<string>(keys[i][2]);
        float val = bp::extract<float>(matrix_dict[keys[i]]);
        coefficientMatrix[ruleName][coefficientMap[target]][coefficientMap[effector]] = val;
        //cout << ruleName << " " << target << " " << effector << " " << coefficientMatrix[ruleName][enumMap[target]][enumMap[effector]] << endl;
    }
    
    addAllRules(); // Update the rules with new matrix.
}


void setMatrix(boost::python::numpy::ndarray matrix)
{
    namespace bp = boost::python;
    
    const int numRules = 10;
    string ruleNames[numRules] = {"EC_activation", "EC_midhealthy","EC_unhealthy", "pmn_primed", "pmn_burst", "mono_function", "mono_activation", "mono_unactivated", "TH1_threshold", "TH2_threshold"};
    
    assert(bp::len(matrix) == numRules*NUM_COEFFICIENTS);
    
    cout << "Setting matrix..." << endl;
    
    for(int i = 0; i < numRules*NUM_COEFFICIENTS; i++) {
        int blockIndex = i/NUM_COEFFICIENTS; // Index of "block" of matrix
        int rowIndex = i%NUM_COEFFICIENTS; // Index of row within a block of matrix
        for(int j = 0; j < NUM_COEFFICIENTS; j++) {
            coefficientMatrix[ruleNames[blockIndex]][rowIndex][j] = bp::extract<float>(matrix[i][j]);
        }
    }
    
    addAllRules(); // Update the rules with the new matrix.
}


float getCoefficient(boost::python::list arglt)
{
    namespace bp = boost::python;

    string ruleName = bp::extract<string>(arglt[0]);
    string target   = bp::extract<string>(arglt[1]);
    string effector = bp::extract<string>(arglt[2]);
    float val = coefficientMatrix[ruleName][coefficientMap[target]][coefficientMap[effector]];
    return val;
}


boost::python::list getMatrixIds()
{
    const int numRules = 10;
    string ruleNames[numRules] = {"EC_activation", "EC_midhealthy","EC_unhealthy", "pmn_primed", "pmn_burst", "mono_function", "mono_activation", "mono_unactivated", "TH1_threshold", "TH2_threshold"};
    boost::python::list matrix_PyList;

    for(int i = 0; i < numRules; i++) 
    {
      for(int j = 0; j < NUM_COEFFICIENTS; j++) 
      for(int k = 0; k < NUM_COEFFICIENTS; k++) {
         boost::python::list Lids;
         Lids.append(ruleNames[i]);
         Lids.append(coefficientNames[j]);
         Lids.append(coefficientNames[k]);
         matrix_PyList.append(Lids);
//       matrix_PyList.append((ruleNames[i],coefficientNames));
       }
      }

    return matrix_PyList;
}

// This is a bit hacky...
boost::python::numpy::ndarray getMatrix()
{
    const int numRules = 10;
    string ruleNames[numRules] = {"EC_activation", "EC_midhealthy","EC_unhealthy", "pmn_primed", "pmn_burst", "mono_function", "mono_activation", "mono_unactivated", "TH1_threshold", "TH2_threshold"};

    // Need to start with list of tuples since this version of Boost.Python does not support multidimensional ndarrays
    boost::python::list matrix_PyList;
    for(int i = 0; i < numRules; i++) {
        for(int j = 0; j < NUM_COEFFICIENTS; j++) {
            matrix_PyList.append(boost::python::make_tuple(
                                                       coefficientMatrix[ruleNames[i]][j][0],
                                                       coefficientMatrix[ruleNames[i]][j][1],
                                                       coefficientMatrix[ruleNames[i]][j][2],
                                                       coefficientMatrix[ruleNames[i]][j][3],
                                                       coefficientMatrix[ruleNames[i]][j][4],
                                                       coefficientMatrix[ruleNames[i]][j][5],
                                                       coefficientMatrix[ruleNames[i]][j][6],
                                                       coefficientMatrix[ruleNames[i]][j][7],
                                                       coefficientMatrix[ruleNames[i]][j][8],
                                                       coefficientMatrix[ruleNames[i]][j][9],
                                                       coefficientMatrix[ruleNames[i]][j][10],
                                                       coefficientMatrix[ruleNames[i]][j][11],
                                                       coefficientMatrix[ruleNames[i]][j][12],
                                                       coefficientMatrix[ruleNames[i]][j][13],
                                                       coefficientMatrix[ruleNames[i]][j][14],
                                                       coefficientMatrix[ruleNames[i]][j][15],
                                                       coefficientMatrix[ruleNames[i]][j][16],
                                                       coefficientMatrix[ruleNames[i]][j][17],
                                                       coefficientMatrix[ruleNames[i]][j][18],
                                                       coefficientMatrix[ruleNames[i]][j][19],
                                                       coefficientMatrix[ruleNames[i]][j][20]));
        }
    }
    
    return boost::python::numpy::array(matrix_PyList); // Convert to 2d ndarray
}

#define addelement(var_map,var_list,python_list,var) \
    for (int j = 0; j < boost::python::len(python_list); ++j) \
    { \
        string var_name = boost::python::extract<string>(python_list[j]); \
        if (var_name == getname(var)) \
        { \
            std::ostringstream oss; \
            oss << &var; \
            std::string address = oss.str(); \
            if(verbosity) cout << "  ==: " << var_name << "   +++ " << var << endl; \
            var_list.append(address.c_str()); \
            var_map[var_name] = &var;\
            break; \
        } \
    }

boost::python::numpy::ndarray run(boost::python::dict parameters)
{
    namespace bp = boost::python;
    
    //bp::numeric::array::set_module_and_type("numpy", "ndarray");

    bp::list Lvars; // TBD: Remove this debug variable
    map<string,float*>float_var_map; // Map from float variable name to pointer to that variable
    map<string,int*>int_var_map; // Map from int variable name to pointer to that variable

    // First handle outputs:
    bp::list outputs; // List of ABM output variables (e.g. total_IL1)
    if(parameters.has_key("outputs"))       outputs = bp::extract<bp::list>             (parameters["outputs"]);

    //-- santiago10
    validation=0;
    bp::list validationSteps; // List of timesteps
    if(parameters.has_key("timesteps")) { validation=1;      validationSteps = bp::extract<bp::list>             (parameters["timesteps"]);}

    // HACK: Manually turn on verbosity first.
    if(parameters.has_key("verbosity") && bp::extract<int>(parameters["verbosity"]) == 1)
        verbosity = 1;
        
    // Map names of outputs to their variables
    addelement(float_var_map,Lvars,outputs,oxyDeficit)
    addelement(float_var_map,Lvars,outputs,totalInfection)
    addelement(float_var_map,Lvars,outputs,total_cytotox)
    addelement(float_var_map,Lvars,outputs,total_endotoxin)
    addelement(float_var_map,Lvars,outputs,total_PAF)
    addelement(float_var_map,Lvars,outputs,total_TNF)
    addelement(float_var_map,Lvars,outputs,total_sTNFr)
    addelement(float_var_map,Lvars,outputs,total_IL1)
    addelement(float_var_map,Lvars,outputs,total_sIL1r)
    addelement(float_var_map,Lvars,outputs,total_IL1ra)
    addelement(float_var_map,Lvars,outputs,total_IFNg)
    addelement(float_var_map,Lvars,outputs,total_IL4)
    addelement(float_var_map,Lvars,outputs,total_IL6)
    addelement(float_var_map,Lvars,outputs,total_IL8)
    addelement(float_var_map,Lvars,outputs,total_IL10)
    addelement(float_var_map,Lvars,outputs,total_IL12)
    addelement(float_var_map,Lvars,outputs,total_GCSF)
    addelement(float_var_map,Lvars,outputs,total_IL6)
    addelement(float_var_map,Lvars,outputs,total_sIL6r)
    addelement(float_var_map,Lvars,outputs,system_oxy)
    addelement(float_var_map,Lvars,outputs,timeperiod)
    
    // Next handle all other keys (i.e. those that are ABM parameter names)
    bp::list keys = bp::list(parameters.keys());
    
    // Map names of parameters to their variables
    addelement(int_var_map,Lvars,keys,xDim)
    addelement(int_var_map,Lvars,keys,yDim)
    addelement(float_var_map,Lvars,keys,antibioticMultiplier)
    addelement(int_var_map,Lvars,keys,numTimeSteps)
    addelement(float_var_map,Lvars,keys,bolusSize)
    addelement(float_var_map,Lvars,keys,infusionSize)
    addelement(float_var_map,Lvars,keys,infusionSizeMultiplier)
    addelement(float_var_map,Lvars,keys,endotoxinInfusionStop)
    addelement(float_var_map,Lvars,keys,meteringRate)
    addelement(int_var_map,Lvars,keys,inj_number)
    addelement(int_var_map,Lvars,keys,numInfectRepeat)
    addelement(float_var_map,Lvars,keys,oxyHeal)
    addelement(int_var_map,Lvars,keys,infectSpread)
    addelement(int_var_map,Lvars,keys,numRecurInj)
    addelement(int_var_map,Lvars,keys,seed)
    addelement(int_var_map,Lvars,keys,injuryStep)
    addelement(int_var_map,Lvars,keys,verbosity)
    addelement(int_var_map,Lvars,keys,validation)
    addelement(float_var_map,Lvars,keys,TNF_scale)
    addelement(int_var_map,Lvars,keys,boundary_mode)
    addelement(int_var_map,Lvars,keys,neighborhood_mode)
    addelement(int_var_map,Lvars,keys,interaction_length_mode)
    
    // Set each parameter with the given value
    for(int i = 0; i < bp::len(keys); i++) {
        string key = bp::extract<string>(keys[i]);
        if(float_var_map.count(key) >= 1)       *float_var_map[key] = bp::extract<float>(parameters[key]);
        else if(int_var_map.count(key) >= 1)    *int_var_map[key]   = bp::extract<int>(parameters[key]);
    }

    // Instantiating the Simulation calls simulationFunctions::initialize().
    //      xDim and yDim must be set before this!
    //      The cytokine update matrix must be set after this!
    Simulation sim = Simulation();

    // Now that the default matrix has been created, either replace it or edit it.
    if(parameters.has_key("matrix"))        setMatrix(bp::extract<bp::numpy::ndarray>   (parameters["matrix"]));
    if(parameters.has_key("matrix_dict"))   editMatrix(bp::extract<bp::dict>             (parameters["matrix_dict"]));

    // Run the model
    bp::list timepoints; // List of 1D ndarrays; converted into 2D ndarray when returning.
 
//    vector<int> validationTimeSteps = {0, 9, 13, 17, 26, 34, 51, 69}; // Ignored if validation==false.
    vector<int> validationTimeSteps;
    if (validation)
    {
       for(int i = 0; i < bp::len(validationSteps); i++) 
         validationTimeSteps.push_back(bp::extract<int>(validationSteps[i]));
         std::vector<int>::iterator result = std::max_element(validationTimeSteps.begin(), validationTimeSteps.end());
         numTimeSteps =  *result +1;
       // cout<<"numTimeSteps = "<<numTimeSteps<<endl;
       }
    
    if(validation) {
        numTimeSteps = validationTimeSteps.back() + 1; // Need to go an extra time step to get the last data point.
    }

    timeperiod=0;
    for(int t = 0; t < numTimeSteps; t++) {
        sim.doStep();
        bool printThisStep = true;
        if(validation) { // If in validation mode, only print according to validation time steps.
            printThisStep = false;
            for(int v = 0; v < validationTimeSteps.size(); v++) {
                if(validationTimeSteps[v] == t) {
                    printThisStep = true;
                    break;
                }
            }
        }
        if(!printThisStep) continue;
        timeperiod = t;
        bp::list cytokine_totals;
        for(int i = 0; i < bp::len(outputs); i++) {
            string name = bp::extract<string>(outputs[i]);
            cytokine_totals.append(*float_var_map[name]);
        }

        //debug
        if (0)
        { 
           cout <<" seed: "<< seed<<" time: "<<t;
           for(int i = 0; i < bp::len(cytokine_totals); i++)
               cout<<" "<<bp::extract<float>(cytokine_totals[i]);
           cout<<endl;
           }

        timepoints.append(boost::python::numpy::array(cytokine_totals));
    }
    
    if(verbosity) cout << "Done!" << endl;
    
    return boost::python::numpy::array(timepoints); // Convert list to 2D ndarray
}

Simulation::Simulation() {
//--- added by santiago10
//--- LPS injection 
    ic.set(endotoxin, meteringRate);

    generator.seed(seed); // Default 32
    initialize();
    i = 0; // Time
    actionIsApplied = false;
}

// Sets the four parameters Chase is using for IL6 mechanism sweep
void setIL6Parameters(float v1, float v2, float v3, float v4)
{   
  coefficientMatrix["mono_activation"][coefficientMap["IL6"]][coefficientMap["infection"]] = v1; // Matrix element: mono_activation_IL6_infection
  coefficientMatrix["mono_activation"][coefficientMap["IL6"]][coefficientMap["constant"]] = v2; // Matrix element: mono_activation_IL6_constant
  addAllRules();
  proTH2_coefficient = v3; // Coefficient outside proTH2 expression
  proTH2_sIL6r_coefficient = v4; // Coefficient of sIL6r inside proTH2 expression
}

BOOST_PYTHON_MODULE(sepsis)
{
    //boost::python::numeric::array::set_module_and_type("numpy", "ndarray");
    Py_Initialize();
    boost::python::numpy::initialize();

    // Initialize global parameters
    initializeCoefficientMap(); // Initialize map from coefficient names to position along rows/columns.
    setDefaultMatrixParameters(); // Start with default matrix parameters.
    addAllRules(); // Actually add the rules to each Cell type.

    // Define module-level constants
    boost::python::scope().attr("NUM_STATES") = NUM_DIFFUSIBLE_CYTOKINES + 2 + 5; // Add [oxyDeficit, totalInfection] and 5 cell counts
    boost::python::scope().attr("NUM_ACTIONS") = NUM_DIFFUSIBLE_CYTOKINES - 2; // Exclude [cytotox, endotoxin]
    boost::python::scope().attr("COEFFICIENT_MAP") = toPythonDict(coefficientMap); // Map from string to column index

    // Simulation class and methods
    class_<Simulation>("Simulation", init<int, float, int, int, int, int, bool, bool, bool>())
    .def("step", &Simulation::doStep)
    .def("getState", &Simulation::getState, getState_overloads(args("channel")))
    .def("getTotals", &Simulation::getTotals)
    .def("getCellCounts", &Simulation::getCellCounts)
    .def("getCellLocations", &Simulation::getCellLocations)
    .def("getPriorityCell", &Simulation::getPriorityCell)
    .def("getNumTimeSteps", &Simulation::getNumTimeSteps)
    .def("receiveAction", &Simulation::receiveAction, args("multipliers"))
    ;
    
    // Module-level methods
    def("run", run, boost::python::arg("dt"));
    def("setMatrix", setMatrix, boost::python::arg("matrix"));
    def("editMatrix", editMatrix, boost::python::arg("matrix_dict"));
    def("getMatrix", getMatrix);
    def("getCoefficient", getCoefficient, boost::python::arg("arglt"));
    def("getMatrixIds", getMatrixIds);
    def("setIL6Parameters", setIL6Parameters, (boost::python::arg("v1"), boost::python::arg("v2"), boost::python::arg("v3"), boost::python::arg("v4")));
}
