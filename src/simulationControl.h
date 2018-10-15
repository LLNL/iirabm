#ifndef SIMULATIONCONTROL_H
#define SIMULATIONCONTROL_H

using namespace std;

class Simulation
{
public:
    Simulation(int inj_number_arg, float oxyHeal_arg, int infectSpread_arg, int numInfectRepeat_arg, int numRecurInj_arg, int seed_arg, bool new_pmn, bool break_ties_arg, bool big_shuffle_arg); // Constructor
    Simulation(); // Constructor
    void doStep(); // Single time step of the simulation
    boost::python::numpy::ndarray getState(int channel = -1); // Returns the spatial state of the simulation at the given channel
    boost::python::numpy::ndarray getTotals(); // Returns the aggregate measures
    boost::python::numpy::ndarray getCellCounts(); // Returns the Cell counts
    boost::python::numpy::ndarray getCellLocations(); // Returns the number of each Cell type at each GridPoint
    boost::python::numpy::ndarray getPriorityCell(); // Returns the priority Cell at each GridPoint (zero means no Cells)
    void receiveAction(boost::python::numpy::ndarray action);
    void clearIntervention(); // Stops the intervention at the end of the time step
    
    int getNumTimeSteps() { return numTimeSteps; }

private:
    int i;
    bool actionIsApplied;
//--- added by santiago10
//--- LPS injection 
    IntroCompartment ic;
};

boost::python::numpy::ndarray run(boost::python::dict dt); // Runs the whole model and returns all results (used for optimization/ABC)
void editMatrix(boost::python::dict matrix_dict); // Edit cytokine matrix given existing matrix and dictionary of new values.
void setMatrix(boost::python::numpy::ndarray matrix); // Set cytokine matrix given.
boost::python::numpy::ndarray getMatrix(); // Get the cytokine matrix as a 2D ndarray.
float getCoefficient(boost::python::list arglt);
boost::python::list getMatrixIds();
void setIL6Parameters(float, float, float, float); // Sets the 4 parameters Chase is using to sweep over IL6 mechanisms.

// Converts a C++ map to a Python dict
template <class K, class V>
boost::python::dict toPythonDict(std::map<K, V> map) {
    typename std::map<K, V>::iterator iter;
    boost::python::dict dictionary;
    for (iter = map.begin(); iter != map.end(); ++iter) {
        dictionary[iter->first] = iter->second;
    }
    return dictionary;
}

#endif
