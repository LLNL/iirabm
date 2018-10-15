
#ifndef SEPSISIOMACROS_H
#define SEPSISIOMACROS_H

#include <string>
#include <map>
#include <list>
#include <vector>
#include <fstream>
#include <iostream>
#include <typeinfo>
#include <typeindex>
#include <unordered_map>
#include <sstream>  
#include <iostream>
#include <iomanip> 

extern std::ostringstream os;

#define EXE "/nfs/apps/mpich2/1.4.1/bin/mpiexec"
#define SEPSISEXE "IIR_Sepsis"
#define updatestr(var,suff) \
        var = vname; \
        var += suff


#define getname(varname) retname(#varname)

#define SETVAL(var) setval(tmap.first,fmap,var,getname(var)) 
#define GETVAL(var) getval(tmap,fmap,var,getname(var))

#define SAVEPARAM(os,tmap,var) \
         SAVEPARAMMAP(os,tmap,var,getname(var))

#define SAVEOUTPARAM(os,tmap,var) \
         SAVEPARAMMAP(os,tmap,var,getname(var))

#define SAVEPARAMMAP(os,tmap,var,namevar) \
         os.str(std::string()); \
         os<< std::fixed<<std::showpoint<<std::setprecision(PRECISION)<<var; \
         tmap[namevar] = string(os.str()) 


#define SCANPARAM(ArgMc,param) \
        ArgMc(param); \
        ArgMc(param ## Start); \
        ArgMc(param ## Increment); \
        ArgMc(param ## Max)

#define PRECISION 3

#define PRNVAR(param,mpout,tmap,Lparams) \
     for (list<string>::iterator it = mpout[Lparams].begin();it != mpout[Lparams].end(); it++) \
      {  if (*it == getname(param)) \
       { SAVEOUTPARAM(os,tmap,param); break; } }


#define SETVAL_0(param) \
        SETVAL(param)

#define SETVAL_1(param) \
        SCANPARAM(SETVAL,param)

#define EXPANDPARAM_0(param,scan,...) \
        SETVAL_ ## scan (param)

#define EXPANDPARAM_1(param,scan,...) \
          PRNVAR(param,__VA_ARGS__)

//-- scan is set to 1 if variable is used to scan parameters, i.e., if variable has "companion" variables with suffixes "Start", "Increment" and "Max", 0 o.w.
#define EXECPARAM(param,scan,func,...) \
       EXPANDPARAM_ ## func (param,scan,__VA_ARGS__)

#define SETPARAM(p,i,...) \
    EXECPARAM(p,i,__VA_ARGS__) 

//-- variables to be output should be added here
#define LISTPARAMS(...)\
    SETPARAM(xDim,                0,__VA_ARGS__); \
    SETPARAM(yDim,                0,__VA_ARGS__); \
    SETPARAM(antibioticMultiplier,0,__VA_ARGS__); \
    SETPARAM(parameterInput,      0,__VA_ARGS__); \
    SETPARAM(numTimeSteps,        0,__VA_ARGS__); \
    SETPARAM(bolusSize,           0,__VA_ARGS__); \
    SETPARAM(infusionSize,        0,__VA_ARGS__); \
    SETPARAM(infusionSizeMultiplier,0,__VA_ARGS__); \
    SETPARAM(meteringRate,        0,__VA_ARGS__); \
    SETPARAM(inj_number,          1,__VA_ARGS__); \
    SETPARAM(numInfectRepeat,     1,__VA_ARGS__); \
    SETPARAM(oxyHeal,             1,__VA_ARGS__); \
    SETPARAM(infectSpread,        1,__VA_ARGS__); \
    SETPARAM(numRecurInj,         1,__VA_ARGS__); \
    SETPARAM(seed,                1,__VA_ARGS__); \
    SETPARAM(mainseed,            0,__VA_ARGS__); \
    SETPARAM(injuryStep,          0,__VA_ARGS__); \
    SETPARAM(oxyInfectFinal,      0,__VA_ARGS__); \
    SETPARAM(dailyOutput,         0,__VA_ARGS__); \
    SETPARAM(writeEverything,     0,__VA_ARGS__); \
    SETPARAM(fileName,            0,__VA_ARGS__); \
    SETPARAM(paramset,            0,__VA_ARGS__); \
    SETPARAM(verbosity,           0,__VA_ARGS__); \
    SETPARAM(validation,          0,__VA_ARGS__); \
    SETPARAM(runID,               0,__VA_ARGS__); \
    SETPARAM(timeperiod,          0,__VA_ARGS__); \
    SETPARAM(oxyDeficit,          0,__VA_ARGS__); \
    SETPARAM(total_cytotox,       0,__VA_ARGS__); \
    SETPARAM(total_endotoxin,     0,__VA_ARGS__); \
    SETPARAM(total_PAF,           0,__VA_ARGS__); \
    SETPARAM(total_TNF,           0,__VA_ARGS__); \
    SETPARAM(total_sTNFr,         0,__VA_ARGS__); \
    SETPARAM(total_IL1,           0,__VA_ARGS__); \
    SETPARAM(total_sIL1r,         0,__VA_ARGS__); \
    SETPARAM(total_IL1ra,         0,__VA_ARGS__); \
    SETPARAM(total_IFNg,          0,__VA_ARGS__); \
    SETPARAM(total_IL4,           0,__VA_ARGS__); \
    SETPARAM(total_IL8,           0,__VA_ARGS__); \
    SETPARAM(total_IL10,          0,__VA_ARGS__); \
    SETPARAM(total_IL12,          0,__VA_ARGS__); \
    SETPARAM(total_GCSF,          0,__VA_ARGS__); \
    SETPARAM(totalInfection,      0,__VA_ARGS__); \
    SETPARAM(system_oxy,          0,__VA_ARGS__); \
    SETPARAM(TNF_scale,           0,__VA_ARGS__)

#define INITOUTPUT(mnames,loutput) \
    mnames[loutput].clear()
   
#define INITREPEATPARAMS(moutmap) \
    INITOUTPUT(moutmap,outputBegin); \
    INITOUTPUT(moutmap,output); \
    INITOUTPUT(moutmap,outputEnd)

#endif





