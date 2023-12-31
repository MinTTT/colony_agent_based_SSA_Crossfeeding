#ifndef CONSTANTS_H_
#define CONSTANTS_H_

// default values for all variables
#define CELL_V 1. // calculate rets is　1023
// Toggle
extern double Green_conc;
extern double Red_conc;
extern double tau_G;
extern double tau_R;
extern double kG;
extern double kR;
extern double nG;
extern double nR;
extern double maxCellVol;
extern int reactionMetric[4][2];

// physical variables
extern double cellRadius;
extern double L_divide;		
extern double k_cc;
extern double k_wc;
extern double varL;
extern double varAngle;
extern double var_pos;
extern double viscosity;		
extern double wall_rough;
extern double gamma_t;
extern double gamma_n;
extern double cell_mu;
extern double wall_mu;
extern double density_threshold;
extern double tension;
extern double DH;

// time
extern double t_max;
extern double initial_dt;
extern double OutputTime;
extern double t0;
extern double UpdateTime;

// non-physical constants
extern int BoxX;
extern int BoxY;
extern int BoxZ;
extern int BoxZAgar;
extern int maxLevels;
extern double BoxLength;
extern int FilterLen;
extern double e;
// directory name
extern char DirName[500];
extern bool FieldOutPut;
// number max cells to simulate (for allocating memory)
extern int maxCells;	// maximum number of cells in the simulation

// nutrient constants
extern double Tortuosity;
extern double KC; 
extern double C_rate;
extern double DiffColony;
extern double DiffAgar;
extern double maxCarbon;
extern double maxGrowthRate;
extern double Maintenance_rate;
extern double Cdt;
extern double ConvCrit;
extern int minIter;
extern int maxIter;
extern int InterfaceCondition;
extern bool NutrientGSI;
extern double Rc;
extern double deltaP;

extern double yield ;
extern double rho0 ;
extern double molCarbon ;

// colony constants
extern int refinementGridHeight;
extern double SigmaMoL;  // Pan's code

#endif /* CONSTANTS_H_ */
