// =========================================================================
// HYSIM V1.0 (November 2017)
// CREATED BY NACHO MOLINA AND SAMUEL ZAMBRANO
// =========================================================================

#include <assert.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <string>
#include <sstream>
#include <time.h>
#include <map>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <sys/time.h>
#include "hysim.h"
using namespace std;


// =========================================================================
// MAIN 
// =========================================================================
int main(int argc, char** argv) {

  // Initialitation: 
  // ----------------------------------------------------

  // Arguments
  if(argc != 9) {cerr << "Usage: <seed> <numcells> <timetotal> <dt> <file_variables> <file_reactions> <fileout_times> <fileout_trace>" << endl;exit(1);}
  long seed = atoi(argv[1]);
  double numc = atof(argv[2]);
  double tlast = atof(argv[3]);
  double tstep = atof(argv[4]);
  char *filevar = argv[5];
  char *filerea = argv[6];
  char *filetimes = argv[7];
  char *filetrace = argv[8];
  numcells = numc;
  tfinal = tlast;
  dt = tstep;

  // Time and seed
  time_t start, finish;
  double deltatime;
  start = time(NULL);
  seed  = time_msec(seed); 
  srandom(seed);
  cerr << "SEED: " << seed << endl; 
  
  // Read variables
  read_variables(filevar);

  // Read reactions
  read_reactions(filerea);

  // Initialize ODE solver
  init_ode_solver();
    

  // Stochatic simulation: 
  // ----------------------------------------------------

  for (int c = 1; c <= numcells; c++) {
    cerr << "CELL: " << c << endl;

    // Initializtaion 
    tperturb = tperturb0;
    init_variables();
    double t0 = 0;
    double t = 0;
    int r = 0;
    
    // open files
    open_files(filetimes,filetrace,c);
    
    // Running simmulation
    while (t < tlast) {
      t0 = t;
      
      // Choosing time of next reaction
      t = get_next_time(t0);
      
      // Choosing next reaction
      r = get_next_reaction();
      
      // Printing
      print_times(t-t0);
            
      // Updating state 
      update_state(r);
      
    }
    cout << c << " " << t << endl;
  }
  
  // Time
  // ----------------------------------------------------
  finish = time(NULL);
  deltatime = (finish - start)/60.0;
  cerr << "TIME: " << deltatime << " min" << endl;
}


// =========================================================================
// Core functions 
// =========================================================================

// setting perturbation
void init_variables() {
  for (int n = 0; n < numvariables; n++) {
      variables[n].value = variables[n].init;
  }
}

// setting perturbation
void set_perturb(double y[]) {
  tperturb = INFINITY;//1./0.;//numeric_limits<double>::infinity();
  for (int n = 0; n < numvariables; n++) {
    if (variables[n].perturb == 1) {
      variables[n].value = variables[n].newvalue;
    }
  }
  for (int n = 0; n < numdetvar; n++) {y[n+1]=variables[detvar[n]].value;}
}

// open out files 
void open_files(char *file1, char *file2, int c) {
  printtimes = 0;
  printtrace = 0;

  // Check if previous files are open and if so close them
  if (out1.is_open()) {out1.close();}
  if (out2.is_open()) {out2.close();}
  
  // New files 
  string N("N");
  string filetag1(file1); 
  string filetag2(file2);
  string fileout1; 
  string fileout2;
  if (numcells == 1) {
    fileout1 = filetag1;
    fileout2 = filetag2;
  }
  else {
    fileout1 = filetag1 + "." + itos(c) + ".out";
    fileout2 = filetag2 + "." + itos(c) + ".out";
  }
  
  // Open files if requiered
  if (!(filetag1.compare(N) == 0)) {
    printtimes = 1;
    out1.open(fileout1.c_str());
  }

  if (!(filetag2.compare(N) == 0)) {
    printtrace = 1;
    out2.open(fileout2.c_str());
  }
}

// printing all variables
void print_trace(double t) {
  out2 << t << " ";
  for (int n = 0; n < numvariables; n++) {
    out2 << variables[n].value << " ";
  }
  out2 << endl;  
}

// printing times
void print_times(double t) {
  if (printtimes == 1) {
    out1 << t << " ";
    for (int n = 0; n < numstovar; n++) {
      out1 << variables[stovar[n]].value << " ";
    }
    out1 << endl;
  }
}

// update variables
void update_state(int r) {
  int n;
  
  for (n = 0; n < numvariables; n++) {
    variables[n].value += reactions[r].stoichvec[n];
  }
}

// get total rate
double get_total_rate() {
  double factor = 0;
  double Q = 0;
  int i,j,n,m,k;
  
  for (i = 0; i < numstorea; i++) {
    j = storea[i];
    
    factor = reactions[j].param;
    for (n = 0; n < reactions[j].numvar; n++) {
      
      m = reactions[j].variables[n];
      if (reactions[j].type != 2) {
        for (k = 0; k < reactions[j].coeffreac[n]; k++) {
          factor *= (variables[m].value - k)/(k+1);
        }
      }
      else {
        if (variables[m].value == 0) {
          factor  = 0;
        }
      }
    }
    Q += factor;
  }
  
  return Q;
}

// get cumulative propensity
double get_cum_prop(double cumprop[]) {
  double factor = 0;
  double Q = 0;
  int i,j,n,m,k;

  for (i = 0; i < numstorea; i++) {
    j = storea[i];

    factor = reactions[j].param;
    for (n = 0; n < reactions[j].numvar; n++) {
      
      m = reactions[j].variables[n];
      if (reactions[j].type != 2) {
	for (k = 0; k < reactions[j].coeffreac[n]; k++) {
	  factor *= (variables[m].value - k)/(k+1);
	}
      }
      else {
	if (variables[m].value == 0) {
	  factor  = 0;
	}
      }
    }
    
    Q += factor;
    cumprop[i] = Q;
  }
  
  return Q;
}

// get next reaction
int get_next_reaction() {
  if (numstorea != 0) {
    double cumprop[numstorea]; 
    double totprop = 0;;
    double R = 0;
    int i = 0;
    
    // getting cumulative propensities
    totprop = get_cum_prop(cumprop);
    
    // choosing reaction
    R = totprop*RNG;
    while (cumprop[i] < R) {i++;}
    
    return storea[i];
  }
  else {
    return 0;
  }
}

// getting next reaction time
double get_next_time(double t) {
 
  if (numstovar == 0) {
    t = get_next_time_determ(t);
  }
  else {
    t = get_next_time_stoch(t);
  }
 
  return(t);
}

// getting next reaction time
double get_next_time_stoch(double t) {
  double R = -log(1-RNG);
  double h = 1e-6;
  double tl,tf,newt;
  double yl[dimsys];
  double y[dimsys];

  // Initialize variables
  for (int n = 0; n < numdetvar; n++) {y[n+1]=variables[detvar[n]].value;}
  y[0] = 0;
  
  // Finding final time
  while (y[0] < R) {

    // Perturbing 
    if (t > tperturb) {
      set_perturb(y);
    }
    
    // Printing 
    if (printtrace == 1 && fmod(t,dt) == 0) {
      for (int n = 0; n < numdetvar; n++) {variables[detvar[n]].value = y[n+1];}
      print_trace(t);
    }
    
    // Recording initial value
    set_ode_var(yl,y);
    tl = t;

    // Calculating next step
    tf = dt*(1.+int(t/dt));
    gsl_odeiv_evolve_apply(evolve_ptr, control_ptr, step_ptr, &my_system, &t, tf, &h, y);
  
  }
  tf = t;

  // Initialization
  set_ode_var(y,yl);

  // Refining search (finding root)
  while (sqrt(pow((y[0]-R),2)) > ZERO) {

    t = tl;
    newt = (tf+tl)/2;
    gsl_odeiv_evolve_apply(evolve_ptr, control_ptr, step_ptr, &my_system, &t, newt, &h, y);
    
    if (y[0] > R) {
      tf = newt;
      set_ode_var(y,yl);
    }
    else {
      tl = newt;
      set_ode_var(yl,y);
    }
  }
  
  // update state
  for (int n = 0; n < numdetvar; n++) {
    variables[detvar[n]].value = y[n+1];
  }
  
  // printing
  //if (fmod(t,dt) == 0 && printtrace == 1) {
  //print_trace(t);
  //}

  return(t);
}

// getting next reaction time
double get_next_time_determ(double t) {
  double h = 1e-6;
  double y[dimsys];
  double tf = 0;

  // Initialize variables
  for (int n = 0; n < numdetvar; n++) {y[n+1]=variables[detvar[n]].value;}
  y[0] = 0;
  
  // Finding final time
  while (t <= tfinal) {
    
    // Perturbing 
    if (t > tperturb) {set_perturb(y);}

    // Printing trace
    if (printtrace == 1 && fmod(t,dt) == 0) {
      for (int n = 0; n < numdetvar; n++) {variables[detvar[n]].value = y[n+1];}
      print_trace(t);
    }
    
    // Calculating next step
    tf = dt*(1.+int(t/dt)); 
    gsl_odeiv_evolve_apply(evolve_ptr, control_ptr, step_ptr, &my_system, &t, tf, &h, y);
  }
  
  // update state
  for (int n = 0; n < numdetvar; n++) {
    variables[detvar[n]].value = y[n+1];
  }
  
  // printing
  //if (printtrace == 1) {
  //print_trace(t);
  //}
  
  return(t);
}

// =========================================================================
// Functions to solve ODE
// =========================================================================

// Set ode variables
void set_ode_var (double a[], double b[]) {
  for (int n = 0; n < dimsys; n++) {
    a[n] = b[n];
  }
}

// System of differential equtions
int diffeq(double t, const double y[], double f[], void *params_ptr) {
  double factor = 0;
  double Q = 0;

  // Updating 
  for (int n = 0; n < numdetvar; n++) {variables[detvar[n]].value = y[n+1];}

  // get total rate
  Q = get_total_rate();

  // Integrating total rate
  f[0] = Q;

  // Deterministic variables
  int n,m,i,j,l,k;
  for (n = 0; n < numdetvar; n++) {
    m = detvar[n];
    f[n+1] = 0;

    for (i = 0; i < variables[m].numreac; i++) {
      j = variables[m].reactions[i];
      factor = reactions[j].stoichvec[m]*reactions[j].param;

      for (l = 0; l < reactions[j].numvar; l++) {
	k = reactions[j].variables[l];
	factor *= pow(variables[k].value,reactions[j].coeffreac[l]);
      }
      f[n+1] += factor;
    }
  }
    
  return GSL_SUCCESS;           /* GSL_SUCCESS defined in gsl/errno.h as 0 */
}

// Some agly variables needed for the ODE solver
void init_ode_solver() {
  
  // gsl rnd
  int seed = time(NULL);
  const gsl_rng_type *gsl_R;
  gsl_rng *gsl_r;
  gsl_rng_env_setup();      
  gsl_R = gsl_rng_default;
  gsl_r = gsl_rng_alloc(gsl_R);
  gsl_rng_set(gsl_r,seed);
  
  // ode solver varibales
  double eps_abs = 1.e-8;       // absolute error requested
  double eps_rel = 1.e-10;      // relative error requested
  
  // allocate/initialize the stepper, the control function, and the evolution function.
  type_ptr    = gsl_odeiv_step_rkf45;
  step_ptr    = gsl_odeiv_step_alloc(type_ptr,dimsys);
  control_ptr = gsl_odeiv_control_y_new(eps_abs,eps_rel);
  evolve_ptr  = gsl_odeiv_evolve_alloc(dimsys);

  // load values into the my_system structure
  my_system.function  = diffeq; // the right-hand-side functions dy[i]/dt
  my_system.dimension = dimsys; // number of diffeq's
  my_system.params    = NULL;   // parameters to pass to rhs and jacobian
}


// =========================================================================
// Functions to parsing files 
// =========================================================================
void read_variables(char *file) {
  double newvalue;
  double value;
  string perturb;
  string type;
  string name;
  string stmp;
  int perturbint;
  int typeint;


  ifstream in(file);
  if(!in) {
    cerr << "ERROR: file could not be opened" << endl;
    exit(1);
  }

  numdetvar = 0;
  numstovar = 0;
  numvariables = 0;
  while (!in.eof()) {
    if (in >> stmp >> numvariables >> tperturb) {
      if (tperturb == 0) {tperturb = INFINITY;}
      tperturb0 = tperturb;

      for (int n = 0; n < numvariables; n++) {
	in >> type >> name >> value >> perturb >> newvalue;

	if (type == "d:") {
	  detvar.push_back(n);
	  numdetvar++;
	  typeint = 0;
	}
	else if (type == "s:") {
	  stovar.push_back(n);
	  numstovar++;
	  typeint = 1;
	}

	if (perturb == "p:") {
	  perturbint = 1;
	}
	else {
	  perturbint = 0;
	}

	variables.push_back(*(new variable));
	variables[n].type  = typeint;
	variables[n].name  = name;
	variables[n].init  = value;
	variables[n].value = value;
	variables[n].newvalue = newvalue;
	variables[n].perturb = perturbint;
	variables[n].numreac = 0;
	variables[n].reactions = *(new vector <int>);
	name2int[name] = n;
      }
    }
  }
  dimsys = numdetvar+1;

  //for (int n = 0; n < numvariables; n++) {
  //cerr << variables[n].name << "\t" << variables[n].perturb << "\t" << variables[n].value << "\t" << variables[n].newvalue << "\t" << tperturb << endl;
  //}
}

void read_reactions(char* file) {
  double param;
  string type;
  string name;
  string stmp;
  int typeint;
  int number;

  ifstream in(file);
  if(!in) {
    cerr << "ERROR: file could not be opened" << endl;
    exit(1);
  }

  numdetrea = 0;
  numstorea = 0;
  numreactions = 0;
  while (!in.eof()) {
    if (in >> stmp >> numreactions) {
      for (int i = 0; i < numreactions; i++) {
	in >> type >> name >> param >> stmp;
	if (type == "d:") {
	  detrea.push_back(i);
	  numdetrea++;
	  typeint = 0;
	}
	else if (type == "s:") {
	  storea.push_back(i);
	  numstorea++;
	  typeint = 1;
	}
	else if (type == "q:") {
	  storea.push_back(i);
          numstorea++;
          typeint = 2;
	}

	reactions.push_back(*(new reaction));
	reactions[i].type  = typeint;
	reactions[i].name  = name;
	reactions[i].param = param;
	reactions[i].numvar = 0;
	reactions[i].variables = *(new vector <int>);
	reactions[i].coeffreac = *(new vector <int>);
	reactions[i].stoichvec = *(new vector <int> (numvariables,0));

	stmp = "";
	int n = 0;
	int sign = -1;
	while (stmp != ":") {
	  in >> number >> name >> stmp;
	  n = name2int[name];
	  reactions[i].stoichvec[n] += sign*number;
	  if (sign == -1) {
	    reactions[i].numvar++;
	    reactions[i].variables.push_back(n);
	    reactions[i].coeffreac.push_back(number);
	  }
	  if (stmp == "->") {sign = 1;}
	}
      }
    }
  }

  for (int i = 0; i < numreactions; i++) {
    for (int n = 0; n < numvariables; n++) {
      if (reactions[i].stoichvec[n] != 0) {
	variables[n].reactions.push_back(i);
	variables[n].numreac++;
      }
    }    
  }
}

// integer to string
string itos(int i) {
  stringstream ss;
  ss << i;
  return(ss.str());
}

long time_msec(long seed) {
  struct timeval start;
  long mtime, seconds, useconds;
  
  if (seed == 0) {
    gettimeofday(&start, NULL);
    seconds  = start.tv_sec;
    useconds = start.tv_usec;
    
    mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
  }
  else {
    mtime = seed;
  }

  return mtime;
}
