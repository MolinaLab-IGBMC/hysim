#include <assert.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <string>
#include <map>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
using namespace std;

// Definitions
#define RNG ((random()+0.1)/2147483648.0)//RAND_MAX
#define ZERO 10e-6

// Structures
struct variable {
  int type;
  int numreac;
  int perturb;
  string name;
  double init;
  double value;
  double newvalue;
  vector <int> reactions;
};

struct reaction {
  int type;
  int numvar;
  string name;
  double param;
  vector <int> variables;
  vector <int> coeffreac;
  vector <int> stoichvec;
};

// ODE solver variables 
gsl_odeiv_step *step_ptr;
gsl_odeiv_control *control_ptr;
gsl_odeiv_evolve *evolve_ptr;
gsl_odeiv_system my_system;
const gsl_odeiv_step_type *type_ptr;

// Global variables
int dimsys;
int numcells;
int numdetvar;
int numdetrea;
int numstovar;
int numstorea;
int numvariables;
int numreactions;
int printtimes;
int printtrace;
double tperturb;
double tperturb0;
double tfinal;
double dt;
ofstream out1;
ofstream out2;
vector <int> detvar;
vector <int> detrea;
vector <int> stovar;
vector <int> storea;
vector <reaction> reactions;
vector <variable> variables;
map <string, int> name2int;

// Functions
void init_variables();
void init_ode_solver();
void read_reactions(char *file);
void read_variables(char *file);
void open_files(char *file1, char *file2, int c);
void set_ode_var (double a[], double b[]);
void print_times(double t);
void print_trace(double t);
void update_state(int r);
void set_perturb(double y[]);
double get_next_time(double t);
double get_next_time_stoch(double t);
double get_next_time_determ(double t);
double get_cum_prop(double cumprop[]);
double get_total_rate();
int get_next_reaction();
int diffeq(double t, const double y[], double f[], void *params_ptr);
string itos(int i);
long time_msec(long seed);
