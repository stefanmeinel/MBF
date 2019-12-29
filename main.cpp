#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <algorithm>
#include <iomanip>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_errno.h>

#include "fitter.h"
#include "parser.h"
#include "abstract_prior.h"
#include "gaussian_prior.h"
#include "abstract_model.h"
#include "chisqr_extra_term.h"

#include "parse_model.h"

#include "multi_alt_exp_Asqr_BC_model.h"
#include "multi_alt_exp_Asqr_expE_BC_model.h"
#include "multi_alt_exp_Asqr_expE_model.h"
#include "multi_alt_exp_Asqr_expE_vec_BC_model.h"
#include "multi_alt_exp_Asqr_expE_vec_model.h"
#include "multi_alt_exp_Asqr_model.h"
#include "multi_alt_exp_Asqr_vec_BC_model.h"
#include "multi_alt_exp_Asqr_vec_model.h"
#include "multi_alt_exp_BC_model.h"
#include "multi_alt_exp_expE_BC_model.h"
#include "multi_alt_exp_expE_mat_model.h"
#include "multi_alt_exp_expE_model.h"
#include "multi_alt_exp_expE_nonsym_mat_model.h"
#include "multi_alt_exp_expE_vec_BC_model.h"
#include "multi_alt_exp_expE_vec_model.h"
#include "multi_alt_exp_mat_model.h"
#include "multi_alt_exp_model.h"
#include "multi_alt_exp_nonsym_mat_model.h"
#include "multi_alt_exp_vec_BC_model.h"
#include "multi_alt_exp_vec_model.h"
#include "multi_exp_Asqr_BC_model.h"
#include "multi_exp_Asqr_expE_BC_model.h"
#include "multi_exp_Asqr_expE_model.h"
#include "multi_exp_Asqr_expE_vec_BC_model.h"
#include "multi_exp_Asqr_expE_vec_model.h"
#include "multi_exp_Asqr_model.h"
#include "multi_exp_Asqr_vec_BC_model.h"
#include "multi_exp_Asqr_vec_model.h"
#include "multi_exp_BC_model.h"
#include "multi_exp_expE_BC_model.h"
#include "multi_exp_expE_mat_model.h"
#include "multi_exp_expE_mat_II_model.h"
#include "multi_exp_expE_mat_upper_model.h"
#include "multi_exp_expE_mat_II_upper_model.h"
#include "multi_exp_expE_model.h"
#include "multi_exp_expE_nonsym_mat_model.h"
#include "multi_exp_expE_vec_BC_model.h"
#include "multi_exp_expE_vec_model.h"
#include "multi_exp_mat_model.h"
#include "multi_exp_mat_II_model.h"
#include "multi_exp_mat_upper_model.h"
#include "multi_exp_mat_II_upper_model.h"
#include "multi_exp_model.h"
#include "multi_exp_nonsym_mat_model.h"
#include "multi_exp_vec_BC_model.h"
#include "multi_exp_vec_model.h"
#include "multi_exp_const_model.h"
#include "multi_exp_expE_const_model.h"
#include "multi_exp_Asqr_const_model.h"
#include "multi_exp_Asqr_expE_const_model.h"
#include "multi_alt_exp_const_model.h"
#include "multi_alt_exp_expE_const_model.h"
#include "multi_alt_exp_Asqr_const_model.h"
#include "multi_alt_exp_Asqr_expE_const_model.h"
#include "multi_exp_vec_const_model.h"
#include "multi_exp_expE_vec_const_model.h"
#include "multi_exp_Asqr_vec_const_model.h"
#include "multi_exp_Asqr_expE_vec_const_model.h"
#include "multi_alt_exp_vec_const_model.h"
#include "multi_alt_exp_expE_vec_const_model.h"
#include "multi_alt_exp_Asqr_vec_const_model.h"
#include "multi_alt_exp_Asqr_expE_vec_const_model.h"
#include "multi_exp_BC_const_model.h"
#include "multi_exp_expE_BC_const_model.h"
#include "multi_exp_Asqr_BC_const_model.h"
#include "multi_exp_Asqr_expE_BC_const_model.h"
#include "multi_alt_exp_BC_const_model.h"
#include "multi_alt_exp_expE_BC_const_model.h"
#include "multi_alt_exp_Asqr_BC_const_model.h"
#include "multi_alt_exp_Asqr_expE_BC_const_model.h"
#include "multi_exp_vec_BC_const_model.h"
#include "multi_exp_expE_vec_BC_const_model.h"
#include "multi_exp_Asqr_vec_BC_const_model.h"
#include "multi_exp_Asqr_expE_vec_BC_const_model.h"
#include "multi_alt_exp_vec_BC_const_model.h"
#include "multi_alt_exp_expE_vec_BC_const_model.h"
#include "multi_alt_exp_Asqr_vec_BC_const_model.h"
#include "multi_alt_exp_Asqr_expE_vec_BC_const_model.h"
#include "threept_multi_alt_exp_expE_model.h"
#include "threept_multi_alt_exp_expE_vec_model.h"
#include "threept_multi_alt_exp_model.h"
#include "threept_multi_alt_exp_vec_model.h"
#include "threept_multi_exp_expE_model.h"
#include "threept_multi_exp_expE_vec_model.h"
#include "threept_multi_exp_model.h"
#include "threept_multi_exp_vec_model.h"

#include "settingsmap.h"

#include "binary_io.h"

#include <limits>

const double empty_double=std::numeric_limits<int>::max()-1;
const int empty_int=std::numeric_limits<int>::max()-1;

const double version=5.63;

const int start_n_functions=1;
const int start_n_variables=1;
const int start_n_parameters=1;
const int start_n_constants=0;
const bool start_bayes=true;
const bool start_num_diff=false;
const bool start_second_deriv_minimization=false;
const bool start_second_deriv_covariance=true;
const double start_num_diff_step=1e-08;
const bool start_use_bse=false;
const bool start_restrict_data=false;
const bool start_bootstrap_normalization=false;
const int start_start_data=1;
const int start_stop_data=1000;
const double start_startlambda=0.001;
const double start_lambdafac=10.0;
const double start_tolerance=0.001;
const int start_svd=0;
const double start_svd_ratio=0.000001;
const double start_svd_value=0.000000000001;
const int start_steps=100;
const int start_bin=1;
const int start_n_exp=1;
const int start_m_exp=1;
const int start_n_vec=2;

const string A_name="A";
const string B_name="B";
const string C_name="C";
const string E_name="E";
const string dE_name="dE";
const string E_initial_name="E";
const string dE_initial_name="dE";
const string E_final_name="F";
const string dE_final_name="dF";
const string t_name="t";
const string T_name="T";

enum model  // keep existing entries for backward compatibility!
{
  PARSE=0,

  MULTIEXP=1,
  MULTIEXPASQR=7,
  MULTIEXPEXPA=8,
  MULTIEXPEXPE=37,
  MULTIEXPEXPAEXPE=38,
  MULTIALTEXP=4,
  MULTIALTEXPASQR=9,
  MULTIALTEXPEXPA=10,
  MULTIALTEXPEXPE=39,
  MULTIALTEXPEXPAEXPE=40,

  MULTIEXPVEC=2,
  MULTIEXPASQRVEC=17,
  MULTIEXPEXPAVEC=18,
  MULTIEXPEXPEVEC=41,
  MULTIEXPEXPAEXPEVEC=42,
  MULTIALTEXPVEC=5,
  MULTIALTEXPASQRVEC=19,
  MULTIALTEXPEXPAVEC=20,
  MULTIALTEXPEXPEVEC=43,
  MULTIALTEXPEXPAEXPEVEC=44,

  MULTIEXPMAT=3,
  MULTIEXPEXPEMAT=27,
  MULTIALTEXPMAT=6,
  MULTIALTEXPEXPEMAT=30,

  MULTIEXPBC=45,
  MULTIEXPASQRBC=46,
  MULTIEXPEXPABC=47,
  MULTIEXPEXPEBC=48,
  MULTIEXPEXPAEXPEBC=49,
  MULTIALTEXPBC=50,
  MULTIALTEXPASQRBC=51,
  MULTIALTEXPEXPABC=52,
  MULTIALTEXPEXPEBC=53,
  MULTIALTEXPEXPAEXPEBC=54,

  MULTIEXPVECBC=55,
  MULTIEXPASQRVECBC=56,
  MULTIEXPEXPAVECBC=57,
  MULTIEXPEXPEVECBC=58,
  MULTIEXPEXPAEXPEVECBC=59,
  MULTIALTEXPVECBC=60,
  MULTIALTEXPASQRVECBC=61,
  MULTIALTEXPEXPAVECBC=62,
  MULTIALTEXPEXPEVECBC=63,
  MULTIALTEXPEXPAEXPEVECBC=64,

  MULTIEXPASQREXPE=101,
  MULTIEXPASQREXPEBC=102,

  MULTIALTEXPASQREXPE=103,
  MULTIALTEXPASQREXPEBC=104,

  MULTIEXPASQREXPEVEC=105,
  MULTIEXPASQREXPEVECBC=106,

  MULTIALTEXPASQREXPEVEC=107,
  MULTIALTEXPASQREXPEVECBC=108,

  THREEPT2VARMULTIEXP=70,
  THREEPT2VARMULTIALTEXP=75,
  THREEPT2VARMULTIEXPVEC=170,
  THREEPT2VARMULTIALTEXPVEC=175,
  THREEPT2VARMULTIEXPEXPE=180,
  THREEPT2VARMULTIALTEXPEXPE=181,
  THREEPT2VARMULTIEXPEXPEVEC=182,
  THREEPT2VARMULTIALTEXPEXPEVEC=183,

  MULTIEXPNONSYMMAT=303,
  MULTIEXPEXPENONSYMMAT=304,
  MULTIALTEXPNONSYMMAT=305,
  MULTIALTEXPEXPENONSYMMAT=306,

  MULTIEXPCONST=400,
  MULTIEXPASQRCONST=401,
  MULTIEXPEXPECONST=402,
  MULTIALTEXPCONST=403,
  MULTIALTEXPASQRCONST=404,
  MULTIALTEXPEXPECONST=405,
  MULTIEXPVECCONST=406,
  MULTIEXPASQRVECCONST=407,
  MULTIEXPEXPEVECCONST=408,
  MULTIALTEXPVECCONST=409,
  MULTIALTEXPASQRVECCONST=410,
  MULTIALTEXPEXPEVECCONST=411,
  MULTIEXPBCCONST=412,
  MULTIEXPASQRBCCONST=413,
  MULTIEXPEXPEBCCONST=414,
  MULTIALTEXPBCCONST=415,
  MULTIALTEXPASQRBCCONST=416,
  MULTIALTEXPEXPEBCCONST=417,
  MULTIEXPVECBCCONST=418,
  MULTIEXPASQRVECBCCONST=419,
  MULTIEXPEXPEVECBCCONST=420,
  MULTIALTEXPVECBCCONST=421,
  MULTIALTEXPASQRVECBCCONST=422,
  MULTIALTEXPEXPEVECBCCONST=423,
  MULTIEXPASQREXPECONST=424,
  MULTIEXPASQREXPEBCCONST=425,
  MULTIALTEXPASQREXPECONST=426,
  MULTIALTEXPASQREXPEBCCONST=427,
  MULTIEXPASQREXPEVECCONST=428,
  MULTIEXPASQREXPEVECBCCONST=429,
  MULTIALTEXPASQREXPEVECCONST=430,
  MULTIALTEXPASQREXPEVECBCCONST=431,

  MULTIEXPMATII=500,
  MULTIEXPEXPEMATII=501,

  MULTIEXPMATIIUPPER=600,
  MULTIEXPEXPEMATIIUPPER=601,
  MULTIEXPMATUPPER=603,
  MULTIEXPEXPEMATUPPER=627

};

bool _in(const string& sub, const string& s)
{
  if(s.find(sub)!=string::npos)
  {
    return true;
  }
  else
  {
    return false;
  }
}


int n_functions;
int n_variables;
int n_parameters;
int n_constants;

model current_model;

fitter* _fitter;
gaussian_prior* _gaussian_prior;
abstract_model* _model;

chisqr_extra_term* _chisqr_extra_term;
bool chisqr_extra_term_enabled;
string chisqr_extra_term_function;
double chisqr_extra_term_num_diff_step;
int chisqr_extra_term_n_constants;
vector< string > chisqr_extra_term_constant_names;
vector< double > chisqr_extra_term_constant_values;

vector< string > functions;
vector< string > variables;
vector< string > parameters;
vector< string > constants;
vector< vector< string > > derivatives;

vector< string > fit_min;
vector< string > fit_max;

vector< double > start_values;
vector< double > priors;
vector< double > sigmas;
vector< double > constant_values;

vector< vector< double > > file_arguments;
vector< vector< vector< double > > > file_data;

vector< vector< double > > fit_arguments;
vector< vector< vector< double > > > fit_data;

bool bayesian;

bool num_diff;
bool second_deriv_minimization;
bool second_deriv_covariance;
double num_diff_step;

bool use_bse;
bool restrict_data;

bool bootstrap_normalization;

int start_n_data;
int stop_n_data;

double start_lambda;
double lambda_factor;
double chisqr_tolerance;
int svd_cut;
double svd_ratio;
double svd_value;
inversion_method inv_method;
int max_steps;
int bin_size;
bool boot_prior;
int bssamples;
string bse_file;
int data_file_type;

string data_file;
string output_dir;

int multiexpdialog_nexp;
bool multiexpdialog_BC;
bool multiexpdialog_constant;
int multiexpvecdialog_nexp;
int multiexpvecdialog_nvec;
bool multiexpvecdialog_BC;
bool multiexpvecdialog_constant;
int multiexpmatdialog_nexp;
int multiexpmatdialog_dim_1;
int multiexpmatdialog_dim_2;
int multiexpmatupperdialog_nexp;
int multiexpmatupperdialog_dim;
int multialtexpdialog_nexp;
int multialtexpdialog_mexp;
bool multialtexpdialog_BC;
bool multialtexpdialog_constant;
int multialtexpvecdialog_nexp;
int multialtexpvecdialog_mexp;
int multialtexpvecdialog_nvec;
bool multialtexpvecdialog_BC;
bool multialtexpvecdialog_constant;
int multialtexpmatdialog_nexp;
int multialtexpmatdialog_mexp;
int multialtexpmatdialog_dim_1;
int multialtexpmatdialog_dim_2;
int threeptmultiexpdialog_nexpinitial;
int threeptmultiexpdialog_nexpfinal;
int threeptmultiexpvecdialog_nexpinitial;
int threeptmultiexpvecdialog_nexpfinal;
int threeptmultiexpvecdialog_nvec;
int threeptmultialtexpdialog_nexpinitial;
int threeptmultialtexpdialog_nexpfinal;
int threeptmultialtexpdialog_mexpinitial;
int threeptmultialtexpdialog_mexpfinal;
int threeptmultialtexpvecdialog_nexpinitial;
int threeptmultialtexpvecdialog_nexpfinal;
int threeptmultialtexpvecdialog_mexpinitial;
int threeptmultialtexpvecdialog_mexpfinal;
int threeptmultialtexpvecdialog_nvec;

int n_parameters_dof;

void init();
bool loadFile(const string& fileName);
bool prepare_fit();
bool fit();
bool bootstrap(string mbf_file_name);
bool reset_fitter();
bool load_data_file();
bool set_fit_data();
string stripfilename(string fullname);
void print_usage();

void init()
{
  current_model=PARSE;
  _fitter=NULL;
  _gaussian_prior=NULL;
  _model=NULL;
  _chisqr_extra_term=NULL;
  chisqr_extra_term_enabled=false;
  chisqr_extra_term_function="";
  chisqr_extra_term_num_diff_step=start_num_diff_step;
  chisqr_extra_term_n_constants=0;
  bayesian=start_bayes;
  num_diff=start_num_diff;
  second_deriv_minimization=start_second_deriv_minimization;
  second_deriv_covariance=start_second_deriv_covariance;
  num_diff_step=start_num_diff_step;
  use_bse=start_use_bse;
  restrict_data=start_restrict_data;
  bootstrap_normalization=start_bootstrap_normalization;
  start_n_data=start_start_data;
  stop_n_data=start_stop_data;
  start_lambda=start_startlambda;
  lambda_factor=start_lambdafac;
  chisqr_tolerance=start_tolerance;
  svd_cut=start_svd;
  svd_ratio=start_svd_ratio;
  svd_value=start_svd_value;
  max_steps=start_steps;
  bin_size=start_bin;
  bse_file="";
  boot_prior=true;
  bssamples=500;
  data_file_type=0;
  n_parameters_dof=0;
  multiexpdialog_nexp=start_n_exp;
  multiexpdialog_BC=false;
  multiexpdialog_constant=false;
  multiexpvecdialog_nexp=start_n_exp;
  multiexpvecdialog_nvec=start_n_vec;
  multiexpvecdialog_BC=false;
  multiexpvecdialog_constant=false;
  multiexpmatdialog_nexp=start_n_exp;
  multiexpmatdialog_dim_1=start_n_vec;
  multiexpmatdialog_dim_2=start_n_vec;
  multiexpmatupperdialog_nexp=start_n_exp;
  multiexpmatupperdialog_dim=start_n_vec;
  multialtexpdialog_nexp=start_n_exp;
  multialtexpdialog_mexp=start_m_exp;
  multialtexpdialog_BC=false;
  multialtexpdialog_constant=false;
  multialtexpvecdialog_nexp=start_n_exp;
  multialtexpvecdialog_mexp=start_m_exp;
  multialtexpvecdialog_nvec=start_n_vec;
  multialtexpvecdialog_BC=false;
  multialtexpvecdialog_constant=false;
  multialtexpmatdialog_nexp=start_n_exp;
  multialtexpmatdialog_mexp=start_m_exp;
  multialtexpmatdialog_dim_1=start_n_vec;
  multialtexpmatdialog_dim_2=start_n_vec;
  threeptmultiexpdialog_nexpinitial=start_n_exp;
  threeptmultiexpdialog_nexpfinal=start_n_exp;
  threeptmultiexpvecdialog_nexpinitial=start_n_exp;
  threeptmultiexpvecdialog_nexpfinal=start_n_exp;
  threeptmultiexpvecdialog_nvec=start_n_vec;
  threeptmultialtexpdialog_nexpinitial=start_n_exp;
  threeptmultialtexpdialog_nexpfinal=start_n_exp;
  threeptmultialtexpdialog_mexpinitial=start_m_exp;
  threeptmultialtexpdialog_mexpfinal=start_m_exp;
  threeptmultialtexpvecdialog_nexpinitial=start_n_exp;
  threeptmultialtexpvecdialog_nexpfinal=start_n_exp;
  threeptmultialtexpvecdialog_mexpinitial=start_m_exp;
  threeptmultialtexpvecdialog_mexpfinal=start_m_exp;
  threeptmultialtexpvecdialog_mexpfinal=start_n_vec;
}


bool loadFile(const string& fileName)
{
  settingsmap m;
  if(!m.load_file(fileName))
  {
    cerr << "cannot read file " << fileName << endl;
    return false;
  }

  if(m.exists("version"))
  {
    double v=m.get_double("version");
    if(v>version)
    {
      cerr << "WARNING: input file was saved with a newer QMBF version: " << v << endl << endl;
    }
  }


  if(!m.exists("nfunctions"))
  {
    cerr << "input file does not contain number of functions" << endl;
    return false;
  }
  if(!m.exists("nvariables"))
  {
    cerr << "input file does not contain number of variables" << endl;
    return false;
  }
  if(!m.exists("nparameters"))
  {
    cerr << "input file does not contain number of parameters" << endl;
    return false;
  }
  if(!m.exists("nconstants"))
  {
    cerr << "input file does not contain number of constants" << endl;
    return false;
  }
  n_functions=m.get_int("nfunctions");
  functions.resize(n_functions);
  n_variables=m.get_int("nvariables");
  variables.resize(n_variables);
  fit_min.resize(n_variables);
  fit_max.resize(n_variables);
  n_parameters=m.get_int("nparameters");
  parameters.resize(n_parameters);
  start_values.resize(n_parameters);
  priors.resize(n_parameters);
  sigmas.resize(n_parameters);
  vector< string > der_temp_vec;
  der_temp_vec.resize(n_parameters);
  derivatives.resize(n_functions, der_temp_vec);
  n_constants=m.get_int("nconstants");
  constants.resize(n_constants);
  constant_values.resize(n_constants);
  if(m.exists("currentmodel"))
  {
    current_model=static_cast<model>(m.get_int("currentmodel"));
  }
  for(int f=0; f<n_functions; ++f)
  {
    functions[f]=m.get_string("function", f);
  }
  for(int v=0; v<n_variables; ++v)
  {
    variables[v]=m.get_string("variable", v);
    fit_min[v]=m.get_string("fitmin", v);
    fit_max[v]=m.get_string("fitmax", v);
  }
  for(int p=0; p<n_parameters; ++p)
  {
    parameters[p]=m.get_string("parameter", p);
    if(m.exists("startval", p))
    {
      double start_val=m.get_double("startval", p);
      if(start_val!=empty_double)
      {
       start_values[p]=start_val;
      }
    }
    if(m.exists("prior", p))
    {
      double prior=m.get_double("prior", p);
      if(prior!=empty_double)
      {
        priors[p]=prior;
      }
    }
    if(m.exists("sigma", p))
    {
      double sigma=m.get_double("sigma", p);
      if(sigma!=empty_double)
      {
        sigmas[p]=sigma;
      }
    }
  }
  for(int c=0; c<n_constants; ++c)
  {
    constants[c]=m.get_string("constant", c);
    if(m.exists("constantvalue", c))
    {
      double val=m.get_double("constantvalue", c);
      if(val!=empty_double)
      {
        constant_values[c]=val;
      }
    }
  }
  for(int f=0; f<n_functions; ++f)
  {
    for(int p=0; p<n_parameters; ++p)
    {
      derivatives[f][p]=m.get_string("derivative", f, p);
    }
  }
  if(m.exists("bayesian"))
  {
    bayesian=m.get_bool("bayesian");
  }
  if(m.exists("numdiff")) num_diff=m.get_bool("numdiff");
  if(m.exists("secondderivminimization")) second_deriv_minimization=m.get_bool("secondderivminimization");
  if(m.exists("secondderivcovariance")) second_deriv_covariance=m.get_bool("secondderivcovariance");

  if(m.exists("numdiffstep"))
  {
    double numdiffstep=m.get_double("numdiffstep");
    if(numdiffstep!=empty_double)
    {
      num_diff_step=numdiffstep;
    }
  }

  if(m.exists("bootprior"))
  {
    boot_prior=m.get_bool("bootprior");
  }
  else
  {
    boot_prior=true;
  }

  if(m.exists("usebse"))
  {
    use_bse=m.get_bool("usebse");
  }
  else
  {
    if(m.exists("bsefile"))
    {
      use_bse=true;
    }
    else
    {
      use_bse=false;
    }
  }

  if(m.exists("bsefile"))
  {
    bse_file=m.get_string("bsefile");
  }

  if(m.exists("bssamples"))
  {
    int _bssamples=m.get_int("bssamples");
    if(_bssamples!=empty_int)
    {
      bssamples=_bssamples;
    }
  }

  if(m.exists("restrictdatarange")) restrict_data=m.get_bool("restrictdatarange");
  if(m.exists("datarangemin")) start_n_data=m.get_int("datarangemin");
  if(m.exists("datarangemax")) stop_n_data=m.get_int("datarangemax");

  if(m.exists("startlambda"))
  {
    double startlambda=m.get_double("startlambda");
    if(startlambda!=empty_double)
    {
      start_lambda=startlambda;
    }
  }
  if(m.exists("lambdafactor"))
  {
    double lambdafactor=m.get_double("lambdafactor");
    if(lambdafactor!=empty_double)
    {
      lambda_factor=lambdafactor;
    }
  }
  if(m.exists("tolerance"))
  {
    double tolerance=m.get_double("tolerance");
    if(tolerance!=empty_double)
    {
      chisqr_tolerance=tolerance;
    }
  }
  if(m.exists("svdcut"))
  {
    int svdcut=m.get_int("svdcut");
    if(svdcut!=empty_int)
    {
      svd_cut=svdcut;
    }
  }
  if(m.exists("svdratio"))
  {
    double ratio=m.get_double("svdratio");
    if(ratio!=empty_double)
    {
      svd_ratio=ratio;
    }
  }
  if(m.exists("svdvalue"))
  {
    double value=m.get_double("svdvalue");
    if(value!=empty_double)
    {
      svd_value=value;
    }
  }
  if(m.exists("inversionmethod"))
  {
    inv_method=static_cast<inversion_method>(m.get_int("inversionmethod"));
  }
  else
  {
    if(svd_cut>0)
    {
      inv_method=simple_cut;
    }
    else
    {
      inv_method=LU_inversion;
    }
  }
  if(m.exists("maxsteps"))
  {
    int maxsteps=m.get_int("maxsteps");
    if(maxsteps!=empty_int)
    {
      max_steps=maxsteps;
    }
  }
  if(m.exists("binsize"))
  {
    int binsize=m.get_int("binsize");
    if(binsize!=empty_int)
    {
      bin_size=binsize;
    }
  }
  if(m.exists("bootstrapnormalization"))
  {
    bootstrap_normalization=m.get_bool("bootstrapnormalization");
    if(bootstrap_normalization)
    {
      bin_size=1;
    }
  }

  if(m.exists("datafiletype"))
  {
    data_file_type=m.get_int("datafiletype");
  }
  data_file=m.get_string("datafile");
  output_dir=m.get_string("outputdir");
  if(*(output_dir.end()-1)!='/')
  {
    output_dir.append("/");
  }
  if(m.exists("multiexpdialog:nexp")) multiexpdialog_nexp=m.get_int("multiexpdialog:nexp");
  if(m.exists("multiexpdialog:BC")) multiexpdialog_BC=m.get_bool("multiexpdialog:BC");
  if(m.exists("multiexpdialog:constant")) multiexpdialog_BC=m.get_bool("multiexpdialog:constant");
  if(m.exists("multiexpvecdialog:nexp"))
  {
    multiexpvecdialog_nexp=m.get_int("multiexpvecdialog:nexp");
    multiexpmatdialog_nexp=m.get_int("multiexpvecdialog:nexp");
  }
  if(m.exists("multiexpvecdialog:nvec"))
  {
    multiexpvecdialog_nvec=m.get_int("multiexpvecdialog:nvec");
    multiexpmatdialog_dim_1=m.get_int("multiexpvecdialog:nvec");
    multiexpmatdialog_dim_2=m.get_int("multiexpvecdialog:nvec");
  }
  if(m.exists("multiexpvecdialog:BC")) multiexpvecdialog_BC=m.get_bool("multiexpvecdialog:BC");
  if(m.exists("multiexpvecdialog:constant")) multiexpvecdialog_BC=m.get_bool("multiexpvecdialog:constant");

  if(m.exists("multiexpmatdialog:nexp")) multiexpmatdialog_nexp=m.get_int("multiexpmatdialog:nexp");
  if(m.exists("multiexpmatdialog:dim1")) multiexpmatdialog_dim_1=m.get_int("multiexpmatdialog:dim1");
  if(m.exists("multiexpmatdialog:dim2")) multiexpmatdialog_dim_2=m.get_int("multiexpmatdialog:dim2");

  if(m.exists("multiexpmatupperdialog:nexp")) multiexpmatupperdialog_nexp=m.get_int("multiexpmatupperdialog:nexp");
  if(m.exists("multiexpmatupperdialog:dim"))  multiexpmatupperdialog_dim=m.get_int("multiexpmatupperdialog:dim");

  if(m.exists("multialtexpdialog:nexp")) multialtexpdialog_nexp=m.get_int("multialtexpdialog:nexp");
  if(m.exists("multialtexpdialog:mexp")) multialtexpdialog_mexp=m.get_int("multialtexpdialog:mexp");
  if(m.exists("multialtexpdialog:BC")) multialtexpdialog_BC=m.get_bool("multialtexpdialog:BC");
  if(m.exists("multialtexpdialog:constant")) multialtexpdialog_BC=m.get_bool("multialtexpdialog:constant");
  if(m.exists("multialtexpvecdialog:nexp"))
  {
    multialtexpvecdialog_nexp=m.get_int("multialtexpvecdialog:nexp");
    multialtexpmatdialog_nexp=m.get_int("multialtexpvecdialog:nexp");
  }
  if(m.exists("multialtexpvecdialog:mexp"))
  {
    multialtexpvecdialog_mexp=m.get_int("multialtexpvecdialog:mexp");
    multialtexpmatdialog_mexp=m.get_int("multialtexpvecdialog:mexp");
  }
  if(m.exists("multialtexpvecdialog:nvec"))
  {
    multialtexpvecdialog_nvec=m.get_int("multialtexpvecdialog:nvec");
    multialtexpmatdialog_dim_1=m.get_int("multialtexpvecdialog:nvec");
    multialtexpmatdialog_dim_2=m.get_int("multialtexpvecdialog:nvec");
  }
  if(m.exists("multialtexpvecdialog:BC")) multialtexpvecdialog_BC=m.get_bool("multialtexpvecdialog:BC");
  if(m.exists("multialtexpvecdialog:constant")) multialtexpvecdialog_BC=m.get_bool("multialtexpvecdialog:constant");

  if(m.exists("multialtexpmatdialog:nexp")) multialtexpmatdialog_nexp=m.get_int("multialtexpmatdialog:nexp");
  if(m.exists("multialtexpmatdialog:mexp")) multialtexpmatdialog_mexp=m.get_int("multialtexpmatdialog:mexp");
  if(m.exists("multialtexpmatdialog:dim1")) multialtexpmatdialog_dim_1=m.get_int("multialtexpmatdialog:dim1");
  if(m.exists("multialtexpmatdialog:dim2")) multialtexpmatdialog_dim_2=m.get_int("multialtexpmatdialog:dim2");

  if(m.exists("threeptmultiexpdialog:nexpinitial")) threeptmultiexpdialog_nexpinitial=m.get_int("threeptmultiexpdialog:nexpinitial");
  if(m.exists("threeptmultiexpdialog:nexpfinal")) threeptmultiexpdialog_nexpfinal=m.get_int("threeptmultiexpdialog:nexpfinal");

  if(m.exists("threeptmultiexpvecdialog:nexpinitial")) threeptmultiexpvecdialog_nexpinitial=m.get_int("threeptmultiexpvecdialog:nexpinitial");
  if(m.exists("threeptmultiexpvecdialog:nexpfinal")) threeptmultiexpvecdialog_nexpfinal=m.get_int("threeptmultiexpvecdialog:nexpfinal");
  if(m.exists("threeptmultiexpvecdialog:nvec")) threeptmultiexpvecdialog_nvec=m.get_int("threeptmultiexpvecdialog:nvec");

  if(m.exists("threeptmultialtexpdialog:nexpinitial")) threeptmultialtexpdialog_nexpinitial=m.get_int("threeptmultialtexpdialog:nexpinitial");
  if(m.exists("threeptmultialtexpdialog:nexpfinal")) threeptmultialtexpdialog_nexpfinal=m.get_int("threeptmultialtexpdialog:nexpfinal");
  if(m.exists("threeptmultialtexpdialog:mexpinitial")) threeptmultialtexpdialog_mexpinitial=m.get_int("threeptmultialtexpdialog:mexpinitial");
  if(m.exists("threeptmultialtexpdialog:mexpfinal")) threeptmultialtexpdialog_mexpfinal=m.get_int("threeptmultialtexpdialog:mexpfinal");

  if(m.exists("threeptmultialtexpvecdialog:nexpinitial")) threeptmultialtexpvecdialog_nexpinitial=m.get_int("threeptmultialtexpvecdialog:nexpinitial");
  if(m.exists("threeptmultialtexpvecdialog:nexpfinal")) threeptmultialtexpvecdialog_nexpfinal=m.get_int("threeptmultialtexpvecdialog:nexpfinal");
  if(m.exists("threeptmultialtexpvecdialog:mexpinitial")) threeptmultialtexpvecdialog_mexpinitial=m.get_int("threeptmultialtexpvecdialog:mexpinitial");
  if(m.exists("threeptmultialtexpvecdialog:mexpfinal")) threeptmultialtexpvecdialog_mexpfinal=m.get_int("threeptmultialtexpvecdialog:mexpfinal");
  if(m.exists("threeptmultialtexpvecdialog:nvec")) threeptmultialtexpvecdialog_nvec=m.get_int("threeptmultialtexpvecdialog:nvec");

  if(m.exists("nparametersdof"))
  {
    n_parameters_dof=m.get_int("nparametersdof");
  }
  else
  {
    if(bayesian)
    {
      n_parameters_dof=0;
    }
    else
    {
      n_parameters_dof=n_parameters;
    }
  }


  if(m.exists("chisqrextraterm:enabled"))
  {
    chisqr_extra_term_enabled=m.get_bool("chisqrextraterm:enabled");
    if(m.exists("chisqrextraterm:function")) chisqr_extra_term_function=m.get_string("chisqrextraterm:function");
    if(m.exists("chisqrextraterm:numdiffstep")) chisqr_extra_term_num_diff_step=m.get_double("chisqrextraterm:numdiffstep");
    if(!m.exists("chisqrextraterm:nconstants"))
    {
      cerr << "input file does not contain number of chi^2 extra term constants" << endl;
      return false;
    }
    chisqr_extra_term_n_constants=m.get_int("chisqrextraterm:nconstants");
    chisqr_extra_term_constant_names.resize(chisqr_extra_term_n_constants);
    chisqr_extra_term_constant_values.resize(chisqr_extra_term_n_constants);
    for(int c=0; c<chisqr_extra_term_n_constants; ++c)
    {
      chisqr_extra_term_constant_names[c]=m.get_string("chisqrextraterm:constant", c);
      if(m.exists("chisqrextraterm:constantvalue", c))
      {
        double val=m.get_double("chisqrextraterm:constantvalue", c);
        if(val!=empty_double)
        {
          chisqr_extra_term_constant_values[c]=val;
        }
      }
    }
  }

  return true;
}



bool prepare_fit()
{
  if(!reset_fitter())
  {
    return false;
  }
  for(int f=0; f<n_functions; ++f)
  {
    functions[f]=_model->get_function_name(f);
  }
  for(int p=0; p<n_parameters; ++p)
  {
    parameters[p]=_model->get_parameter_name(p);
  }
  for(int v=0; v<n_variables; ++v)
  {
    variables[v]=_model->get_variable_name(v);
  }
  for(int c=0; c<n_constants; ++c)
  {
    constants[c]=_model->get_constant_name(c);
  }
  _gaussian_prior->set_enabled(bayesian);
  _fitter->set_num_diff(num_diff);
  _fitter->set_second_deriv_minimization(second_deriv_minimization);
  _fitter->set_second_deriv_covariance(second_deriv_covariance);
  _fitter->set_num_diff_step(num_diff_step);

  _fitter->set_initial_lambda(start_lambda);
  _fitter->set_lambda_factor(lambda_factor);
  _fitter->set_tolerance(chisqr_tolerance);
  _model->set_constants(constant_values);
  for(int p=0; p<n_parameters; ++p)
  {
    _fitter->set_starting_value(p, start_values[p]);
  }
  if(bayesian)
  {
    for(int p=0; p<n_parameters; ++p)
    {
      _gaussian_prior->set_prior(p, priors[p]);
      _gaussian_prior->set_sigma(p, sigmas[p]);
    }
  }

  _chisqr_extra_term->set_enabled(chisqr_extra_term_enabled);
  if(chisqr_extra_term_enabled)
  {
    _chisqr_extra_term->set_num_diff_step(chisqr_extra_term_num_diff_step);
    _chisqr_extra_term->set_constants(chisqr_extra_term_constant_values);
    _chisqr_extra_term->set_parameters(start_values);
    _chisqr_extra_term->chi_sqr();
    if(_chisqr_extra_term->no_of_errors()!=0)
    {
      cerr << "Error: Error while evaluating chi^2 extra term function" << endl;
      return false;
    }
  }

  if(!load_data_file())
  {
    return false;
  }
  if(!set_fit_data())
  {
    return false;
  }
  if(!(_fitter->test_model()))
  {
    cerr << "Error in fit model" << endl;
    return false;
  }
  return true;
}


bool fit(string mbf_file_name)
{
  if(!prepare_fit())
  {
    return false;
  }


// begin output report
  cout.setf(ios::left,ios::adjustfield);
  cout << endl << "Data file = " << data_file << endl << endl;

  if(restrict_data)
  {
    cout << "Data range min = " << start_n_data << endl;
    cout << "Data range max = " << stop_n_data << endl << endl;
  }

  if(chisqr_extra_term_enabled)
  {
    cout << "enabled chi^2 additional term = " << chisqr_extra_term_function << endl;
    for(int c=0; c<chisqr_extra_term_n_constants; ++c)
    {
      cout << "   constant  " << chisqr_extra_term_constant_names[c]
           << " = " << chisqr_extra_term_constant_values[c] << endl;
    }
    cout << endl;
  }

  for(int f=0; f<n_functions; ++f)
  {
    cout << "Function " << f+1 << ":" << endl << functions[f] << endl << endl;
  }
  for(int v=0; v<n_variables; ++v)
  {
    cout << "Variable  " << variables[v]
         << "  Fitting min = " << fit_min[v]
         << "  Fitting max = " << fit_max[v] << endl << endl;
  }
  for(int c=0; c<n_constants; ++c)
  {
    cout << "Constant  " << constants[c]
         << " = " << constant_values[c] << endl << endl;
  }
  cout << "Fit algorithm settings:" << endl << endl;
  cout << "Start lambda                 = " << start_lambda << endl;
  cout << "Lambda factor                = " << lambda_factor << endl;
  cout << "Delta(chi^2) tolerance       = " << chisqr_tolerance << endl;

  cout << "Normalization of correlation matrix: ";
  if(bootstrap_normalization)
  {
    cout << "1/(N-1)" << endl;
  }
  else
  {
    cout << "1/(N*(N-1))" << endl;
  }


  cout << "(Pseudo-)Inversion method for correlation matrix: ";
  switch(inv_method)
  {
    case diagonal:
      cout << "Diagonal only (uncorrelated fit)" << endl;
      break;
    case LU_inversion:
      cout << "LU decomposition" << endl;
      break;
    case simple_cut:
      cout << "SVD with fixed cut" << endl;
      cout << "Fixed SVD cut: " << svd_cut << endl;
      break;
    case ratio_cut:
      cout << "SVD with EV ratio cut" << endl;
      cout << "SVD EV ratio cut: " <<svd_ratio << endl;
      break;
    case absolute_cut:
      cout << "SVD with EV value cut" << endl;
      cout << "SVD EV value cut: " << svd_value << endl;
      break;
    default:
      break;
  }
  cout << "Maximum number of iterations = " << max_steps << endl;
  cout << "Bin size                     = " << bin_size << endl;
  if(boot_prior && bayesian)
  {
    cout << "Gaussian random priors for bootstrap/multifit activated." << endl << endl;
  }
// end output report


  cout << endl << "Fit started." << endl << endl;


  if(num_diff)
  {
    cout << "Using numerical differentiation for first derivatives." << endl;
  }
  else
  {
    cout << "Using symbolic differentiation for first derivatives." << endl;
  }

  if(second_deriv_minimization)
  {
    cout << "Using second derivatives for minimization." << endl;
  }

  if(second_deriv_covariance)
  {
    cout << "Using second derivatives for parameter covariance matrix." << endl;
  }

  if(num_diff || second_deriv_minimization || second_deriv_covariance)
  {
    cout << "Numerical differentiation step size = " << num_diff_step << endl << endl;
  }

  string fit_message="";
  int steps_needed=_fitter->fit(max_steps, fit_message);
  if(fit_message!="")
  {
    cout << endl << fit_message << endl << endl;
  }


  if(steps_needed==max_steps+1)
  {
    cout << "Fit did not converge after " << max_steps << " iterations" << endl << endl << endl;
    cout << "Parameter   Start value     Prior value     Prior sigma" << endl;
    cout << "-------------------------------------------------------" << endl;
    for(int p=0; p<n_parameters; ++p)
    {
      cout << setw(12) << parameters[p];
      cout << setw(16) << start_values[p];
      if(bayesian)
      {
        cout << setw(16) << priors[p]
             << setw(16) << sigmas[p];
      }
      else
      {
        cout << "                " << "                ";
      }
      cout << endl;
    }
    return false;
  }
  else
  {
    cout << "Fit converged after " << steps_needed << " iterations." << endl << endl << endl;
    cout << "Parameter   Start value     Prior value     Prior sigma     Fit result      sqrt(cov)" << endl;
    cout << "-------------------------------------------------------------------------------------" << endl;
    for(int p=0; p<n_parameters; ++p)
    {
      cout << setw(12) << parameters[p];
      cout << setw(16) << start_values[p];
      if(bayesian)
      {
        cout << setw(16) << priors[p]
             << setw(16) << sigmas[p];
      }
      else
      {
        cout << "                " << "                ";
      }
      cout << setw(16) << _fitter->get_parameter(p);
      cout << setw(16) << sqrt(_fitter->get_covariance(p, p));
      cout << endl;
    }
    cout << endl;
    double chisqr=_fitter->get_chi_sqr();
    double dof=_fitter->get_dof()-_fitter->get_cut()-n_parameters_dof;
    if(dof<0.0)
    {
      dof=0.0;
    }
    cout << "chi^2/dof        = " << chisqr/dof << endl;
    cout << "Q(dof/2,chi^2/2) = " << gsl_sf_gamma_inc_Q(dof/2.0, chisqr/2.0) << endl;
    cout << "Effective dof    = " << dof << endl << endl;

    if( !(_fitter->covariance_positive_definite()) )
    {
      cout << endl << "WARNING: Hessian of chi^2 is not positive definite" << endl;
    }

// write covariance to file
    string cov_filename=output_dir+mbf_file_name+"_covariance.dat";
    ofstream cov_output(cov_filename.c_str());
    if(!cov_output)
    {
      cerr << "Error: could not write to file \"" << cov_filename << "\"" << endl << endl;
      return true;
    }
    cov_output.precision(14);
    cov_output.setf(ios::left);
    int max_length=1;
    for(int p1=0; p1<n_parameters; ++p1)
    {
      if(parameters[p1].length()>max_length)
      {
        max_length=parameters[p1].length();
      }
    }
    for(int p1=0; p1<n_parameters; ++p1)
    {
      for(int p2=0; p2<n_parameters; ++p2)
      {
        cov_output << setw(max_length+4) << parameters[p1] << setw(max_length+4) << parameters[p2] << setw(0) << _fitter->get_covariance(p1, p2) << endl;
      }
    }
    cov_output.close();

// write results to file
    string results_filename=output_dir+mbf_file_name+"_results.dat";;
    ofstream results_output(results_filename.c_str());
    if(!results_output)
    {
      cerr << "Error: could not write to file \"" << results_filename << "\"" << endl << endl;
      return true;
    }
    results_output.precision(14);
    results_output.setf(ios::left);
    max_length=1;
    for(int p1=0; p1<n_parameters; ++p1)
    {
      if(parameters[p1].length()>max_length)
      {
        max_length=parameters[p1].length();
      }
    }
    for(int p1=0; p1<n_parameters; ++p1)
    {
      results_output << setw(max_length+4) << parameters[p1] << setw(0) << _fitter->get_parameter(p1) << endl;
    }
    results_output.close();
  }
  return true;
}

bool bootstrap(string mbf_file_name)
{
  if(!fit(mbf_file_name))
  {
    return false;
  }

  vector< double > boot_data_temp_1(n_functions, 0.0);
  vector< vector< double > > boot_data_temp_2(fit_arguments.size(), boot_data_temp_1);
  vector< vector< vector< double > > > boot_data(fit_data.size(), boot_data_temp_2);
  gsl_rng_default_seed=0;
  gsl_rng_env_setup();
  const gsl_rng_type* Trng=gsl_rng_mt19937;
  gsl_rng* rng=gsl_rng_alloc(Trng);

  vector< int > bootconfig_temp(fit_data.size(), 0);
  vector< vector< int > > bootconfig;


  if(use_bse)
  {
    if(bse_file.empty())
    {
      cerr << "Error: no bootstrap ensemble file selected" << endl;
      return false;
    }
    ifstream bse(bse_file.c_str());
    if(!bse)
    {
      cerr << "Error: could not open bootstrap ensemble file" << endl;
      return false;
    }
    bse >> bssamples;
    bootconfig.resize(bssamples, bootconfig_temp);
    int bsdata=0;
    bse >> bsdata;
    if(bsdata!=fit_data.size())
    {
      cerr << "Error: bootstrap ensemble file has wrong number of data sets" << endl;
      return false;
    }
    for(int boot=0; boot<bssamples; ++boot)
    {
      for(unsigned int b=0; b<fit_data.size(); ++b)
      {
        bse >> bootconfig[boot][b];
        bootconfig[boot][b]-=1;
      }
    }
  }
  else
  {
    bootconfig.resize(bssamples, bootconfig_temp);
    for(int boot=0; boot<bssamples; ++boot)
    {
      for(unsigned int b=0; b<fit_data.size(); ++b)
      {
        bootconfig[boot][b]=static_cast<int>(gsl_rng_uniform(rng)*fit_data.size());
      }
    }
  }

  cout << endl << "Performing bootstrap with " << bssamples << " ensembles..." << endl;
  vector< ofstream* > all_output;
  all_output.resize(n_parameters);

  for(int p=0; p<n_parameters; ++p)
  {
    string filename;
    filename=output_dir+mbf_file_name+"_bootstrap_"+parameters[p]+".dat";
    ofstream* output = new ofstream(filename.c_str());
    if(!(*output))
    {
      cerr << "Cannot write file "  << filename << endl;
      return false;
    }
    output -> precision(15);
    all_output[p]=output;
  }


  for(int boot=0; boot<bssamples; ++boot)
  {
    for(unsigned int b=0; b<fit_data.size(); ++b)
    {
      int r=bootconfig[boot][b];
      boot_data[b]=fit_data[r];
    }
    _fitter->set_inversion_method(inv_method);
    _fitter->set_bootstrap_normalization(bootstrap_normalization);
    _fitter->set_svd_cut(svd_cut);
    _fitter->set_svd_cut_ratio(svd_ratio);
    _fitter->set_svd_cut_absolute(svd_value);
    string set_data_message;
    _fitter->set_data(fit_arguments, boot_data, set_data_message);
    if(boot_prior && bayesian)
    {
      for(int p=0; p<n_parameters; ++p)
      {
        _gaussian_prior->set_prior(p, priors[p]+gsl_ran_gaussian(rng, sigmas[p]));
      }
    }

    string fit_message;
    int steps_needed=_fitter->fit(max_steps, fit_message);
    double chisqr=_fitter->get_chi_sqr();
    double dof=_fitter->get_dof()-_fitter->get_cut()-n_parameters_dof;
    if(dof<0.0)
    {
      dof=0.0;
    }
    if(steps_needed==max_steps+1)
    {
      cerr << "Bootstrap warning: Fit #" << boot+1 << " did not converge after " << max_steps << " iterations. " << "chi^2/dof = " << chisqr/dof << endl;
    }
    else
    {
      cout << "Bootstrap fit #" << boot+1 << " converged after " << steps_needed << " iterations. " << "chi^2/dof = " << chisqr/dof << endl;
    }
    for(int p=0; p<n_parameters; ++p)
    {
      *all_output[p] << _fitter->get_parameter(p) << endl;
    }
  }
  gsl_rng_free(rng);
  cout << "Bootstrap with " << bssamples << " ensembles completed." << endl << endl;
  for(int p=0; p<n_parameters; ++p)
  {
    all_output[p] -> close();
    delete all_output[p];
  }
  if(bootstrap_normalization)
  {
    cout << "Warning: data correlation matrix normalization is set to 1/(N-1). Maybe you should use multifit instad of bootstrap." << endl << endl;
  }

  return true;
}




bool multifit(string mbf_file_name)
{
  if(!fit(mbf_file_name))
  {
    return false;
  }

  gsl_rng_default_seed=0;
  gsl_rng_env_setup();
  const gsl_rng_type* Trng=gsl_rng_mt19937;
  gsl_rng* rng=gsl_rng_alloc(Trng);

  cout << endl << "Performing multifit with " << fit_data.size() << " samples..." << endl;
  vector< ofstream* > all_output;
  all_output.resize(n_parameters);

  for(int p=0; p<n_parameters; ++p)
  {
    string filename;
    filename=output_dir+mbf_file_name+"_multifit_"+parameters[p]+".dat";
    ofstream* output = new ofstream(filename.c_str());
    if(!(*output))
    {
      cerr << "Cannot write file "  << filename << endl;
      return false;
    }
    output -> precision(15);
    all_output[p]=output;
  }


  for(int n=0; n<fit_data.size(); ++n)
  {
    _fitter->set_average_data(fit_data[n]);

    if(boot_prior && bayesian)
    {
      for(int p=0; p<n_parameters; ++p)
      {
        _gaussian_prior->set_prior(p, priors[p]+gsl_ran_gaussian(rng, sigmas[p]));
      }
    }

    string fit_message;
    int steps_needed=_fitter->fit(max_steps, fit_message);
    double chisqr=_fitter->get_chi_sqr();
    double dof=_fitter->get_dof()-_fitter->get_cut()-n_parameters_dof;
    if(dof<0.0)
    {
      dof=0.0;
    }
    if(steps_needed==max_steps+1)
    {
      cerr << "Multifit warning: Fit #" << n+1 << " did not converge after " << max_steps << " iterations. " << "chi^2/dof = " << chisqr/dof << endl;
    }
    else
    {
      cout << "Multifit fit #" << n+1 << " converged after " << steps_needed << " iterations. " << "chi^2/dof = " << chisqr/dof << endl;
    }
    for(int p=0; p<n_parameters; ++p)
    {
      *all_output[p] << _fitter->get_parameter(p) << endl;
    }
  }
  gsl_rng_free(rng);
  cout << "Multifit with " << fit_data.size() << "samples completed." << endl << endl;
  for(int p=0; p<n_parameters; ++p)
  {
    all_output[p] -> close();
    delete all_output[p];
  }

  return true;
}





bool load_data_file()
{
  if (data_file_type==0)
  {
    ifstream input(data_file.c_str());
    if(!input)
    {
      cerr << "Error: Could not open data file" << endl;
      return false;
    }
    cout << "Reading data file (text format)..." << endl;
    int n_functions_temp;
    int n_variables_temp;
    int n_data_points;
    int n_data_sets;
    input >> n_functions_temp;
    input >> n_variables_temp;
    input >> n_data_points;
    input >> n_data_sets;
    if(n_functions_temp!=n_functions)
    {
      cerr << "Error: data file has wrong number of functions" << endl;
      return false;
    }
    if(n_variables_temp!=n_variables)
    {
      cerr << "Error: data file has wrong number of variables" << endl;
      return false;
    }
    if( (n_data_sets<=0) || (n_data_points<=0) )
    {
      cerr << "Error: file contains no data" << endl;
      return false;
    }
    file_arguments.clear();
    file_data.clear();
    vector< double > file_arg_temp(n_variables, 0.0);
    file_arguments.resize(n_data_points, file_arg_temp);
    vector< double > file_data_temp_1(n_functions, 0.0);
    vector< vector< double > > file_data_temp_2(n_data_points, file_data_temp_1);
    file_data.resize(n_data_sets, file_data_temp_2);
    int m_buf;
    for(int m=0; m<n_data_points; ++m)
    {
      input >> m_buf;
      if(m_buf!=m+1)
      {
        cerr << "Error while reading data file" << endl;
        return false;
      }
      for(int v=0; v<n_variables; ++v)
      {
        input >> file_arguments[m][v];
      }
    }
    int n_buf;
    for(int n=0; n<n_data_sets; ++n)
    {
      for(int m=0; m<n_data_points; ++m)
      {
        input >> n_buf >> m_buf;
        if( (n_buf!=n+1) || (m_buf!=m+1) )
        {
          cerr << "Error while reading data file" << endl;
          return false;
        }
        for(int f=0; f<n_functions; ++f)
        {
          input >> file_data[n][m][f];
        }
      }
    }
    input.close();
  }
  else if (data_file_type==1)
  {
    ifstream input(data_file.c_str(), ios::in | ios::binary);
    if(!input)
    {
      cerr << "Error: Could not open data file" << endl;
      return false;
    }
    cout << "Reading data file (binary 32 bit little endian format)..." << endl;
    int n_functions_temp;
    int n_variables_temp;
    int n_data_points;
    int n_data_sets;
    float temp;
    temp=read_float(input);
    n_functions_temp=static_cast<int>(temp);
    temp=read_float(input);
    n_variables_temp=static_cast<int>(temp);
    temp=read_float(input);
    n_data_points=static_cast<int>(temp);
    temp=read_float(input);
    n_data_sets=static_cast<int>(temp);
    if(n_functions_temp!=n_functions)
    {
      cerr << "Error: data file has wrong number of functions" << endl;
      return false;
    }
    if(n_variables_temp!=n_variables)
    {
      cerr << "Error: data file has wrong number of variables" << endl;
      return false;
    }
    if( (n_data_sets<=0) || (n_data_points<=0) )
    {
      cerr << "Error: file contains no data" << endl;
      return false;
    }
    file_arguments.clear();
    file_data.clear();
    vector< double > file_arg_temp(n_variables, 0.0);
    file_arguments.resize(n_data_points, file_arg_temp);
    vector< double > file_data_temp_1(n_functions, 0.0);
    vector< vector< double > > file_data_temp_2(n_data_points, file_data_temp_1);
    file_data.resize(n_data_sets, file_data_temp_2);
    for(int m=0; m<n_data_points; ++m)
    {
      for(int v=0; v<n_variables; ++v)
      {
        file_arguments[m][v]=read_float(input);
      }
    }
    for(int n=0; n<n_data_sets; ++n)
    {
      for(int m=0; m<n_data_points; ++m)
      {
        for(int f=0; f<n_functions; ++f)
        {
          file_data[n][m][f]=read_float(input);
        }
      }
    }
    input.close();
  }
  return true;
}


bool set_fit_data()
{
  int n_data_points=file_arguments.size();
  int n_data_sets=file_data.size();

  int fit_n_data_sets;
  if(restrict_data)
  {
    fit_n_data_sets=stop_n_data-start_n_data+1;
    if (stop_n_data>n_data_sets)
    {
      cerr << "Error: range of data sets exceeds data file";
      return false;
    }
    if (fit_n_data_sets<5)
    {
      cerr << "Error: range of data sets too small";
      return false;
    }
  }
  else
  {
    start_n_data=1;
    stop_n_data=n_data_sets;
    fit_n_data_sets=n_data_sets;
  }


// select data points in fitting range
  vector< parser* > fit_min_parsers;
  vector< parser* > fit_max_parsers;
  for(int v=0; v<n_variables; ++v)
  {
    fit_min_parsers.push_back(new parser(fit_min[v]));
    fit_max_parsers.push_back(new parser(fit_max[v]));
  }
  map< string, double > table;
  for(int c=0; c<n_constants; ++c)
  {
    table[constants[c]]=constant_values[c];
  }
  fit_arguments.clear();
  vector< int > indices;
  for(int m=0; m<n_data_points; ++m)
  {
    for(int v=0; v<n_variables; ++v)
    {
      table[variables[v]]=file_arguments[m][v];
    }
    bool in_range=true;
    for(int v=0; v<n_variables; ++v)
    {
      double min=fit_min_parsers[v]->parse(table);
      double max=fit_max_parsers[v]->parse(table);
      if( (fit_min_parsers[v]->get_no_of_errors()>0) || (fit_max_parsers[v]->get_no_of_errors()>0))
      {
        for(vector< parser* >::iterator it=fit_min_parsers.begin(); it!=fit_min_parsers.end(); ++it)
        {
          delete *it;
        }
        for(vector< parser* >::iterator it=fit_max_parsers.begin(); it!=fit_max_parsers.end(); ++it)
        {
          delete *it;
        }
        cerr << "Error while parsing fitting range" << endl;
        return false;
      }
      double arg=file_arguments[m][v];
      if( !( (min<=arg) && (arg<=max) ) )
      {
        in_range=false;
      }
    }
    if(in_range)
    {
      fit_arguments.push_back(file_arguments[m]);
      indices.push_back(m);
    }
  }
  for(vector< parser* >::iterator it=fit_min_parsers.begin(); it!=fit_min_parsers.end(); ++it)
  {
    delete *it;
  }
  for(vector< parser* >::iterator it=fit_max_parsers.begin(); it!=fit_max_parsers.end(); ++it)
  {
    delete *it;
  }
// end select data points in fitting range
  int fit_n_data_points=fit_arguments.size();
  if(fit_n_data_points==0)
  {
    cerr << "Error: No data points in fitting range";
    return false;
  }
  if(fit_n_data_sets/bin_size<5)
  {
    cerr << "Error: Bin size too large";
    return false;
  }
  if(fit_n_data_points==0)
  {
    cerr << "Error: No data points in fitting range" << endl;
    return false;
  }

  vector< double > fit_data_temp_1(n_functions, 0.0);
  vector< vector< double > > fit_data_temp_2(fit_n_data_points, fit_data_temp_1);
  fit_data.clear();
  fit_data.resize(fit_n_data_sets/bin_size, fit_data_temp_2);

  for(int m=0; m<fit_n_data_points; ++m)
  {
    for(int f=0; f<n_functions; ++f)
    {
      for(int b=0; b<fit_n_data_sets/bin_size; ++b)
      {
        // average over entries e in bin b
        fit_data[b][m][f]=0.0;
        for(int e=0; e<bin_size; ++e)
        {
          fit_data[b][m][f]+=file_data[start_n_data-1+b*bin_size+e][indices[m]][f];
        }
        fit_data[b][m][f]/=bin_size;
      }
    }
  }
  _fitter->set_inversion_method(inv_method);
  _fitter->set_bootstrap_normalization(bootstrap_normalization);
  _fitter->set_svd_cut(svd_cut);
  _fitter->set_svd_cut_ratio(svd_ratio);
  _fitter->set_svd_cut_absolute(svd_value);
  string set_data_message;
  _fitter->set_data(fit_arguments, fit_data, set_data_message);
  if( (inv_method!=LU_inversion) || (_fitter->get_cut()!=0) )
  {
    cout << "Pseudoinverse of data correlation matrix computed using SVD. Kept "
            << _fitter->get_dof()-_fitter->get_cut() << " out of "
            << _fitter->get_dof() << " singular values." << endl;
  }
  if(set_data_message!="")
  {
    cout << endl << set_data_message << endl << endl;
  }
  return true;
}




bool reset_fitter()
{
  if (current_model==MULTIEXP)
  {
    delete _model;
    _model=new multi_exp_model(E_name, dE_name, A_name, B_name, t_name, multiexpdialog_nexp);
  }
  else if (current_model==MULTIEXPASQR)
  {
    delete _model;
    _model=new multi_exp_Asqr_model(E_name, dE_name, A_name, B_name, t_name, multiexpdialog_nexp);
  }
  else if (current_model==MULTIEXPASQREXPE)
  {
    delete _model;
    _model=new multi_exp_Asqr_expE_model(E_name, dE_name, A_name, B_name, t_name, multiexpdialog_nexp);
  }
//  else if (current_model==MULTIEXPEXPA)
//  {
//    delete _model;
//    _model=new multi_exp_expA_model(E_name, dE_name, A_name, B_name, t_name, multiexpdialog_nexp);
//  }
  else if (current_model==MULTIEXPEXPE)
  {
    delete _model;
    _model=new multi_exp_expE_model(E_name, dE_name, A_name, B_name, t_name, multiexpdialog_nexp);
  }
//  else if (current_model==MULTIEXPEXPAEXPE)
//  {
//    delete _model;
//    _model=new multi_exp_expA_expE_model(E_name, dE_name, A_name, B_name, t_name, multiexpdialog_nexp);
//  }
  else if (current_model==MULTIALTEXP)
  {
    delete _model;
    _model=new multi_alt_exp_model(E_name, dE_name, A_name, B_name, t_name, multialtexpdialog_nexp, multialtexpdialog_mexp);
  }
  else if (current_model==MULTIALTEXPASQR)
  {
    delete _model;
    _model=new multi_alt_exp_Asqr_model(E_name, dE_name, A_name, B_name, t_name, multialtexpdialog_nexp, multialtexpdialog_mexp);
  }
  else if (current_model==MULTIALTEXPASQREXPE)
  {
    delete _model;
    _model=new multi_alt_exp_Asqr_expE_model(E_name, dE_name, A_name, B_name, t_name, multialtexpdialog_nexp, multialtexpdialog_mexp);
  }
//  else if (current_model==MULTIALTEXPEXPA)
//  {
//    delete _model;
//    _model=new multi_alt_exp_expA_model(E_name, dE_name, A_name, B_name, t_name, multialtexpdialog_nexp, multialtexpdialog_mexp);
//  }
  else if (current_model==MULTIALTEXPEXPE)
  {
    delete _model;
    _model=new multi_alt_exp_expE_model(E_name, dE_name, A_name, B_name, t_name, multialtexpdialog_nexp, multialtexpdialog_mexp);
  }
//  else if (current_model==MULTIALTEXPEXPAEXPE)
//  {
//    delete _model;
//    _model=new multi_alt_exp_expA_expE_model(E_name, dE_name, A_name, B_name, t_name, multialtexpdialog_nexp, multialtexpdialog_mexp);
//  }

  else if (current_model==MULTIEXPVEC)
  {
    delete _model;
    _model=new multi_exp_vec_model(E_name, dE_name, A_name, B_name, t_name, multiexpvecdialog_nexp, multiexpvecdialog_nvec);
  }
  else if (current_model==MULTIEXPASQRVEC)
  {
    delete _model;
    _model=new multi_exp_Asqr_vec_model(E_name, dE_name, A_name, B_name, t_name, multiexpvecdialog_nexp, multiexpvecdialog_nvec);
  }
  else if (current_model==MULTIEXPASQREXPEVEC)
  {
    delete _model;
    _model=new multi_exp_Asqr_expE_vec_model(E_name, dE_name, A_name, B_name, t_name, multiexpvecdialog_nexp, multiexpvecdialog_nvec);
  }
//  else if (current_model==MULTIEXPEXPAVEC)
//  {
//    delete _model;
//    _model=new multi_exp_expA_vec_model(E_name, dE_name, A_name, B_name, t_name, multiexpvecdialog_nexp, multiexpvecdialog_nvec);
//  }
  else if (current_model==MULTIEXPEXPEVEC)
  {
    delete _model;
    _model=new multi_exp_expE_vec_model(E_name, dE_name, A_name, B_name, t_name, multiexpvecdialog_nexp, multiexpvecdialog_nvec);
  }
//  else if (current_model==MULTIEXPEXPAEXPEVEC)
//  {
//    delete _model;
//    _model=new multi_exp_expA_expE_vec_model(E_name, dE_name, A_name, B_name, t_name, multiexpvecdialog_nexp, multiexpvecdialog_nvec);
//  }
  else if (current_model==MULTIALTEXPVEC)
  {
    delete _model;
    _model=new multi_alt_exp_vec_model(E_name, dE_name, A_name, B_name, t_name, multialtexpvecdialog_nexp, multialtexpvecdialog_mexp, multialtexpvecdialog_nvec);
  }
  else if (current_model==MULTIALTEXPASQRVEC)
  {
    delete _model;
    _model=new multi_alt_exp_Asqr_vec_model(E_name, dE_name, A_name, B_name, t_name, multialtexpvecdialog_nexp, multialtexpvecdialog_mexp, multialtexpvecdialog_nvec);
  }
  else if (current_model==MULTIALTEXPASQREXPEVEC)
  {
    delete _model;
    _model=new multi_alt_exp_Asqr_expE_vec_model(E_name, dE_name, A_name, B_name, t_name, multialtexpvecdialog_nexp, multialtexpvecdialog_mexp, multialtexpvecdialog_nvec);
  }
//  else if (current_model==MULTIALTEXPEXPAVEC)
//  {
//    delete _model;
//    _model=new multi_alt_exp_expA_vec_model(E_name, dE_name, A_name, B_name, t_name, multialtexpvecdialog_nexp, multialtexpvecdialog_mexp, multialtexpvecdialog_nvec);
//  }
  else if (current_model==MULTIALTEXPEXPEVEC)
  {
    delete _model;
    _model=new multi_alt_exp_expE_vec_model(E_name, dE_name, A_name, B_name, t_name, multialtexpvecdialog_nexp, multialtexpvecdialog_mexp, multialtexpvecdialog_nvec);
  }
//  else if (current_model==MULTIALTEXPEXPAEXPEVEC)
//  {
//    delete _model;
//    _model=new multi_alt_exp_expA_expE_vec_model(E_name, dE_name, A_name, B_name, t_name, multialtexpvecdialog_nexp, multialtexpvecdialog_mexp, multialtexpvecdialog_nvec);
//  }

  else if (current_model==MULTIEXPMAT)
  {
    delete _model;
    _model=new multi_exp_mat_model(E_name, dE_name, A_name, B_name, t_name, multiexpmatdialog_nexp, multiexpmatdialog_dim_1, multiexpmatdialog_dim_2);
  }
  else if (current_model==MULTIEXPEXPEMAT)
  {
    delete _model;
    _model=new multi_exp_expE_mat_model(E_name, dE_name, A_name, B_name, t_name, multiexpmatdialog_nexp, multiexpmatdialog_dim_1, multiexpmatdialog_dim_2);
  }
  else if (current_model==MULTIEXPMATII)
  {
    delete _model;
    _model=new multi_exp_mat_II_model(E_name, dE_name, A_name, B_name, t_name, multiexpmatdialog_nexp, multiexpmatdialog_dim_1, multiexpmatdialog_dim_2);
  }
  else if (current_model==MULTIEXPEXPEMATII)
  {
    delete _model;
    _model=new multi_exp_expE_mat_II_model(E_name, dE_name, A_name, B_name, t_name, multiexpmatdialog_nexp, multiexpmatdialog_dim_1, multiexpmatdialog_dim_2);
  }

  else if (current_model==MULTIEXPMATUPPER)
  {
    delete _model;
    _model=new multi_exp_mat_upper_model(E_name, dE_name, A_name, B_name, t_name, multiexpmatupperdialog_nexp, multiexpmatupperdialog_dim);
  }
  else if (current_model==MULTIEXPEXPEMATUPPER)
  {
    delete _model;
    _model=new multi_exp_expE_mat_upper_model(E_name, dE_name, A_name, B_name, t_name, multiexpmatupperdialog_nexp, multiexpmatupperdialog_dim);
  }
  else if (current_model==MULTIEXPMATIIUPPER)
  {
    delete _model;
    _model=new multi_exp_mat_II_upper_model(E_name, dE_name, A_name, B_name, t_name, multiexpmatupperdialog_nexp, multiexpmatupperdialog_dim);
  }
  else if (current_model==MULTIEXPEXPEMATIIUPPER)
  {
    delete _model;
    _model=new multi_exp_expE_mat_II_upper_model(E_name, dE_name, A_name, B_name, t_name, multiexpmatupperdialog_nexp, multiexpmatupperdialog_dim);
  }


  else if (current_model==MULTIALTEXPMAT)
  {
    delete _model;
    _model=new multi_alt_exp_mat_model(E_name, dE_name, A_name, B_name, t_name, multialtexpmatdialog_nexp, multialtexpmatdialog_mexp, multialtexpmatdialog_dim_1, multialtexpmatdialog_dim_2);
  }
  else if (current_model==MULTIALTEXPEXPEMAT)
  {
    delete _model;
    _model=new multi_alt_exp_expE_mat_model(E_name, dE_name, A_name, B_name, t_name, multialtexpmatdialog_nexp, multialtexpmatdialog_mexp, multialtexpmatdialog_dim_1, multialtexpmatdialog_dim_2);
  }


  else if (current_model==MULTIEXPNONSYMMAT)
  {
    delete _model;
    _model=new multi_exp_nonsym_mat_model(E_name, dE_name, A_name, B_name, t_name, multiexpmatdialog_nexp, multiexpmatdialog_dim_1, multiexpmatdialog_dim_2);
  }
  else if (current_model==MULTIEXPEXPENONSYMMAT)
  {
    delete _model;
    _model=new multi_exp_expE_nonsym_mat_model(E_name, dE_name, A_name, B_name, t_name, multiexpmatdialog_nexp, multiexpmatdialog_dim_1, multiexpmatdialog_dim_2);
  }
  else if (current_model==MULTIALTEXPNONSYMMAT)
  {
    delete _model;
    _model=new multi_alt_exp_nonsym_mat_model(E_name, dE_name, A_name, B_name, t_name, multialtexpmatdialog_nexp, multialtexpmatdialog_mexp, multialtexpmatdialog_dim_1, multialtexpmatdialog_dim_2);
  }
  else if (current_model==MULTIALTEXPEXPENONSYMMAT)
  {
    delete _model;
    _model=new multi_alt_exp_expE_nonsym_mat_model(E_name, dE_name, A_name, B_name, t_name, multialtexpmatdialog_nexp, multialtexpmatdialog_mexp, multialtexpmatdialog_dim_1, multialtexpmatdialog_dim_2);
  }


  else if (current_model==MULTIEXPBC)
  {
    delete _model;
    _model=new multi_exp_BC_model(E_name, dE_name, A_name, B_name, t_name, T_name, multiexpdialog_nexp);
  }
  else if (current_model==MULTIEXPASQRBC)
  {
    delete _model;
    _model=new multi_exp_Asqr_BC_model(E_name, dE_name, A_name, B_name, t_name, T_name, multiexpdialog_nexp);
  }
  else if (current_model==MULTIEXPASQREXPEBC)
  {
    delete _model;
    _model=new multi_exp_Asqr_expE_BC_model(E_name, dE_name, A_name, B_name, t_name, T_name, multiexpdialog_nexp);
  }
//  else if (current_model==MULTIEXPEXPABC)
//  {
//    delete _model;
//    _model=new multi_exp_expA_BC_model(E_name, dE_name, A_name, B_name, t_name, T_name, multiexpdialog_nexp);
//  }
  else if (current_model==MULTIEXPEXPEBC)
  {
    delete _model;
    _model=new multi_exp_expE_BC_model(E_name, dE_name, A_name, B_name, t_name, T_name, multiexpdialog_nexp);
  }
//  else if (current_model==MULTIEXPEXPAEXPEBC)
//  {
//    delete _model;
//    _model=new multi_exp_expA_expE_BC_model(E_name, dE_name, A_name, B_name, t_name, T_name, multiexpdialog_nexp);
//  }
  else if (current_model==MULTIALTEXPBC)
  {
    delete _model;
    _model=new multi_alt_exp_BC_model(E_name, dE_name, A_name, B_name, t_name, T_name, multialtexpdialog_nexp, multialtexpdialog_mexp);
  }
  else if (current_model==MULTIALTEXPASQRBC)
  {
    delete _model;
    _model=new multi_alt_exp_Asqr_BC_model(E_name, dE_name, A_name, B_name, t_name, T_name, multialtexpdialog_nexp, multialtexpdialog_mexp);
  }
  else if (current_model==MULTIALTEXPASQREXPEBC)
  {
    delete _model;
    _model=new multi_alt_exp_Asqr_expE_BC_model(E_name, dE_name, A_name, B_name, t_name, T_name, multialtexpdialog_nexp, multialtexpdialog_mexp);
  }
//  else if (current_model==MULTIALTEXPEXPABC)
//  {
//    delete _model;
//    _model=new multi_alt_exp_expA_BC_model(E_name, dE_name, A_name, B_name, t_name, T_name, multialtexpdialog_nexp, multialtexpdialog_mexp);
//  }
  else if (current_model==MULTIALTEXPEXPEBC)
  {
    delete _model;
    _model=new multi_alt_exp_expE_BC_model(E_name, dE_name, A_name, B_name, t_name, T_name, multialtexpdialog_nexp, multialtexpdialog_mexp);
  }
//  else if (current_model==MULTIALTEXPEXPAEXPEBC)
//  {
//    delete _model;
//    _model=new multi_alt_exp_expA_expE_BC_model(E_name, dE_name, A_name, B_name, t_name, T_name, multialtexpdialog_nexp, multialtexpdialog_mexp);
//  }

  else if (current_model==MULTIEXPVECBC)
  {
    delete _model;
    _model=new multi_exp_vec_BC_model(E_name, dE_name, A_name, B_name, t_name, T_name, multiexpvecdialog_nexp, multiexpvecdialog_nvec);
  }
  else if (current_model==MULTIEXPASQRVECBC)
  {
    delete _model;
    _model=new multi_exp_Asqr_vec_BC_model(E_name, dE_name, A_name, B_name, t_name, T_name, multiexpvecdialog_nexp, multiexpvecdialog_nvec);
  }
  else if (current_model==MULTIEXPASQREXPEVECBC)
  {
    delete _model;
    _model=new multi_exp_Asqr_expE_vec_BC_model(E_name, dE_name, A_name, B_name, t_name, T_name, multiexpvecdialog_nexp, multiexpvecdialog_nvec);
  }
//  else if (current_model==MULTIEXPEXPAVECBC)
//  {
//    delete _model;
//    _model=new multi_exp_expA_vec_BC_model(E_name, dE_name, A_name, B_name, t_name, T_name, multiexpvecdialog_nexp, multiexpvecdialog_nvec);
//  }
  else if (current_model==MULTIEXPEXPEVECBC)
  {
    delete _model;
    _model=new multi_exp_expE_vec_BC_model(E_name, dE_name, A_name, B_name, t_name, T_name, multiexpvecdialog_nexp, multiexpvecdialog_nvec);
  }
//  else if (current_model==MULTIEXPEXPAEXPEVECBC)
//  {
//    delete _model;
//    _model=new multi_exp_expA_expE_vec_BC_model(E_name, dE_name, A_name, B_name, t_name, T_name, multiexpvecdialog_nexp, multiexpvecdialog_nvec);
//  }
  else if (current_model==MULTIALTEXPVECBC)
  {
    delete _model;
    _model=new multi_alt_exp_vec_BC_model(E_name, dE_name, A_name, B_name, t_name, T_name, multialtexpvecdialog_nexp, multialtexpvecdialog_mexp, multialtexpvecdialog_nvec);
  }
  else if (current_model==MULTIALTEXPASQRVECBC)
  {
    delete _model;
    _model=new multi_alt_exp_Asqr_vec_BC_model(E_name, dE_name, A_name, B_name, t_name, T_name, multialtexpvecdialog_nexp, multialtexpvecdialog_mexp, multialtexpvecdialog_nvec);
  }
  else if (current_model==MULTIALTEXPASQREXPEVECBC)
  {
    delete _model;
    _model=new multi_alt_exp_Asqr_expE_vec_BC_model(E_name, dE_name, A_name, B_name, t_name, T_name, multialtexpvecdialog_nexp, multialtexpvecdialog_mexp, multialtexpvecdialog_nvec);
  }
//  else if (current_model==MULTIALTEXPEXPAVECBC)
//  {
//    delete _model;
//    _model=new multi_alt_exp_expA_vec_BC_model(E_name, dE_name, A_name, B_name, t_name, T_name, multialtexpvecdialog_nexp, multialtexpvecdialog_mexp, multialtexpvecdialog_nvec);
//  }
  else if (current_model==MULTIALTEXPEXPEVECBC)
  {
    delete _model;
    _model=new multi_alt_exp_expE_vec_BC_model(E_name, dE_name, A_name, B_name, t_name, T_name, multialtexpvecdialog_nexp, multialtexpvecdialog_mexp, multialtexpvecdialog_nvec);
  }
//  else if (current_model==MULTIALTEXPEXPAEXPEVECBC)
//  {
//    delete _model;
//    _model=new multi_alt_exp_expA_expE_vec_BC_model(E_name, dE_name, A_name, B_name, t_name, T_name, multialtexpvecdialog_nexp, multialtexpvecdialog_mexp, multialtexpvecdialog_nvec);
//  }
  else if (current_model==THREEPT2VARMULTIEXP)
  {
    delete _model;
    _model=new threept_multi_exp_model(E_initial_name, dE_initial_name, E_final_name, dE_final_name, A_name, B_name, t_name, T_name, threeptmultiexpdialog_nexpinitial, threeptmultiexpdialog_nexpfinal);
  }
  else if (current_model==THREEPT2VARMULTIEXPEXPE)
  {
    delete _model;
    _model=new threept_multi_exp_expE_model(E_initial_name, dE_initial_name, E_final_name, dE_final_name, A_name, B_name, t_name, T_name, threeptmultiexpdialog_nexpinitial, threeptmultiexpdialog_nexpfinal);
  }
  else if (current_model==THREEPT2VARMULTIEXPVEC)
  {
    delete _model;
    _model=new threept_multi_exp_vec_model(E_initial_name, dE_initial_name, E_final_name, dE_final_name, A_name, B_name, t_name, T_name, threeptmultiexpvecdialog_nexpinitial, threeptmultiexpvecdialog_nexpfinal, threeptmultiexpvecdialog_nvec);
  }
  else if (current_model==THREEPT2VARMULTIEXPEXPEVEC)
  {
    delete _model;
    _model=new threept_multi_exp_expE_vec_model(E_initial_name, dE_initial_name, E_final_name, dE_final_name, A_name, B_name, t_name, T_name, threeptmultiexpvecdialog_nexpinitial, threeptmultiexpvecdialog_nexpfinal, threeptmultiexpvecdialog_nvec);
  }
  else if (current_model==THREEPT2VARMULTIALTEXP)
  {
    delete _model;
    _model=new threept_multi_alt_exp_model(E_initial_name, dE_initial_name, E_final_name, dE_final_name, A_name, B_name, t_name, T_name, threeptmultialtexpdialog_nexpinitial, threeptmultialtexpdialog_mexpinitial, threeptmultialtexpdialog_nexpfinal, threeptmultialtexpdialog_mexpfinal);
  }
  else if (current_model==THREEPT2VARMULTIALTEXPEXPE)
  {
    delete _model;
    _model=new threept_multi_alt_exp_expE_model(E_initial_name, dE_initial_name, E_final_name, dE_final_name, A_name, B_name, t_name, T_name, threeptmultialtexpdialog_nexpinitial, threeptmultialtexpdialog_mexpinitial, threeptmultialtexpdialog_nexpfinal, threeptmultialtexpdialog_mexpfinal);
  }
  else if (current_model==THREEPT2VARMULTIALTEXPVEC)
  {
    delete _model;
    _model=new threept_multi_alt_exp_vec_model(E_initial_name, dE_initial_name, E_final_name, dE_final_name, A_name, B_name, t_name, T_name, threeptmultialtexpvecdialog_nexpinitial, threeptmultialtexpvecdialog_mexpinitial, threeptmultialtexpvecdialog_nexpfinal, threeptmultialtexpvecdialog_mexpfinal, threeptmultialtexpvecdialog_nvec);
  }
  else if (current_model==THREEPT2VARMULTIALTEXPEXPEVEC)
  {
    delete _model;
    _model=new threept_multi_alt_exp_expE_vec_model(E_initial_name, dE_initial_name, E_final_name, dE_final_name, A_name, B_name, t_name, T_name, threeptmultialtexpvecdialog_nexpinitial, threeptmultialtexpvecdialog_mexpinitial, threeptmultialtexpvecdialog_nexpfinal, threeptmultialtexpvecdialog_mexpfinal, threeptmultialtexpvecdialog_nvec);
  }
  else if (current_model==MULTIEXPCONST)
  {
    delete _model;
    _model=new multi_exp_const_model(E_name, dE_name, A_name, B_name, t_name, C_name, multiexpdialog_nexp);
  }
  else if (current_model==MULTIEXPASQRCONST)
  {
    delete _model;
    _model=new multi_exp_Asqr_const_model(E_name, dE_name, A_name, B_name, t_name, C_name, multiexpdialog_nexp);
  }
  else if (current_model==MULTIEXPASQREXPECONST)
  {
    delete _model;
    _model=new multi_exp_Asqr_expE_const_model(E_name, dE_name, A_name, B_name, t_name, C_name, multiexpdialog_nexp);
  }
  else if (current_model==MULTIEXPEXPECONST)
  {
    delete _model;
    _model=new multi_exp_expE_const_model(E_name, dE_name, A_name, B_name, t_name, C_name, multiexpdialog_nexp);
  }
  else if (current_model==MULTIALTEXPCONST)
  {
    delete _model;
    _model=new multi_alt_exp_const_model(E_name, dE_name, A_name, B_name, t_name, C_name, multialtexpdialog_nexp, multialtexpdialog_mexp);
  }
  else if (current_model==MULTIALTEXPASQRCONST)
  {
    delete _model;
    _model=new multi_alt_exp_Asqr_const_model(E_name, dE_name, A_name, B_name, t_name, C_name, multialtexpdialog_nexp, multialtexpdialog_mexp);
  }
  else if (current_model==MULTIALTEXPASQREXPECONST)
  {
    delete _model;
    _model=new multi_alt_exp_Asqr_expE_const_model(E_name, dE_name, A_name, B_name, t_name, C_name, multialtexpdialog_nexp, multialtexpdialog_mexp);
  }
  else if (current_model==MULTIALTEXPEXPECONST)
  {
    delete _model;
    _model=new multi_alt_exp_expE_const_model(E_name, dE_name, A_name, B_name, t_name, C_name, multialtexpdialog_nexp, multialtexpdialog_mexp);
  }
  else if (current_model==MULTIEXPVECCONST)
  {
    delete _model;
    _model=new multi_exp_vec_const_model(E_name, dE_name, A_name, B_name, t_name, C_name, multiexpvecdialog_nexp, multiexpvecdialog_nvec);
  }
  else if (current_model==MULTIEXPASQRVECCONST)
  {
    delete _model;
    _model=new multi_exp_Asqr_vec_const_model(E_name, dE_name, A_name, B_name, t_name, C_name, multiexpvecdialog_nexp, multiexpvecdialog_nvec);
  }
  else if (current_model==MULTIEXPASQREXPEVECCONST)
  {
    delete _model;
    _model=new multi_exp_Asqr_expE_vec_const_model(E_name, dE_name, A_name, B_name, t_name, C_name, multiexpvecdialog_nexp, multiexpvecdialog_nvec);
  }
  else if (current_model==MULTIEXPEXPEVECCONST)
  {
    delete _model;
    _model=new multi_exp_expE_vec_const_model(E_name, dE_name, A_name, B_name, t_name, C_name, multiexpvecdialog_nexp, multiexpvecdialog_nvec);
  }
  else if (current_model==MULTIALTEXPVECCONST)
  {
    delete _model;
    _model=new multi_alt_exp_vec_const_model(E_name, dE_name, A_name, B_name, t_name, C_name, multialtexpvecdialog_nexp, multialtexpvecdialog_mexp, multialtexpvecdialog_nvec);
  }
  else if (current_model==MULTIALTEXPASQRVECCONST)
  {
    delete _model;
    _model=new multi_alt_exp_Asqr_vec_const_model(E_name, dE_name, A_name, B_name, t_name, C_name, multialtexpvecdialog_nexp, multialtexpvecdialog_mexp, multialtexpvecdialog_nvec);
  }
  else if (current_model==MULTIALTEXPASQREXPEVECCONST)
  {
    delete _model;
    _model=new multi_alt_exp_Asqr_expE_vec_const_model(E_name, dE_name, A_name, B_name, t_name, C_name, multialtexpvecdialog_nexp, multialtexpvecdialog_mexp, multialtexpvecdialog_nvec);
  }
  else if (current_model==MULTIALTEXPEXPEVECCONST)
  {
    delete _model;
    _model=new multi_alt_exp_expE_vec_const_model(E_name, dE_name, A_name, B_name, t_name, C_name, multialtexpvecdialog_nexp, multialtexpvecdialog_mexp, multialtexpvecdialog_nvec);
  }
  else if (current_model==MULTIEXPBCCONST)
  {
    delete _model;
    _model=new multi_exp_BC_const_model(E_name, dE_name, A_name, B_name, t_name, T_name, C_name, multiexpdialog_nexp);
  }
  else if (current_model==MULTIEXPASQRBCCONST)
  {
    delete _model;
    _model=new multi_exp_Asqr_BC_const_model(E_name, dE_name, A_name, B_name, t_name, T_name, C_name, multiexpdialog_nexp);
  }
  else if (current_model==MULTIEXPASQREXPEBCCONST)
  {
    delete _model;
    _model=new multi_exp_Asqr_expE_BC_const_model(E_name, dE_name, A_name, B_name, t_name, T_name, C_name, multiexpdialog_nexp);
  }
  else if (current_model==MULTIEXPEXPEBCCONST)
  {
    delete _model;
    _model=new multi_exp_expE_BC_const_model(E_name, dE_name, A_name, B_name, t_name, T_name, C_name, multiexpdialog_nexp);
  }
  else if (current_model==MULTIALTEXPBCCONST)
  {
    delete _model;
    _model=new multi_alt_exp_BC_const_model(E_name, dE_name, A_name, B_name, t_name, T_name, C_name, multialtexpdialog_nexp, multialtexpdialog_mexp);
  }
  else if (current_model==MULTIALTEXPASQRBCCONST)
  {
    delete _model;
    _model=new multi_alt_exp_Asqr_BC_const_model(E_name, dE_name, A_name, B_name, t_name, T_name, C_name, multialtexpdialog_nexp, multialtexpdialog_mexp);
  }
  else if (current_model==MULTIALTEXPASQREXPEBCCONST)
  {
    delete _model;
    _model=new multi_alt_exp_Asqr_expE_BC_const_model(E_name, dE_name, A_name, B_name, t_name, T_name, C_name, multialtexpdialog_nexp, multialtexpdialog_mexp);
  }
  else if (current_model==MULTIALTEXPEXPEBCCONST)
  {
    delete _model;
    _model=new multi_alt_exp_expE_BC_const_model(E_name, dE_name, A_name, B_name, t_name, T_name, C_name, multialtexpdialog_nexp, multialtexpdialog_mexp);
  }
  else if (current_model==MULTIEXPVECBCCONST)
  {
    delete _model;
    _model=new multi_exp_vec_BC_const_model(E_name, dE_name, A_name, B_name, t_name, T_name, C_name, multiexpvecdialog_nexp, multiexpvecdialog_nvec);
  }
  else if (current_model==MULTIEXPASQRVECBCCONST)
  {
    delete _model;
    _model=new multi_exp_Asqr_vec_BC_const_model(E_name, dE_name, A_name, B_name, t_name, T_name, C_name, multiexpvecdialog_nexp, multiexpvecdialog_nvec);
  }
  else if (current_model==MULTIEXPASQREXPEVECBCCONST)
  {
    delete _model;
    _model=new multi_exp_Asqr_expE_vec_BC_const_model(E_name, dE_name, A_name, B_name, t_name, T_name, C_name, multiexpvecdialog_nexp, multiexpvecdialog_nvec);
  }
  else if (current_model==MULTIEXPEXPEVECBCCONST)
  {
    delete _model;
    _model=new multi_exp_expE_vec_BC_const_model(E_name, dE_name, A_name, B_name, t_name, T_name, C_name, multiexpvecdialog_nexp, multiexpvecdialog_nvec);
  }
  else if (current_model==MULTIALTEXPVECBCCONST)
  {
    delete _model;
    _model=new multi_alt_exp_vec_BC_const_model(E_name, dE_name, A_name, B_name, t_name, T_name, C_name, multialtexpvecdialog_nexp, multialtexpvecdialog_mexp, multialtexpvecdialog_nvec);
  }
  else if (current_model==MULTIALTEXPASQRVECBCCONST)
  {
    delete _model;
    _model=new multi_alt_exp_Asqr_vec_BC_const_model(E_name, dE_name, A_name, B_name, t_name, T_name, C_name, multialtexpvecdialog_nexp, multialtexpvecdialog_mexp, multialtexpvecdialog_nvec);
  }
  else if (current_model==MULTIALTEXPASQREXPEVECBCCONST)
  {
    delete _model;
    _model=new multi_alt_exp_Asqr_expE_vec_BC_const_model(E_name, dE_name, A_name, B_name, t_name, T_name, C_name, multialtexpvecdialog_nexp, multialtexpvecdialog_mexp, multialtexpvecdialog_nvec);
  }
  else if (current_model==MULTIALTEXPEXPEVECBCCONST)
  {
    delete _model;
    _model=new multi_alt_exp_expE_vec_BC_const_model(E_name, dE_name, A_name, B_name, t_name, T_name, C_name, multialtexpvecdialog_nexp, multialtexpvecdialog_mexp, multialtexpvecdialog_nvec);
  }
  else // use PARSE
  {
    for(int f=0; f<n_functions; ++f)
    {
      if(functions[f].empty())
      {
        cerr << "Error: function missing" << endl;
        return false;
      }
      if( _in("{", functions[f]) || _in("&", functions[f]) || _in("&", functions[f]) || _in("!", functions[f]) || _in("}", functions[f])
       || _in("%", functions[f]) || _in("?", functions[f]) || _in("[", functions[f]) || _in(":", functions[f])
       || _in(">", functions[f]) || _in("#", functions[f]) || _in("~", functions[f]) || _in("$", functions[f])
       || _in("]", functions[f]) || _in("<", functions[f]) )
      {
        cerr << "Error: forbidden character in function " << f+1 << endl;
        return false;
      }
    }
    for(int v=0; v<n_variables; ++v)
    {
      if(variables[v].empty())
      {
        cerr << "Error: variable name missing" << endl;
        return false;
      }
      if(count(variables.begin(), variables.end(), variables[v])>1)
      {
        cerr << "Error: variable name repeated" << endl;
        return false;
      }
      if( _in("{", variables[v]) || _in("&", variables[v]) || _in("&", variables[v]) || _in("!", variables[v]) || _in("}", variables[v])
       || _in("%", variables[v]) || _in("?", variables[v]) || _in("[", variables[v]) || _in(":", variables[v])
       || _in(">", variables[v]) || _in("#", variables[v]) || _in("~", variables[v]) || _in("$", variables[v])
       || _in("]", variables[v]) || _in("<", variables[v]) ||
                     _in("arcsin", variables[v]) ||  _in("sinh", variables[v]) ||
                     _in("sin", variables[v]) ||  _in("arccos", variables[v]) ||
                     _in("cosh", variables[v]) ||  _in("cos", variables[v]) ||
                     _in("arctan", variables[v]) ||  _in("tanh", variables[v]) ||
                     _in("tan", variables[v]) ||  _in("exp", variables[v]) ||
                     _in("log", variables[v]) ||  _in("sqrt", variables[v]) ||
                     _in("sqr", variables[v]) ||  _in("alt", variables[v]) )
      {
        cerr << "Error: forbidden character or substring in variable " << v+1 << endl;
        return false;
      }
    }
    for(int p=0; p<n_parameters; ++p)
    {
      if(parameters[p].empty())
      {
        cerr << "Error: parameter name missing" << endl;
        return false;
      }
      if(count(parameters.begin(), parameters.end(), parameters[p])>1)
      {
        cerr << "Error: parameter name repeated" << endl;
        return false;
      }
      if( _in("{", parameters[p]) || _in("&", parameters[p]) || _in("&", parameters[p]) || _in("!", parameters[p]) || _in("}", parameters[p])
       || _in("%", parameters[p]) || _in("?", parameters[p]) || _in("[", parameters[p]) || _in(":", parameters[p])
       || _in(">", parameters[p]) || _in("#", parameters[p]) || _in("~", parameters[p]) || _in("$", parameters[p])
       || _in("]", parameters[p]) || _in("<", parameters[p]) ||
                     _in("arcsin", parameters[p]) ||  _in("sinh", parameters[p]) ||
                     _in("sin", parameters[p]) ||  _in("arccos", parameters[p]) ||
                     _in("cosh", parameters[p]) ||  _in("cos", parameters[p]) ||
                     _in("arctan", parameters[p]) ||  _in("tanh", parameters[p]) ||
                     _in("tan", parameters[p]) ||  _in("exp", parameters[p]) ||
                     _in("log", parameters[p]) ||  _in("sqrt", parameters[p]) ||
                     _in("sqr", parameters[p]) ||  _in("alt", parameters[p]) )
      {
        cerr << "Error: forbidden character or substring in parameter " << p+1 << endl;
        return false;
      }
    }
    for(int c=0; c<n_constants; ++c)
    {
      if(constants[c].empty())
      {
        cerr << "Error: constant name missing" << endl;
        return false;
      }
      if(count(constants.begin(), constants.end(), constants[c])>1)
      {
        cerr << "Error: constant name repeated" << endl;
        return false;
      }
      if( _in("{", constants[c]) || _in("&", constants[c]) || _in("&", constants[c]) || _in("!", constants[c]) || _in("}", constants[c])
       || _in("%", constants[c]) || _in("?", constants[c]) || _in("[", constants[c]) || _in(":", constants[c])
       || _in(">", constants[c]) || _in("#", constants[c]) || _in("~", constants[c]) || _in("$", constants[c])
       || _in("]", constants[c]) || _in("<", constants[c]) ||
                     _in("arcsin", constants[c]) ||  _in("sinh", constants[c]) ||
                     _in("sin", constants[c]) ||  _in("arccos", constants[c]) ||
                     _in("cosh", constants[c]) ||  _in("cos", constants[c]) ||
                     _in("arctan", constants[c]) ||  _in("tanh", constants[c]) ||
                     _in("tan", constants[c]) ||  _in("exp", constants[c]) ||
                     _in("log", constants[c]) ||  _in("sqrt", constants[c]) ||
                     _in("sqr", constants[c]) ||  _in("alt", constants[c]) )
      {
        cerr << "Error: forbidden character or substring in constant " << c+1 << endl;
        return false;
      }
    }
    if(!num_diff)
    {
      for(int f=0; f<n_functions; ++f)
      {
        for(int p=0; p<n_parameters; ++p)
        {
          if(derivatives[f][p].empty())
          {
            cerr << "Error: derivative" << endl;
            cerr << "of function " << f+1 << endl;
            cerr << "by parameter " << p+1 << endl;
            cerr << "missing" << endl;
            return false;
          }
          if( _in("{", derivatives[f][p]) || _in("&", derivatives[f][p]) || _in("&", derivatives[f][p]) || _in("!", derivatives[f][p]) || _in("}", derivatives[f][p])
           || _in("%", derivatives[f][p]) || _in("?", derivatives[f][p]) || _in("[", derivatives[f][p]) || _in(":", derivatives[f][p])
           || _in(">", derivatives[f][p]) || _in("#", derivatives[f][p]) || _in("~", derivatives[f][p]) || _in("$", derivatives[f][p])
           || _in("]", derivatives[f][p]) || _in("<", derivatives[f][p]) )
          {
            cerr << "Error: forbidden character in derivative" << endl;
            cerr << "of function " << f+1 << endl;
            cerr << "by parameter " << p+1 << endl;
            return false;
          }
        }
      }
    }
    delete _model;
    _model=new parse_model(functions, variables, parameters, constants, derivatives);
  }
  delete _gaussian_prior;
  _gaussian_prior=new gaussian_prior(n_parameters);


  if(chisqr_extra_term_enabled)
  {
    if(chisqr_extra_term_function.empty())
    {
      cerr << "Error: chi^2 extra term function missing" << endl;
      return false;
    }
    if( _in("{", chisqr_extra_term_function) || _in("&", chisqr_extra_term_function) || _in("&", chisqr_extra_term_function) || _in("!", chisqr_extra_term_function) || _in("}", chisqr_extra_term_function)
      || _in("%", chisqr_extra_term_function) || _in("?", chisqr_extra_term_function) || _in("[", chisqr_extra_term_function) || _in(":", chisqr_extra_term_function)
      || _in(">", chisqr_extra_term_function) || _in("#", chisqr_extra_term_function) || _in("~", chisqr_extra_term_function) || _in("$", chisqr_extra_term_function)
      || _in("]", chisqr_extra_term_function) || _in("<", chisqr_extra_term_function) )
    {
      cerr << "Error: forbidden character in chi^2 extra term function" << endl;
      return false;
    }
    for(int c=0; c<chisqr_extra_term_n_constants; ++c)
    {
      if(chisqr_extra_term_constant_names[c].empty())
      {
        cerr << "Error: constant name missing" << endl;
        return false;
      }
      if(count(chisqr_extra_term_constant_names.begin(), chisqr_extra_term_constant_names.end(), chisqr_extra_term_constant_names[c])>1)
      {
        cerr << "Error: constant name repeated" << endl;
        return false;
      }
      if( _in("{", chisqr_extra_term_constant_names[c]) || _in("&", chisqr_extra_term_constant_names[c]) || _in("&", chisqr_extra_term_constant_names[c]) || _in("!", chisqr_extra_term_constant_names[c]) || _in("}", chisqr_extra_term_constant_names[c])
       || _in("%", chisqr_extra_term_constant_names[c]) || _in("?", chisqr_extra_term_constant_names[c]) || _in("[", chisqr_extra_term_constant_names[c]) || _in(":", chisqr_extra_term_constant_names[c])
       || _in(">", chisqr_extra_term_constant_names[c]) || _in("#", chisqr_extra_term_constant_names[c]) || _in("~", chisqr_extra_term_constant_names[c]) || _in("$", chisqr_extra_term_constant_names[c])
       || _in("]", chisqr_extra_term_constant_names[c]) || _in("<", chisqr_extra_term_constant_names[c]) ||
                     _in("arcsin", chisqr_extra_term_constant_names[c]) ||  _in("sinh", chisqr_extra_term_constant_names[c]) ||
                     _in("sin", chisqr_extra_term_constant_names[c]) ||  _in("arccos", chisqr_extra_term_constant_names[c]) ||
                     _in("cosh", chisqr_extra_term_constant_names[c]) ||  _in("cos", chisqr_extra_term_constant_names[c]) ||
                     _in("arctan", chisqr_extra_term_constant_names[c]) ||  _in("tanh", chisqr_extra_term_constant_names[c]) ||
                     _in("tan", chisqr_extra_term_constant_names[c]) ||  _in("exp", chisqr_extra_term_constant_names[c]) ||
                     _in("log", chisqr_extra_term_constant_names[c]) ||  _in("sqrt", chisqr_extra_term_constant_names[c]) ||
                     _in("sqr", chisqr_extra_term_constant_names[c]) ||  _in("alt", chisqr_extra_term_constant_names[c]) )
      {
        cerr << "Error: forbidden character or substring in chi^2 extra term constant " << c+1 << endl;
        return false;
      }
    }
    for(int p=0; p<_model->get_n_parameters(); ++p)
    {
      if(std::find(chisqr_extra_term_constant_names.begin(), chisqr_extra_term_constant_names.end(), _model->get_parameter_name(p))!=chisqr_extra_term_constant_names.end())
      {
        cerr << "Error: constant name \"" << _model->get_parameter_name(p) << "\" in chi^2 extra term function is equal to a fit parameter name " << std::endl << std::endl;
        return false;
      }
    }
    for(int c=0; c<_model->get_n_constants(); ++c)
    {
      if(std::find(chisqr_extra_term_constant_names.begin(), chisqr_extra_term_constant_names.end(), _model->get_constant_name(c))!=chisqr_extra_term_constant_names.end())
      {
        cerr << "Error: constant name \"" << _model->get_constant_name(c) << "\" in chi^2 extra term function is equal to a fit constant name " << std::endl << std::endl;
        return false;
      }
    }
  }
  delete _chisqr_extra_term;
  _chisqr_extra_term = new chisqr_extra_term(_model->get_n_parameters(), chisqr_extra_term_constant_names.size(), chisqr_extra_term_function);
  _chisqr_extra_term->set_enabled(chisqr_extra_term_enabled);
  for(int p=0; p<_model->get_n_parameters(); ++p)
  {
    _chisqr_extra_term->set_parameter_name(p, _model->get_parameter_name(p));
  }
  _chisqr_extra_term->init_parameter_list();
  for(int c=0; c<chisqr_extra_term_constant_names.size(); ++c)
  {
    _chisqr_extra_term->set_constant_name(c, chisqr_extra_term_constant_names[c]);
  }

  delete _fitter;
  _fitter=new fitter(_model, _gaussian_prior, _chisqr_extra_term);
  return true;
}

string stripfilename(string fullname)
{
  int pos=fullname.rfind("/");
  fullname.erase(0, pos+1);
  return fullname;
}

void print_usage()
{
  cerr << endl << "usage: MBF [options] file" << endl << endl;
  cerr << "Options:" << endl;
  cerr << "  -b          perform bootstrap" << endl;
  cerr << "  -m          perform multifit" << endl << endl;
  cerr << "MBF always performs a fit and writes fit results to output directory" << endl << endl;
}

int main(int argc, char *argv[])
{
  gsl_set_error_handler_off();
  stringstream version_st;
  version_st.precision(2);
  version_st.setf(ios::fixed);
  version_st << version;
  cout << endl << "This is MBF version " << version_st.str() << endl << endl;
  if(argc==2)
  {
    init();
    if(!loadFile(argv[1]))
    {
      return 1;
    }
    if(fit(stripfilename(argv[1])))
    {
      return 0;
    }
    else
    {
      return 1;
    }
  }
  else if(argc==3)
  {
    if(string(argv[1])==string("-b"))
    {
      init();
      if(!loadFile(argv[2]))
      {
        return 1;
      }
      if(bootstrap(stripfilename(argv[2])))
      {
        return 0;
      }
      else
      {
        return 1;
      }
    }
    else if(string(argv[1])==string("-m"))
    {
      init();
      if(!loadFile(argv[2]))
      {
        return 1;
      }
      if(multifit(stripfilename(argv[2])))
      {
        return 0;
      }
      else
      {
        return 1;
      }
    }
    else
    {
      print_usage();
      return 1;
    }
  }
  else
  {
    print_usage();
    return 1;
  }
}
