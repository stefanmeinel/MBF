#include "multi_exp_expE_mat_model.h"

multi_exp_expE_mat_model::multi_exp_expE_mat_model(string _E_str, string _dE_str, string _A_str, string _B_str, string _t_str, int n_exp, int dim_1, int dim_2)
{
  E_str=_E_str;
  dE_str=_dE_str;
  A_str=_A_str;
  B_str=_B_str;
  t_str=_t_str;
  nexp=n_exp;
  dim1=dim_1;
  dim2=dim_2;
  if(dim1>=dim2)
  {
    dim=dim1;
  }
  else
  {
    dim=dim2;
  }
  params.resize(nexp+dim*nexp);
}


int multi_exp_expE_mat_model::get_n_functions()
{
  return dim1*dim2;
}


int multi_exp_expE_mat_model::get_n_variables()
{
  return 1;
}


int multi_exp_expE_mat_model::get_n_parameters()
{
  return nexp+dim*nexp;
}

int multi_exp_expE_mat_model::get_n_constants()
{
  return 0;
}


void multi_exp_expE_mat_model::set_parameters(const vector< double >& parameters)
{
  params=parameters;
}


void multi_exp_expE_mat_model::set_constants(const vector< double >& constants)
{
}


double multi_exp_expE_mat_model::eval_function(int function, const vector< double >& arguments)
{
  int v1=function/dim2;
  int v2=function%dim2;
  double t=arguments[0];
  double sum=0.0;
  double E=exp(params[0]);
  sum+=exp(-E*t);
  for(int n=1; n<nexp; ++n)
  {
    E+=exp(params[n]);
    sum+=params[nexp+v1*nexp+n]*params[nexp+v2*nexp+n]*exp(-E*t);
  }
  sum*=params[nexp+v1*nexp]*params[nexp+v2*nexp];
  return sum;
}


double multi_exp_expE_mat_model::eval_derivative(int function, int parameter, const vector< double >& arguments)
{
  int v1=function/dim2;
  int v2=function%dim2;
  int pv=(parameter-nexp)/nexp;
  int pe=(parameter-nexp)%nexp;
  double t=arguments[0];
  if(parameter==0)
  {
    return (-exp(params[0])*t)*eval_function(function, arguments);
  }
  else if( (parameter>0) && (parameter<nexp) )
  {
    double sum=0.0;
    double E=exp(params[0]);
    for(int n=1; n<parameter; ++n)
    {
      E+=exp(params[n]);
    }
    for(int n=parameter; n<nexp; ++n)
    {
      E+=exp(params[n]);
      sum+=params[nexp+v1*nexp+n]*params[nexp+v2*nexp+n]*exp(-E*t);
    }
    sum*=(-exp(params[parameter])*t)*params[nexp+v1*nexp]*params[nexp+v2*nexp];
    return sum;
  }
  else if( (parameter>=nexp) && (parameter<nexp+dim*nexp) )
  {
    if( (pv==v1) && (pv==v2) )
    {
      if(pe==0)
      {
        double sum=0.0;
        double E=exp(params[0]);
        sum+=exp(-E*t);
        for(int n=1; n<nexp; ++n)
        {
          E+=exp(params[n]);
          sum+=params[nexp+v1*nexp+n]*params[nexp+v2*nexp+n]*exp(-E*t);
        }
        sum*=2.0*params[nexp+v1*nexp];
        return sum;
      }
      else if( (pe>0) && (pe<nexp) )
      {
        double E=exp(params[0]);
        for(int n=1; n<=pe; ++n)
        {
          E+=exp(params[n]);
        }
        return params[nexp+v1*nexp]*params[nexp+v2*nexp]*2.0*params[nexp+v1*nexp+pe]*exp(-E*t);
      }
      else
      {
        return 0.0;
      }
    }
    else if( (pv==v1) && (pv!=v2) )
    {
      if(pe==0)
      {
        double sum=0.0;
        double E=exp(params[0]);
        sum+=exp(-E*t);
        for(int n=1; n<nexp; ++n)
        {
          E+=exp(params[n]);
          sum+=params[nexp+v1*nexp+n]*params[nexp+v2*nexp+n]*exp(-E*t);
        }
        sum*=params[nexp+v2*nexp];
        return sum;
      }
      else if( (pe>0) && (pe<nexp) )
      {
        double E=exp(params[0]);
        for(int n=1; n<=pe; ++n)
        {
          E+=exp(params[n]);
        }
        return params[nexp+v1*nexp]*params[nexp+v2*nexp]*params[nexp+v2*nexp+pe]*exp(-E*t);
      }
      else
      {
        return 0.0;
      }
    }
    else if( (pv!=v1) && (pv==v2) )
    {
      if(pe==0)
      {
        double sum=0.0;
        double E=exp(params[0]);
        sum+=exp(-E*t);
        for(int n=1; n<nexp; ++n)
        {
          E+=exp(params[n]);
          sum+=params[nexp+v1*nexp+n]*params[nexp+v2*nexp+n]*exp(-E*t);
        }
        sum*=params[nexp+v1*nexp];
        return sum;
      }
      else if( (pe>0) && (pe<nexp) )
      {
        double E=exp(params[0]);
        for(int n=1; n<=pe; ++n)
        {
          E+=exp(params[n]);
        }
        return params[nexp+v1*nexp]*params[nexp+v2*nexp]*params[nexp+v1*nexp+pe]*exp(-E*t);
      }
      else
      {
        return 0.0;
      }
    }
    else
    {
      return 0.0;
    }
  }
  else
  {
    return 0.0;
  }
}


int multi_exp_expE_mat_model::no_of_errors()
{
  return 0;
}


string multi_exp_expE_mat_model::get_function_name(int function, string variable)
{
  int v1=function/dim2;
  int v2=function%dim2;
  stringstream function_st;
  function_st << A_str << "__" << v1+1 << "*" << A_str << "__" << v2+1 << "*(exp(-exp(" << E_str << ")*" << variable << ")";
  for(int n=1; n<nexp; ++n)
  {
    function_st << "+" << B_str << "_" << n << "__" << v1+1 << "*" << B_str << "_" << n << "__" << v2+1 << "*exp(-(exp(" << E_str << ")";
    for(int k=1; k<=n; ++k)
    {
      function_st << "+exp(" << dE_str << "_" << k << ")";
    }
    function_st << ")*" << variable << ")";
  }
  function_st << ")";
  return function_st.str();
}


string multi_exp_expE_mat_model::get_function_name(int function)
{
  return get_function_name(function, t_str);
}


string multi_exp_expE_mat_model::get_parameter_name(int parameter)
{
  vector<string> names(nexp+dim*nexp, "");
  names[0]=E_str;
  for(int n=1; n<nexp; ++n)
  {
    stringstream param_st;
    param_st << dE_str << "_" << n;
    names[n]=param_st.str();
  }
  for(int v=0; v<dim; ++v)
  {
    stringstream a1st;
    a1st << A_str << "__" << v+1;
    names[nexp+v*nexp]=a1st.str();
    for(int n=1; n<nexp; ++n)
    {
      stringstream b_st;
      b_st << B_str << "_" << n << "__" << v+1;
      names[nexp+v*nexp+n]=b_st.str();
    }
  }
  return names[parameter];
}


string multi_exp_expE_mat_model::get_variable_name(int variable)
{
  if(variable!=0)
  {
    return "";
  }
  return t_str;
}

string multi_exp_expE_mat_model::get_constant_name(int constant)
{
  return "";
}
