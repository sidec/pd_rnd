#include "m_pd.h"
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

static t_class * tabplot_class;

typedef struct _tabplot {
  t_object x_obj;
  t_symbol *x_arrayname;
  t_outlet * l_out;
  t_outlet * b_out;
} t_tabplot;


void tabplot_gaussian(t_tabplot * x,
                      t_floatarg sigma,
                      t_floatarg x_a,
                      t_floatarg x_b){
  int i,vecsize;
  t_garray *a;
  t_word *vec;
  double delta = abs(x_b - x_a);

  if (!(a = (t_garray *)pd_findbyclass(x->x_arrayname, garray_class))){
    pd_error(x, "%s: no such array", x->x_arrayname->s_name);
  }
  else if (!garray_getfloatwords(a, &vecsize, &vec)){
    pd_error(x, "%s: bad template for tabwrite", x->x_arrayname->s_name);
  } else {
    double step = delta/vecsize;
    double x_n = (float)x_a;
    int n;
    for(n=0; n < vecsize; n++){
      x_n += step;
      double pdf = gsl_ran_gaussian_pdf(x_n, sigma);;
      vec[n].w_float = (float)pdf;

      t_atom  to_out[2];
      SETFLOAT(to_out, (float)x_n);
      SETFLOAT(to_out+1, (float)pdf);
      outlet_list(x->l_out, gensym("list"), 2, to_out );
    }
    outlet_bang(x->b_out);
    garray_redraw(a);
  }
}

void tabplot_exponential(t_tabplot * x,
                         t_floatarg mu,
                         t_floatarg x_a,
                         t_floatarg x_b){
  int vecsize;
  t_garray *a;
  t_word *vec;
  double delta = abs(x_b - x_a);

  if (!(a = (t_garray *)pd_findbyclass(x->x_arrayname, garray_class))){
    pd_error(x, "%s: no such array", x->x_arrayname->s_name);
  }
  else if (!garray_getfloatwords(a, &vecsize, &vec)){
    pd_error(x, "%s: bad template for tabwrite", x->x_arrayname->s_name);
  } else {
    double step = delta/vecsize;
    double x_n = (float)x_a;
    int n;
    for(n=0; n < vecsize; n++){
      x_n += step;
      double pdf = gsl_ran_exponential_pdf(x_n, mu);;
      vec[n].w_float = (float)pdf;

      t_atom  to_out[2];
      SETFLOAT(to_out, (float)x_n);
      SETFLOAT(to_out+1, (float)pdf);
      outlet_list(x->l_out, gensym("list"), 2, to_out );
    }
    outlet_bang(x->b_out);
    garray_redraw(a);
  }
}

void tabplot_gamma(t_tabplot * x,
                   t_floatarg alpha,
                   t_floatarg theta,
                   t_floatarg x_a,
                   t_floatarg x_b){
  int i,vecsize;
  t_garray *a;
  t_word *vec;
  double delta = abs(x_b - x_a);

  if (!(a = (t_garray *)pd_findbyclass(x->x_arrayname, garray_class))){
    pd_error(x, "%s: no such array", x->x_arrayname->s_name);
  }
  else if (!garray_getfloatwords(a, &vecsize, &vec)){
    pd_error(x, "%s: bad template for tabwrite", x->x_arrayname->s_name);
  } else {
    double step = delta/vecsize;
    double x_n = (float)x_a;
    int n;
    for(n=0; n < vecsize; n++){
      x_n += step;
      double pdf = gsl_ran_gamma_pdf(x_n, alpha, theta);;
      vec[n].w_float = (float)pdf;

      t_atom  to_out[2];
      SETFLOAT(to_out, (float)x_n);
      SETFLOAT(to_out+1, (float)pdf);
      outlet_list(x->l_out, gensym("list"), 2, to_out );
    }
    outlet_bang(x->b_out);
    garray_redraw(a);
  }
}



void tabplot_set(t_tabplot *x, t_symbol *s)
{
    x->x_arrayname = s;
}


// constructor
void * tabplot_new(t_symbol * s){
  t_tabplot * x = (t_tabplot *) pd_new(tabplot_class);
  x->x_arrayname = s;
  x->l_out = outlet_new(&x->x_obj, &s_list);
  x->b_out = outlet_new(&x->x_obj, &s_bang);
  return (void *) x;
}

void tabplot_setup(void){
  tabplot_class = class_new(gensym("tabplot"),
                            (t_newmethod) tabplot_new,
                            0,
                            sizeof(t_tabplot),
                            CLASS_DEFAULT, A_DEFSYM, 0);
  class_addmethod(tabplot_class,
                  (t_method)tabplot_gaussian,
                  gensym("gaussian"),
                  A_DEFFLOAT, A_DEFFLOAT, A_DEFFLOAT,0);
  class_addmethod(tabplot_class,
                  (t_method)tabplot_exponential,
                  gensym("exponential"),
                  A_DEFFLOAT, A_DEFFLOAT, A_DEFFLOAT,0);
  class_addmethod(tabplot_class,
                  (t_method)tabplot_gamma,
                  gensym("gamma"),
                  A_DEFFLOAT, A_DEFFLOAT, A_DEFFLOAT, A_DEFFLOAT,0);


    // help file
  class_sethelpsymbol(tabplot_class, gensym("tabplot-help"));

}
