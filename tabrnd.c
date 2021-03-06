#include "m_pd.h"
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

static t_class * tabrnd_class;

typedef struct _tabrnd
{
    t_object x_obj;
    t_symbol *x_arrayname;
    gsl_rng * r;
    /* t_outlet * l_out; */
    t_outlet * b_out;
} t_tabrnd;


void tabrnd_gaussian(t_tabrnd * x,
                     t_floatarg sigma,
                     t_floatarg mu)
{
    int vecsize;
    t_garray *a;
    t_word *vec;

    if (!(a = (t_garray *)pd_findbyclass(x->x_arrayname, garray_class)))
    {
        pd_error(x, "%s: no such array", x->x_arrayname->s_name);
    }
    else if (!garray_getfloatwords(a, &vecsize, &vec))
    {
        pd_error(x, "%s: bad template for tabwrite", x->x_arrayname->s_name);
    }
    else
    {
        double sigma_ = (sigma==0)? 1.0 : (double)sigma;
        int n;
        for(n=0; n < vecsize; n++)
        {

            double rnd_x = gsl_ran_gaussian(x->r, sigma_);
            rnd_x += mu;
            vec[n].w_float = (float)rnd_x;

            /* t_atom  to_out[2]; */
            /* SETFLOAT(to_out, (float)n); */
            /* SETFLOAT(to_out+1, (float)rnd_x); */
            /* outlet_list(x->l_out, gensym("list"), 2, to_out ); */
        }
        outlet_bang(x->b_out);
        garray_redraw(a);
    }
}

void tabrnd_exponential(t_tabrnd * x,
                        t_floatarg lambda)
{
    int vecsize;
    t_garray *a;
    t_word *vec;

    if (!(a = (t_garray *)pd_findbyclass(x->x_arrayname, garray_class)))
    {
        pd_error(x, "%s: no such array", x->x_arrayname->s_name);
    }
    else if (!garray_getfloatwords(a, &vecsize, &vec))
    {
        pd_error(x, "%s: bad template for tabwrite", x->x_arrayname->s_name);
    }
    else
    {
        double l = (lambda==0)? 1.0 : (double)lambda;
        int n;
        for(n=0; n < vecsize; n++)
        {
            double rnd_x = gsl_ran_exponential(x->r, l);

            /* double pdf = gsl_ran_exponential_pdf(x_n, mu);; */
            vec[n].w_float = (float)rnd_x;

            /* t_atom  to_out[2]; */
            /* SETFLOAT(to_out, (float)x_n); */
            /* SETFLOAT(to_out+1, (float)pdf); */
            /* outlet_list(x->l_out, gensym("list"), 2, to_out ); */
        }
        outlet_bang(x->b_out);
        garray_redraw(a);
    }
}

void tabrnd_gamma(t_tabrnd * x,
                  t_floatarg alpha,
                  t_floatarg beta)
{
    int vecsize;
    t_garray *a;
    t_word *vec;

    if (!(a = (t_garray *)pd_findbyclass(x->x_arrayname, garray_class)))
    {
        pd_error(x, "%s: no such array", x->x_arrayname->s_name);
    }
    else if (!garray_getfloatwords(a, &vecsize, &vec))
    {
        pd_error(x, "%s: bad template for tabwrite", x->x_arrayname->s_name);
    }
    else
    {
        double kappa = (alpha==0)? 2.0 : (double)alpha;
        double theta = (beta==0)? 2.0 : (double)beta;
        int n;
        for(n=0; n < vecsize; n++)
        {
            double rnd_x = gsl_ran_gamma(x->r, kappa, theta);
            vec[n].w_float = (float)rnd_x;

            /* t_atom  to_out[2]; */
            /* SETFLOAT(to_out, (float)x_n); */
            /* SETFLOAT(to_out+1, (float)pdf); */
            /* outlet_list(x->l_out, gensym("list"), 2, to_out ); */
        }
        outlet_bang(x->b_out);
        garray_redraw(a);
    }
}



void tabrnd_set(t_tabrnd *x, t_symbol *s)
{
    x->x_arrayname = s;
}


// constructor
void * tabrnd_new(t_symbol * s, int argc, t_atom * argv)
{
    t_tabrnd * x = (t_tabrnd *) pd_new(tabrnd_class);
    //  t_symbol array_name = atom_getsymbol(argv);
    x->x_arrayname = atom_getsymbol(argv);// array_name->s_name;
    long seed = (argc==3)? (long)atom_getint(argv+2):abs((long)(rand()*LONG_MAX));
    /* x->l_out = outlet_new(&x->x_obj, &s_list); */
    x->b_out = outlet_new(&x->x_obj, &s_bang);
    x->r = gsl_rng_alloc(gsl_rng_ranlxs2);
    gsl_rng_set(x->r,seed);

    return (void *) x;
}
// destroyer
void tabrnd_free(t_tabrnd  * x)
{
    char info[80];
    sprintf(info,"%s generator has been freed.", gsl_rng_name(x->r));
    post(info);
    gsl_rng_free(x->r);
}

void tabrnd_setup(void)
{
    tabrnd_class = class_new(gensym("tabrnd"),
                             (t_newmethod) tabrnd_new,
                             (t_method) tabrnd_free,
                             sizeof(t_tabrnd),
                             CLASS_DEFAULT, A_GIMME, 0);
    class_addmethod(tabrnd_class,
                    (t_method)tabrnd_gaussian,
                    gensym("gaussian"),
                    A_DEFFLOAT, A_DEFFLOAT,0);
    class_addmethod(tabrnd_class,
                    (t_method)tabrnd_exponential,
                    gensym("exponential"),
                    A_DEFFLOAT,0);
    class_addmethod(tabrnd_class,
                    (t_method)tabrnd_gamma,
                    gensym("gamma"),
                    A_DEFFLOAT, A_DEFFLOAT,0);
    class_addmethod(tabrnd_class,
                    (t_method)tabrnd_set,
                    gensym("set"),
                    A_DEFSYMBOL,0);


    // help file
    class_sethelpsymbol(tabrnd_class, gensym("tabrnd-help"));

}
