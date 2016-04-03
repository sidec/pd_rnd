#include "m_pd.h"
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

static t_class * tabperm_class;

typedef struct _tabperm
{
    t_object x_obj;
    t_symbol *x_arrayname;
    gsl_rng * r;
    /* t_outlet * l_out; */
    t_outlet * b_out;
} t_tabperm;


void tabperm_bang(t_tabperm * x)
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
        float * table = malloc(sizeof(float)*vecsize);


        int n;
        for(n=0; n < vecsize; n++)
        {
            table[n] = vec[n].w_float;
        }
        gsl_ran_shuffle(x->r, table, vecsize, sizeof (float));
        garray_redraw(a);

        for(n=0; n < vecsize; n++)
        {
            vec[n].w_float = table[n];
        }
        free(table);
    }
}


void tabperm_set(t_tabperm *x, t_symbol *s)
{
    x->x_arrayname = s;
}


// constructor
void * tabperm_new(t_symbol * s)
{
    t_tabperm * x = (t_tabperm *) pd_new(tabperm_class);
    //  t_symbol array_name = atom_getsymbol(argv);
    x->x_arrayname = s; //atom_getsymbol(argv);// array_name->s_name;
    x->b_out = outlet_new(&x->x_obj, &s_bang);
    x->r = gsl_rng_alloc(gsl_rng_ranlxs2);
    long seed = abs((long)(rand()*LONG_MAX));
    gsl_rng_set(x->r,seed);

    return (void *) x;
}
// destroyer
void tabperm_free(t_tabperm  * x)
{
    char info[80];
    sprintf(info,"%s generator has been freed.", gsl_rng_name(x->r));
    post(info);
    gsl_rng_free(x->r);
}

void tabperm_setup(void)
{
    tabperm_class = class_new(gensym("tabperm"),
                              (t_newmethod) tabperm_new,
                              (t_method) tabperm_free,
                              sizeof(t_tabperm),
                              CLASS_DEFAULT, A_DEFSYMBOL, 0);
    class_addbang(tabperm_class, (t_method)tabperm_bang);
    class_addmethod(tabperm_class,
                    (t_method)tabperm_set,
                    gensym("set"),
                    A_DEFSYMBOL,0);



    // help file
    class_sethelpsymbol(tabperm_class, gensym("tabperm-help"));

}
