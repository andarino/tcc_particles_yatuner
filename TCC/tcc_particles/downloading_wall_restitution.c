//para compilar: gcc dowloading_wall_resitution.c -o tt -lgsl -lgslcblas -lm
//               ./tt (tempo total) (graos na vertical) (graos na horizontal) (parâmetro de afilamento)

/*
* =====================================================================================
*
*       Filename:  dowloading_wall_resitution.c
*
*    Description:
*
*        Version:  0.2
*        Created:  19-07-2017 13:11:41 BRT
*       Revision:  none
*       Compiler:  gcc
*
*         Author:  Luis Paulo Machado
*        Company:  Universidade Federal da Pará
*
* =====================================================================================
*/
#include <stdio.h>
#include <stdlib.h>		//também trata de alocação de memória
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <unistd.h>		//para rodar o gnuplot
#include <time.h>

int _num_layer,_grain_layer;	//número de camadas e número de grãos em cada camada
double _q,RESTITUTION,DOWLOADING_WALL;

//N=1+(número de camadas x grãos por camada)+((número de camadas-1)x(grãos por camada+1))
//1=grão incidente
//número de grãos principais=(número de camadas x grãos por camada)
//número de grãos intersticiais=(número de camadas-1)x(grãos por camada+1)

#define N      ((2*_num_layer*_grain_layer)+_num_layer-_grain_layer)
#define STANDARD_G 9.80665
#define W_INITIAL 0.999

void novos_vetores ();

void interaction ();
void interaction_inter ();
void top_wall ();
void bottom_wall ();
void temporary_wall();

void energy_inter ();
void energy_temporary_wall();
void energy_between_grains ();
void energyDOWLOADING_WALL ();
void energy_wall_bottom ();
void KE_GE_MOMENTUM ();

double _r,_d,_delta;


typedef struct
{  double *alpha,*beta,*temporary_wall,*alpha_inter,*a,*b,*radii,*mass;
  double *max_alpha,*max_beta,*max_temporary_wall,*max_alpha_inter;
  double *RESTITUTION;//,*dowloading_wall;
} data;

void novos_vetores (double *x,const double *y,void *params,double *fx,double *fy)
{
 int n;  data edo_params = *(data *) params;
 for (n = 0; n < N; n++)
   {
     x[n]=edo_params.a[n]+y[n];	x[N+n]=edo_params.b[n]+y[N+n];
     x[2*N+n]=y[2*N+n];		x[3*N+n]=y[3*N+n];
     fx[n]=0.;				fy[n]=0.;
   }
}

void interaction (const double *x, void *params,int n1, int n2,
     double *fx, double *fy, int alfa)
{
 data edo_params = *(data *) params;
 double d_rel_x, d_rel_y,v_rel;

 d_rel_x = x[n2] - x[n1];  d_rel_y = x[N + n2] - x[N + n1];
 _r = hypot (d_rel_x, d_rel_y);

 _delta = edo_params.radii[n1]+edo_params.radii[n2] - _r;

 if(_delta<=0.)_d=edo_params.max_alpha[alfa]=0.;
 else
   {
     if(_delta>edo_params.max_alpha[alfa])
 {
   edo_params.max_alpha[alfa]=_delta;
   _d = edo_params.alpha[alfa]*pow(_delta,1.5)/_r;
 }
     else _d=(1.-RESTITUTION)*edo_params.alpha[alfa]
                                                 *pow(_delta,1.5)/_r;
   }//printf("%g\t%d\n",RESTITUTION,alfa);

 fx[n2]+= _d * d_rel_x;   fy[n2]+=_d * d_rel_y;
 fx[n1]-=_d * d_rel_x;	   fy[n1]-=_d * d_rel_y;

}

void interaction_inter (const double *x, void *params,int n1, int n2,
                       double *fx, double *fy, int alfa)
{
 data edo_params = *(data *) params;
 double d_rel_x, d_rel_y,v_rel;

 d_rel_x = x[n2] - x[n1];  d_rel_y = x[N + n2] - x[N + n1];
 _r = hypot (d_rel_x, d_rel_y);

 _delta = edo_params.radii[n1]+edo_params.radii[n2] - _r;

 if(_delta<=0.)_d=edo_params.max_alpha_inter[alfa]=0.;
 else
   {
     if(_delta>edo_params.max_alpha_inter[alfa])
 {
   edo_params.max_alpha_inter[alfa]=_delta;
   _d = edo_params.alpha_inter[alfa]*pow(_delta,1.5)/_r;
 }
     else _d=(1.-RESTITUTION)*edo_params.alpha[alfa]*
                                                  pow(_delta,1.5)/_r;
   }

 fx[n2]+= _d * d_rel_x;   fy[n2]+=_d * d_rel_y;
 fx[n1]-=_d * d_rel_x;	   fy[n1]-=_d * d_rel_y;
}

void top_wall (const double *x,void *params,int n,double *fy,int i)
{

 data edo_params = *(data *) params;	//parede de cima
 _delta= (x[N+n]+edo_params.radii[n])-(edo_params.b[1]+edo_params.radii[1]);

 if(_delta<=0.) {fy[n]+=0.; edo_params.max_beta[i]=0.;}
 else
   {
     if(_delta>edo_params.max_beta[i])
 {
   edo_params.max_beta[i]=_delta;
   fy[n] -= edo_params.beta[i]*pow(_delta,1.5);
 }
     else fy[n]-=(1.-RESTITUTION)*edo_params.beta[i]
                                                    *pow(_delta,1.5);
   }
}

void bottom_wall (const double *x,void *params,int n,double *fy,int i)
{
 data edo_params = *(data *) params;	//parede de baixo

 _delta= (x[N+n]-edo_params.radii[n])-
   (edo_params.b[_grain_layer]-edo_params.radii[_grain_layer]);

 if(_delta>=0.){fy[n]+=0.; edo_params.max_beta[i]=0.;}
 else
   {
     if(_delta<edo_params.max_beta[i])
 {
   edo_params.max_beta[i]=_delta;
   fy[n] += edo_params.beta[i]*pow(-_delta,1.5);
 }
     else fy[n]+=(1.-RESTITUTION)*edo_params.beta[i]
                                                   *pow(-_delta,1.5);
   }
}

void temporary_wall(const double *y,void *params,int n_bottom,
       int n_inter,int n_top,double *fx,int alfa_inter,int alfa_top)
{
 data edo_params = *(data *) params;	//parede de baixo
 int k=n_bottom+1; double wall;

 //--------------------------------------------------------------------
 if(y[k]>=0.){fx[k]+=0.; edo_params.max_temporary_wall[n_bottom]=0.;}
 else
   {
     if(y[k]<edo_params.max_temporary_wall[n_bottom])
 {
   edo_params.max_temporary_wall[n_bottom]=y[k];
   fx[k] += edo_params.temporary_wall[n_bottom]*pow(-y[k],1.5);
 }
     else fx[k]+=(1.-RESTITUTION)
                 *edo_params.temporary_wall[n_bottom]*pow(-y[k],1.5);
   }
 //--------------------------------------------------------------------
 wall=edo_params.a[n_inter]+edo_params.radii[n_inter]+y[n_inter]
                                                    -DOWLOADING_WALL;
 if(wall<=0.)
   {fx[n_inter]+=0.; edo_params.max_temporary_wall[alfa_inter]=0.;}
 else
   {
     if(wall>edo_params.max_temporary_wall[alfa_inter])
 {
   edo_params.max_temporary_wall[alfa_inter]=wall;
   fx[n_inter]-=edo_params.temporary_wall[alfa_inter]*pow(wall,1.5);
 }
     else fx[n_inter]-=(1.-RESTITUTION)
                *edo_params.temporary_wall[alfa_inter]*pow(wall,1.5);
   }
 //--------------------------------------------------------------------
 wall=edo_params.a[n_top]+edo_params.radii[n_top]+y[n_top]
                                                    -DOWLOADING_WALL;
 if(wall<=0.)
   {fx[n_top]+=0.; edo_params.max_temporary_wall[alfa_top]=0.;}
 else
   {
     if(wall>edo_params.max_temporary_wall[alfa_top])
 {
   edo_params.max_temporary_wall[alfa_top]=wall;
   fx[n_top]-=edo_params.temporary_wall[alfa_top]*pow(wall,1.5);
 }
     else fx[n_top]-=(1.-RESTITUTION)*
                   edo_params.temporary_wall[alfa_top]*pow(wall,1.5);
   }
 //--------------------------------------------------------------------
 if(alfa_top==2*_grain_layer+1)
   {
     wall=edo_params.a[n_inter-1]+edo_params.radii[n_inter-1]+
                                  y[n_inter-1]-DOWLOADING_WALL;

     if(wall<=0.)
 {fx[n_inter-1]+=0.;edo_params.max_temporary_wall[alfa_inter-1]=0.;}
     else
 {
   if(wall>edo_params.max_temporary_wall[alfa_inter-1])
     {
       edo_params.max_temporary_wall[alfa_inter-1]=wall;
       fx[n_inter-1]-=edo_params.temporary_wall[alfa_inter-1]*pow(wall,1.5);
     }
   else fx[n_inter-1]-=(1.-RESTITUTION)*
               edo_params.temporary_wall[alfa_inter-1]*pow(wall,1.5);
 }
   }
}

int func (double t,const double y[],double f[],void *params)
{
 double x[4*N],full_forceX[N],full_forceY[N];
 int l, p, k, i, j, n;

 data edo_params = *(data *) params;
 novos_vetores (x, y, &edo_params,full_forceX,full_forceY);

 //----------------------Eixo 'X'-----------------------------------------------
 //ENTRE OS GRÃOS PRINCIPAIS EM TODAS AS CADEIAS
 //  printf("\nGRAOS PRINCIPAIS EM TODAS AS CADEIAS\n");

 i=3*_grain_layer; //QUANTIDADE DE GRAOS INICIAIS
 for (k = 1; k < _num_layer; k++)	//camada dos grãos principais
   {
     p=(k-1)*(2*_grain_layer+1)+1;
     l=k*(2*_grain_layer+1)+1;
     for (j = 0; j < _grain_layer; j++)
 {
   interaction(x,&edo_params,p,l,full_forceX,full_forceY,i); //0=alfa
   i++;p++;l++;
   }//printf("\n");
     i+=5*_grain_layer-1;
   }//exit(1);
 //printf("-------------------------------------------\n");

 //ENTRE OS GRÃOS PRINCIPAIS E DECORADORES DA FRENTE
 //  printf("\nGRAOS PRINCIPAIS E DECORADORES DA FRENTE\n");

 i=_grain_layer;   n = _grain_layer + 1;
 for (k = 1; k < _num_layer; k++)	//camada dos grãos principais
   {
     p=(k-1)*(2*_grain_layer+1)+1;
     for (j = 0; j < _grain_layer; j++)
 {
   for (l = 1; l < 3; l++)	//camada dos grãos principais
     {
       interaction(x,&edo_params,p,n,full_forceX,full_forceY,i);
       i++; n++;
     }//printf("\n");
         p++; n--;
 }//printf("------------------------------------------------\n");

     n += _grain_layer + 1;
     i+=4*_grain_layer-1;
   }//exit(1);
 //printf("***********************************************\n");//exit


 //ENTRE OS GRÃOS INTERTICIAIS E OS PRINCIPAIS DA FRENTE
 //printf("\nGRAOS INTERTICIAIS E OS PRINCIPAIS DA FRENTE\n");

 i=4*_grain_layer;   n = _grain_layer + 1;
 for (k = 1; k < _num_layer; k++)	//camada dos grãos principais
   {
     l=k*(2*_grain_layer+1)+1;
     for (j = 0; j < _grain_layer; j++)
 {
   for (p = 1; p < 3; p++)	//camada dos grãos principais
     {
       interaction(x,&edo_params,n,l,full_forceX,full_forceY,i);
       n++; i++;
     }//printf("\n");
         n--; l++;
 }//printf("\n");
     n += _grain_layer + 1;
     i+=4*_grain_layer-1;
   }//exit(1);

 ////////////// ENTRE OS GRÃOS INTERTICIAIS /////////////////////////
 //printf("\nGRAOS INTERTICIAIS\n");

 n = _grain_layer + 1;//primeiro grão intersticial
 i = 0;
 //  printf("grao principal e o grao decorador da frente\n");

 for (k = 1; k < _num_layer-1; k++)	//camada dos grãos principais
   {
     for (j = 0; j < _grain_layer; j++)
 {
   interaction_inter(x,&edo_params,n,n+1,full_forceX,full_forceY,i+j);
   interaction_inter(x,&edo_params,n,n+2*_grain_layer+1,
                           full_forceX,full_forceY,i+j+_grain_layer);
   n++;
 }
     interaction_inter(x,&edo_params,n,n+2*_grain_layer+1,
     full_forceX,full_forceY,i+j+_grain_layer);
     n += _grain_layer + 1;
     i += 2*_grain_layer + 1;
   }

 for (j = 0; j < _grain_layer; j++)
   {interaction_inter(x,&edo_params,n,n+1,full_forceX,full_forceY,i+j);
     n++;}//exit(1);


 //----------------------Eixo 'Y'-----------------------------------------------
 //ENTRE OS GRÃOS PRINCIPAIS EM TODAS AS CAMADAS
 //printf("\nGRAOS PRINCIPAIS EM TODAS AS CAMADAS\n");

 l = p= 1;
 //  printf("\ngraos principais em uma MESMA camada\n");
 for (k = 1; k < (_num_layer + 1); k++)	//k é o número de camadas
   {
     for (j = 0; j < _grain_layer-1; j++)
 {interaction(x,&edo_params,l+j,l+j+1,full_forceX,full_forceY,p); //0=alfa
   p++;}//printf("\n");
     p+= 5*_grain_layer;
     l=k*(2*_grain_layer+1)+1;
   }//exit(1);
 //printf("-------------------------------------------\n");

 //PARA A PAREDE DE CIMA
 //printf("Paredes de CIMA\n");

 i=0;
 for(k=1;k<_num_layer;k++)
   { p=(k-1)*(2*_grain_layer+1)+1;

     top_wall(x,&edo_params,p,full_forceY,i); //grãos principais
     top_wall(x,&edo_params,p+_grain_layer,full_forceY,i+2); //grãos decoradores
     i+=4;
   }

 p=(k-1)*(2*_grain_layer+1)+1;
 top_wall(x,&edo_params,p,full_forceY,i); //grão principal
 //exit(1);

 //PARA A PAREDE DE BAIXO
 //printf("Paredes de BAIXO\n");

 i=1;
 for(k=1;k<_num_layer;k++)
   {p=_grain_layer*(2*k-1)+k-1;

     bottom_wall(x,&edo_params,p,full_forceY,i); //grãos principais
     bottom_wall(x,&edo_params,p+_grain_layer+1,full_forceY,i+2); //grãos intersticicias
     i+=4;
   }

 p=_grain_layer*(2*k-1)+k-1;
 bottom_wall(x,&edo_params,p,full_forceY,i); //grãos principais
 //exit(1);

 //PAREDE TEMPORÁRIA DA PRIMEIRA CAMADA
 for(k=0;k<_grain_layer;k++)
   temporary_wall(y,&edo_params,k,N-1-_grain_layer-k,N-1-k,
      full_forceX,2*_grain_layer-k,(3*_grain_layer)-k);//exit(1);
 /*----------------------------------------------------------------------------
   ------------ATUALIZAÇÃO DAS COMPONENTES NO EIXO 'X' e 'Y'---------------------
   -----------------------------------------------------------------------------*/

 for(n=1;n<N;n++)
   {
     f[n] = x[2*N+n];
     f[2*N+n] = (full_forceX[n]/edo_params.mass[n])-STANDARD_G;

     f[N+n] = x[3*N+n];
     f[3*N+n] = full_forceY[n] / edo_params.mass[n];

   } f[0]=f[2*N]=f[N]=f[3*N]=0.;//printf("-------------------\n");
 //  exit(1);   //printf("%g\t%g\n",full_forceX[0],full_forceX[1]);

 return GSL_SUCCESS;
}

void energy_between_grains (const double *x,void *params,int n1,int n2,double *PE,int alfa)
{
 data edo_params = *(data *) params;
 double d_rel_x, d_rel_y; //n1 e n2 são os números dos grãos

 d_rel_x= x[n2]-x[n1];		d_rel_y=x[N+n2]-x[N+n1];
 _r = hypot (d_rel_x, d_rel_y);

 _delta = edo_params.radii[n1]+edo_params.radii[n2] - _r;

 if(_delta<=0.)_d=0.;
 else
   {
     if(_delta>edo_params.max_alpha[alfa])
 _d = 0.2*edo_params.alpha[alfa]*pow(_delta,2.5);
     else _d=(1.-RESTITUTION)*0.2*edo_params.alpha[alfa]
                                                    *pow(_delta,2.5);
   }

 PE[n1]+= _d;     PE[n2]+= _d;
}

void energy_inter(const double *x,void *params,int n1,int n2,double *PE,int alfa)
{
 data edo_params = *(data *) params;
 double d_rel_x, d_rel_y;

 d_rel_x = x[n2] - x[n1];  d_rel_y = x[N + n2] - x[N + n1];
 _r = hypot (d_rel_x, d_rel_y);

 _delta = edo_params.radii[n1]+edo_params.radii[n2] - _r;

 if(_delta<=0.)_d=0.;
 else
   {
     if(_delta>edo_params.max_alpha_inter[alfa])
 _d = 0.2*edo_params.alpha[alfa]*pow(_delta,2.5);
     else _d=(1.-RESTITUTION)*0.2*edo_params.alpha[alfa]
                                                    *pow(_delta,2.5);
   }

 PE[n1]+= _d;     PE[n2]+= _d;
}

void energyDOWLOADING_WALL(const double *y,void *params,int n,double *PE,int i)
{
 data edo_params = *(data *) params;	//parede de cima
 _delta= (y[N+n]+edo_params.radii[n])-(edo_params.b[1]+edo_params.radii[1]);

 if(_delta<=0.) PE[n]+=0.;
 else
   {
     if(_delta>edo_params.max_beta[i])
   PE[n] += 0.4*edo_params.beta[i]*pow(_delta,2.5);

     else PE[n] +=(1.-RESTITUTION)*0.4*edo_params.beta[i]
                                                    *pow(_delta,2.5);
   }

 //  printf("PE_%d=%g\tbeta%d=%g\ty_%d=%g\n",n,PE[n],i,edo_params.beta[i],n,y[N+n]);
}

void energy_wall_bottom(const double *y,void *params,int n,double *PE,int i)
{
 data edo_params = *(data *) params;

 _delta= (y[N+n]-edo_params.radii[n])-
         (edo_params.b[_grain_layer]-edo_params.radii[_grain_layer]);

 if(_delta>=0.)PE[n]+=0.;
 else
   {
     if(_delta<edo_params.max_beta[i])
   PE[n] += 0.4*edo_params.beta[i]*pow(-_delta,2.5);
     else PE[n]+=(1.-RESTITUTION)*0.4*edo_params.beta[i]
                                                   *pow(-_delta,2.5);
   }
}

void KE_GE_MOMENTUM(const double *y,void *params,int n,double *KE,double *GE,double *P)
{
 data edo_params = *(data *) params;	//parede de cima
 KE[n] = 0.5*edo_params.mass[n]*(y[2*N+n]*y[2*N+n]+y[3*N+n]*y[3*N+n]);
 GE[n] = edo_params.mass[n]*STANDARD_G*y[n];
 P[n]  = edo_params.mass[n]*y[2*N+n];
 P[N+n]= edo_params.mass[n]*y[3*N+n];
}

void energy_temporary_wall(const double *y,void *params,int n_bottom,
       int n_inter,int n_top,double *PE_wall,int alfa_inter,int alfa_top)
{
 data edo_params = *(data *) params;	//parede de baixo
 int k=n_bottom+1; double wall;

 //--------------------------------------------------------------------
 if(y[k]>=0.)PE_wall[k]+=0.;
 else
   {
     if(y[k]<edo_params.max_temporary_wall[n_bottom])
       PE_wall[k]+=0.4*edo_params.temporary_wall[n_bottom]*pow(-y[k],2.5);

     else PE_wall[k]+=(1.-RESTITUTION)
             *0.4*edo_params.temporary_wall[n_bottom]*pow(-y[k],2.5);
   }
 //--------------------------------------------------------------------
 wall=edo_params.a[n_inter]+edo_params.radii[n_inter]+y[n_inter]
      -DOWLOADING_WALL;
 if(wall<=0.) PE_wall[n_inter]+=0.;
 else
   {
     if(wall>edo_params.max_temporary_wall[alfa_inter])
   PE_wall[n_inter] += 0.4*edo_params.temporary_wall[alfa_inter]
                                                      *pow(wall,2.5);
     else PE_wall[n_inter]+=(1.-RESTITUTION)
                *0.4*edo_params.temporary_wall[alfa_inter]*pow(wall,2.5);
   }
 //--------------------------------------------------------------------
 wall=edo_params.a[n_top]+edo_params.radii[n_top]+y[n_top]
      -DOWLOADING_WALL;

 if(wall<=0.) PE_wall[n_top]+=0.;
 else
   {
     if(wall>edo_params.max_temporary_wall[alfa_top])
   PE_wall[n_top]+=    0.4*edo_params.temporary_wall[alfa_top]
                                                  *pow(wall,2.5);
     else PE_wall[n_top]+=(1.-RESTITUTION)*0.4*
               edo_params.temporary_wall[alfa_top]*pow(wall,2.5);
   }
 //--------------------------------------------------------------------

 if(alfa_top==2*_grain_layer+1)
   {
     wall=edo_params.a[n_inter-1]+edo_params.radii[n_inter-1]+
                                              y[n_inter-1]-DOWLOADING_WALL;
     if(wall<=0.) PE_wall[n_inter-1]+=0.;
     else
 {
   if(wall>edo_params.max_temporary_wall[alfa_inter-1])
     PE_wall[n_inter-1]+=0.4*edo_params.temporary_wall[alfa_inter-1]
                                                      *pow(wall,2.5);
   else PE_wall[n_inter-1]+=(1.-RESTITUTION)*
           0.4*edo_params.temporary_wall[alfa_inter-1]*pow(wall,2.5);
 }
   }

}

int main (int argc, char **argv)
{

int n, k, i, j, l, p;
 double t, h, Tmax;
 double *x, *x_err, *dydt_in, *dydt_out, *f, *swap, *a, *b;
 double *alpha, *beta,*temporary_wall,*alpha_inter;
 double *radii, *mass, *several_dates;
 data edo_params;
 const gsl_odeiv_step_type *T = gsl_odeiv_step_rk4;
 double *max_alpha,*max_beta,*max_temporary_wall,*max_alpha_inter;

 //////////////// Coeficiente da Dissipação /////////////////////////
 RESTITUTION=W_INITIAL;

 if (argc < 1)
   {printf ("sintax: Em %s determine o valor de S,_lambda e Tmax!!!\n",
      argv[0]);exit (1);}

 Tmax = atof (argv[1]);
 if (Tmax < 0){printf("O tempo deve ser maior\n");exit (1);}

 _num_layer = atoi (argv[2]);
 if (_num_layer < 1)
   { printf ("Deve ter + de 2 camadas. Cadeias na vertical.\n");
     exit (1);}

 _grain_layer = atoi (argv[3]);
 if (_grain_layer < 2)
   { printf
 ("Deve ter + de 2 graos em cada camada. Cadeias na horizontal.\n");
     exit (1);}

 _q = atof (argv[4]);
 if (_q<0){printf("O parametro de afilamento deve ser maior que zero.\n");
   exit (1);}


 printf ("\nNumero de graos: %d\n\n", N);	//exit(1);

 // alocação dinâmica de memória
 if ((radii = (double *) malloc (N * sizeof (double)))==NULL)
   {printf ("cannot allocate memory\n");exit (1);}
 if ((mass = (double *) malloc (N * sizeof (double)))==NULL)
   {printf ("cannot allocate memory\n");exit (1);}
 if ((a = (double *) malloc ((N) * sizeof (double)))==NULL)
   {printf ("cannot allocate memory\n");exit (1);}
 if ((b = (double *) malloc ((N) * sizeof (double)))==NULL)
   {printf ("cannot allocate memory\n");exit (1);}
 if ((x = (double *) malloc ((4 * N) * sizeof (double)))==NULL)
   {printf ("cannot allocate memory\n");exit (1);}
 if ((several_dates=(double *)malloc((4*N)*sizeof(double)))==NULL)
   {printf ("cannot allocate memory\n");exit (1);}
 if ((x_err = (double *) malloc ((4 * N) * sizeof (double)))==NULL)
   {printf ("cannot allocate memory\n");exit (1);}
 if ((dydt_in = (double *) malloc ((4 * N) * sizeof (double)))==NULL)
   {printf ("cannot allocate memory\n");exit (1);}
 if ((dydt_out = (double *) malloc ((4 * N) * sizeof (double)))==NULL)
   {printf ("cannot allocate memory\n");exit (1);}
 if ((f = (double *) calloc ((4 * N), sizeof (double)))==NULL)
   {printf ("cannot allocate memory\n");exit (1);}

//ELASTIC CONSTANTS
 if((alpha=(double*)malloc(
   ((_num_layer-1)*(6*_grain_layer-1)+_grain_layer)*sizeof(double)))==NULL)
   {printf ("cannot allocate memory\n");exit (1);}
 if ((beta=(double*)malloc((4*_num_layer-2)*sizeof(double)))==NULL)
   {printf ("cannot allocate memory\n");exit (1);}
 if ((temporary_wall=(double*)malloc((3*_grain_layer+1)*sizeof(double)))==NULL)
   {printf ("cannot allocate memory\n");exit (1);}
 if ((alpha_inter=(double*)malloc(
    (_grain_layer*(_num_layer-1)+(_grain_layer+1)*(_num_layer+1))
     *sizeof(double)))==NULL){printf("cannot allocate memory\n");exit(1);}

//MAXIMUM COMPRESSIONS
 if((max_alpha=(double*)malloc(
   ((_num_layer-1)*(6*_grain_layer-1)+_grain_layer)*sizeof(double)))==NULL)
   {printf ("cannot allocate memory\n");exit (1);}
 if ((max_beta=(double*)malloc((4*_num_layer-2)*sizeof(double)))==NULL)
   {printf ("cannot allocate memory\n");exit (1);}
 if ((max_temporary_wall=(double*)malloc((3*_grain_layer+1)*sizeof(double)))==NULL)
   {printf ("cannot allocate memory\n");exit (1);}
 if ((max_alpha_inter=(double*)malloc(
    (_grain_layer*(_num_layer-1)+(_grain_layer+1)*(_num_layer+1))
     *sizeof(double)))==NULL){printf("cannot allocate memory\n");exit(1);}

 /******************************************************************/
 /************* MASSAS, RAIOS E POSIÇÕES INICIAIS ******************/
 /******************************************************************/
 FILE *ler;         char name1[100];
 double random_dates;  int ok;

 sprintf(name1, "initial_random_dates_square_q%g_%dx%d.dat",
   _q,_num_layer,_grain_layer);

 if (( ler = fopen(name1,"r")) == NULL )
   {printf("O arquivo da massa, raio e posicoes"
     "\nNAO pode ser aberto!!!\n");exit(1);}

 //posições iniciais
 n=0;
 while (1)
   { //printf("aqui %d\n",n);
     ok = fscanf (ler, "%lg",&random_dates);
     if( (n==(4*N)) || (ok != 1) )  break;
     several_dates[n]=random_dates;    n ++;
   } fclose(ler); //for(n=0;n<(4*N);n++) printf("%d\t%g\n",n,several_dates[n]); exit(1);

 for(n=0;n<N;n++)
   {mass[n]=several_dates[n];     	radii[n]=several_dates[N+n];
     a[n]=several_dates[(2*N)+n];	        b[n]=several_dates[(3*N)+n];
     x[n]=x[N+n]=x[2*N+n]=x[3*N+n]=0.;}   free(several_dates);

 /*   for (n = 0; n < N; n++)
      printf("a%d=%g\t\tb%d=%g\t\tr%d=%g\t\tm%d=%g\n",
      n,a[n],n,b[n],n,radii[n],n,mass[n]); exit(1);*/

 /******************************************************************/
 /*************CÁLCULO DAS CONSTANTES ELÁSTICAS*********************/
 /******************************************************************/
 FILE *reading;         char name20[100];
 double *several_constants; double INITIAL_WALL;

 sprintf(name20, "random_constants_square_q%g_%dx%d.dat",
   _q,_num_layer,_grain_layer);

 if (( reading = fopen(name20,"r")) == NULL )
   {printf("O arquivo das constantes""\nNAO pode ser aberto!!!\n");
    exit(1);}

 //Lê as posições iniciais de um arquivo
 i=(_num_layer-1)*(6*_grain_layer-1)+_grain_layer;
 l=4*_num_layer-2;
 j=3*_grain_layer+1;
 k=_grain_layer*(_num_layer-1)+(_grain_layer+1)*(_num_layer-2);
 p=i+l+j+k+1; //'2' is due the dissipative constants and '1' is for 'DOWLOADING_WALL'

 if ((several_constants=(double *)malloc((p)*sizeof(double)))==NULL)
   {printf ("cannot allocate memory\n");exit (1);}


 n=0;
 while (1)
   { //printf("aqui %d\n",n);
     ok = fscanf (reading, "%lg",&random_dates);
     if( (n==p) || (ok != 1) )  break;
     several_constants[n]=random_dates;    n++;
   } fclose(reading); //for(n=0;n<p;n++) printf("%d\t%g\n",n,several_constants[n]); exit(1);

 INITIAL_WALL=DOWLOADING_WALL=several_constants[0]; i++;
 //Elastic Constants
 for(n=1;n<i;n++)
  {alpha[n-1]=several_constants[n];max_alpha[n-1]=0.;
                  /*printf("alfa%d=%g\n",n-1,alpha[n-1]);*/} i+=l;

 p=0; for(n;n<i;n++)
  {beta[p]=several_constants[n]; max_beta[p]=0.;
                      /*printf("beta%d=%g\n",p,beta[p]);*/p++;} i+=j;

 p=0;  for(n;n<i;n++)
  {temporary_wall[p]=several_constants[n];max_temporary_wall[p]=0.;
                  /*printf("tem%d=%g\n",p,temporary_wall[p]);*/p++;}i+=k;

 p=0;   for(n;n<i;n++)
  {alpha_inter[p]=several_constants[n];max_alpha_inter[p]=0.;
                  /*printf("inter%d=%g\n",p,alpha_inter[p]);*/p++;}
i+=(_num_layer-1)*(6*_grain_layer-1)+_grain_layer;
free(several_constants);// exit(1);

 /*------------------ Cinemática Inicial --------------------------*/
 for (n = 0; n < N; n++) x[n] = x[N+n]=x[2*N+n]=x[3*N+n]=0.;


 // inializa parâmetros para sistema de EDOs
 edo_params.alpha = alpha;  edo_params.beta = beta;
 edo_params.temporary_wall = temporary_wall;
 edo_params.alpha_inter = alpha_inter;

 edo_params.max_alpha = max_alpha;  edo_params.max_beta = max_beta;
 edo_params.max_temporary_wall = max_temporary_wall;
 edo_params.max_alpha_inter = max_alpha_inter;

 edo_params.a = a;        edo_params.b = b;
 edo_params.mass = mass;  edo_params.radii = radii;

 t = 0.;  h = 1e-8;			//tempo inicial; passo do RK;

 // inicializa variáveis da biblioteca GSL para solução de EDOs
 gsl_odeiv_step *s = gsl_odeiv_step_alloc (T, 4 * N);
 gsl_odeiv_system sys = { func, NULL, 4 * N, &edo_params };
 GSL_ODEIV_FN_EVAL (&sys, t, x, dydt_in);

 int grain1, grain2, alfa, alfa1, alfa2;
 double elastic[N], gravity[N], kinetic[N],momentum[2*N],y[4*N];
 double KINETIC,ELASTIC,GRAVITY,MOMENTUM_X,MOMENTUM_Y;
 double E_INITIAL,TOTAL;

 double energy_before=100.,tolerance=5.e-4,relative_erro;
//  double inter_t=0.01,print=(inter_t/h);
 double inter_t=1e-2,print=(inter_t/h);
 double wall_t=1e-4,verify_wall=(wall_t/h);
 int ciclo=0,num=1000,n_inicial,flag_break=1,max_flag_break=2;
 n_inicial = num;

 FILE *sequence1;                 char name2[100];
 FILE *position_dissipated;       char name7[100];
 FILE *energy;                    char name8[100];
 FILE *sequence2;                 char name4[100];

 sprintf(name7, "step0_position_after_restitution_%dx%d.dat",_num_layer,_grain_layer);
 sprintf(name8, "energy_behaviour_%dx%d.dat",_num_layer,_grain_layer);

 if ((energy=fopen(name8,"w"))==NULL) {
 printf("erro no print de 'energy_behaviour.dat'!!!\n");exit(1);}
 fprintf(energy,"#D_e\tE\terro\tKE\tEE\tGE\tt\tw\tWALL\n");//exit(1);
 fclose(energy);

 sprintf(name4, "sequencia_%dx%d.plot", _num_layer,_grain_layer);
 if ((sequence2 = fopen(name4,"w"))==NULL) {
 printf("erro no print de 'sequencia' do GIF!!!\n");exit(1);}
 fclose(sequence2);

 double WALL_TOP;
 while (t <= Tmax)
   {
     /*********************** DOWLOADING WALL *********************/
     if ((ciclo%(int) verify_wall)==0)
 {
   novos_vetores (y, x, &edo_params,elastic,kinetic); WALL_TOP=0.;
   for (n = 0; n < _grain_layer; n++)
     {
       k=N-1-n; p=N-1-_grain_layer-n;
       if(y[k]+radii[k]>WALL_TOP){WALL_TOP=y[k]+radii[k];}
       if(y[p]+radii[p]>WALL_TOP){WALL_TOP=y[p]+radii[p];}
     }
   p=N-1-_grain_layer-n;
   if(y[p]+radii[p]>WALL_TOP){WALL_TOP=y[p]+radii[p];}
   if(WALL_TOP<DOWLOADING_WALL)DOWLOADING_WALL=WALL_TOP;
   //printf("wall=%g\tt=%g\n",DOWLOADING_WALL,t);
       }

     if ((ciclo%(int) print)==0)
 {
         /*-------------- Dados para fazer o filme ----------------*/
         novos_vetores (y, x, &edo_params,elastic,kinetic);

   sprintf (name2, "%d_random_square_q%g_%dx%d.dat", num,_q,_num_layer,_grain_layer);
   if ((sequence1 = fopen (name2, "w+")) == NULL)
     {printf ("erro no print da saída!!!\n");exit (1);} num++;

         for (n = 0; n < N; n++)
     fprintf (sequence1, "%g \t %g \t %g\t %d\n",
        y[n]/radii[0],y[N+n]/radii[0],radii[n]/radii[0],n);
         fclose(sequence1);
         /*--------------------------------------------------------*/

   //////////////////////////////////////////////////////////////////////
   /////////// ENTRE OS GRÃOS PRINCIPAIS EM TODAS AS CAMADAS ////////////
   //////////////////////////////////////////////////////////////////////
   //printf("\ngraos principais em uma MESMA camada\n");

   l = p= 1;
   for (k = 1; k < (_num_layer + 1); k++)	//k é o número de camadas
     {
       for (j = 0; j < _grain_layer-1; j++)
   {
     energy_between_grains(y,&edo_params,l+j,l+j+1,elastic,p); //0=alfa
     KE_GE_MOMENTUM(y,&edo_params,l+j,kinetic,gravity,momentum);
     p++;
   }
       KE_GE_MOMENTUM(y,&edo_params,l+j,kinetic,gravity,momentum);
       //printf("----------------------------------------------\n");

       p+=5*_grain_layer;	      l=k*(2*_grain_layer+1)+1;
     }//exit(1);
   //printf("***************************************************\n");


   //ENTRE OS GRÃOS PRINCIPAIS EM TODAS AS CADEIAS
   //printf("graos principais entre DUAS camadas\n");
   i=3*_grain_layer;

   for (k = 1; k < _num_layer; k++)	//camada dos grãos principais
     {
       p=(k-1)*(2*_grain_layer+1)+1;
       l=k*(2*_grain_layer+1)+1;
       for (j = 0; j < _grain_layer; j++)
   {
     energy_between_grains(y,&edo_params,p,l,elastic,i);
     i++;p++;l++;
   }//printf("--------------------------------------------\n");
       i+=5*_grain_layer-1;
     }//exit(1);
   //printf("***************************************************\n");


   //ENTRE OS GRÃOS PRINCIPAIS E DECORADORES DA FRENTE
   //printf("\nGRAOS PRINCIPAIS E DECORADORES DA FRENTE\n");

   i=_grain_layer;   n = _grain_layer + 1;//número do grão decorador
   for (k = 1; k < _num_layer; k++)	//camada dos grãos principais
     {
       p=(k-1)*(2*_grain_layer+1)+1;
       KE_GE_MOMENTUM(y,&edo_params,n,kinetic,gravity,momentum);

       //printf("KE%d=%g\n",n,kinetic[n]);

       for (j = 0; j < _grain_layer; j++)
   {
     for (l = 1; l < 3; l++)	//camada dos grãos principais
       { energy_between_grains(y,&edo_params,p,n,elastic,i);
         i++; n++;
       }//printf("\n");
     p++; n--;
     KE_GE_MOMENTUM(y,&edo_params,n,kinetic,gravity,momentum);

     //printf("KE%d=%g\n",n,kinetic[n]);

   }//printf("------------------------------------------------\n");
       n+=_grain_layer+1;	      i+=4*_grain_layer-1;
     }//exit(1);
   //printf("***********************************************\n");//exit(1);


   //ENTRE OS GRÃOS DECORADORES E OS PRINCIPAIS DA FRENTE
   //printf("\nGRAOS DECORADORES E OS PRINCIPAIS DA FRENTE\n");

   i=4*_grain_layer;   n = _grain_layer + 1;//número do grão decorador

   for (k = 1; k < _num_layer; k++)	//camada dos grãos principais
     {
       l=k*(2*_grain_layer+1)+1;
       for (j = 0; j < _grain_layer; j++)
   {
     for (p = 1; p < 3; p++)	//camada dos grãos principais
       {
         energy_between_grains(y,&edo_params,n,l,elastic,i);
         i++; n++;
       }//printf("\n");
     n--; l++;
   }//printf("------------------------------------------------\n");
       n+=_grain_layer+1;  i+=4*_grain_layer-1;
     }//printf("***********************************************\n");//exit(1);

   ////////////////////////////////////////////////////////////////////
   ////////////// ENTRE OS GRÃOS INTERTICIAIS /////////////////////////
   ////////////////////////////////////////////////////////////////////
   //printf("\nGRAOS INTERTICIAIS\n");

   n = _grain_layer + 1;//primeiro grão intersticial
   i = 0;
   //  printf("grao principal e o grao decorador da frente\n");

   for (k = 1; k < _num_layer-1; k++)	//camada dos grãos principais
     {
       for (j = 0; j < _grain_layer; j++)
   {
     energy_inter(y,&edo_params,n,n+1,elastic,i+j);
     energy_inter(y,&edo_params,n,n+2*_grain_layer+1,
            elastic,i+j+_grain_layer);
     n++;
   }//printf("----------------------------------------\n");
       energy_inter(y,&edo_params,n,n+2*_grain_layer+1,
        elastic,i+j+_grain_layer);
       n += _grain_layer + 1;
       i += 2*_grain_layer + 1;
     }

   for (j = 0; j < _grain_layer; j++)
     {energy_inter(y,&edo_params,n,n+1,elastic,i+j);n++;}
         //exit(1);
   ////////////////////////////////////////////////////////////////////


   //PARA A PAREDE DE CIMA
   //printf("Paredes de CIMA\n");
   energyDOWLOADING_WALL (y,&edo_params,1,elastic,0);

   i=2;
   for(k=1;k<_num_layer;k++)
     {
       p=(k-1)*(2*_grain_layer+1)+1;

       energyDOWLOADING_WALL (y,&edo_params,p+_grain_layer,elastic,i);
       //printf("grain=%d\n",p+_grain_layer);

       i+=2;
     }//printf("\n***********************************************\n");//exit(1);

   //PARA A PAREDE DE BAIXO
   //printf("\nParedes de BAIXO\n");
   energy_wall_bottom (y,&edo_params,_grain_layer,elastic,1);

   i=3;
   for(k=1;k<_num_layer;k++)
     {
       p=(k-1)*(2*_grain_layer+1)+1;
       energy_wall_bottom (y,&edo_params,p+(2*_grain_layer),elastic,i);
       i+=2;
     }//exit(1);

   //          printf("x_0=%g\t",x[0]);
   for(n=0;n<_grain_layer;n++)
     energy_temporary_wall(x,&edo_params,n,N-1-_grain_layer-n,
         N-1-n,elastic,2*_grain_layer-k,(3*_grain_layer)-n);

   if ((position_dissipated=fopen(name7,"w"))==NULL) {
     printf("erro no print da 'POSITION_DISSIPATED'!!!\n");exit(1);}

   //CÁLCULO DE ENERGIA TOTAL
   KINETIC=ELASTIC=MOMENTUM_X=MOMENTUM_Y=GRAVITY=0.;
   fprintf(position_dissipated,"%g\n",x[0]);
   for(n=1;n<N;n++)
     {
       fprintf(position_dissipated,"%g\n",x[n]);
       ELASTIC+=elastic[n];
       KINETIC+=kinetic[n];
       GRAVITY+=gravity[n];
       MOMENTUM_X+=momentum[n];
       MOMENTUM_Y+=momentum[N+n];
     }//exit(1);

   if(ciclo==0){E_INITIAL=KINETIC+ELASTIC+GRAVITY;}
   TOTAL = KINETIC+ELASTIC+GRAVITY;
   relative_erro=(E_INITIAL-TOTAL)/E_INITIAL;


   printf("D_e=%g\tE=%g\terro=%g\tKE=%g\tEE=%g\tGE=%g\tt=%g\tw=%g\twall=%g\n",
    energy_before-TOTAL,TOTAL,relative_erro,KINETIC,
    ELASTIC,GRAVITY,t,RESTITUTION,DOWLOADING_WALL);

   energy=fopen(name8,"a");
   fprintf(energy,"%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",energy_before-TOTAL,TOTAL,relative_erro,KINETIC,ELASTIC,GRAVITY,t,RESTITUTION,DOWLOADING_WALL);
         fclose(energy);

   if((energy_before-TOTAL)<tolerance) RESTITUTION-=0.033;

         for(n=N;n<2*N;n++)fprintf(position_dissipated,"%g\n",x[n]);
         fprintf(position_dissipated,"#t=%g\tEi-Ef=%g\tEf=%g",
     t,E_INITIAL-TOTAL,TOTAL);fclose(position_dissipated);
   energy_before = TOTAL;//exit(1);

   /********************** SEQUÊNCIA PARA O GIF ******************/
         sequence2 = fopen(name4,"a");
   fprintf(sequence2,
     "plot \"./%d_random_square_q%g_%dx%d.dat\"using 1:2:3:3 with circles lc palette lw 2 lt 1\n",num-1,_q,_num_layer,_grain_layer);
   fprintf(sequence2,
     "set arrow from %g,%g to %g,%g nohead\n\n",
     1.,(b[1]+radii[1])/radii[0],1.,(b[_grain_layer]-radii[_grain_layer])/radii[0]);
   fprintf(sequence2,
     "set arrow from %g,%g to %g,%g nohead\n\n",
     1.,(b[1]+radii[1])/radii[0],INITIAL_WALL/radii[0],(b[1]+radii[1])/radii[0]);
   fprintf(sequence2,
     "set arrow from %g,%g to %g,%g nohead\n\n",
     1.,(b[_grain_layer]-radii[_grain_layer])/radii[0],INITIAL_WALL/radii[0],
     (b[_grain_layer]-radii[_grain_layer])/radii[0]);
   fprintf(sequence2,
     "set arrow from %g,%g to %g,%g nohead\n\n",
     INITIAL_WALL/radii[0],(b[_grain_layer]-radii[_grain_layer])/radii[0],
     INITIAL_WALL/radii[0],(b[1]+radii[1])/radii[0]);
   fprintf(sequence2,
     "set arrow from %g,%g to %g,%g nohead lc 1\n\n",
     DOWLOADING_WALL/radii[0],(b[_grain_layer]-radii[_grain_layer])/radii[0],
     DOWLOADING_WALL/radii[0],(b[1]+radii[1])/radii[0]);
   fprintf(sequence2,"pause 1\n\n");
         fclose(sequence2);
   /**************************************************************/

         if(RESTITUTION<=0.01)break;
 }

     int status
 = gsl_odeiv_step_apply (s, t, h, x, x_err, dydt_in, dydt_out, &sys);
     // passo do rk

     if (status != GSL_SUCCESS)// verifica se o passo foi bem sucedido
 {printf ("problema\n");	  break;}
     // troca dydt_in e dydt_out
     swap = dydt_in;dydt_in = dydt_out;dydt_out = swap;

     t += h; ciclo++;
   }


 ////////////////////////////////////////////////////////////////////
 /////////////////// GERANDO VIDEO PELO GNUPLOT /////////////////////
 ////////////////////////////////////////////////////////////////////

 FILE *video;       char name3[100];
//  FILE *sequence2;       char name4[100];

 sprintf(name3, "evolution_desorder_q%g_%dx%d.gp",_q,_num_layer,_grain_layer);
 video = fopen(name3,"w+");

 fprintf(video,
 "set palette defined (0.25 'black', 0.5 'green',1 'blue')\n");
 //Acima esta a definição da escala de cores
 fprintf(video,"set cbrange [0.2:1]\n"); //tamanho da escala de cores
 fprintf(video,"set style fill transparent solid 0.9\n");
 //Acima: 90% das cores são preenchidas e permite transparência

 fprintf(video,"unset key\n");

 fprintf(video,"set yrange[-12:12]\n");
 fprintf(video,"set xrange[-2:27]\n\n");
 fprintf(video,"set ylabel \"y\"\n");
 fprintf(video,"set xlabel \"x\"\n");
 fprintf(video,"set cblabel \"normalized radii\"\n\n");

 fprintf(video,"load \'sequencia_%dx%d.plot\'\n\n",_num_layer,_grain_layer);

/*  sprintf(name4, "sequencia.plot");
 sequence2 = fopen(name4,"w+");
 for(n=n_inicial;n<num;n++)
   {
     fprintf(sequence2,
       "plot \"./%d_random_square_q%g_%dx%d.dat\"using 1:2:3:3 with circles lc palette lw 2 lt 1\n",n,_q,_num_layer,_grain_layer);

     fprintf(sequence2,
       "set arrow from %g,%g to %g,%g nohead\n\n",
    1.,(b[1]+radii[1])/radii[0],1.,(b[_grain_layer]-radii[_grain_layer])/radii[0]);

     fprintf(sequence2,
       "set arrow from %g,%g to %g,%g nohead\n\n",
    1.,(b[1]+radii[1])/radii[0],INITIAL_WALL/radii[0],(b[1]+radii[1])/radii[0]);

     fprintf(sequence2,
       "set arrow from %g,%g to %g,%g nohead\n\n",
    1.,(b[_grain_layer]-radii[_grain_layer])/radii[0],INITIAL_WALL/radii[0],
    (b[_grain_layer]-radii[_grain_layer])/radii[0]);

     fprintf(sequence2,
       "set arrow from %g,%g to %g,%g nohead\n\n",
    INITIAL_WALL/radii[0],(b[_grain_layer]-radii[_grain_layer])/radii[0],
    INITIAL_WALL/radii[0],(b[1]+radii[1])/radii[0]);

     fprintf(sequence2,
       "set arrow from %g,%g to %g,%g nohead lc 1\n\n",
    DOWLOADING_WALL/radii[0],(b[_grain_layer]-radii[_grain_layer])/radii[0],
    DOWLOADING_WALL/radii[0],(b[1]+radii[1])/radii[0]);

     fprintf(sequence2,"pause 1\n\n");
   }*/

 ////////////////////////////////////////////////////////////////////
 /////////////////// GERANDO 'GIF' PELO GNUPLOT /////////////////////
 ////////////////////////////////////////////////////////////////////
 FILE *gif;       char name5[100];

 sprintf(name5, "animate_q%g_%dx%d.gp",_q,_num_layer,_grain_layer);
 gif = fopen(name5,"w+");

 fprintf(gif,"reset \n set term gif animate\n");
 fprintf(gif,"set output \"animate_q%g_%dx%d.gif\"\n",_q,_num_layer,_grain_layer);
 fprintf(gif,
 "set palette defined (0.25 'black', 0.5 'green',1 'blue')\n");
 //Acima esta a definição da escala de cores
 fprintf(gif,"set cbrange [0.2:1]\n"); //tamanho da escala de cores
 fprintf(gif,"set style fill transparent solid 0.9\n");
 //Acima: 90% das cores são preenchidas e permite transparência

 fprintf(gif,"unset key\n");

 fprintf(gif,"set yrange[-12:12]\n");
 fprintf(gif,"set xrange[-2:27]\n\n");
 fprintf(gif,"set ylabel \"y\"\n");
 fprintf(gif,"set xlabel \"x\"\n");
 fprintf(gif,"set cblabel \"normalized radii\"\n\n");

 fprintf(gif,"load \'sequencia.plot\'\n\n");

//  gsl_odeiv_step_free (s);
 free (mass);
 free (alpha);
 free (beta);
 free (temporary_wall);
 free (f);
 free (x_err);
 free (dydt_in);
 free (dydt_out);

 free (radii);
 free (x);
 free (a);
 free (b);
 
 return 0;
}
