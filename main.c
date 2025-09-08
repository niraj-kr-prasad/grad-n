 #include <complex.h>

// Total viscoelastic viscosity 
#define MU0 1.  
#define MUR 1.        
// Ratio of the solvent viscosity to the total viscosity
#define BETA (0.11)
// Polymer viscosity
#define MUP ((1. - BETA)*MU0)
// Solvent viscosity 
#define MUS (BETA*MU0)
// Relaxation viscoelastic time
#define LAM 1.
// Average velocity steady flow
#include "tecplot.h"
#include "navier-stokes/centered.h"
#include "log-conform.h"
#include "fene-p.h"
#include "two-phase.h"
#include "tension.h"
#include "tag.h"
//#include "axi.h"


#define L2 5.
#define M 1.0  
#define Ca 0.1    // Capillary number
#define Deb 1   // Deborah number
#define Re 2. // Reynold number
#define We (Ca*Re) 
#define l1_dim 40.0 //dimensional length of the entry region 
//#define l2_dim 34.0 //dimensional length of the entry region - tapering
#define l3_dim 40.0  //dimensional lenth of the central region
//#define l4_dim 34.0 //dimensional lenth of the exit region - tapering
#define l5_dim 40.0 //dimensional length of the exit region 
#define D_dim 6.0 //dimensional dia of the HeLa cell
#define D1_dim 24.0  //dimensional dia of the entry region
#define D2_dim 4.0  //dimensional dia of the exit region
#define L_dim (l1_dim+l3_dim+l5_dim)  //total length of the domain in the x direction
//all lengths are defined in micrometers
//now we non dimensionalise them
//#define L 1.0  //total length of the domain in the x direction
#define l1 (l1_dim*1.0/L_dim) //non dimensional length of the entry region 
//#define l2 (l2_dim/L_dim) //non dimensional length of the entry region - tapering
#define l3 (l3_dim/L_dim)  //non dimensional lenth of the central region
//#define l4 (l4_dim/L_dim) //non dimensional lenth of the exit region - tapering
#define l5 (l5_dim/L_dim) //non dimensional length of the exit region 
#define D  (D_dim/L_dim)  
#define D1 (D1_dim/L_dim ) //non dimensional dia of the entry region
#define D2 (D2_dim/L_dim)  //non dimensional dia of the exit region
//non dimensional dia of the HeLa cell i.e. this is taken as the characteristic length scale


//#define m1 ((D1-D2)/(2.0*l2)) //slope of entry tapering region(- for top)
//#define m2 ((D1-D2)/(2.0*l4)) //slope of exit tapering region(+ for top)


#define DOMAINSIZE 1


int lev;
//int MAXLEVEL = 12;
scalar lam[], mupc[], SDXX[], SYY[], SXY[], TR[];


int main()
{
 size (DOMAINSIZE);
   origin (0.,-0.5);
//p[left] = dirichlet(1.0);
   const scalar lam[] = LAM;
  lambda = lam;
  const scalar mupc[] = MUP;
  mup = mupc;
  const face vector mus[] = {MUS,MUS};
  mu = mus;
  mu1 = MUR*BETA/Re;
  mu2 = 1./Re;
  rho1 = M;
  rho2 = 1.;
  f.sigma = 1./We;
  p[right] = dirichlet(0.);
   u.n[left] = dirichlet(1.);
  //u.n[right] = neumann(0.0);
  u.t[left] = dirichlet(0.0);
   //u.t[right] = dirichlet(0.0);
  //f[bottom] = 0.0;
  //f[top] = 0.0;
  u.t[top] = dirichlet(0.0);
  u.t[bottom] = dirichlet(0.0);
  //u.n[top]= dirichlet(0.0);
  //u.n[bottom] = dirichlet(0.0);
    //p[right] = dirichlet(0.0);
    SDXX[left] = neumann(0.);
    SDXX[top] = neumann(0.);
    SDXX[bottom] = neumann(0.);
      SXY[left] = neumann(0.);
    SXY[top] = neumann(0.);
    SXY[bottom] = neumann(0.);
      SYY[left] = neumann(0.);
    SYY[top] = neumann(0.);
    SYY[bottom] = neumann(0.);
  //SD[right] = neumann(0.);
  lev = 10;
  init_grid (1024);
  run();
}
   
event init (i = 0) 

{  
 //refine(y < D1/2. && y > -D1/2.  && level<9);
  //refine((sq(D/2.0) - (sq(x- l1/2.0)+ sq(y))) && level<9);
   mask (y > D1/2 && x<l1 ? top : none);
   mask(y < -D1/2 && x<l1 ? bottom : none);
   mask (y>D2/2 && x>l1 && x<=l1+l3 ? top : none);
   mask (y<-D2/2 && x>=l1 && x<=l1+l3 ? bottom : none);
   mask(y>D1/2 && x>l1+l3 ? top : none);
  mask(y<-D1/2 && x>l1+ l3 ? bottom:none);
   fraction(f, sq(D/2.0) - (sq(x- l1/2.0)+ sq(y)));

  foreach()
   u.x[] = 0;
 boundary ({f});
}
event TAU(t += 0.001)
{
    FILE * fp;
foreach () 
  {
    SDXX[] = (tau_p.x.x[])*clamp(f[],0,1);
    SYY[] = (tau_p.y.y[])*clamp(f[],0,1);
     SXY[] = (tau_p.y.y[])*clamp(f[],0,1);
    TR[] = (trA[])*clamp(f[],0,1);
  }
  boundary ({SDXX,SYY,SXY,TR});
     if (i== 0)
  fp= fopen("Tau.dat","w");
else
  fp= fopen("Tau.dat","a");

  fprintf (fp, "%g %g %g %g \n", SDXX, SYY, SXY, TR);
  fflush(fp);
  fclose(fp);

}
event properties (i++) 
{
  foreach () 
  {
    lam[] = Deb*clamp(f[],0,1);
    mupc[] = MUR*MUP*clamp(f[],0,1)/Re;
  }
  boundary ({lam,mupc});

}


event deformation (t += 0.001; t< 1)
{
  FILE * fp;

  double maxx = -HUGE, maxy = -HUGE, minx = HUGE, miny = HUGE ;

  foreach (reduction(max:maxx) reduction(min:minx) reduction(max:maxy) reduction(min:miny)) 
  {
    if (f[] > 1e-6 && f[] < 1.0 - 1e-6 )
      {
        coord p;
        coord n = mycs (point, f);
        double alp = plane_alpha (f[], n);
        plane_area_center (n, alp, &p);
        double xval  = x + Delta*p.x;
        double yval  = y + Delta*p.y;
        
        if (xval > maxx) maxx = xval;
        if (yval > maxy) maxy = yval;
        if (xval < minx) minx = xval;
        if (yval < miny) miny = yval;
      }
  }
  
  double Le = maxx - minx;
  double B = maxy - miny;
  double De = (Le - B) / (Le + B);
  double zd = (((maxx - l1/2.) - Le/2)/2.);
   if (i== 0)
  fp= fopen("Deformation.dat","w");
else
  fp= fopen("Deformation.dat","a");

  fprintf (fp, "%g %g %g %g %g %g \n", t, dt, zd, De);
  fflush(fp);
  fclose(fp);

  /*fprintf (fp, "%g %g %g %g \n", t, dt, zd, De);
 fflush (fout);
 fprintf (stderr, "%g %g %g\n", t, zd, De)*/

}
event droplets (t += 0.001)
{
  FILE *fp;
  scalar m[];
  foreach()
    m[] = f[] > 1e-3;
  int n = tag (m);
  double v[n];
  coord b[n];
  double vel[n];
  for (int j = 0; j < n; j++)
    v[j] = b[j].x = b[j].y = b[j].z = vel[j]=0.;
  foreach_leaf()
    if (m[] > 0) {
      int j = m[] - 1;
      v[j] += dv()*f[];
      vel[j] += dv()*f[]*u.x[];
      coord p = {x,y,z};
      foreach_dimension()
  b[j].x += dv()*f[]*p.x;
    }

  if (i == 0)
    fp = fopen("Position.dat","w");
  else
    fp = fopen("Position.dat","a");

  for (int j = 0; j < n; j++)
    fprintf (fp, "%g %g %g %g \n", t, v[j], b[j].x/v[j], vel[j]/v[j]);
  fflush (fp);
  fclose (fp);

  for (int j = 0; j < n; j++)
    fprintf (ferr, "%g %g %g %g %g\n", t, dt, v[j], b[j].x/v[j], vel[j]/v[j]);

  fflush (fout);
}

event adapt (i++) {
  adapt_wavelet ({f, u}, (double[]){1e-5, 1e-5, 1e-5}, lev, lev - 1);
}
event tecplot_output(t += 0.01;t< 0.5)
{
  struct OutputTec tec;
  scalar * list_Q = {u,p,f, SDXX, SYY, SXY, TR};
  tec.tec_cc = list_copy(list_Q);
  sprintf(tec.varname,"X Y U V P VoF SDXX SDYY SXY TR");
  output_tec(tec,i,t, 1.);
  free(tec.tec_cc);
}

/*{
  struct OutputTec tec;
  scalar * list_Q = {u,p,f};
  tec.tec_cc = list_copy(list_Q);
  sprintf(tec.varname,"X Y U V P VoF");
  output_tec(tec,i,t, 1.);
  free(tec.tec_cc);
}*/

event movie (t += 0.02)
{
  char name[100];
  sprintf (name, "IMG_%.2f.ppm", t);

  FILE * fp = fopen (name, "w");

  output_ppm (f, fp,1024, linear = true);
  fclose (fp); 
}
