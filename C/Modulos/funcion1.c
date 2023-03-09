#include "allvars.h"

extern Particulas *part;
extern io_encabezado encabezado;
extern double G;
extern double etha;
extern double eps;
extern int nTotal;

// evolve the system using a leap-frog integrator
int evolve(double totalTime, double dt)
{

  int i, counter=0, counter2=0;
  double t =0.0;
  char nombreArchivo[1000];

  sprintf(nombreArchivo,"./output/sim_eps_galaxia_%.4d",counter);
  imprimirSnapshot(nombreArchivo, t);
  
  // calcula aceleracion sobre todas las particulas
  acceleration();

  while( t<totalTime )
    {

      tiempoAdactativo(&dt);

      printf("tiempo = %.10lf  dt = %.10lf\n",t,dt);

      // drift
      for( i=0; i<nTotal; i++ )
	{
	  part[i].pos[X] = part[i].pos[X] + 0.5*dt*part[i].vel[X];
	  part[i].pos[Y] = part[i].pos[Y] + 0.5*dt*part[i].vel[Y];
	  part[i].pos[Z] = part[i].pos[Z] + 0.5*dt*part[i].vel[Z];
	}
      
      // calcula aceleracion sobre todas las particulas
      acceleration();
      
      // kick
      for( i=0; i<nTotal; i++ )
	{
	  part[i].vel[X] = part[i].vel[X] + dt*part[i].accel[X];
	  part[i].vel[Y] = part[i].vel[Y] + dt*part[i].accel[Y];
	  part[i].vel[Z] = part[i].vel[Z] + dt*part[i].accel[Z];
	}
      
      // drift
      for( i=0; i<nTotal; i++ )
	{
	  part[i].pos[X] = part[i].pos[X] + 0.5*dt*part[i].vel[X];
	  part[i].pos[Y] = part[i].pos[Y] + 0.5*dt*part[i].vel[Y];
	  part[i].pos[Z] = part[i].pos[Z] + 0.5*dt*part[i].vel[Z];
	}
      
      t += dt;
      counter++;
      
      //if( counter%100 == 0 )
      //{
      //  counter2++;
      // sprintf(nombreArchivo,"sim_eps_galaxia_%.3d",counter2);
	  sprintf(nombreArchivo,"./output/sim_eps_galaxia_%.4d",counter);
	  imprimirSnapshot(nombreArchivo, t);
	  //}

    }

  return 0;
}

int acceleration(void)
{

  int i, j;
  double dr[3]; 
  double r, r3;
  double eps2;

  eps2 = eps*eps;

  for( i=0; i<nTotal; i++ )
    {

      part[i].accel[X] = 0.0;
      part[i].accel[Y] = 0.0;
      part[i].accel[Z] = 0.0;

      for( j=0; j<nTotal; j++ )
      {
	if( j != i )
	  {

	    dr[X] = part[i].pos[X] - part[j].pos[X];
	    dr[Y] = part[i].pos[Y] - part[j].pos[Y];
	    dr[Z] = part[i].pos[Z] - part[j].pos[Z];

	    r = sqrt( dr[X]*dr[X] + dr[Y]*dr[Y] + dr[Z]*dr[Z] + eps2 );
	    r3 = r*r*r;

	    part[i].accel[X] = part[i].accel[X] - G*part[j].masa*( dr[X]/r3 );
	    part[i].accel[Y] = part[i].accel[Y] - G*part[j].masa*( dr[Y]/r3 );
	    part[i].accel[Z] = part[i].accel[Z] - G*part[j].masa*( dr[Z]/r3 );

	  }
      }

      part[i].accelMag = sqrt( part[i].accel[X]*part[i].accel[X] +
			       part[i].accel[Y]*part[i].accel[Y]+
			       part[i].accel[Z]*part[i].accel[Z] );

    }

  return 0;
}

int imprimirSnapshot(char *nombreArchivo, double t)
{
  int i, dummy;
  FILE *fSnaps=fopen(nombreArchivo,"w");

  encabezado.time = t;

  dummy = 256;
  
  // imprimo encabezado
  fwrite(&dummy, sizeof(int), 1, fSnaps);
  fwrite(&encabezado, sizeof(io_encabezado), 1, fSnaps);
  fwrite(&dummy, sizeof(int), 1, fSnaps);

  // imprimo posiciones
  dummy = nTotal*3*sizeof(float);
  fwrite(&dummy, sizeof(int), 1, fSnaps);
  for( i=0; i<nTotal; i++ )
    fwrite(&part[i].pos[0], sizeof(float), 3, fSnaps);
  fwrite(&dummy, sizeof(int), 1, fSnaps);

  // imprimo velocidades
  dummy = nTotal*3*sizeof(float);
  fwrite(&dummy, sizeof(int), 1, fSnaps);
  for( i=0; i<nTotal; i++ )
    fwrite(&part[i].vel[0], sizeof(float), 3, fSnaps);
  fwrite(&dummy, sizeof(int), 1, fSnaps);
  
  // imprimo ids
  dummy = nTotal*sizeof(int);
  fwrite(&dummy, sizeof(int), 1, fSnaps);
  for( i=0; i<nTotal; i++ )
    fwrite(&part[i].id, sizeof(int), 1, fSnaps);
  fwrite(&dummy, sizeof(int), 1, fSnaps);

  // imprimo masas
  dummy = nTotal*sizeof(float);
  fwrite(&dummy, sizeof(int), 1, fSnaps);
  for( i=0; i<nTotal; i++ )
    fwrite(&part[i].masa, sizeof(float), 1, fSnaps);
  fwrite(&dummy, sizeof(int), 1, fSnaps);
  
  fclose(fSnaps);
  
  return 0;
}

int tiempoAdactativo(double *dt)
{

  int i;
  double dtAux, accel;

  *dt = etha/sqrt( part[0].accelMag );

  for( i=1; i<nTotal; i++ )
    {

      dtAux =  etha/sqrt( part[i].accelMag );
     
      if( dtAux<*dt )
	{
	  *dt = dtAux;
	  accel = part[i].accelMag;
	}

    } 

  printf("accel = %.10lf  dt = %.10lf\n",accel,*dt);

  return 0;
}


