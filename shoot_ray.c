#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "ray.h"
#define ARR_SIZE 60000
#define deg2rad M_PI/180.
#define rad2deg 180./M_PI

/*
   Written by: samhaug 
   2021 Jan 18 01:54:07 PM
*/

int jdebug = 0;

int usage();
int shoot(struct Ray *ray, double new_rad, double new_vel);
int prepare_model(FILE *fp, double rad[], double vp[], double vs[], int *num);

int main(int argc, char **argv){
   FILE *fp; 
   fp = fopen("smoothprem_nodensity_flip.txt","r");
   //fp = fopen("interp_prem.txt","r");
   double rad[ARR_SIZE];
   double vp[ARR_SIZE];
   double vs[ARR_SIZE];
   //initialize ray for P velocity
   struct Ray p_ray;
   //initialize ray for S velocity
   //struct Ray s_ray;
   int num=0;

   if (!prepare_model(fp, rad, vp, vs, &num)){
      fprintf(stderr,"Problem with prepare_model\n");
      exit(1);
   } 
   
   if ( argc != 7 ){
      usage();
      exit(1);
   }

   // Starting incidence angle
   double angle = atof(argv[1]);
   // Starting radius
   double start_rad = atof(argv[2]);
   // Ending radius
   double end_rad = atof(argv[3]);
   // Starting latitude
   double lat_1 = atof(argv[4]);
   // Starting longitude
   double lon_1 = atof(argv[5]);
   // bearing CW from north
   double bearing = atof(argv[6]);

   double lat_2, lon_2, alpha21;

   //Find index of start_rad
   int i=0;
   double min = 100;

   while ( min > 2 ){
     min = fabs(start_rad - rad[i]);
     i++;
   }
   if (jdebug) printf("Closest rad: %5d %8.4lf\n", i, rad[i]);

   //initialize ray
   p_ray.radius = rad[i];
   p_ray.angle = angle;
   p_ray.angle_rad = deg2rad*p_ray.angle;
   p_ray.vel = 1./vp[i];
   p_ray.p_sph = p_ray.radius * p_ray.vel * sin(p_ray.angle_rad);
   p_ray.time = 0;
   p_ray.dist = 0;

   while(p_ray.angle < 88.5){
      if (!shoot(&p_ray, rad[i], 1./vp[i])){
         fprintf(stderr,"Problem with ray shooting\n");
         exit(1);
      }
      //vincenty_direct_sphere(lat_1, lon_1 , bearing, 111195.*rad2deg*p_ray.dist, 
      //           1000*.p_ray.radius, &lat_2, &lon_2, &alpha21);
      printf("%8.2lf %8.2lf %8.2lf\n", rad2deg*p_ray.dist, p_ray.radius, p_ray.time);
      //printf("%8.2lf %8.2lf %8.2lf\n", rad2deg*p_ray.dist, p_ray.radius, p_ray.time);
      //lat_1 = lat_2;
      //lon_1 = lon_2;
      i++;
   }

   p_ray.angle = 180.-p_ray.angle;
   i=i-1;

   // Upward branch
   while(p_ray.radius < end_rad){
      if (!shoot(&p_ray, rad[i], 1./vp[i])){
         fprintf(stderr,"Problem with ray shooting\n");
         exit(1);
      }
      //vincenty_direct_sphere(lat_1, lon_1 , bearing, 111195.*rad2deg*p_ray.dist, 
      //         1000.*p_ray.radius, &lat_2, &lon_2, &alpha21);
      //printf("%8.2lf %8.2lf %8.2lf %8.2lf\n",lat_2, lon_2, p_ray.radius, p_ray.time);
      //lat_1 = lat_2;
      //lon_1 = lon_2;
      printf("%8.2lf %8.2lf %8.2lf\n", rad2deg*p_ray.dist, p_ray.radius, p_ray.time);
      i=i-1;
   }

}

int prepare_model(FILE *fp, double rad[], double vp[], double vs[], int *num){
  *num = 0;
  // Read in three column Earth model. 
  // Radius (km)  Vp (km/s) Vs (km/s)
  // Largest radius should be first line. Center of Earth should be last.
  while (!feof(fp)){
    if (!fscanf(fp,"%lf %lf %lf",rad+*num, vp+*num, vs+*num)){
      fprintf(stderr,"Error reading input model file\n");
    }
    if (jdebug) fprintf(stderr,"%3d %lf %lf %lf \n",*num,rad[*num],vp[*num],vs[*num]);
    (*num)+=1;
  }
  return 1;
}

int shoot(struct Ray *ray, double new_rad, double new_vel){
    // Update ray location, angle, and total time.
    double d_r;
    double time_inc,dist_inc;
    d_r = fabs(ray->radius - new_rad);
    //Discrete addition to time (Shearer Eqn. 4.43)
    time_inc = ( d_r * pow(new_vel*new_rad,2) / 
                 ( sqrt(pow(new_vel*new_rad,2) - pow(ray->p_sph,2)) * new_rad ));
    ray->time += time_inc;
    // Discrete addition to distance (Shearer Eqn. 4.43)
    dist_inc = ( d_r * ray->p_sph /
                 ( sqrt( pow(new_vel*new_rad,2) - pow(ray->p_sph,2)) * new_rad ));
    ray->dist += dist_inc;
    ray->radius = new_rad;
    ray->vel = new_vel;
    ray->angle_rad = asin(ray->p_sph / (ray->radius * ray->vel) );
    ray->angle = rad2deg*(ray->angle_rad);
    return 1;
}


int usage(){
    fprintf(stderr," \n");
    fprintf(stderr,"USAGE : shoot_ray takeoff_angle r1 r2 lat lon az\n");
    fprintf(stderr," \n");
    fprintf(stderr,"takeoff_angle: Takeoff angle (degrees) from vertical \n");
    fprintf(stderr,"lat/lon: starting coordiantes \n");
    fprintf(stderr,"az: Azimuth (cw from north) \n");
    fprintf(stderr,"r1: Starting radius (km) \n");
    fprintf(stderr,"r2: Ending radius (km) \n");
    fprintf(stderr," \n");
    exit(1);
}





