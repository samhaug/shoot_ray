#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "ray.h"
#define ARR_SIZE 60000
#define M_PI 3.141592654
#define deg2rad M_PI/180.
#define rad2deg 180./M_PI

/*
   Written by: samhaug 
   2021 Jan 18 01:54:07 PM
*/

int jdebug = 0;

int shoot(struct Ray *ray, double new_rad, double new_vel);
int prepare_model(FILE *fp, double rad[], double vp[], double vs[], int *num);

int vincenty_direct(double,double,double,double,double *,double *,double *);
int shoot_ray(double rad[], double vp[], double vs[],  double angle, double start_rad, 
              double end_rad, double lat_1, double lon_1,
              double bearing, int *ray_count, double lat_list[], double lon_list[], 
              double dep_list[], double time_list[]);

int shoot_ray(double rad[], double vp[], double vs[], double angle, double start_rad, 
              double end_rad, double lat_1, double lon_1,
              double bearing, int *ray_count, double lat_list[], double lon_list[], 
              double dep_list[], double time_list[]){


   //initialize ray for P velocity
   struct Ray p_ray;
   //initialize ray for S velocity
   //struct Ray s_ray;
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
      // end program if ray reaches CMB
      if (p_ray.radius < 3400.) return 0;
      vincenty_direct(lat_1, lon_1 , bearing, 111195.*rad2deg*p_ray.dist, 
                 &lat_2, &lon_2, &alpha21);
      lat_list[*ray_count] = lat_2;
      lon_list[*ray_count] = lon_2;
      dep_list[*ray_count] = p_ray.radius;
      time_list[*ray_count] = p_ray.time;
      *ray_count = *ray_count+1;
      i++;

   }

   if (jdebug) printf("Ray turns\n");

   p_ray.angle = 180.-p_ray.angle;
   i=i-1;

   // Upward branch
   while(p_ray.radius < end_rad){
      if (!shoot(&p_ray, rad[i], 1./vp[i])){
         fprintf(stderr,"Problem with ray shooting\n");
         exit(1);
      }
      vincenty_direct(lat_1, lon_1 , bearing, 111195.*rad2deg*p_ray.dist, 
                &lat_2, &lon_2, &alpha21);
      lat_list[*ray_count] = lat_2;
      lon_list[*ray_count] = lon_2;
      dep_list[*ray_count] = p_ray.radius;
      *ray_count = *ray_count+1;
      i=i-1;
   }

return 0;

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
    //if (jdebug) fprintf(stderr,"%3d %lf %lf %lf \n",*num,rad[*num],vp[*num],vs[*num]);
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


