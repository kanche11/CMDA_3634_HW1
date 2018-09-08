#include <math.h>
#include "intersections.h"
#include <stdio.h>

//if there is no intersection, this function should return 0
//Otherwise, populate 'intersection' with the point of intersection and return 1
int ray_sphere_intersection(ray_t observer, sphere_t obj, vector_t *intersection) {
  
  //linear term inside the norm
  vector_t h_term = observer.dir;
  //constant term inside the norm
  vector_t const_term= difference(observer.start,obj.center);

  double a = dot_product(h_term,h_term);
  
  double b = 2*dot_product(h_term,const_term);
  
  double c = dot_product(const_term,const_term) - obj.radius;

  double discriminant = b*b - 4*a*c;

  //make sure the equation has real roots
  if (discriminant < 0.) return 0;
  //find the smallest positive h
  double h;
  if ((-1*b - sqrt(discriminant))/(2*a) > 0) h = (-1*b - sqrt(discriminant))/(2*a); //sphere in front of us
  else if ((-1*b + sqrt(discriminant))/(2*a) > 0) h = (-1*b + sqrt(discriminant))/(2*a); //we are inside the sphere
  else return 0; //no positive roots, so the sphere is behind us

  vector_t solution = scaled_sum(1.,observer.start,h,observer.dir);

  copy_vector(solution,intersection);
  return 1;  
}

//if there is no intersection, this function should return 0
//Otherwise, populate 'intersection' with the point of intersection and return 1
int ray_disk_intersection(ray_t observer, disk_t obj, vector_t *intersection) {
  //Question 3: Modify this function to compute an intersection
  double t = dot_product(difference(obj.center,observer.start),obj.normal)/dot_product(observer.dir,obj.normal);
  vector_t new = scaled_sum(1.,observer.start,t,observer.dir);
  if(distance(obj.center, new)<= obj.radius && t >= 0){
    copy_vector(new,intersection);
    return 1;
  }
  return 0;
}

//if there is no intersection, this function should return 0
//Otherwise, populate 'intersection' with the point of intersection and return 1
int ray_cylinder_intersection(ray_t observer, cylinder_t obj, vector_t *intersection) {
  //Question 5: Modify this function to compute an intersection

  
  vector_t d = observer.dir;
  vector_t diff = difference(observer.start,obj.center);
  vector_t const1 = sum(d,scalar_product(dot_product(obj.axis,d),obj.axis));
  vector_t const2 = difference(diff, scalar_product(dot_product(obj.axis,diff),obj.axis));
  
  double a = dot_product(const1,const1);
  double b = 2*(dot_product(const1,const2));
  double c = dot_product(const2,const2)-obj.radius;
  double discriminant = b*b - 4*a*c;

   if (discriminant < 0.) return 0;
  
  double h;
  if ((-1*b - sqrt(discriminant))/(2*a) > 0) h = (-1*b - sqrt(discriminant))/(2*a); 
  else if ((-1*b + sqrt(discriminant))/(2*a) > 0) h = (-1*b + sqrt(discriminant))/(2*a); 
  else return 0;

  double lowBound = dot_product(obj.axis, difference(sum(observer.start,scalar_product(h,d)),obj.center));

  double upBound = dot_product(obj.axis, difference(sum(observer.start,scalar_product(h,d)),sum(obj.center,scalar_product(obj.height,obj.axis))));
  
  if (lowBound <= 0 || upBound >= 0) return 0;

  vector_t solution = scaled_sum(1.,observer.start,h,observer.dir);

  copy_vector(solution,intersection);
  return 1; 
}

int ray_cone_intersection(ray_t observer, cone_t obj, vector_t *intersection) {
  //Question 7: Modify this function to compute an intersection
  return 0;
}
