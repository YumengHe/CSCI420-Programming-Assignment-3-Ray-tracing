/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: Yumeng He
 * *************************
*/

#ifdef WIN32
  #include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
  #include <GL/gl.h>
  #include <GL/glut.h>
#elif defined(__APPLE__)
  #include <OpenGL/gl.h>
  #include <GLUT/glut.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef WIN32
  #define strcasecmp _stricmp
#endif

#include <imageIO.h>

// helper file
#include "vec3.h"
#include "color.h"
#include "ray.h"
#include <cstdlib>
#include <ctime>
#include <random>

#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

char * filename = NULL;

// The different display modes.
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;

// While solving the homework, it is useful to make the below values smaller for debugging purposes.
// The still images that you need to submit with the homework should be at the below resolution (640x480).
// However, for your own purposes, after you have solved the homework, you can increase those values to obtain higher-resolution images.
// #define WIDTH 640
// #define HEIGHT 480

#define WIDTH 640
#define HEIGHT 480

double ratio = WIDTH / HEIGHT;

// The field of view of the camera, in degrees.
#define fov 60.0

// Buffer to store the image when saving it to a JPEG.
unsigned char buffer[HEIGHT][WIDTH][3];

#define MAX_DEPTH 5

char pushed_token[100] = "";

struct Vertex {
  // variables
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
  double reflectivity;
  // constructor
  Vertex() : reflectivity(0.0) {}
};

struct Triangle
{
  Vertex v[3];
};

struct Sphere {
  // variables
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double reflectivity;
  double radius;
  // constructor
  Sphere() : reflectivity(0.0) {}
};

struct Light
{
  // variables
  double position[3];
  double color[3];
  bool is_area_light;
  vec3 area_u;     // one egde
  vec3 area_v;     // another edge
  int samples;     // number of samples for soft shadows
  // constructor
  Light() : is_area_light(false), area_u(vec3(0, 0, 0)), area_v(vec3(0, 0, 0)), samples(1) {}

};

struct Material {
  color diffuse;
  color specular;
  double shininess;
  double reflectivity;
};

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void parse_check(const char *expected, const char *found);

// ---------- new functions -----------

vec3 interpolate_normal(const Triangle& triangle, double u, double v) {
  // vertices normal
  vec3 n0(triangle.v[0].normal[0], triangle.v[0].normal[1], triangle.v[0].normal[2]);
  vec3 n1(triangle.v[1].normal[0], triangle.v[1].normal[1], triangle.v[1].normal[2]);
  vec3 n2(triangle.v[2].normal[0], triangle.v[2].normal[1], triangle.v[2].normal[2]);

  // interpolate normals
  // u + v + w = 1
  double w = 1.0 - u - v;
  vec3 normal = w * n0 + u * n1 + v * n2;

  // normalize
  return unit_vector(normal);
}

color interpolate_color(const Triangle& triangle, double u, double v, bool diffuse) {
  double c0[3], c1[3], c2[3];

  if (diffuse) {
    memcpy(c0, triangle.v[0].color_diffuse, 3 * sizeof(double));
    memcpy(c1, triangle.v[1].color_diffuse, 3 * sizeof(double));
    memcpy(c2, triangle.v[2].color_diffuse, 3 * sizeof(double));
  } else {
    memcpy(c0, triangle.v[0].color_specular, 3 * sizeof(double));
    memcpy(c1, triangle.v[1].color_specular, 3 * sizeof(double));
    memcpy(c2, triangle.v[2].color_specular, 3 * sizeof(double));
  }

  double w = 1.0 - u - v;
  double r = w * c0[0] + u * c1[0] + v * c2[0];
  double g = w * c0[1] + u * c1[1] + v * c2[1];
  double b = w * c0[2] + u * c1[2] + v * c2[2];

  return color(r, g, b);
}

double interpolate_shininess(const Triangle& triangle, double u, double v) {
  double s0 = triangle.v[0].shininess;
  double s1 = triangle.v[1].shininess;
  double s2 = triangle.v[2].shininess;

  double w = 1.0 - u - v;
  return w * s0 + u * s1 + v * s2;
}

color clamp_color(const color& c) {
  return color(std::min(1.0, c.x), std::min(1.0, c.y), std::min(1.0, c.z));
}

double hit_sphere(vec3 center, double radius, const ray& r){
  vec3 oc = center - r.origin();
  // a = d * d
  double a = r.direction().length_squared();
  double h = dot(r.direction(), oc);
  double c = oc.length_squared() - radius * radius;
  double discriminant = h * h - a * c;

  // find normal
  if (discriminant < 0) {
    return -1.0;
  } else {
    return (h - std::sqrt(discriminant)) / a;
  }
}

vec3 reflect(const vec3& v, const vec3& n) {
  return v - 2 * dot(v, n) * n;
}

bool hit_triangle(const Triangle& triangle, const ray& r, double& t, double& u, double& v) {
  const double EPSILON = 1e-8;
  // triangle vertices
  vec3 v0(triangle.v[0].position[0], triangle.v[0].position[1], triangle.v[0].position[2]);
  vec3 v1(triangle.v[1].position[0], triangle.v[1].position[1], triangle.v[1].position[2]);
  vec3 v2(triangle.v[2].position[0], triangle.v[2].position[1], triangle.v[2].position[2]);

  // two edges
  vec3 edge1 = v1 - v0;
  vec3 edge2 = v2 - v0;

  // determinant
  vec3 h = cross(r.direction(), edge2);
  double a = dot(edge1, h);

  // the ray lies in the plane of the triangle
  if (fabs(a) < EPSILON)
      return false; // the ray is parallel to triangle

  double f = 1.0 / a;
  vec3 s = r.origin() - v0;
  u = f * dot(s, h);

  // check intersection
  if (u < 0.0 || u > 1.0)
    return false;

  vec3 q = cross(s, edge1);
  v = f * dot(r.direction(), q);

  // check intersection
  if (v < 0.0 || u + v > 1.0)
    return false;

  t = f * dot(edge2, q);

  if (t > EPSILON) { // ray intersection
    return true;
  } else {
    return false; // line intersection but not a ray intersection
  }
}

double interpolate_property(const Triangle& triangle, double u, double v, const char* property) {
  double p0, p1, p2;
  if (strcmp(property, "shininess") == 0) {
    p0 = triangle.v[0].shininess;
    p1 = triangle.v[1].shininess;
    p2 = triangle.v[2].shininess;
  } else if (strcmp(property, "reflectivity") == 0) {
    p0 = triangle.v[0].reflectivity;
    p1 = triangle.v[1].reflectivity;
    p2 = triangle.v[2].reflectivity;
  } else { // error or default cases
    return 0.0;
  }
  double w = 1.0 - u - v;
  return w * p0 + u * p1 + v * p2;
}

color ray_color(const ray& r, int depth)
{ 
  // background color
  if (depth > MAX_DEPTH) {
    return color(0, 0, 0);
  }

  double closest_t = std::numeric_limits<double>::infinity();
  color final_color(0, 0, 0);
  bool hit_anything = false;

  vec3 intersection_point;
  vec3 normal_at_intersection;
  Material material_at_intersection;

  // sphere intersections
  for (int i = 0; i < num_spheres; i++) {
    vec3 center(spheres[i].position[0], spheres[i].position[1], spheres[i].position[2]);
    double t = hit_sphere(center, spheres[i].radius, r);
    if (t > 0.0 && t < closest_t) {
      closest_t = t;
      hit_anything = true;

      // intersection point
      intersection_point = r.at(t);

      // normal at intersection point
      normal_at_intersection = unit_vector(intersection_point - center);

      // material properties
      material_at_intersection.diffuse = color(spheres[i].color_diffuse[0], spheres[i].color_diffuse[1], spheres[i].color_diffuse[2]);
      material_at_intersection.specular = color(spheres[i].color_specular[0], spheres[i].color_specular[1], spheres[i].color_specular[2]);
      material_at_intersection.shininess = spheres[i].shininess;
      material_at_intersection.reflectivity = spheres[i].reflectivity;
    }
  }

  // triangle intersections
  for (int i = 0; i < num_triangles; i++) {
    double t, u, v;
    if (hit_triangle(triangles[i], r, t, u, v) && t < closest_t) {
      closest_t = t;
      hit_anything = true;

      // intersection point
      intersection_point = r.at(t);

      // normal using barycentric coordinates (u, v)
      normal_at_intersection = interpolate_normal(triangles[i], u, v);

      // material properties
      material_at_intersection.diffuse = interpolate_color(triangles[i], u, v, true);  // true for diffuse
      material_at_intersection.specular = interpolate_color(triangles[i], u, v, false); // false for specular
      material_at_intersection.shininess = interpolate_property(triangles[i], u, v, "shininess");
      material_at_intersection.reflectivity = interpolate_property(triangles[i], u, v, "reflectivity");
    }
  }

  if (hit_anything) {
    // ambient light contribution
    color ambient(ambient_light[0], ambient_light[1], ambient_light[2]);
    final_color = ambient * material_at_intersection.diffuse;

    for (int i = 0; i < num_lights; i++) {
      Light& light = lights[i];

      // light contribution accumulator
      double light_contrib = 0.0;

      // determine the number of samples
      int num_samples = light.is_area_light ? light.samples : 1;

      // initialize random number generator
      static bool rand_initialized = false;
      if (!rand_initialized) {
        srand((unsigned int)time(NULL));
        rand_initialized = true;
      }

      for (int s = 0; s < num_samples; ++s) {
        vec3 light_sample_pos;
        
        // area lights
        if (light.is_area_light) 
        {
          // Sample a point on the area light using stratified sampling
          int sqrt_samples = (int)sqrt((double)light.samples);
          int u_index = s % sqrt_samples;
          int v_index = s / sqrt_samples;
          double rand_u = (u_index + ((double)rand() / RAND_MAX)) / sqrt_samples;
          double rand_v = (v_index + ((double)rand() / RAND_MAX)) / sqrt_samples;

          // Adjust to [-0.5, 0.5] range
          rand_u -= 0.5;
          rand_v -= 0.5;

          // Sample point on area light
          light_sample_pos = vec3(light.position[0], light.position[1], light.position[2]) + rand_u * light.area_u + rand_v * light.area_v;
        } 
        // point light
        else 
        {
          light_sample_pos = vec3(light.position[0], light.position[1], light.position[2]);
        }

        // Compute direction to the light sample
        vec3 L = unit_vector(light_sample_pos - intersection_point);

        // Shadow ray
        vec3 shadow_ray_origin = intersection_point + 1e-4 * L;
        ray shadow_ray(shadow_ray_origin, L);

        // Distance to light sample
        double light_distance = (light_sample_pos - intersection_point).length();

        // Shadow check
        bool in_shadow = false;

        // Check for occlusion with spheres
        for (int j = 0; j < num_spheres; ++j) {
          vec3 center(spheres[j].position[0], spheres[j].position[1], spheres[j].position[2]);
          double t_shadow = hit_sphere(center, spheres[j].radius, shadow_ray);
          if (t_shadow > 0.0 && t_shadow < light_distance) {
            in_shadow = true;
            break;
          }
        }

        // Check for occlusion with triangles
        if (!in_shadow) 
        {
          for (int j = 0; j < num_triangles; ++j) {
            double t_shadow, u_shadow, v_shadow;
            if (hit_triangle(triangles[j], shadow_ray, t_shadow, u_shadow, v_shadow) &&
              t_shadow < light_distance) {
              in_shadow = true;
              break;
            }
          }
        }

        // accumulate light contribution
        if (!in_shadow) 
        {
          light_contrib += 1.0;
        }
      }

      // Average the light contribution
      light_contrib /= num_samples;

      // If there's any light contribution, compute shading
      if (light_contrib > 0.0) {
        // Use the direction to the center of the light for shading calculations
        vec3 light_dir = unit_vector(vec3(light.position[0], light.position[1], light.position[2]) - intersection_point);
        double LdotN = std::max(0.0, dot(light_dir, normal_at_intersection));

        // Compute reflection vector
        vec3 R = reflect(-light_dir, normal_at_intersection);

        // Direction to viewer (camera is at origin)
        vec3 V = unit_vector(-r.direction());

        // Compute specular component
        double RdotV = std::max(0.0, dot(R, V));
        double specular_factor = pow(RdotV, material_at_intersection.shininess);

        // Light color
        color light_color(light.color[0], light.color[1], light.color[2]);

        // Accumulate color contribution scaled by light_contrib
        final_color = final_color + light_contrib * light_color * (material_at_intersection.diffuse * LdotN + material_at_intersection.specular * specular_factor);
      }
    }

    // Reflection
    if (material_at_intersection.reflectivity > 0.0) {
      // Compute reflection ray
      vec3 reflect_dir = reflect(unit_vector(r.direction()), normal_at_intersection);
      ray reflect_ray(intersection_point + 1e-4 * reflect_dir, reflect_dir);

      // Recursively trace the reflection ray
      color reflect_color = ray_color(reflect_ray, depth + 1);

      // Blend the reflection color with the local color
      final_color = (1.0 - material_at_intersection.reflectivity) * final_color + material_at_intersection.reflectivity * reflect_color;
    }

    // Clamp final color to [0, 1]
    final_color = clamp_color(final_color);

    return final_color;
  } 
  else // background color
  {
    return color(1.0, 1.0, 1.0);
  }
}

// detect key input
void keyboardFunc(unsigned char key, int x, int y)
{
  switch (key)
  {
    case 27: // ESC key
      exit(0); // exit the program
    break;
  }
}

// render
void draw_scene()
{
  // Camera setup
  double focal_length = 1.0;
  double viewport_height = 2.0;
  double viewport_width = viewport_height * ((double)WIDTH / (double)HEIGHT); // Maintain aspect ratio
  vec3 camera_center = vec3(0, 0, 0);

  // Viewport vectors
  vec3 viewport_u = vec3(viewport_width, 0, 0); // x-direction
  vec3 viewport_v = vec3(0, viewport_height, 0); // y-direction (positive)

  // Delta increments for each pixel
  vec3 delta_u = viewport_u / WIDTH;
  vec3 delta_v = viewport_v / HEIGHT;

  // Upper-left corner of the viewport
  vec3 viewport_upper_left = camera_center - vec3(0, 0, focal_length) - viewport_u / 2 + viewport_v / 2;
  vec3 pixel00 = viewport_upper_left + 0.5 * (delta_u - delta_v);

  // Render the scene

  // Number of samples per pixel (e.g., 4 for 2x2 grid)
  const int samples_per_pixel = 16; // Adjust as needed for quality vs. performance

  // Random number generator setup
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(0.0, 1.0);

  // Render the scene
  for (int y = 0; y < HEIGHT; y++) {
    glPointSize(2.0);
    glBegin(GL_POINTS);
    for (int x = 0; x < WIDTH; x++) 
    {
      color pixel_color(0, 0, 0);

      // Supersampling
      for (int s = 0; s < samples_per_pixel; ++s)
      {
        // Generate random offsets within the pixel
        double rand_u = distribution(generator);
        double rand_v = distribution(generator);

        vec3 pixel_sample = pixel00 + (x + rand_u) * delta_u - (y + rand_v) * delta_v;
        vec3 ray_direction = unit_vector(pixel_sample - camera_center);
        ray r(camera_center, ray_direction);

        pixel_color = pixel_color + ray_color(r, 0);
      }

      // Average the color
      pixel_color = pixel_color / samples_per_pixel;

      // Clamp color components to [0,1]
      pixel_color = clamp_color(pixel_color);

      // Convert color components from [0,1] to [0,255]
      unsigned char red = static_cast<unsigned char>(255.999 * pixel_color.x);
      unsigned char green = static_cast<unsigned char>(255.999 * pixel_color.y);
      unsigned char blue = static_cast<unsigned char>(255.999 * pixel_color.z);

      plot_pixel(x, y, red, green, blue);
    }
    glEnd();
    glFlush();
  }
  printf("Ray tracing completed.\n");
  fflush(stdout);
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  buffer[y][x][0] = r;
  buffer[y][x][1] = g;
  buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  int adjusted_y = HEIGHT - y - 1; // flip the y coordinate
  plot_pixel_display(x, adjusted_y, r, g, b);
  if(mode == MODE_JPEG)
    plot_pixel_jpeg(x, adjusted_y, r, g, b); // use adjusted_y here
}

void save_jpg()
{
  printf("Saving JPEG file: %s\n", filename);

  ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
  if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
    printf("Error in saving\n");
  else 
    printf("File saved successfully\n");
}

void read_token(FILE* file, char* token)
{
  if (pushed_token[0] != '\0') {
    strcpy(token, pushed_token);
    pushed_token[0] = '\0';  // clear token
    return;
  }

  // read the next token
  int ch;
  int index = 0;

  while (1) {
    ch = fgetc(file);

    // Handle EOF
    if (ch == EOF) {
      token[0] = '\0';
      return;
    }

    // whitespace
    if (isspace(ch))
      continue;

    // comments
    if (ch == '#')
    {
      while ((ch = fgetc(file)) != EOF && ch != '\n');
      continue;
    }
    if (ch == '/' && (ch = fgetc(file)) == '/') 
    {
      while ((ch = fgetc(file)) != EOF && ch != '\n');
      continue;
    } 
    else if (ch == '/') 
    {
      ungetc(ch, file);
      ch = '/';
    }
    break;
  }

  // read the token
  do {
    if (ch == '#' || (ch == '/' && (ch = fgetc(file)) == '/')) {
      if (ch == '/') {
        ungetc(ch, file);
      }
      ungetc('/', file);
      break;
    }
    token[index++] = ch;
    ch = fgetc(file);
  } while (ch != EOF && !isspace(ch));

  // null-terminate
  token[index] = '\0';

  // skip rest of the line if comment
  if (ch == '#') {
    while ((ch = fgetc(file)) != EOF && ch != '\n');
  } else if (ch == '/' && (ch = fgetc(file)) == '/') {
    while ((ch = fgetc(file)) != EOF && ch != '\n');
  } else if (ch != EOF) {
    ungetc(ch, file);
  }
}

void push_token(const char* token)
{
  strcpy(pushed_token, token);
}

void parse_check(const char *expected, const char *found)
{
  if (strcasecmp(expected, found))
  {
    printf("Expected '%s' found '%s'\n", expected, found);
    printf("Parsing error; abnormal program abortion.\n");
    exit(0);
  }
}

void parse_bool(FILE* file, const char* check, bool* value) {
  char str[100];
  read_token(file, str);
  if (strcasecmp(str, "true") == 0) {
      *value = true;
  } else if (strcasecmp(str, "false") == 0) {
      *value = false;
  } else {
      printf("Error parsing boolean value for '%s'\n", check);
      exit(0);
  }
  printf("%s %s\n", check, *value ? "true" : "false");
}

void parse_int(FILE* file, const char* check, int* value) {
  char str[100];
  read_token(file, str);
  if (str[0] == '\0') {
      printf("Error parsing integer value for '%s'\n", check);
      exit(0);
  }
  *value = atoi(str);
  printf("%s %d\n", check, *value);
}

void parse_doubles(FILE* file, const char *check, double p[3], bool check_property = true)
{
  char str[100];
  if (check_property) {
    read_token(file, str);
    parse_check(check, str);
  }
  for (int i = 0; i < 3; ++i) {
    read_token(file, str);
    if (str[0] == '\0') {
      printf("Error parsing doubles for '%s'\n", check);
      exit(0);
    }
    p[i] = atof(str);
  }
  printf("%s %lf %lf %lf\n", check, p[0], p[1], p[2]);
}

void parse_rad(FILE* file, double* r)
{
  char str[100];
  read_token(file, str);
  parse_check("rad:", str);
  read_token(file, str);
  if (str[0] == '\0') {
    printf("Error parsing radius\n");
    exit(0);
  }
  *r = atof(str);
  printf("rad: %f\n", *r);
}

void parse_shi(FILE* file, double* shi)
{
  char s[100];
  read_token(file, s);
  parse_check("shi:", s);
  read_token(file, s);
  if (s[0] == '\0') {
      printf("Error parsing shininess value\n");
      exit(0);
  }
  *shi = atof(s);
  printf("shi: %f\n", *shi);
}

void parse_reflectivity(FILE* file, double* reflectivity)
{
  char s[100];
  read_token(file, s);
  if (strcasecmp("ref:", s) == 0)
  {
    read_token(file, s);
    if (s[0] == '\0') {
        printf("Error parsing reflectivity value\n");
        exit(0);
    }
    *reflectivity = atof(s);
    printf("ref: %f\n", *reflectivity);
  }
  else
  {
    *reflectivity = 0.0;
    push_token(s);
  }
}

int loadScene(char *argv)
{
  FILE * file = fopen(argv,"r");
  if (!file)
  {
    printf("Unable to open input file %s. Program exiting.\n", argv);
    exit(0);
  }

  int number_of_objects;
  char type[50];
  Triangle t;
  Sphere s;
  Light l;
  char str[100];

  // Read number of objects
  read_token(file, str);
  if (str[0] == '\0') {
    printf("Error reading number of objects\n");
    exit(0);
  }
  number_of_objects = atoi(str);
  printf("number of objects: %i\n", number_of_objects);

  // Read ambient light
  parse_doubles(file, "amb:", ambient_light);

  for (int i = 0; i < number_of_objects; i++)
  {
    read_token(file, type);
    if (type[0] == '\0') {
      printf("Error reading object type\n");
      exit(0);
    }
    printf("%s\n", type);
    if (strcasecmp(type, "triangle") == 0)
    {
      printf("found triangle\n");
      for (int j = 0; j < 3; j++)
      {
        parse_doubles(file, "pos:", t.v[j].position);
        parse_doubles(file, "nor:", t.v[j].normal);
        parse_doubles(file, "dif:", t.v[j].color_diffuse);
        parse_doubles(file, "spe:", t.v[j].color_specular);
        parse_shi(file, &t.v[j].shininess);
        parse_reflectivity(file, &t.v[j].reflectivity);
      }

      if (num_triangles == MAX_TRIANGLES)
      {
        printf("too many triangles, you should increase MAX_TRIANGLES!\n");
        exit(0);
      }
      triangles[num_triangles++] = t;
    }
    else if (strcasecmp(type, "sphere") == 0)
    {
      printf("found sphere\n");

      parse_doubles(file, "pos:", s.position);
      parse_rad(file, &s.radius);
      parse_doubles(file, "dif:", s.color_diffuse);
      parse_doubles(file, "spe:", s.color_specular);
      parse_shi(file, &s.shininess);
      parse_reflectivity(file, &s.reflectivity);

      if (num_spheres == MAX_SPHERES)
      {
        printf("too many spheres, you should increase MAX_SPHERES!\n");
        exit(0);
      }
      spheres[num_spheres++] = s;
  }
    else if (strcasecmp(type, "light") == 0) {
    printf("found light\n");
    parse_doubles(file, "pos:", l.position);
    parse_doubles(file, "col:", l.color);

    // area light properties
    char next_token[100];
    while (true) {
      read_token(file, next_token);
      if (strcasecmp(next_token, "area_light:") == 0) {
        parse_bool(file, "area_light:", &l.is_area_light);
      } else if (strcasecmp(next_token, "area_u:") == 0) {
        double area_u[3];
        parse_doubles(file, "area_u:", area_u, false);
        l.area_u = vec3(area_u[0], area_u[1], area_u[2]);
      } else if (strcasecmp(next_token, "area_v:") == 0) {
        double area_v[3];
        parse_doubles(file, "area_v:", area_v, false);
        l.area_v = vec3(area_v[0], area_v[1], area_v[2]);
      } else if (strcasecmp(next_token, "samples:") == 0) {
        parse_int(file, "samples:", &l.samples);
      } else {
        push_token(next_token);
        break;
      }
    }

    if (num_lights == MAX_LIGHTS) {
      printf("too many lights, you should increase MAX_LIGHTS!\n");
      exit(0);
    }
    lights[num_lights++] = l;
    }
  else
    {
      printf("unknown type in scene description:\n%s\n", type);
      exit(0);
    }
  }
  return 0;
}

void display()
{
}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  // Hack to make it only draw once.
  static int once=0;
  if(!once)
  {
    draw_scene();
    if(mode == MODE_JPEG)
      save_jpg();
  }
  once=1;
}

int main(int argc, char ** argv)
{
  if ((argc < 2) || (argc > 3))
  {  
    printf ("Usage: %s <input scenefile> [output jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
  {
    mode = MODE_JPEG;
    filename = argv[2];
  }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  #ifdef __APPLE__
    // This is needed on recent Mac OS X versions to correctly display the window.
    glutReshapeWindow(WIDTH - 1, HEIGHT - 1);
  #endif
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  // callback for pressing the keys on the keyboard
  glutKeyboardFunc(keyboardFunc);
  glutMainLoop();
}

