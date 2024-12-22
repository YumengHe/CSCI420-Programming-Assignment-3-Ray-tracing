# CSCI420 Programming Assignment 3: Ray tracing

**Assignment Link:** [Ray tracing](https://odedstein.com/teaching/hs-2024-csci-420/assign3/)

**Final Grade:** 120 out of 100

## Introduction 
In this assignment, we built a ray tracer. The ray tracer is be able to handle opaque surfaces with lighting and shadows.

## MANDATORY FEATURES
1) Ray tracing triangles 
2) Ray tracing sphere 
3) Triangle Phong Shading 
4) Sphere Phong Shading 
5) Shadows rays 
6) Still images 

## Usage
To run the program, use the follwing command:
```bash
./hw3 <scene file> <image name>
./hw3 <scene file>
```

Examples:
```bash
./hw3 test1.scene test1.jpg
./hw3 test2.scene test2.jpg
./hw3 test2edit.scene test2edit.jpg
./hw3 new.scene new.jpg
./hw3 snow.scene snow.jpg
./hw3 SIGGRAPH.scene SIGGRAPH.jpg
```

## Keyboard Commands
- **esc:** Exit the render.

## Additional Features Implemented

### Recursive Reflection and Extended Material Properties

- **Recursive Reflection**: Implemented recursive ray tracing to handle reflections, allowing surfaces to reflect other objects in the scene. When a ray hits a reflective surface, a reflection ray is calculated using the incident ray and the surface normal. The color from the reflected ray is blended with the local shading based on the material's reflectivity. This enables realistic rendering of mirrors and shiny objects.

- **Extended Scene Description with "ref" and "shi" Properties**: Enhanced the scene file format to include two new properties for materials:

  - **"ref"**: Reflectivity coefficient, specifying how much a surface reflects other objects. This property is parsed from the scene file and used in the recursive reflection calculations to control the amount of reflection for each material.

  - **"shi"**: Shininess coefficient, specifying the specular exponent in the Phong reflection model. By parsing this property from the scene file, we allow materials to have varying levels of shininess, affecting the sharpness of specular highlights.

  These extensions provide greater flexibility in defining material properties directly within the scene files, enabling more detailed and realistic rendering of surfaces.

### Good Antialiasing

- **Supersampling Antialiasing**: Implemented supersampling by casting multiple rays per pixel at slightly different offsets within the pixel area. The colors from these rays are averaged to compute the final pixel color. This reduces aliasing artifacts such as jagged edges, resulting in smoother and higher-quality images. The number of samples per pixel can be adjusted to balance between image quality and rendering time.

### Soft Shadows

- **Area Lights and Multiple Shadow Rays**: Introduced area lights by allowing light sources to have spatial extent instead of being point sources. Multiple shadow rays are cast towards different sample points on the area light source. The proportion of unblocked shadow rays determines the softness of the shadow, creating more realistic soft shadow effects with penumbra regions. This adds depth and realism to the scene by simulating how light behaves in the real world.

### Möller-Trumbore Algorithm for Triangle Intersection

- **Efficient Ray-Triangle Intersection**: Implemented the Möller-Trumbore algorithm for ray-triangle intersection testing. This algorithm is efficient and numerically robust, computing intersections using barycentric coordinates without explicitly calculating the plane equation of the triangle. It improves performance and accuracy when rendering scenes with many triangles, especially for complex models and meshes.

### Comments in Scene Files

- **Support for Comments**: Enhanced the scene file parser to allow comments within the scene description files. Lines starting with `#` or `//` are treated as comments and ignored during parsing. This allows for better organization and documentation within scene files without affecting the parsing process. Users can add explanations, disable parts of the scene, or provide any notes directly in the scene file.

## Helper Files

- **vec3.h**: Defines a 3D vector class `vec3` with common vector operations such as addition, subtraction, dot product, cross product, and normalization. Used extensively for points, directions, normals, and color calculations.

- **ray.h**: Defines a `ray` class representing a ray with an origin and a direction. This is fundamental for ray casting, as rays are traced from the camera through each pixel into the scene to detect intersections with objects.

- **color.h**: Defines a `color` class for RGB color representation, supporting operations like addition, multiplication, scaling, and clamping to ensure color values remain within valid ranges. Used to compute and store the color at each pixel based on lighting and material properties.

## Sample Image
- output of new.scene
![new.scene](<hw3-starterCode/new.jpg>)

- output of test2.scene
![test2.scene](<hw3-starterCode/test2.jpg>)

- output of test2edit.scene
![test2edit.scene](<hw3-starterCode/test2edit.jpg>)

- output of SIGGRAPH.scene
![SIGGRAPH.scene](<hw3-starterCode/SIGGRAPH.jpg>)