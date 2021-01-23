# Converting a triangle mesh to a tetrahedral OpenVolumeMesh

Requires
* [OpenFlipper](www.openflipper.org) for mesh operations, visualisation, and crucially 
  OpenVolumeMesh .ovm format support.
* [pygalmesh](https://github.com/nschloe/pygalmesh) to generate tetrahedral meshes from triangle
  meshes.
* [meshio](https://github.com/nschloe/meshio) to massage the file format into something 
  OpenFlipper can understand (this is a dependency of `pygalmesh` so does not 
  need installing separately).

1. Open triangle mesh in OpenFlipper and use "Hole Filler" to ensure a 
   closed mesh.
2. Then
```sh
# Read a surface mesh (.obj works, but other formats should work too) and generate a MEDIT .mesh
# tetrahedron mesh (.vtk may work too).
pygalmesh-volume-from-surface triangles.obj tetrahedrons.mesh \
    --max-radius-surface-delaunay-ball 0.001 --max-circumradius-edge-ratio 1.5 
# OpenFlipper only supports ASCII .vtk files (pygalmesh vtk output would have been binary),
# so convert.
meshio-convert --ascii --input-format medit --output-format vtk tetrahedrons.mesh tetrahedrons.vtk
```
3. Use OpenFlipper GUI to load .vtk file and save out as .ovm.