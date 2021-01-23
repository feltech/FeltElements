# Stanford Bunny

Bunny triangle mesh "zipper" reconstruction .ply file retrieved from 
http://graphics.stanford.edu/data/3Dscanrep on 01/2021.

Holes filled using [OpenFlipper](www.openflipper.org) and the output saved as .obj. 

Converted to volume tetrahedral MEDIT .mesh via
[pygalmesh](https://github.com/nschloe/pygalmesh/#volume-meshes-from-surface-meshes), converted 
from .mesh to ASCII .vtk via [meshio](https://github.com/nschloe/meshio), then finally converted 
from .vtk to .ovm via [OpenFlipper](www.openflipper.org).

Possibly other combinations of filetypes would work, but this is the first combo that worked for
me.



