## Compile
Prerequisites
`cmake, hdf5, gcc-9, gmp, mpfr`

Compile, 
```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make breve_bin -j 4
```

## Thingi10k scripts
* Download Thingi10k models without self intersection and manifold. (5172 results, and we provide the list in `fullname_10k_good.txt`)
`https://ten-thousand-models.appspot.com/results.html?q=with+not+self-intersection%2C+is+manifold`
* Use your favorite job scheduler to loop through the files. An example bash script is as follows.
```bash
binary=PATH_TO_BUILD/breve_bin
for i in {1..5172}
do
   ( let "fid = $i"
    filename=$(sed ${fid}'q;d' fullname_10k_good.txt)
    $binary -i $INPUT_DIR$filename -o $OUTPUT_DIR -l $LOG_DIR --loglevel 2 --serial shell --simplify-shell --init-thick 0.0001 --angle 0.001 --thick 0.05)
done
```

## Reproduce Results

* Fig. 1: Teaser Wheel (One componet from `1517923.stl`, `h=0.1`)
* Fig. 3: Example
`V3_solid_virus.stl`  `https://www.thingiverse.com/thing:4166787`
* Fig. 9: Initial bevel pattern and optimized result (`screwdriver.off`, `h=0.05`)
* Fig. 11: Varied target thickness (`78481.stl`, `h=0.5,0.1,0.01,0.001`)
* Fig. 13: Pinching singularity (`1146172.stl`, `h=0.05`)
* Fig. 15: `Thingi10k` Gallery. 
`91006.stl` chains, `56630.stl` earth.
`63785.stl` mouse skull.  `69930` octopus.
`76540.stl` cheese.
`53750.stl` space filling curve
`100506.stl` a long tube


* Fig. 16: Postprocess for cleaner shell extraction. (`camel_b.obj`)

* Fig. 17: Heat Geodesic Comparision (`121868.stl` from `Thingi10k`)
`h=0.005, e=0.02`

* Fig. 18: Making Tetrahedralization Easier. (`knot100K.obj` from `mpz14`)
`knot100k.obj --thick 0.005`
`e=0.02`

* Fig. 19: Post-process Boolean Operation. (Two meshes from `Thingi10k`, processed to `bird_engine.obj`)
```python
import pymesh
mesh = pymesh.load_mesh('../tests/data/bird_engine.obj')
bird, engine = pymesh.separate_mesh(mesh)
birdengine = pymesh.boolean(bird, engine, operation = 'union')
pymesh.save_mesh('birdengine_boolean.ply', birdengine, 'source', 'source_face')
```
`birdengine_boolean.ply --thick 0.01`
`--edge 0.05`

* Fig. 20: Detail Map for Bunny. (`bunny.obj` from `https://github.com/libigl/libigl-tests-data`)

* Fig. 21: Shell Texture. (`animal.obj`)
`h=0.05`
* Fig. 22: Armadillo Multi-Layer Cage. (`armadillo.obj`)

* Fig. 23: Extension to self intersecting or non orientable surfaces.
`leg-intersect.obj`

## Visualize Results