# cabaret-triangle

Trying to implement triangle computational grid to our perfect Cabaret Scheme.

## Grid generation 

We use Triangle for 2-D grid generation. It requires .poly file on input to generate mesh files.

By definition, a .poly file is just a list of vertices and segments. Here is a template for it:

```
First line: <# of vertices> <dimension (must be 2)> <# of attributes> <# of boundary markers (0 or 1)>
Following lines: <vertex #> <x> <y> [attributes] [boundary marker]
One line: <# of segments> <# of boundary markers (0 or 1)>
Following lines: <segment #> <endpoint> <endpoint> [boundary marker]
One line: <# of holes>
Following lines: <hole #> <x> <y>
Optional line: <# of regional attributes and/or area constraints>
Optional following lines: <region #> <x> <y> <attribute> <maximum area>
```

A .poly file can also contain information about holes and concavities, as well as regional attributes and constraints on the areas of triangles.

Steps for grid generation:

Edit PolyGenerator.cpp to set up (boundary) and (amount of cells) values.
Then compile and run it via:
```
g++ PolyGenerator.cpp && ./a.out
```

To generate mesh files by Triangle use command(you should make triangle before):
```
./triangle -pqeD Mesh.poly
```

Now you have all mesh files to run main programm.