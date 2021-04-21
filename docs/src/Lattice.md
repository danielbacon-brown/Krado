# Lattice


## Defining a lattice

A 1D lattice object can be created using:
    lattice = Lattice( 300*nm )
where the lattice vector is the X direction.
    
The periodicity of the lattice vector can be altered by calling with a vector:
```julia
    lattice = Lattice( [200*nm, 200*nm] ) 
```
where L₁ is a real-valued 2-vector.
    
A 2D lattice can be created using two lattice vectors:
```julia
    lattice = Lattice( [100*nm, 0], [-100/2*nm, sqrt(3)/2 * 100nm])
```
where L₁ and L₂ are both real-valued 2-vectors.  The line above creates a hexagonal lattice.
        
You may wish to change the location of the unit cell.  By default, the (0,0) UV coordinate is at the (0,0) XY coordinate.  If both lattice vectors are positive, the bottom left corner of the unit cell origin will be at the origin.
To put the center of the unit cell at the origin, call the function with the optional parameter originOffsetUV as follows:
```julia
    Lattice([100*nm, 0], [-100/2*nm, sqrt(3)/2 * 100nm]; originOffsetUV = [-0.5, -0.5] )
```
which shifts the unit by half of the lattice vectors.
To shift the position of the unit cell in XY coordinates, use originOffsetXY:
```julia
    Lattice([100*nm, 0], [-100/2*nm, sqrt(3)/2 * 100nm]; originOffsetUV = [300*nm, 500nm] )
```

```@docs
Lattice
Lattice(L₁,L₂;originOffsetUV=[],originOffsetXY=[])
Lattice(L₁; originOffsetUV = [], originOffsetXY = [])
Lattice(l₁::Real; originOffsetUV = [], originOffsetXY = [])
```

## Useful lattice functions
Sometimes you want to convert between an X,Y coordinate system with a U,V coordinate system.  Use these functions to perform the conversion for any given lattice.
```@docs
convertUVtoXY(lattice::Lattice, uv)
convertXYtoUV(lattice::Lattice, xu)
getLatticeCenterCoordinates(lattice::Lattice)
```


# Advanced functions

## Calculating reciprocal vectors
These functions calculate the reciprocal lattice vectors based on the input real-space vectors.
```@docs
calcReciprocals(L₁,L₂)
calcReciprocals(L)
is1D(lattice::Lattice)
```
## Lattice boundary
```@docs
calcLatticeBoundaryLine(lattice::Lattice)
calcReciprocalLatticeBoundaryLine(lattice::Lattice)
```
