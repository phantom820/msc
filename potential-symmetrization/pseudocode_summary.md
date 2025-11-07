# Summary 
### We've converted 2,003,909 Markdown files to PDF and counting!

Note assume we have a system with 8 atoms, and total of 6 displacements possible for an atom ( 2 for the +- and 3 directions)

### Step 1:
Build a lookup table of how the atoms map to each other under different displacements + rotation. So we would have something like 

```
0 x + R 1 -> 1 y - 
0 x + R2 -> 2 z +
0 x - R2 -> 1 x +
......
```

This will be an exhaustive mapping that we will use further down as the lookup say its $M(N, d, s, R)$ where $N, d, s , R$ represent the atom number, displacement axis $x|y|z$ , the sign + | - and the rotation matrix

### Step 2 Generate All configurations
This would be picking an atom and then applying all possible displacements to it. So we have that 1 atom generates

```
1 Atom -> 6 Configurations 
```

In total we would have $8 \times 6 = 48$ configurations. We save all of these somewhere say a dictionary. These are all "unique"

$$ C = \{ C_1, C_2, ....., C_{N}\}$$


### Step 3 Pick a configuration and find its symmetries 
From all set of configurations remove configuration say $C_i$. Now take $C_i$ apply a rotation matrix to it and get $C_i^*$. For $C_i^*$ you know  $N, d, s , R$  and from the lookup table if you plug these values in you get a set of mappings that were holding true initially and you then need to check if these are still true i.e the following equivalence holds

$$ C_i^*(N, r) =  C_j([M(N, d, s, R)]_{C_i^*})$$

The above is intended to read for an alternate configuration that was initially mapped/symmetric to $C_i^*$ check if the equivalence still holds if it does, then note this symetry. 

You need to do this till the initial set of all configurations is empty, at which point you will end up with map of the form


```
C_1 -> { R_1 -> C_2, R_4 -> {C_5}}
C_2 -> {R_3 -> C_1 , R_12 -> C_6}
```

This is my understanding and wanted to confirm this correct before doing cleanup and implementation