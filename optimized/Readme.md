### Atom Mapping
#### Equilibrium System 
The represantion for the atom to atom mapping needs to be able to easily support answering the following question. 

<em>Which atoms are symmetrical to each other under particular operators ?</em>

It needs to capture the relation $i \to j$ under a set of operators $\{O_1, O_2, .....,O_k\} $
```json
{
  "Atom Index" : {
    "Atom Index" : "List of operator indices"
  }
}

```
i.e snippet for atom 1
```json
 "1": {
    "1": [
      0,
      1,
      6,
      7,
      10,
      11
    ],
    "5": [
      2,
      3,
      8,
      9,
      12,
      13
    ]
 }

```
In the above represention we have atom $1 \to 1$ under operators $\{0, 1, 6, 7, 10, 11\}$ and $1 \to 5$ under operators $\{2, 3, 8, 9, 12, 13\}$. Note $5 \to 1$ under $\{2, 3, 8, 9, 12, 13\}^{-1}$.

#### Displaced System
For the displaced system, the atom to atom mapping needs to incorporate the displacement axis $x,y,z$ and direction $+, -$ under particular operators. It needs to capture the relation $(i, axis, direction) \to (j, axis, direction) $ under a set of operators $\{O_1, O_2, .....,O_k\} $. So we can use 

```json
{
  "Atom Index|Axis|Direction" : {
    "Atom Index|Axis|Direction" : [
      "List of operator indices"
    ]
  }
}

```
i.e
```json
"1|x|-": {
    "1|x|-": [
      0,
      1
    ],
    "1|y|-": [
      6,
      7
    ],
    "1|z|-": [
      10,
      11
    ],
    "5|y|+": [
      2,
      3
    ],
    "5|x|+": [
      8,
      9
    ],
    "5|z|-": [
      12,
      13
    ]
}
```
In the above represention we have atom $(1, x, -) \to (5,y, +)$ under operators $\{2, 3\}$ and so forth. 


### Potential Mapping

#### Equilibrium potential
The potential is represented in real space grid and want to represent how an index $i \to j$ under an operator $O$. The representation needs to answer 

<em>What is the symmetric potential under operator $O$ ?</em>

```json
"Operator Index": "List of alternate indices"

```

i.e
```json
"1": [
    0,
    10,
    20,
    30,
    40,
    50,
    .....
]
```

the above captures that operator $0$ induces the indice mapping $0 \to 0, 1 \to 10, 2 \to 20, 3 \to 30, ...$. This index mapping means that the potential can be reordered using the new index and the following relation holds under an operaror $O$
$$ V_{eq}[i] \approx  V_{eq}[j]$$

For a set of operators we have 

$$ V_{eq}[i] \approx  \frac{1}{N} \sum_{n = 1} ^ NV_{eq}[j^{O_n}]$$

From the above we would be interested in how well symmetry can 

#### Displaced Potential
This is the potential for a displaced system. Similarily we would want to know the symmetrized potential. However we do not need to figure out which symmetries still hold under the displacement $(i, axis, direction)$ as this would already be determined and stored by the mapping of atoms in displaced system say $M(i, axis, direction)$. Note $M$ is mapping function and kicks out all $(j, axis, direction, operator)$. Therefore for a potential under displacement say $V_{1,x,+}$ with say $M(1,x, +) = \{(2, z, -, 3), (4, y, +, 2)\}$ 

$$ V_{1,x,+} \approx \frac{1}{2}[V_{2,z, -}^{O_3} + V_{4,y, +}^{O_2}]$$


<em>the right hand side of the above reads take the avergage of (the potential under displacement of atom 2 on z axis in - direction reordered by operator 3 + the potential under displacement of atom 4 on y axis in + direction reordered under operator 2) </em>

Note how an operator reorders the potential indices needs to only be done for the equilibrium potential of that system and can be used as a lookup.