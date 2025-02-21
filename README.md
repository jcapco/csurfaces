# Toolkit for Cubic Surfaces

This is a toolkit for computing many things regarding cubic surfaces.

## Installation

Currently I am providing `Maple` codes as files with `mpl` extension (text-based). As things progress, I plan to include ports of the code to `Python` (with `Giacpy` or `Sympy`) and/or `Singular` and/or `Macaulay2`. For visualization I would also include some `Mathematica` codes. This is work in progress

## Usage

We say that a quaternary cubic form $F$ (or its zero set $X$) in the variables $w,x,y,z$ has a *normal form* if: 
 - $X$ has no lines parallel to the $x$, $y$ or $z$ axes.
 - $X$ has no lines at the plane of infinity $w=0$.

Throughout when we write *cubic surface*, we mean one that is defined over $\mathbb Q$, normal (i.e. finitely many singular points) and has only finitely many number of lines.

### *Cubic surfaces not in normal form*
For a cubic surface $X$ not in normal form, there is a very easy way to compute the lines on the surface. One views the cubic surface as an affine surface, parametrizes the lines with respect to one of the variables, in our code we use $x$ as the parameter. The other variables $w,y,z$ are linear equations in $x$ with coefficients that need to be computed. These coefficients are computed by substitution in the cubic equation (via $w=1$ dehomogenization of the cubic form) and solving for zeros of the coefficients of $1,x,x^2,x^3$. This final part can be done via Groebner basis and is encoded in the `find_gb_lines` procedure of the library.

Here is an example how this is done
  ```maple
  read("lib_cubic_surface.mpl"):
  F := randpoly([w,x,y,z],dense, degree=3):  
  gb := find_gb_lines(subs(w=1,F)):
  ```
In the example above since $F$ was chosen randomly, we would expect a smooth cubic surface not in normal form. The affine lines are parametrized via $(x,a+bx,c+dx)$ and the Groebner basis `gb` solves for $d$. So that `gb[1]` must be of univariate in $d$ of degree 27 and the remaining elements of `gb` are linear in $a,b,c$. In this way we get the 27 lines of the smooth cubic surface $X$. The lines can be solved, by solving for all the roots of $d$ in `gb[1]` and solving the unique values of $a,b,c$ for each of the linear equations in `gb[i]`. where $i>1$.

The above procedure also works for normal cubic surfaces (not necessarily smooth) as long as they are not in normal form and as long as it contains only finitely many lines!

### *Cubic surfaces in normal form*
In many applications, it is convenient to study cubic surfaces in normal form. However, we cannot compute the lines using the method described previously. 

Given is a normal cubic surface. One way to check if the previous procedure would work, is to try look at the polynomial `gb[1]`. If this polynomial is not of degree 27 or if it is not univariate in $d$, then we need to find another way to compute the number of lines. 

An easy way (not necessarily optimal) would be to projectively transform the cubic surface into a general position (i.e. make it into a non-normal form) and then use the above method. This is encoded in the procedure `generic_cubic_transform`

We illustrate this with a singular cubic in the Sylvester normal form :

---
**`normal_form`**
<a id="normal_form"></a>
  ```
  read("lib_cubic_surface.mpl"):
  #cubic surface singluar at (2:2:2:-3)
  F := w^3 + x^3 + y^3 + 4*z^3/9 - 4*(w+x+y+z)^3/9:
  M:=linalg:-randmatrix(4,4,entries=rand(-4..4)):
  linalg:-det(M);
  g,G,gb := generic_cubic_transform(subs(w=1,F), M):

  ```
---
We compute the determinant of $M$ to make sure that it is a non-singular matrix (thus a valid projective transformation). The new transformed cubic form is $G(w,x,y,z)$ and it is dehomogenized as $g(x,y,z)$. We could check for the degrees of polynomials in `gb`
  ```
  seq([degree(gb[i],a),degree(gb[i],b),degree(gb[i],c),degree(gb[i],d)],i=1..nops(gb)); 
  ```
Note that if `gb[1]` is not univariate and of degree 27 (in our case the roots come with multiplicities, as we have a singular cubic surface), then we were not successful in our transformation and we need to find a better $M$. In this case,  we could repeat this procedure with random $M$ but small random entries until we are successful (or increase the bits of the random entries, but this would mean that the computation would take longer).

### Number of Lines

Although the procedure discussed here will work for non-normal cubic forms, it is always recommended to use normal forms (such as those expressed by Sylvester normal forms). With low bit coefficients, the procedure for computing lines is much faster for these normal forms.

From `gb[1]`, we can get the number of lines contained in the cubic surface by executing `no_lines(gb)`. However, we can get more information than that. Provided that the cubic surface is such that the roots of the univariate polynomial in `gb[1]` are algebraic numbers in low bits. We can explicitly compute for the lines. This is the case for Sylvester normal form provided in [`normal_form`](#normal_form). In this particular case we can make use of the procedure `create_lines` that takes $M$ and $gb$ ( [`normal_form`](#normal_form)) see as input:

---
  **`compute_lines`**
  <a id="compute_lines"></a>  
  ```
  Lines, mlines := create_lines(gb,M):
  ```
---
You can give the identity matrix for $M$ if a non-normal form was used in the first place.

In our particular example, the outputs are 
  ```
  Lines := [[w, 0, -w, z], [0, -y, y, z], [w, -w, 0, z], [1/3*(-5+2*6^(1/2))*(4*z*6^(1/2)+3*y+12*z), 1/3*(-2+6^(1/2))*(-z*6^(1/2)+3*y), y, z], [w, x, -w+2*x, 1/4*6^(1/2)*(-x*6^(1/2)+2*w-2*x)], [-1/75*(16+6^(1/2))*(z*6^(1/2)+6*y), 1/46875*(131+16*6^(1/2))*(128*z*6^(1/2)+375*y-48*z), y, z], [w, x, -8/5*w-x, -1/20*6^(1/2)*(-6^(1/2)*w+16*w+20*x)], [w, -5/8*w-5/8*y, y, 1/2000*(-3+8*6^(1/2))*(-16*y*6^(1/2)+125*w-131*y)], [-1/60*6^(1/2)*(8*x*6^(1/2)-3*x+10*z), x, 1/60*6^(1/2)*(-8*x*6^(1/2)-3*x+10*z), z], [w, -y, y, -1/6*(3+2*6^(1/2))*w], [2/5*(3+2*6^(1/2))*z, x, -x, z], [w, -w, y, -1/6*(3+2*6^(1/2))*y], [-x, x, 2/5*(3+2*6^(1/2))*z, z], [1/3*(-2+6^(1/2))*(-z*6^(1/2)+3*x), x, 1/3*(-5+2*6^(1/2))*(4*z*6^(1/2)+3*x+12*z), z], [w, 2*w-y, y, -1/2*(6^(1/2)+3)*(-y*6^(1/2)+w+2*y)], [-y, 2/5*(3+2*6^(1/2))*z, y, z], [w, x, -w, -1/6*(3+2*6^(1/2))*x], [1/3*(-5+2*6^(1/2))*(4*z*6^(1/2)+3*x+12*z), x, 1/3*(-2+6^(1/2))*(-z*6^(1/2)+3*x), z], [w, -w+2*y, y, 1/4*6^(1/2)*(-y*6^(1/2)+2*w-2*y)], [-1/46875*(-131+16*6^(1/2))*(-128*z*6^(1/2)+375*x-48*z), x, 1/75*(-16+6^(1/2))*(-z*6^(1/2)+6*x), z], [w, -8/5*y-w, y, 1/20*6^(1/2)*(y*6^(1/2)+20*w+16*y)]]:
  mlines := [0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0]:
  ```
As can be seen `Lines` are the lines parametrized by two of the coordinates say $(w:z)$. In our example, the lines can be parametrized over a quadratic number field with low bits. Thus, the procedure `create_lines` can terminate in short amount of time. For more complicated cubic surfaces defined over higher bits, one has no other choice than to use floating point approximations to compute the lines and the procedure `no_lines` (which is exact and terminate almost immediately, if `gb` was already solved). The second output `mlines` is a list with entries 0 or 1 which describe multiplicity of the lines in `Lines`. If the line is non-simple and has multiplicity $>1$ then the entry is 1 otherwise it is 0. Clearly if the cubic surface is smooth all entries of `mlines` are 0 and there are 27 of them.

If we use 
  ```
  add(dlines[i],i=1..nops(dlines))
  ```
We get the number of multiple lines.

### Incidence of Lines

If we are able to compute `Lines` (see [`compute_lines`](#compute_lines)), then we can use the output to compute the incidence relations. This terminates if the lines are defined over algebraic numbers of low degree and low bits. The procedure that allows this is called `compute_incidence` that takes `Lines` as input. 
  ```
  incs, poses := compute_incidence(Lines):
  ```
The output `incs` is a list of lists of integer indices, such that any element in $j$ in `inc[i]` satisfy $j>i$ and indicates that `Lines[j]` intersects with `Lines[i]`. The output `poses` is a list of lists of points in $\mathbb P^3$ (represented as 4-tuples of algebraic numbers). Such that the $k$-th point in `poses[i]` is the point of intersection of `Lines[incs[k]]` and `Lines[i]`. This therefore encodes incidence relation of all the lines.
  
The Eckardt points are those points on the cubic surface for which 3 or more lines are concurrent. If the surface is smooth, then 3 of the lines can only be concurrent and they are in fact coplanar. For singular cubics we can have more such lines and the concurrent lines need not be coplanar. 

## Acknowledgement

Many thanks to the Research Institute for Symbolic Computation, Johannes Kepler University, Linz, Austria.

## Citing

Once there is sufficient content I would give instruction on how to cite. Currently, you can use the code as is and just name me or this site if you insist on citing.

