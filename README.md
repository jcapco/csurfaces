# Toolkit for Cubic Surfaces

This is a toolkit for computing many things regarding cubic surfaces.

## Installation

Currently I am providing `Maple` codes in `mpl` (text-based). This is work in progress

## Usage

We say that a quaternary cubic form $F$ (or its zero set $X$) in the variables $w,x,y,z$ has a *normal form* if: 
 - $X$ has no lines parallel to the $x$, $y$ or $z$ axes.
 - $X$ has no lines at the plane of infinity $w=0$.
Throughout when we write *cubic surface*, we mean one that is normal (i.e. finitely many singular points) and has only finitely many number of lines.

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
  ```maple
  read("lib_cubic_surface.mpl"):
  #cubic surface singluar at (2:2:2:-3)
  F := w^3 + x^3 + y^3 + 4*z^3/9 - 4*(w+x+y+z)^3/9:
  M:=linalg:-randmatrix(4,4,entries=rand(-4..4)):
  linalg:-det(M);
  g,G,gb := generic_cubic_transform(subs(w=1,F), M):

  ```
We compute the determinant of $M$ to make sure that it is a non-singular matrix (thus a valid projective transformation). The new transformed cubic form is $G(w,x,y,z)$ and it is dehomogenized as $g(x,y,z)$. We could check for the degrees of polynomials in `gb`
  ```
  seq([degree(gb[i],a),degree(gb[i],b),degree(gb[i],c),degree(gb[i],d)],i=1..nops(gb)); 
  ```
Note that if `gb[1]` is not univariate and of degree 27 (in our case the roots come with multiplicities, as we have a singular cubic form), then we were not successful in our transformation and we need to find a better $M$. In this case,  we could repeat this procedure with random $M$ but small random entries until we are successful (or increase the bits of the random entries, but this would mean that the computation would take longer).

## Acknowledgement

Many thanks to the Research Institute for Symbolic Computation, Johannes Kepler University, Linz, Austria.

## Citing

Once there is sufficient content I will give instruction on how to cite. Currently, you can use the code as is and just name me if you insist on citing.

