Linear programming in computational geometry
============================================

This is adapted from the lecture notes for CS558, Computational Geometry which was taught at the University of Wisconsin-Madison in Fall 2013.

## What is a linear program?

Contrary to what the name might suggest, [linear programming](http://en.wikipedia.org/wiki/Linear_programming) isn't really about "programming" per-se, at least in the sense which we would usually think of it.  Rather, like [dynamic programming](http://en.wikipedia.org/wiki/Dynamic_programming), the term programming here is archaic and was originally meant to connotate connections with optimization theory.

Linear programming is really about optimizing linear functions which are defined on convex polytopes. Specifically, a linear program is a system of (closed) linear inequalities together with some linear function called the *objective* which we seek to minimize (or equivalently maximize).  For example, the following is an instance of a linear program:

```
Maximize:    2 x + 3 y

Such that:

            x + y <= 1
          2 x - y >= 2
```

Linear programming is ubiquitious with many applications in signal processing, machine learning, operations research, economics and industrial engineering.  It is in some sense a fundamental algorithmic building block, generalizing common algorithms like network, in addition to being the basis for many approaches to approximately solving NP-hard optimization problems.

In the context of computational geometry, linear programming serves many purposes like finding intersections between convex sets or constructing minimal bounding volumes. However, one thing that is somewhat more specialized in the use of linear programming in computational geometry is that the dimension of the linear spaces under consideration is usually quite small (say 2 or 3) compared to the number of constraints. This can be exploited in various ways through algorithmic tradeoffs that push more of the solving complexity onto the dimension of the space.  For example, we will later discuss an algorithm which has the expected complexity of `O(d! n)`, which in high dimensions is completely impractical.  Yet for low dimensional problems as commonly arise in computational geometry, this cost would be completely reasonable and holding `d` constant would even be optimal!

### Geometric intuition

But before going into more detail on how to solve a linear program, we must first say more about what they are and how to translate geometric problems into linear programs. In doing this, it will be useful to build up a geometric intuition of what a linear program actually means.

One way to think about a linear program is that each inequality, or *constraint* cuts out some linear halfspace.  The intersection of these halfspaces, or the *feasible region*, is the set of all possible vectors which satisfies the inequalities. It is always a convex polytope, though it may not be bounded.

(Insert picture of convex polytope + inequalities)

The objective function can be thought of as a vector field which is defined at every point in space. Physically, we can imagine that this vector field "pushes" each point in the convex set along the direction determined by the slope of the objective function. Following this flow field while staying within the feasible region, one eventually will arrive at one of the vertices the polytope (if it exists that is). This vertex is the solution of the linear program.

(Insert picture of flow field)

Each vertex of the feasible region is formed by the intersection of exactly d+1 halfspaces, which suggests a naive way to solve a linear program by just iterating over all `d+1` tuples of constraints. This basic idea turns out to be a useful starting point for many approaches to linear programming as we will see shortly.


### Matrix notation

To simplify writing linear programs, it is helpful to introduce matrix notation. The main idea is that we rewrite the system of constraints using the following new variables:

* `A` which is an `n`-by-`d` matrix of constraints
* `b` which is an `n` dimensional vector representing the right hand side of the constraints
* `c` which is a `d` dimensional vector encoding the objective function
* `x` which is a `d` dimensional vector encoding the solution

Then the linear program can be concisely written as follows:

```
Minimize:   x^T c

Such that:
          A x >= b
```

Where the `>=` here is applied component-wise.

### Variations of linear programming

It will also be useful to introduce a few closely related variations of linear programming. In general, these formulations change the type of objective function by replacing it with a polynomial of greater or lesser degree. In particular, the two special cases we will care about are feasibility (which replaces the objective with a constant function) and quadratic programming (which replaces the objective with a degree 2 polynomial).

#### Feasibility

The first of these problem, or feasibility, is to test if a system of linear equations is solvable - that is if the feasible set is empty or not. Feasibility is strictly a subset of linear programming, and can be solved with slightly less work. Feasibility is important in problems related to intersections of convex polytopes, like testing if a ray/triangle intersect for example.

#### Quadratic programming

Quadratic programming generalizes linear programming by replacing the linear objective function with a general positive semidefinite quadratic function.  That is given a `d`-by-`d` matrix `Q` where for any test vector `x`, we require that:

```
x^T Q x >= 0
```

This requirement is essential for quadratic programming to be efficiently solvable.  Specifically if `Q` has any negative eigenvalues, then finding an optimal solution is NP-hard. Given such a matrix `Q`, a quadratic program can be written as follows:

```
Minimize:
      x^T Q x + x^T c

Such that:
        A x >= b
```

It is easy to see that quadratic programming generalizes linear programming, since we can just take `Q = 0` and recover a linear program. Visually, the main difference in a quadratic program is that the flow field generated by the quadratic program is now some sequence of ellipses flowing toward a single point.  

(Insert picture of quadratic program flow field)

Again, it is intuitive that these fields will push any point in the feasible region towards a unique point, however unlike in linear programs this point is likely to be a vertex. In the context of computational geometry, it is not hard to adapt many algorithms to solve quadratic programs, and so we shall switch between these systems freely.


### KKT conditions

In optimization problems with equality constraints, Lagrange multipliers can be used to convert a constrained optimization problem into an unconstrained problem. A similar technique can be applied to optimization problems with inequality constraints, via a result known as the KKT theorem (short for Karush-Kuhn-Tucker).We will now state this result in the context of *quadratic programming* since it generalizes the case of linear programming/feasibility:

#### KKT Theorem

A feasible solution `x'` to a quadratic program,

```
Minimize:
      x^T Q x + x^T c

Such that:
        A x = b
         0 <= x
```

is optimal if and only if there exists an n-dimensional vector `l` and a d-dimensional vector `u >= 0` such that:

```
c^T + 2 x'^T Q = -l^T A + u^T
        u^T x' = 0
```

## Closest point and intersection problems

One of the most frequently encountered low level tasks in computational geometry is to detect if two convex polytopes intersect, or if they are separated then to compute the clearance distance between them.

These sorts of problems are ubiquitous in physics simulations (where it is usually lumped in with "collision detection"), computer graphics (for ray casting and clipping), and GIS systems. Because they are used so frequently, a lot of effort has been spent on optimizing these problems for special cases, like in this [table for instance](http://www.realtimerendering.com/intersections.html). Older texts like [Graphics Gems](http://tog.acm.org/resources/GraphicsGems/) have whole chapters devoted to solving special instances of these intersection problems. Taken independently, these techniques can seem a bit ad-hoc, however at a deeper level they are all closely related -- specifically they are all instances of quadratic programming.

### Affine sets

To understand how this works, we first study a version of this problem for [affine subspaces](http://en.wikipedia.org/wiki/Affine_space#Affine_subspaces), like lines points, planes and so on. In general there are two ways that an affine set can be defined:

* Primal form: As an affine combination of points

```
S = { p_0 + a_1 p_1 + .... a_n p_n : a_i is a real number }
```

Where `p_0` is a point and each `p_i` is a vector, and each `a_i` is some affine coefficient

* Dual form: As an intersection of hyperplanes

```
    {     p_1^T (x - p_0) = 0   }       
S = { x : p_2^T (x - p_0) = 0   }
    {          ...              }
    {     p_n^T (x - p_0) = 0   }
```

Again `p_0` is a point and each `p_i` is a vector

#### Primal formulation

In the primal formulation, an affine set is written as an [affine combination](http://en.wikipedia.org/wiki/Affine_combination) of points, where we take one of the points, p_0, to be fixed.  Given a pair of affine spaces:

```
S = p_0 + a_1 p_1 + ... + a_n p_n
T = q_0 + b_1 q_1 + ... + b_m q_m
```

The squared distance between them is defined as follows:

```
d(S,T)^2 = min_{x in S, y in T}  |x - y|^2
```

Algebraically, this is equivalent to solving the following linear least squares problem:

```
min_{a_i,b_j} (p_0 + a_1 p_1 + ... + a_n p_n - q_0 - b_1 q_1 + ... b_m q_m)^2
```

To simplify notation, let us define the auxiliary matrices:

```
     ( p_1^T )
P =  ( p_2^T )
     (  ...  )
     ( p_n^T )


     ( q_1^T )
Q =  ( q_2^T )
     (  ...  )
     ( q_m^T )

M  = ( P P^T   P Q^T )
     ( Q P^T   Q Q^T )
```

And the vectors:


```
    ( a_1 )
a = ( a_2 )
    ( ... )
    ( a_n )


    ( b_1 )
b = ( b_2 )
    ( ... )
    ( b_m )

x = ( a )
    ( b )

c = ( P (p_0 - q_0) )
    ( Q (p_0 - q_0) )
```

Then the above minimization problem can be rewritten as:

```
min   x^T M x - x^T c
 x
```

Which is a linear least squares problem, and so the solution to `x` is:

```
M x = c
```

The interpretation of `x` is that it is the coordinates of the pair of closest points in `S` and `T` in barycentric coordinates.  If this solution is not unique, then the above equation will be rank deficient but still solvable, and the solution space will give the collection of pairs of closest points.  The two sets will intersect if and only if:

```
P^T a + p_0 = Q^T b + q_0
```

#### Dual formulation

The story for the dual formulation is quite similar, though instead of solving for the coordinates of the closest pair of points in barycentric coordinates we will just directly compute the closest pair. In this version the two sets S and T will be determined by the *intersection* of a collection of hyperplanes:


```
    {     p_1^T (x - p_0) = 0   }       
S = { x : p_2^T (x - p_0) = 0   }
    {          ...              }
    {     p_n^T (x - p_0) = 0   }


    {     q_1^T (y - q_0) = 0   }       
T = { y : q_2^T (y - q_0) = 0   }
    {          ...              }
    {     q_m^T (y - q_0) = 0   }
```

The closest pair of points is the solution to the following equation:

```
  min   |x - y|^2
x in S
y in T
```

Again, we will reintroduce the following matrices as we did in the primal case:

```
     ( p_1^T )
P =  ( p_2^T )
     (  ...  )
     ( p_n^T )


     ( q_1^T )
Q =  ( q_2^T )
     (  ...  )
     ( q_m^T )
```

With this notation, we can define the vectors:

```
d = P p_0
e = Q q_0
```

So the constraints on `x` and `y` can be written as:

```
P x - d = 0
Q y - e = 0
```

Introducing Lagrange multipliers `s` and `u` for these constraints, the system can be written equivalently as:

```
min s^T (P x - d) + u^T (Q y - e) + (x - y)^T (x - y)
```

This gives the following system of constraints:

```
s^T P - u^T Q = 0
      P x - d = 0
      Q y - e = 0
```

### Convex sets

The situation for convex sets is roughly analogous.  Again, the formulation of the problem is highly dependent on whether the polytopes are represented as intersections or combinations.  Specifically, for convex polytopes there are two basic ways they are encoded:

* V-Polytopes: As a convex combination of vertices
* H-Polytopes:  As an intersection of halfspaces

Again these two representations are dual to one another, though the precise way the problem is formulated in each case is subtly different.

#### Vertex formulation

A V-polytope representation of a convex polytope is given by enumerating the vertices of the convex set.  Specifically, the polytope `S` is determined by a convex combination of vertices:

```
S = { a_0 p_0 + a_1 p_1 + ... + a_n p_n }
```

Where:

```
      0 <= a_i <= 1
a_0 + a_1 + ... + a_n = 1
```

Following the conventions from the affine/primal case, we can again define `T` as:

```
T = { b_0 q_0 + b_1 q_1 + ... + b_m q_m }
```

So the distance between `S` and `T` is:

```
min  (a_0 p_0 + a_1 p_1 + ... + a_n p_n - b_0 q_0 - ... - b_m q_m)^2

s.t.
      0 <= a_i
      0 <= b_j
      a_0 + ... + a_n = 1
      b_0 + ... + b_m = 1
```

Which is a quadratic program.

#### Halfspace formulation

```
min   (x - y)^2

s.t.
      P (x - p_0) = 0
      Q (y - q_0) = 0
```

## Minimal enclosing ball

Another common geometry problem is to find the smallest enclosing ball containing a collection of points.  That is, given a point set, `S = { p_1, p_2, ..., p_n }` in Euclidean space, we want to find a circle with center `p` and smallest radius `r` such that for each `p_i` in S,

```
(p_i - p)^2 <= r^2
```

Or equivalently, 

```
min  max   (p - p_k)^2
 p   p_k
```

Now observe that the center of the minimal enclosing ball must reamin within the convex hull of the point set, so without loss of generality:

```
p = a_1 p_1 + ... + a_n p_n

        0 <= a_i
a_1 + a_2 + ... + a_n = 1
```

Define the matrix/vector pair:

```
     ( p_1^T )
P =  ( p_2^T )
     (  ...  )
     ( p_n^T )

    ( a_1 )
a = ( a_2 )
    ( ... )
    ( a_n )
```

So,

```
p = P^T a
```

Then the following quadratic program can be used to solve the minimal enclosing ball problem:

```
Minimize:   a^T P P^T a - sum_i p_i^T p_i a_i

Such that:

        0 <= a_i
        a_1 + a_2 + ... + a_n = 1
```

## Algorithms

### Simplex method

### Seidel's algorithm

### Clarkson's algorithm

### Warm restarting
