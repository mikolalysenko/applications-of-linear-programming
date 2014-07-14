Linear programming computational geometry
=========================================

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

### Variations of linear programming

It will also be useful to introduce a few closely related variations of linear programming. In general, these formulations change the type of objective function by replacing it with a polynomial of greater or lesser degree. In particular, the two special cases we will care about are feasibility (which replaces the objective with a constant function) and quadratic programming (which replaces the objective with a degree 2 polynomial).

#### Feasibility

The first of these problem, or feasibility, is to test if a system of linear equations is solvable - that is if the feasible set is empty or not. Feasibility is strictly a subset of linear programming, and can be solved with slightly less work. Feasibility is important in problems related to intersections of convex polytopes, like testing if a ray/triangle intersect for example.

#### Quadratic programming

Quadratic programming generalizes linear programming by replacing the linear objective function with a general positive semidefinite quadratic function.


### Duality

## Solving problems with linear programming

### Closest point and intersection problems

### Minimal enclosing sphere


## Algorithms

### Simplex method

### Seidel's algorithm

### Clarkson's algorithm

