# poisson-solver
Node Implementation of a finite diference poisson solver

##Usage

```javascript
var Poisson = require('.');

//Set conditions
var conditions = {
  w: 1,		//width or x dimension
  h: 1,		//height or y dimension
  n: 50, 	//width or x dimension
  m: 50
};

//Initialize solver
var poisson = new Poisson(conditions);

var N = [];
var S = [];
var E = [];
var W = [];

for (var i = 0; i < conditions.n; i++) {
  N[i] = S[i] = E[i] = W[i] = 0;
}

//Set Dirichelt boundary conditions to zero
poisson.setBoundaryConditions(N, S, E, W);

//Maximum number of iterations allowed
var maxItterations = 100000000;

//Maximum residue allowed
var maxResidue = 1E-9 ;

//Run solver 
var output = poisson.solver( maxItterations, maxResidue);

console.log('Solver converged with', output.iterations, 'iterations and', output.residue, 'residue.');
//Print solution
poisson.print('./field.txt', poisson.u.old);

//Print analitical solution
poisson.analitical();
poisson.print('./analitical.txt', poisson.u.analitical);

```

##Browser version

If you have [browserify installed]() you can create a browser version of the poisson solver. Just run

`npm run browserify`

and then include a script tag on your html file.

```html
<script type='text/javascript', src='http://latex.codecogs.com/latexit.js'></script>
```

Now you can require the solver in the browser like you would do with node

```javascript
var Poisson = require('poisson');
```