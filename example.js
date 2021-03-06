var Poisson = require('.');

//Set conditions
var conditions = {
  w: 1,	 
  h: 1,   
  n: 50, 	
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

//Set boundary conditions to zero
poisson.setBoundaryConditions(N, S, E, W);

//Maximum number of iterations allowed
var maxItterations = 100000000;

//Maximum residue allowed
var maxResidue = 1E-9 ;

//Run solver 
var output = poisson.solver( maxItterations, maxResidue);

console.log('Solver converged with', output.iterations, 'iterations and', output.residue, 'residue.');
//Print 
poisson.print('./field.txt', poisson.u.old);

poisson.analytical();
poisson.print('./analytical.txt', poisson.u.analytical);








