var Poisson = require('.poisson');

//Set conditions
var conditions = {
  w: 1,		//width or x dimension
  h: 1,		//height or y dimension
  n: 1000, 	//width or x dimension
  m: 1000
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
poisson.solver( maxItterations, maxResidue);

//Print 
poisson.print('./field.txt', poisson.u.old);

poisson.analitical();
poisson.print('./analitical.txt', poisson.u.analitical);







