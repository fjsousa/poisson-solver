var Poisson = require('../poisson');
var assert  = require('assert');

module.exports = function () {

  var opt = {
    w: 0.5,
    h: 0.5,
    n: 5, 
    m: 5 
  };

  var poisson = new Poisson(opt);

  N = [1,2,3,4,5];
  S = [7,8,9,10,11];
  E = [11,12,13,14,15];
  W = [17,18,19,20,21];

  var bcErr = poisson.setBoundaryConditions(N, S, E, W);

  assert(bcErr === null, bcErr);
  assert(poisson.u.old[1] === N[1], 'North Boundary set');
  assert(poisson.u.old[(poisson.n-1)*poisson.n + 1 ] === S[1], 'South Boundary set');
  assert(poisson.u.old[poisson.n] === W[1], 'West Boundary set');
  assert(poisson.u.old[poisson.n*2-1 ] === E[1], 'West Boundary set');

  // wrong bc array size
  N = S = E = W = [0,0,0,0];
  bcErr = poisson.setBoundaryConditions(N, S, E, W);
  assert(bcErr !== null, 'Error on boundary conditions is trigered by incorrect array size');	
};
