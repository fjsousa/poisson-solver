var Poisson = require('../poisson');
var assert  = require('assert');

module.exports = function () {

  var opt = {
    w: 1,
    h: 1,
    n: 4, //y # of nodes
    m: 5//x # of nodes
  };

  var poisson = new Poisson(opt);
  
  var d_y = opt.h / (opt.n - 1);
  
  assert(poisson.d_y === d_y, 'Parameter d_y not set.');
  assert(poisson.m === opt.m, 'Parameter m not set.');
  assert(poisson.h === opt.h, 'Condition h not set.');
  assert(poisson.n === opt.n, 'Condition n not set.');

  var N = [0,0,0,0,0];
  var S = [0,0,0,0,0];
  var E = [0,0,0,0];
  var W = [0,0,0,0];

  var err = poisson.setBoundaryConditions(N, S, E, W);

  assert(!err, 'Boundaries not set!');

  var maxItterations = 100000;
  var maxResidue = 1E-9;
  var output = poisson.solver( maxItterations, maxResidue);

  assert(!isNaN(output.residue), 'Residue !== NaN');
  assert(output.residue !== undefined, 'Residue !== undefined');
  assert(output.residue !== null, 'Residue !== undefined');
  assert(!isNaN(output.iterations), 'Iterations !== NaN');
  assert(output.iterations !== undefined, 'Iterations !== undefined');
  assert(output.iterations !== null, 'Iterations !== undefined');

  assert(output.iterations >= 1, 'Iterations larger than one');


  if (output.residue <= maxResidue) {
    poisson.analytical();

    var res = 0;
    for (var i = 0; i < opt.n*opt.m; i++) {
      res += Math.abs(poisson.u.new[i] - poisson.u.analytical[i]);
    }

    res = res/poisson.n/poisson.m;

    assert(res < 1E-1,'Numerical equals analytical.'); 

  }
};