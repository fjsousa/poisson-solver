var Poisson = require('../poisson');
var assert  = require('assert');

module.exports = function () {

  var getRow = function (row) {

    var array = new Array(poisson.bm);
    for (var i = 0; i < array.length; i++) {
      array[i] = poisson.u.analytical[row*poisson.bm + i];
    }
    return array;
  };

  var getCol = function (col) {

    var array = new Array(poisson.bn);
    for (var i = 0; i < array.length; i++) {
      array[i] = poisson.u.analytical[i*poisson.bn + col];
    }
    return array;
  };

  var opt = {
    w: 1,
    h: 1,
    n: 5, 
    m: 5,
    bRow: 0,
    bCol: 2,
    bRows: 3,
    bCols: 3 
  };

  var poisson = new Poisson(opt);

  var d_x = opt.w / (opt.bCols*opt.m - 1);
  var d_y = opt.h / (opt.bRows*opt.n - 1); 
  var bm = opt.m + 1;
  var bw = d_x * bm;

  assert(poisson.d_x === d_x, 'Property d_x.');
  assert(poisson.d_y === d_y, 'Property d_y.');
  assert(poisson.bw === bw, 'Property bw.');
  assert(poisson.bm === bm, 'Property bm.');

  poisson.analytical();

  var nw = u(d_x*(opt.m * opt.bCol - 1), 0);
  var sw = u(d_x*(opt.m * opt.bCol - 1), d_y*opt.n);

  assert(nw === poisson.u.analytical[0], 'NW analytical ok.');
  assert(sw === poisson.u.analytical[opt.n*poisson.bm], 'SW analytical ok.');

  var errb = poisson.setBoundaryConditions(getRow(0), getRow(opt.n), getCol(opt.m), getCol(0));
  assert(!errb, 'Boundaries not set ');

  if (!errb) {

    var maxItterations = 100000;
    var maxResidue = 1E-10 ;
    var output = poisson.solver( maxItterations, maxResidue);

    assert(!isNaN(output.residue), 'Residue !== NaN');
    assert(output.residue !== undefined, 'Residue !== undefined');
    assert(output.residue !== null, 'Residue !== undefined');
    assert(!isNaN(output.iterations), 'Iterations !== NaN');
    assert(output.iterations !== undefined, 'Iterations !== undefined');
    assert(output.iterations !== null, 'Iterations !== undefined');

    assert(output.iterations >= 1, 'Iterations larger than one');


    if (output.residue <= maxResidue) {

      var res = 0;
      for (var i = 0; i < poisson.bn*poisson.bm; i++) {
        res += Math.abs(poisson.u.new[i]-poisson.u.analytical[i]);
      }

      res = res/poisson.bn/poisson.bm;

      assert(res < 1E-4, 'Numerical equals analytical.');

    } 
  }

  function u(x,y) {
    return Math.cos(Math.PI*(0.5-x/poisson.w)) * Math.sin(Math.PI*y/poisson.h);
  }

};