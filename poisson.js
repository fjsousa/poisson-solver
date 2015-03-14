var fs = require('fs');

var Poisson = function (opt) {

  this.h = this.bh = opt.h;
  this.w = this.bw = opt.w;
  this.n = this.bn = opt.n; //y # of nodes
  this.m = this.bm = opt.m; //x # of nodes

  this.d_x = this.bw / (opt.m - 1);
  this.d_y = this.bh / (opt.n - 1);

  this.bRow = opt.bRow || 0;
  this.bCol = opt.bCol || 0;
  this.bRows = opt.bRows || 0;
  this.bCols = opt.bCols || 0;
  
  if (!!this.bCols && !!this.bRows) {

    this.d_x = opt.w / (opt.bCols*opt.m -1);
    this.d_y = opt.h / (opt.bRows*opt.n -1);

    if (isBoundaryBlock(opt.bCol, opt.bCols)) {
      this.bm = this.m + 1; 
      this.bw = this.d_x * this.bm;
    } else {
      this.mCore = this.m + 2; 
      this.bw = this.d_x * (this.mCore + 1);
    }

    if (isBoundaryBlock(opt.bRow, opt.bRows)) {
      this.bn = this.n + 1; 
      this.bh = this.d_y * this.bn;
    } else {
      this.nCore = this.n + 2; 
      this.bh = this.d_y * (this.nCore + 1);
    }

  }

  this.AP = new Array(this.bn*this.bm);
  this.AE = new Array(this.bn*this.bm);
  this.AW = new Array(this.bn*this.bm);
  this.AN = new Array(this.bn*this.bm);
  this.AS = new Array(this.bn*this.bm);

  for (var i = 0; i < this.AP.length; i++ ) {
    this.AP[i] = -2/(this.d_x*this.d_x) -2/(this.d_y*this.d_y);
    this.AE[i] = 1/(this.d_x*this.d_x);
    this.AW[i] = 1/(this.d_x*this.d_x);
    this.AN[i] = 1/(this.d_y*this.d_y);
    this.AS[i] = 1/(this.d_y*this.d_y);
  }

  this.u = {
    'old': new Array(this.bn*this.bm),
    'new': new Array(this.bn*this.bm)
  };

  this.initialize();

  function isBoundaryBlock(b, bTotal ) {
    return b === 0 || b === bTotal - 1;
  }
};

//Needs to return used iterations
//returns null if solver doesn't converge
Poisson.prototype.solver = function (maxItt, maxRes) {

  var itt = 0;
  var res = 1000000;
  while ( res > maxRes && itt < maxItt ) {

    this.copy();

    for (var i = 1; i < this.bn - 1; i++) { //y
      for (var j = 1; j < this.bm - 1; j++) { //x

        var idx = i*(this.bm) + j;

        var sum = this.AW[idx]*this.u.new[i * this.bm + j - 1] +
         this.AE[idx]*this.u.new[i * this.bm + j + 1] +
         this.AS[idx]*this.u.new[(i + 1) * this.bm + j] +
         this.AN[idx]*this.u.new[(i - 1) * this.bm + j];

        var x = this.xx(j);
        var y = this.yy(i);

        var pi = Math.PI;
        var b = -pi*pi/this.h/this.h*Math.sin(pi*y/this.h)*Math.cos(pi*(0.5-x/this.w)) -
          Math.sin(pi*y/this.h)*(pi*pi/this.w/this.w)*Math.cos(pi*(0.5-x/this.w));

        this.u.new[idx] = (-sum + b)/this.AP[idx];

      }
    }

    res = this.residue();
    if (++itt%100 === 0) {
      // console.log('Solver: iteration #', itt, 'with residue of', res);
    }

    // this.swap();

  }

  return {residue: res, iterations: itt};
};

Poisson.prototype.setBoundaryConditions = function(N, S, E, W) {

  if ( N.length !== this.bm || S.length !== this.bm || E.length !== this.bn || W.length !== this.bn ) {
    return 'Boundary conditions size don\'t match';
  }

  for (var col = 0; col < N.length; col ++) {
    this.u.new[col] = this.u.old[col] = N[col];
    this.u.new[(this.bn -1) * this.bm + col] = this.u.old[(this.bn -1) * this.bm + col] = S[col];
  }

  for (var row = 0; row < E.length; row ++ ){
    this.u.new[row*this.bm] = this.u.old[row*this.bm] = W[row];
    this.u.new[row*this.bm + this.bm - 1] = this.u.old[row*this.bm + this.bm - 1] = E[row];

  }

  return null;

};

Poisson.prototype.residue = function () {

  res = 0;
  for (var i = 0; i < this.bn*this.bm; i++) {
    res += Math.abs(this.u.new[i]-this.u.old[i]);
  }
  return res/this.bm/this.bn;

};

Poisson.prototype.print = function (filename, dataSet) {

  if (fs.existsSync(filename)) {
    fs.unlinkSync(filename);
}
  for (var i = 0; i < this.bn; i++) { //y
    for (var j = 0; j < this.bm; j++) { //x
      fs.appendFileSync(filename, '' + dataSet[i*this.bm + j] + ' ');
    }
    fs.appendFileSync(filename, '\n');
  }
};

Poisson.prototype.copy = function () {

  for (var i = 0; i < this.u.new.length; i++) {
    this.u.old[i] = this.u.new[i];
  }

};

Poisson.prototype.swap = function () {

  var temp = this.u.old;
  this.u.old = this.u.new;
  this.u.new = temp;

};

Poisson.prototype.initialize = function () {

  for (var i = 0; i < this.bn*this.bm; i++) {
    this.u.new[i] = this.u.old[i] = 0;
  }

};

Poisson.prototype.xx = function (j) {
  if (this.bCol === 0) {
    return j * this.d_x;
  } else {
    return (this.bCol * this.m  - 1 + j)*this.d_x;  
  }
};

Poisson.prototype.yy = function (i) {
  if (this.bRow === 0) {
    return  i * this.d_y;
  } else {
    return (this.bRow * this.n  - 1 + i)*this.d_y;
  }
};

Poisson.prototype.analytical = function () {

  this.u.analytical = new Array(this.bm*this.bn);
  for (var i = 0; i < this.bn; i++) { //y
    for (var j = 0; j < this.bm; j++) { //x  

      var x = this.xx(j);
      var y = this.yy(i);  

      this.u.analytical[i*this.bm + j] = Math.cos(Math.PI*(0.5-x/this.w)) * Math.sin(Math.PI*y/this.h);
    }
  }
};
module.exports = Poisson;