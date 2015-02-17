var tests = [
  require('./basic')
  ,require('./boundary-conditions')
  ,require('./parallel')
];

tests.forEach( function(test) {
  test();
});


console.log('Test ok');