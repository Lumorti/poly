var x0 >= -0.57735 <= 0.57735;
var x1 >= -0.57735 <= 0.57735;
var x2 >= -0.57735 <= 0.57735;
var x3 >= -0.57735 <= 0.57735;
var x4 >= -0.57735 <= 0.57735;
var x5 >= -0.57735 <= 0.57735;
var x6 >= -0.57735 <= 0.57735;
var x7 >= -0.57735 <= 0.57735;
var x8 >= -0.57735 <= 0.57735;
var x9 >= -0.57735 <= 0.57735;
minimize dummy_obj: 0;
subject to eqn0: -0.333333+x9*x9+x8*x8 = 0;
subject to eqn1: x4*x6+x0*x2+x1*x3+0.333333+x5*x7+-1.000000*x8 = 0;
subject to eqn2: x0*x6+-1.000000*x2*x4+-1.000000*x3*x5+x1*x7+-1.000000*x9 = 0;
subject to eqn3: x4*x4+-0.333333+x0*x0 = 0;
subject to eqn4: -0.333333+x1*x1+x5*x5 = 0;
subject to eqn5: 0.666667*x0*x1+0.333333*x0*x0+-0.222222+0.333333*x1*x1+0.666667*x4*x5+0.384900*x0+0.333333*x4*x4+0.384900*x1+0.333333*x5*x5 = 0;
subject to eqn6: -0.333333+x6*x6+x2*x2 = 0;
subject to eqn7: -0.333333+x7*x7+x3*x3 = 0;
subject to eqn8: 0.333333*x3*x3+0.666667*x2*x3+0.666667*x6*x7+-0.222222+0.333333*x2*x2+0.333333*x7*x7+0.384900*x3+0.384900*x2+0.333333*x6*x6 = 0;
subject to pos0: 4.000000*x4+-8.000000*x2+-2.000000*x3+x5+-1.000000*x7+8.000000*x0+-4.000000*x6+2.000000*x1 >= 0;
subject to pos1: 8.000000*x0+-4.000000*x5+-2.000000*x3+4.000000*x4+2.000000*x2+-1.000000*x7+-8.000000*x1+x6 >= 0;
subject to pos2: x7 >= 0;
solve;

