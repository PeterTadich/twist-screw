//npm install https://github.com/PeterTadich/pseudo-inverse https://github.com/PeterTadich/matrix-computations

// mcts = twist-screw

// ECMAScript module

import * as hlao from 'matrix-computations';
import * as pinv from 'pseudo-inverse';

//examples ref:
//Modern Robotics, Chapter 3.3.2: Twists (Part 1 of 2)
//https://www.youtube.com/watch?v=mvGZtO_ruj0&list=PLggLP4f-rq02vX0OQQ5vrCxbJrzamYDfx&index=17

/*
//example 1
var w = [[0.0],[0.0],[1.0]]; //spatial rotational velocity
var v = [[0.0],[-2.0],[0.0]]; //spatial linear velocity
//var q = [[2.0],[0.0],[0.0]];
*/
/*
//example 2
var w = [[0.0],[0.0],[1.0]]; //spatial rotational velocity
var v = [[-2.0/Math.sqrt(2.0)],[2.0/Math.sqrt(2.0)],[0.0]]; //spatial linear velocity
//var q = [[-2.0/Math.sqrt(2.0)],[-2.0/Math.sqrt(2.0)],[0.0]];
*/
/*
//example 3
var w = [[0.0],[-1.0],[0.0]]; //spatial rotational velocity
var v = [[0.0],[0.0],[0.0]]; //spatial linear velocity
//var q = [[0.0],[0.0],[0.0]];
*/
/*
//example 3.23 ref: Lynch, Modern Robotics
var w = [[0.0],[0.0],[2.0]]; //spatial rotational velocity
var v = [[-2.0],[-4.0],[0.0]]; //spatial linear velocity
//var q = [[2.0],[-1.0],[0.0]];
*/

//fig 6.8 ref: Kajita, Humanoid Robotics
var w = [[1.0],[0.0],[0.0]]; //spatial rotational velocity
var v = [[0.3],[0.0],[1.0]]; //spatial linear velocity
//var q = [-0.0500,0.2500,0.0500]; //not confirmed

/*
//fig 6.9 ref: Kajita, Humanoid Robotics
var w = [[1.0],[0.0],[1.0]]; //spatial rotational velocity
var v = [[0.5],[0.1],[0.0]]; //spatial linear velocity
//var q = [[0.0],[0.0],[0.0]];
*/

function screwAxisFromTwist(w,v){
    var w_mag = magnitude(w); //Lynch where 'w' is s_hat*theta_dot. Murray does not assume 'w' is a unit rotation vector
    var w_u = hlao.matrix_multiplication_scalar(w,1.0/w_mag); //unit vector

    if(w_mag != 0.0){
        //ref: Murray page 45
        //Pitch
        var h = hlao.vector_dot(w,v)/(w_mag*w_mag); //Murray, dot(w,v)/(w_mag*w_mag) where 'w' is not a unit vector
        //Axis
        var l = hlao.matrix_arithmetic(hlao.matrix_multiplication_scalar(hlao.vector_cross(w,v),1.0/(w_mag*w_mag)),w,'+'); //where 'w' is not a unit vector
        //Magnitude
        var M = w_mag;
        
        //ref: Lynch
        var s_hat = w_u;
        var theta_dot = w_mag;
        var h = hlao.vector_dot(w_u,v)/theta_dot; //IMPORTANT: Lynch 'w_u' is a unit vector
        
        //solve for q
        var frag1 = pinv.pinv(hlao.matrix_multiplication_scalar(S(s_hat),theta_dot));
        var frag2 = hlao.matrix_arithmetic(v,hlao.matrix_multiplication_scalar(s_hat,h*theta_dot),'-');
        var q = hlao.matrix_multiplication(frag1,frag2);
    } else { //w === 0.0
        var v_mag = magnitude(v);
        var v_u = hlao.matrix_multiplication_scalar(v,1.0/v_mag); //unit vector
        
        //ref: Murray page 45
        //Pitch
        var h = NaN; //infinite pitch
        //Axis
        var l = hlao.matrix_arithmetic([[0.0],[0.0],[0.0]],v,'+'); //where 'v' is not a unit vector ('l' is a line going through the origin)
        //Magnitude
        var M = magnitude(v);
        
        //ref: Lynch
        var s_hat = v_u;
        var theta_dot = v_mag;
        var h = NaN; //Lynch = Murray
        
        //solve for q
        var frag1 = pinv.pinv(hlao.matrix_multiplication_scalar(S(s_hat),theta_dot));
        var frag2 = hlao.matrix_arithmetic(v,hlao.matrix_multiplication_scalar(s_hat,1.0*theta_dot),'-'); //h = 1
        var q = hlao.matrix_multiplication(frag1,frag2);
    }

    //screw
    console.log("screw S {q,s_hat,h}:");
    console.log('   - q: [' + q[0][0].toFixed(4) + ',' + q[1][0].toFixed(4) + ',' + q[2][0].toFixed(4) + ']^T');
    console.log('   - s_hat: [' + s_hat[0][0].toFixed(4) + ',' + s_hat[1][0].toFixed(4) + ',' + s_hat[2][0].toFixed(4) + ']^T');
    console.log('   - h: ' + h.toFixed(4));
    console.log("   where theta_dot: " + theta_dot.toFixed(4));
    console.log('   Murray axis "l" (a line): [' + l[0][0].toFixed(4) + ',' + l[1][0].toFixed(4) + ',' + l[2][0].toFixed(4) + ']^T');

    //twist
    var tv_1 = hlao.vector_cross(hlao.matrix_multiplication_scalar(s_hat,-1.0),q);
    var tv_2 = hlao.matrix_multiplication_scalar(s_hat,h);
    var tv = hlao.matrix_multiplication_scalar(hlao.matrix_arithmetic(tv_1,tv_2,'+'),theta_dot);
    var tw = s_hat;
    var twist = [[tw[0][0]],[tw[1][0]],[tw[2][0]],[tv[0][0]],[tv[1][0]],[tv[2][0]]];
    console.log("twist V = [w,v]^T:");
    console.log('   - [' + w[0][0].toFixed(4) + ',' + w[1][0].toFixed(4) + ',' + w[2][0].toFixed(4) + ',' + tv[0][0].toFixed(4) + ',' + tv[1][0].toFixed(4) + ',' + tv[2][0].toFixed(4) + ']^T');
    
    return([q,s_hat,h]);
}

function magnitude(a){
    return(Math.sqrt(a[0][0]*a[0][0] + a[1][0]*a[1][0] + a[2][0]*a[2][0]));
}

function S(a){
    var ax = a[0][0]; var ay = a[1][0]; var az = a[2][0];
    return(
        [
            [0.0, -az,  ay],
            [ az, 0.0, -ax],
            [-ay,  ax, 0.0]
        ]
    );
}

export {
    screwAxisFromTwist
};