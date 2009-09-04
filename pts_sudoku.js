#! /usr/bin/smjs
// by pts@fazekas.hu at Fri Sep  4 11:01:47 CEST 2009

// Instantiate variables in the equation (list of variables) eq. Returns
// null on contradiction, or the new insti.
function inst(eq, vars, insts, insti) {
  var eql = eq.length;
  var i, j, v, eqj;
  var nulli1 = null, nulli2 = null;
  for (i = 0; i < eql; ++i) {
    v = vars[eq[i]];
    if (v == true) {
      for (j = 0; j < eql; ++j) {
        // more than one ``true''
        if (vars[eq[j]] == true && i != j) return null;
      }
      for (j = 0; j < eql; ++j) {
        eqj = eq[j];
        v = vars[eq[j]];
        if (v == null) {
          vars[eqj] = false;
          insts[insti++] = eqj;
        }
      }
      return insti;
    } else if (v == null) {
      if (nulli1 == null) {
        nulli1 = i;
      } else {
        nulli2 = i;
      }
    }
  }
  if (nulli1 == null) return null;  // all ``false''
  if (nulli2 == null) {
    eqj = eq[nulli1];
    vars[eqj] = true;
    insts[insti++] = eqj;
  }
  return insti;
}

// The callback will be called with vars for each different solution.
function solve(eqs, vars, callback, insts, insti) {
  //print('S ' + vars);
  var eqsl = eqs.length;
  var varl = vars.length;
  var insti0 = insti;
  var i, insti2;
  while (true) {
    //print('Q ' + vars)
    for (i = 0; i < eqsl; ++i) {
      insti2 = inst(eqs[i], vars, insts, insti);
      //print('T' + i + insti2);
      if (insti2 == null) {
        i = -1;
        break;
      }
      insti = insti2;
    }
    if (i < 0) break;
    //print('I ' + vars);
    i = 0;  // TODO(pts): Pass this as well.
    while (i < varl && vars[i] != null) ++i;
    // !! varl can be less than vars.length
    if (i == varl) {  // all variables are known
      callback(vars);
      break;
    } else {
      vars[i] = false;
      insts[insti++] = i;
      // TODO(pts): Pass i as an argument.
      solve(eqs, vars, callback, insts, insti);
      vars[i] = true;
      // Continue the while loop with ``true'', just as with:
      // solve(eqs, vars, callback, insts, insti); break;
    }
  }
  while (insti > insti0) {
    vars[insts[--insti]] = null;
  }
}

var vars = [null, null, null, null];
var eqs  = [[1, 2, 3], [0, 2]];
var insts = Array(4);
solve(eqs, vars, print, insts, 0);

//var vars = [false, true, true, false];
//var eqs  = [[1, 2, 3], [0, 2]];
//var insts = Array(4);
//print(inst([1, 2, 3], vars, insts, 0));
