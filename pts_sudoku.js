#! /usr/bin/smjs
// by pts@fazekas.hu at Fri Sep  4 11:01:47 CEST 2009

// Instantiate variables in the equation (list of variables) eq. Returns
// a bool indicating if there is a solution.
function inst(eq, vars) {
  var eql = eq.length;
  var i, j, v;
  var nulli1 = null, nulli2 = null;
  for (i = 0; i < eql; ++i) {
    v = vars[eq[i]];
    if (v == true) {
      for (j = 0; j < eql; ++j) {
        v = vars[eq[j]];
        if (v == true) {
          if (j != i) return false;  // more than one ``true''
        } else if (v == null) {
          vars[eq[j]] = false;
        }
      }
      return true;
    } else if (v == null) {
      if (nulli1 == null) {
        nulli1 = i;
      } else {
        nulli2 = i;
      }
    }
  }
  if (nulli1 == null) return false;  // all ``false''
  if (nulli2 == null) vars[eq[nulli1]] = true;
  return true;
}

// Copies an array of integers.
//function copyIntArray(ary) {
//  return ary.length ? ary.join(',').split(',') : [];
//}

// The callback will be called with vars for each different solution.
function solve(eqs, vars, callback) {
  print('S' + vars);
  var eqsl = eqs.length;
  var varl = vars.length;
  var i;
  for (i = 0; i < eqsl; ++i) {
    if (!inst(eqs[i], vars)) return;
  }
  print('I' + vars);
  var nulli = 0;
  while (nulli < varl && vars[nulli] != null) ++nulli;
  if (nulli == varl) {  // all variables are known
    callback(vars);
  } else {
    vars[nulli] = false;
    // TODO(pts): Pass nulli as an argument.
    solve(eqs, vars, callback);
    vars[nulli] = true;
    solve(eqs, vars, callback);  // TODO(pts): Maybe tail recursion?
    vars[nulli] = null;
  }
}

var vars = [null, null, null, null];
var eqs  = [[1, 2, 3]];
solve(eqs, vars, print);
