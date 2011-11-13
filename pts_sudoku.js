#! /usr/bin/smjs
//
// pts_sudoku.js: a sudoku-solver in JavaScript
// by pts@fazekas.hu at Fri Sep  4 11:01:47 CEST 2009
//
// TODO(pts): Micro-optimize.
//
// More advanced sudoku solvering methods:
// http://www.sudoku.com/boards/viewtopic.php?t=6556&start=0
//
// A solver on the web: http://zfacts.com/giz/sudoku.js
// http://zfacts.com/p/sudoku_solver.html
//
// The solver on http://www.omegajunior.net/code/sudokusolver.html
// seems to be using doMiddleIntersection.

// Instantiate variables in the equation (list of variables) eq. Returns
// null on contradiction, or the new insti.
function inst(eq, vars, insts, insti) {
  var eql = eq.length;
  var i, j, v, eqj;
  var nulli1 = null, nulli2 = null;
  // TODO(pts): Micro-optimize this function.
  // TODO(pts): Move an equation to the end if all variables found.
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
// Returns the number of true--false choice points encountered.
function solve(eqs, vars, callback, insts, insti) {
  //print('S ' + vars);
  var eqsl = eqs.length;
  var varl = vars.length;
  var insti0 = insti;
  var i, insti3, insti2;
  var numChoices = 0;
  while (true) {
    //print('Q ' + vars)
    insti3 = insti;
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
    if (insti != insti3)
      continue;  // Try to instantiate more.
    //print('I ' + vars);
    i = 0;  // TODO(pts): Pass this as well.
    while (i < varl && vars[i] != null) ++i;
    if (i == varl) {  // all variables are known
      callback(vars);
      break;
    } else {
      ++numChoices;
      if (0) { // !! debug
        print('CH' + i);
        var nn = 9;
        var nnnn = nn * nn;
        var solution = Array(nnnn);
        var l, m;
        for (l = 0; l < nnnn; ++l) {
          solution[l] = '?';
          for (m = 0; m < nn; ++m) {
            if (vars[l * nn + m])
              solution[l] = m + 1;
          }
        }
        for (l = 0; l < nn; ++l) {
          var s = '';
          for (h = 0; h < nn; ++h) {
            s += solution[l * nn + h] + ',';
          }
          print(s);
        }
        print('.');
      }

      vars[i] = false;
      insts[insti++] = i;
      // TODO(pts): Pass i as an argument.
      numChoices += solve(eqs, vars, callback, insts, insti);
      vars[i] = true;
      // Continue the while loop with ``true'', just as with:
      // numChoices += solve(eqs, vars, callback, insts, insti); break;
    }
  }
  while (insti > insti0) {
    vars[insts[--insti]] = null;
  }
  return numChoices;
}

// Test:
//var vars = Array(4);  // [null, null, null, null];
//var eqs  = [[1, 2, 3], [0, 2]];
//var insts = Array(4);
//solve(eqs, vars, print, insts, 0);

// Test:
//var vars = [false, true, true, false];
//var eqs  = [[1, 2, 3], [0, 2]];
//var insts = Array(4);
//print(inst([1, 2, 3], vars, insts, 0));

// Example problem:
//   var problem = [1,0,0,0,
//                  0,0,0,2,
//                  0,0,2,0,
//                  0,3,0,0];
// For the example problem, the callback gets called with:
//   [1,2,3,4,
//    3,4,1,2,
//    4,1,2,3,
//    2,3,4,1]
function sudoku(problem, callback, doMiddleIntersection) {
  var nnnn = problem.length;
  var nn = Math.sqrt(nnnn) ^ 0;
  if (nn * nn != nnnn) throw 'n^4 expected (not a perfect square)'
  var n = Math.sqrt(nn) ^ 0;
  if (n * n != nn) throw 'n^4 expected';
  var solution = Array(nnnn);
  if (doMiddleIntersection == null)
    doMiddleIntersection = true;

  function lowCallback(vars) {
    var l, m;
    for (l = 0; l < nnnn; ++l) {
      for (m = 0; m < nn; ++m) {
        if (vars[l * nn + m])
          solution[l] = m + 1;
      }
    }
    callback(solution);
  }

  var vars = Array(nn * nnnn + 2 * nnnn * n * ! !doMiddleIntersection);
  var vari = nn * nnnn;
  var eqs = [];
  var i, j, ii, jj, iii, jjj, k, eq, eqa, eqb, eqi;
  for (i = 0; i < nnnn; ++i) {
    if (problem[i] > 0) {
      if (problem[i] > nn) throw 'cell too large';
      vars[problem[i] - 1 + i * nn] = true;
    }
  }
  
  // Values in a single cell.
  for (i = 0; i < nnnn; ++i) {
    eq = Array(j);
    for (j = 0; j < nn; ++j) {
      eq[j] = j + i * nn;  // TODO(pts): Get rid of the `* nn', everywhere.
    }
    eqs.push(eq);
  }

  // Values in a row. (row i, column j, value k)
  for (k = 0; k < nn; ++k) {
    for (i = 0; i < nn; ++i) {
      eq = Array(nn);
      for (j = 0; j < nn; ++j) {
        eq[j] = k + j * nn + i * nnnn;
      }
      eqs.push(eq);
    }
  }

  // Values in a column. (row i, column j, value k)
  for (k = 0; k < nn; ++k) {
    for (j = 0; j < nn; ++j) {
      eq = Array(nn);
      for (i = 0; i < nn; ++i) {
        eq[i] = k + j * nn + i * nnnn;
      }
      eqs.push(eq);
    }
  }

  // Values in a middle square. (row (ii, i), column (jj, j), value k)
  for (k = 0; k < nn; ++k) {
    for (ii = 0; ii < n; ++ii) {
      for (jj = 0; jj < n; ++jj) {
        eq = Array(nn);
        for (i = 0; i < n; ++i) {
          for (j = 0; j < n; ++j) {
            eq[i * n + j] = k + (jj * n + j) * nn + (ii * n + i) * nnnn;
          }
        }
        eqs.push(eq);
      }
    }
  }

  if (doMiddleIntersection) {
    // TODO(pts): Investigate why enabling doMiddleIntersection makes the
    // algorithm slower for nn=9 (even though it decreases the number of
    // choices). Maybe because the so many extra equations?
    // TODO(pts): Use doMiddleIntersection lazily.
    
    // Intersection of a middle square and a column. (Middle row ii, middle
    // column jj, column j, full row i, value k).
    for (k = 0; k < nn; ++k) {
      for (ii = 0; ii < n; ++ii) {
        for (jj = 0; jj < n; ++jj) {
          for (j = 0; j < n; ++j) {
            eqa = Array(nn - n + 1);
            eqb = Array(nn - n + 1);
            eqi = 0;
            for (i = 0; i < ii * n; ++i) {
              eqa[eqi]   = k + (jj * n + j) * nn + i * nnnn;
              eqb[eqi++] = k + (jj * n + j) * nnnn + i * nn;
            }
            for (i += n; i < nn; ++i) {
              eqa[eqi  ] = k + (jj * n + j) * nn + i * nnnn;
              eqb[eqi++] = k + (jj * n + j) * nnnn + i * nn;
            }
            if (eqi != nn - n)
              throw 'bad equation length';
            eqa[eqi] = vari++;
            eqb[eqi] = vari++;
            eqs.push(eqa);
            eqs.push(eqb);

            eqa = Array(nn - n + 1);
            eqb = Array(nn - n + 1);
            eqi = 0;
            for (iii = 0; iii < n; ++iii) {
              for (jjj = 0; jjj < n; ++jjj) {
                if (jjj != j) {
                  eqa[eqi]   = k + (jj * n + jjj) * nn + (ii * n + iii) * nnnn;
                  eqb[eqi++] = k + (jj * n + jjj) * nnnn + (ii * n + iii) * nn;
                }
              }
            }
            if (eqi != nn - n)
              throw 'bad equation length';
            // Share the variable with the previous eqa, eqb equations.
            eqa[eqi] = vari - 2;
            eqb[eqi] = vari - 1;
            eqs.push(eqa);
            eqs.push(eqb);
          }
        }
      }
    }
    if (vari != nn * nnnn + 2 * nnnn * n)
      throw 'variable count mismatch';
  }

  // print(eqs.join('\n'));
  return solve(eqs, vars, lowCallback, Array(nn * nnnn), 0);
}

print('--- 1:')
var problem1 = [1,0,0,0,
                0,0,0,2,
                0,0,2,0,
                0,3,0,0];
// 1,2,3,4,3,4,1,2,4,1,2,3,2,3,4,1
print('SUMCH' + sudoku(problem1, print));

// TODO(pts): Investigate why we need 5 choices here (even with
// doMiddleIntersection).
//
// Gnome Sudoku very hard.
print('--- 2:')
var problem2 = [8,0,0,1,0,3,0,5,0,
                3,1,0,4,0,0,8,0,0,
                0,0,5,0,0,0,0,0,0,
                0,9,0,2,0,0,0,0,5,
                0,5,8,7,0,9,3,6,0,
                4,0,0,0,0,5,0,8,0,
                0,0,0,0,0,0,2,0,0,
                0,0,2,0,0,1,0,3,7,
                0,3,0,8,0,2,0,0,1];
// 8,6,4,1,2,3,7,5,9,3,1,9,4,5,7,8,2,6,7,2,5,6,9,8,4,1,3,6,9,3,2,8,4,1,7,5,2,5,8,7,1,9,3,6,4,4,7,1,3,6,5,9,8,2,1,4,7,5,3,6,2,9,8,5,8,2,9,4,1,6,3,7,9,3,6,8,7,2,5,4,1
print('SUMCH' + sudoku(problem2, print));

// TODO(pts): Investigate why we need 2 choices here (with or without
// doMiddleIntersection).
//
// Gnome Sudoku very hard.
print('--- 3:')
var problem3 = [0,8,0,0,0,2,7,0,0,
                5,0,3,7,0,0,4,0,0,
                0,6,7,0,5,0,0,0,3,
                6,9,0,0,7,4,0,0,5,
                0,0,5,0,0,0,1,0,0,
                2,0,0,9,1,0,0,7,4,
                3,0,0,0,2,0,8,6,0,
                0,0,2,0,0,6,5,0,7,
                0,0,6,5,0,0,0,3,0];
// 1,8,9,4,3,2,7,5,6,5,2,3,7,6,9,4,1,8,4,6,7,8,5,1,9,2,3,6,9,1,2,7,4,3,8,5,7,4,5,6,8,3,1,9,2,2,3,8,9,1,5,6,7,4,3,5,4,1,2,7,8,6,9,8,1,2,3,9,6,5,4,7,9,7,6,5,4,8,2,3,1
print('SUMCH' + sudoku(problem3, print, true));

// http://www.forbeginners.info/sudoku-puzzles/extreme-20.htm
print('--- 4:')
var problem4 = [6,0,0,0,0,0,0,4,0,
                0,0,5,0,0,2,0,0,7,
                7,2,9,0,0,0,0,0,3,
                0,9,0,0,4,0,0,0,1,
                0,0,0,0,6,0,0,0,0,
                4,0,0,0,8,0,0,7,0,
                3,0,0,0,0,0,1,6,5,
                2,0,0,4,0,0,8,0,0,
                0,5,0,0,0,0,0,0,4];
print('SUMCH' + sudoku(problem4, print, true));

// Hardest sudoku on 2006-07-11:
// http://www.usatoday.com/news/offbeat/2006-11-06-sudoku_x.htm
// Seems to be very easy for the software (0 choices with doMiddleIntersection,
// 7 choices without it).
print('--- 5:')
var problem5 = [8,5,0,0,0,2,4,0,0,
                7,2,0,0,0,0,0,0,9,
                0,0,4,0,0,0,0,0,0,
                0,0,0,1,0,7,0,0,2,
                3,0,5,0,0,0,9,0,0,
                0,4,0,0,0,0,0,0,0,
                0,0,0,0,8,0,0,7,0,
                0,1,7,0,0,0,0,0,0,
                0,0,0,0,3,6,0,4,0];
print('SUMCH' + sudoku(problem5, print, true));

// http://www.sudokuoftheday.co.uk/sudoku.php?ACTION=sotd&RES=1024&LEVEL=A
print('--- 6:')
var problem6 = [6,4,0,9,0,0,0,1,5,
                0,8,9,0,1,0,0,0,2,
                0,0,7,0,0,8,0,0,0,
                0,6,0,0,4,0,3,0,0,
                3,0,1,6,0,0,0,0,0,
                0,0,0,0,0,0,0,2,0,
                0,0,0,3,0,0,0,0,0,
                9,0,0,0,0,5,0,3,0,
                0,0,0,0,0,0,0,4,0];
print('SUMCH' + sudoku(problem6, print, true));

// http://www.sudoku.com/boards/viewtopic.php?t=5336&postdays=0&postorder=asc&start=43
// Lots of more, similar puzzles on the site above.
// This is really hard (117 choices with doMiddleIntersection, 153 choices
// without it).
print('--- 7:')
var problem7 = [5,0,0,0,0,0,0,0,9,
                0,2,0,1,0,0,0,7,0,
                0,0,8,0,0,0,3,0,0,
                0,4,0,6,0,0,0,0,0,
                0,0,0,0,5,0,0,0,0,
                0,0,0,2,0,7,0,1,0,
                0,0,3,0,0,0,8,0,0,
                0,6,0,0,0,4,0,2,0,
                9,0,0,0,0,0,0,0,5];
// 5,3,6,4,7,2,1,8,9,4,2,9,1,8,3,5,7,6,1,7,8,9,6,5,3,4,2,7,4,1,6,9,8,2,5,3,6,8,2,3,5,1,4,9,7,3,9,5,2,4,7,6,1,8,2,5,3,7,1,9,8,6,4,8,6,7,5,3,4,9,2,1,9,1,4,8,2,6,7,3,5
print('SUMCH' + sudoku(problem7, print, true));

// http://www.sudoku.com/boards/viewtopic.php?t=6539&start=0
// Lots of more, similar puzzles on the site above.
// This is really hard (580 choices with doMiddleIntersection, 669 choices
// without it).
print('--- 8:')
var problem8 = [0,0,0,0,0,0,0,1,2,
                0,0,0,0,0,0,0,0,3,
                0,0,2,3,0,0,4,0,0,
                0,0,1,8,0,0,0,0,5,
                0,6,0,0,7,0,8,0,0,
                0,0,0,0,0,9,0,0,0,
                0,0,8,5,0,0,0,0,0,
                9,0,0,0,4,0,5,0,0,
                4,7,0,0,0,6,0,0,0];
// 8,3,9,4,6,5,7,1,2,1,4,6,7,8,2,9,5,3,7,5,2,3,9,1,4,8,6,3,9,1,8,2,4,6,7,5,5,6,4,1,7,3,8,2,9,2,8,7,6,5,9,3,4,1,6,2,8,5,3,7,1,9,4,9,1,3,2,4,8,5,6,7,4,7,5,9,1,6,2,3,8
print('SUMCH' + sudoku(problem8, print, false));
