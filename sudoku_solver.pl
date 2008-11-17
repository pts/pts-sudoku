%
% sudoku_solver.pl -- Sudoku solver in standard Prolog.
% by pts@fazekas.hu at Sun Nov 16 21:29:57 CET 2008
% --- Mon Nov 17 02:10:17 CET 2008
%
% Use the sudoku/2 and the solve/3 predicates to solve a sudoku puzzle. For
% example:
%
% ?- sudoku(s(2,[[[2],[],[4],[]],[[],[4],[],[2]],
%                [[],[2],[1],[3]],[[3],[],[],[4]]]),s(2,X)).
% X = [[2,3,4,1],[1,4,3,2],[4,2,1,3],[3,1,2,4]] ? ;
% no
%
% See more examples in test*/0 below.
%

% gs_valid(GS0) binds all G in GS except if GS has no t and at least 2 vars.
gs_valid([G|GS]) :-
  ( G == t -> gs_all_f(GS)
  ; var(G) -> 
    ( gs_first_nonf(GS, G2, GS2) ->
      ( var(G2) ->
        ( gs_has_t(GS2) -> G = f, G2 = f, gs_maxone_t(GS2)
        ; true  % We don't choose here.
        )
      ; G = f, gs_all_f(GS2)
      )
    ; G = t
    )
  ; gs_valid(GS)
  ).

gs_all_f([]).
gs_all_f([f|GS]) :-
  gs_all_f(GS).

gs_first_nonf([G|GS], G2, GS2) :-
  ( G == f -> gs_first_nonf(GS, G2, GS2)
  ; G2 = G, GS2 = GS
  ).

gs_has_t([G|GS]) :-
  ( G == t -> true
  ; gs_has_t(GS)
  ).

gs_maxone_t([]).
gs_maxone_t([G|GS]) :-
  ( G == t -> gs_all_f(GS)
  ; G = f, gs_maxone_t(GS)
  ).

% This may bind.
gss_valid([]).
gss_valid([GS|GSS]) :-
  gs_valid(GS),
  gss_valid(GSS).

gs_varcount([], C, C).
gs_varcount([G|GS], C1, C3) :-
  ( var(G) -> C2 is C1 + 1
  ; C2 = C1
  ),
  gs_varcount(GS, C2, C3).

gss_varcount([], C, C).
gss_varcount([GS|GSS], C1, C3) :-
  gs_varcount(GS, C1, C2),
  gss_varcount(GSS, C2, C3).

% Useful.
gss_label(GSS) :-
  gss_varcount(GSS, 0, C1),
  gss_label_next(GSS, C1).

gss_label_next(GSS, C1) :-
  gss_valid(GSS),
  gss_varcount(GSS, 0, C2),
  ( C1 == C2 -> gss_choice_without_ground_gs(GSS)
  ; gss_label_next(GSS, C2)
  ).

gss_choice_without_ground_gs(GSS) :-
  ( GSS = [] -> true
  ; GSS = [GS|GSS2],
    ( ground(GS) -> gss_choice_without_ground_gs(GSS2)
    ; gs_fvar(GS, G2, GS2),
      ( G2 = t, gs_all_f(GS2), gss_label(GSS2)
      ; G2 = f, gss_label(GSS)
      )
    )
  ).

gs_fvar([G|GS], G2, GS2) :-
  ( var(G) -> G2 = G, GS2 = GS
  ; G \== t, gs_fvar(GS, G2, GS2)
  ).

% ---

% Useful. Convert list to bitmap.
list_bitmap(L, RR, BM) :-
  ( L = [] -> fromto_any(1, RR, BM)
  ; sort(L, L1),
    L1 = [B|L2],
    sorted_bitmap(1, B, L2, RR, BM)
  ).

fromto_any(A, B, BM) :-
  ( A =< B -> BM = [_|BM1], A1 is A + 1, fromto_any(A1, B, BM1)
  ; BM = []
  ).

fromto_f(A, B, BM) :-
  ( A =< B -> BM = [f|BM1], A1 is A + 1, fromto_f(A1, B, BM1)
  ; BM = []
  ).

sorted_bitmap(A, B, L, RR, BM) :-
  ( A > RR -> BM = []
  ; A < B -> A1 is A + 1, BM=[f|BM1], sorted_bitmap(A1, B, L, RR, BM1)
  ; A == B -> A1 is A + 1, BM=[_|BM1], sorted_bitmap_low(L, A1, RR, BM1)
  ).

sorted_bitmap_low([], A, RR, BM) :-
  fromto_f(A, RR, BM).
sorted_bitmap_low([B|L], A, RR, BM) :-
  ( B >= A -> sorted_bitmap(A, B, L, RR, BM)
  ; sorted_bitmap_low(L, A, RR, BM)
  ).

bitmap_one([G|GS], A, J) :-
  ( G == t -> J = A
  ; A1 is A + 1, bitmap_one(GS, A1, J)
  ).   

hss_append_levels(HSS, A, RR, GSS1, GSS3) :-
  ( A > RR -> GSS3 = GSS1
  ; hss_append_level(HSS, A, GSS1, GSS2),
    A1 is A + 1, hss_append_levels(HSS, A1, RR, GSS2, GSS3)
  ).

hss_append_level([], _A, GSS, GSS).
hss_append_level([HS|HSS], A, GSS1, [GS|GSS2]) :-
  hs_level(HS, A, GS),
  hss_append_level(HSS, A, GSS1, GSS2).

hs_level([], _A, []).
hs_level([H|HS], A, [G|GS]) :-
  my_nth1(A, H, G),  % TODO: speed up to build multiple GSSs at once.
  hs_level(HS, A, GS).

my_nth1(A, [G|GS], G2) :-
  ( A =< 1 -> A == 1, G2 = G
  ; A1 is A-1, my_nth1(A1, GS, G2)
  ).

xss_bitmap_one([], []).
xss_bitmap_one([XS|XSS], [VS|VSS]) :-
  xs_bitmap_one(XS, VS),
  xss_bitmap_one(XSS, VSS).

xs_bitmap_one([], []).
xs_bitmap_one([X|XS], [V|VS]) :-
  bitmap_one(X, 1, V),
  xs_bitmap_one(XS, VS).

xss_list_bitmap([], _RR, []).
xss_list_bitmap([XS|XSS], RR, [VS|VSS]) :-
  xs_list_bitmap(XS, RR, VS),
  xss_list_bitmap(XSS, RR, VSS).

xs_list_bitmap([], _RR, []).
xs_list_bitmap([X|XS], RR, [V|VS]) :-
  list_bitmap(X, RR, V),
  xs_list_bitmap(XS, RR, VS).

xss_structure([], []).
xss_structure([XS|XSS], [VS|VSS]) :-
  xs_structure(XS, VS),
  xss_structure(XSS, VSS).

xs_structure([], []).
xs_structure([_|XS], [_|VS]) :-
  xs_structure(XS, VS).

slow_sqrt(RR, R) :-
  slow_sqrt_incr(0, RR, R).

slow_sqrt_incr(A, RR, R) :-
  AA is A * A,
  ( AA < RR -> A1 is A + 1, slow_sqrt_incr(A1, RR, R)
  ; AA == RR -> R = A
  ).

all_length([], _N).
all_length([L|LS], N) :-
  length(L, N),
  all_length(LS, N).

flatten([], GSS, GSS).
flatten([GSS|GSSS], GSS1, GSS3) :-
  my_append(GSS, GSS2, GSS3),
  flatten(GSSS, GSS1, GSS2).

my_append([], X, X).
my_append([H|T], L, [H|M]) :-
  my_append(T, L, M).

transpose(L, M) :-
  ( L = [K|_],
    ( K = [] -> M = []
    ; column_separate(L, HS, TS),
      M = [HS|MS],
      transpose(TS, MS)
    )
  ).

column_separate([], [], []).
column_separate([[H|T]|Q], [H|HS], [T|TS]) :-
  column_separate(Q, HS, TS).

transpose_r(L, R, M) :-
  ( L = [K|_],
    ( K = [] -> M = []
    ; column_separate_r(L, R, HS, TS),
      M = [HS|MS],
      transpose_r(TS, R, MS)
    )
  ).

column_separate_r([], _R, [], []).
column_separate_r([P|Q], R, HX, [T|TS]) :-
  length(H, R),
  my_append(H, T, P),
  my_append(H, HS, HX),
  column_separate_r(Q, R, HS, TS).

process_r(L1, R, L2) :-
  ( L1 = [] -> L2 = []
  ; length(H1, R),
    my_append(H1, T1, L1),
    transpose_r(H1, R, H2),
    my_append(H2, T2, L2),
    process_r(T1, R, T2)
  ). 


solve(Input, Solution, R) :-
  length(Input, RR),
  slow_sqrt(RR, R),
  all_length(Input, RR),
  xss_structure(Input, Solution),
  xss_list_bitmap(Input, RR, BitmapSS),
  flatten(BitmapSS, [], GSS1),  % cell constraints
  process_r(BitmapSS, R, HSS1),  % small square constraints
  my_append(BitmapSS, HSS1, HSS2),  % row contraints
  transpose(BitmapSS, HSS3),  % column contraints
  my_append(HSS3, HSS2, HSS4),
  hss_append_levels(HSS4, 1, RR, GSS1, GSS2),
  gss_label(GSS2),  % This is the only nondeterministic call.
  xss_bitmap_one(BitmapSS, Solution).

sudoku(s(R, Input), s(R, Solution)) :-
  solve(Input, Solution, R).

test0 :-
  sudoku(s(2,[[[2],[],[4],[]],[[],[4],[],[2]],
              [[],[2],[1],[3]],[[3],[],[],[4]]]),s(2,Solution)),
  print(Solution), nl, fail.

% [2,3,4,1]
% [1,4,3,2]
% [4,2,1,3]
% [3,1,2,4]
%
test1 :-
  solve([[[2],[],[4],[]],[[],[4],[],[2]],[[],[2],[1],[3]],[[3],[],[],[4]]],
      [[A, B, C, D], [E, F, G, H], [I, J, K, L], [M, N, O, P]], _R),
  print([A, B, C, D]), nl,
  print([E, F, G, H]), nl,
  print([I, J, K, L]), nl,
  print([M, N, O, P]), nl,
  nl, fail.

test2 :-
  solve([[[],[],[],[]],[[],[],[],[]],[[],[],[],[]],[[],[],[],[]]],
      [[A, B, C, D], [E, F, G, H], [I, J, K, L], [M, N, O, P]], _R),
  print([A, B, C, D]), nl,
  print([E, F, G, H]), nl,
  print([I, J, K, L]), nl,
  print([M, N, O, P]), nl,
  nl, fail.

test3 :-
  solve([[[ ],[2],[1],[4],[ ],[ ],[8],[5],[ ]],
         [[ ],[ ],[5],[8],[ ],[ ],[ ],[ ],[6]],
         [[ ],[7],[ ],[6],[9],[5],[ ],[2],[ ]],
         [[9],[ ],[ ],[ ],[ ],[8],[5],[ ],[ ]],
         [[ ],[3],[4],[ ],[ ],[ ],[7],[9],[ ]],
         [[ ],[ ],[6],[7],[ ],[ ],[ ],[ ],[1]],
         [[ ],[4],[ ],[5],[3],[7],[ ],[8],[ ]],
         [[7],[ ],[ ],[ ],[ ],[4],[9],[ ],[ ]],
         [[ ],[6],[2],[ ],[ ],[1],[3],[4],[ ]]], S, _R),
  print(S), fail.
