(* ::Package:: *)

(* ::Section:: *)
(*RGE for normal hierarchy*)


(* ::Text:: *)
(*Be careful that 1307.3536 uses (d X)/(d Log[\[Mu]^2]) as the LHS of the RGE. So, for example, SMRGE`RGE[g1sq] == d (g1^2)/(d Log[\[Mu]^2]) == g1 \[Mu] dg1/d\[Mu], and SMRGE`RGE[lam] == d (lam)/(d Log[\[Mu]^2]) == (\[Mu]/2) d\[Lambda]/d\[Mu].*)
(*So we have wrapped the original RGE by factor two:*)
(*RGE[x] == d x/(d Log[\[Mu]]) == \[Mu] d x/d\[Mu] = 2*SMRGE`RGE[x].*)


SetDirectory[NotebookDirectory[]];
Needs["SMRGE`", "1307_3536.wl"];
$ContextPath = Cases[$ContextPath, Except["SMRGE`"]];
(* rewrite RGEs *)
RGE[a_] := 2 SMRGE`RGE[a] /. {SMRGE`FourPi -> (4\[Pi])}
(* rename variables *)
(Evaluate[Symbol["Global`"<>SymbolName[#]]] = #) &/@ SMRGE`couplings;


(* ::Input:: *)
(*(* to reload *)*)
(*<<1307_3536.wl*)
(*$ContextPath = Cases[$ContextPath, Except["SMRGE`"]];*)


$Assumptions = Join[Evaluate[# \[Element] Reals & /@ {delta, c1, c2, c3, s1, s2, s3}], {s1^2+c1^2 == 1, s2^2+c2^2 == 1, s3^2+c3^2 == 1}];


RGE[m\[Nu][k_]] := 1/(16\[Pi]^2) * (-m\[Nu][k](3yasq Abs[V[3,k]]^2 + 3g2sq - 4lam - 6ytsq - 6ybsq - 2yasq));
RGE[\[Theta]1] = 1/(16\[Pi]^2) * (3/2)(yasq/c2)Re[s3 V[3,1] Vc[3,3] + c3 T32];
RGE[\[Theta]2] = 1/(16\[Pi]^2) * (3yasq/2) Re[Exp[-I delta](-c3 V[3,1]Vc[3,3] + s3 T32)];
RGE[\[Theta]3] = 1/(16\[Pi]^2) * (-3yasq/2)(s2/c2)Re[(c2/s2)V[3,1]Vc[3,2] + Exp[-I delta](s3 V[3,1]Vc[3,3] - c3 T32)];
RGE[delta] = (3yasq/2) Im[Total[{
  (V[3,1] Vc[3,2])/(c3 s3),
  -s3/(s1 c2 c3) * V[3,1]Vc[2,2]Vc[3,3],
  -Exp[-I delta]/(s2 c1 c2) * V[3,1]V[2,2]Vc[3,3],
  -T32 Exp[-I delta] V[2,1]/(c1 c2 s2),
  +T32 c3 Vc[2,1]/(s1 s3 c2)}]];
RGE[phi] = 3yasq Im[Total[{
  c3 V[3,1]Vc[3,2]/s3,
  Vc[2,1]Vc[3,3]V[3,1] / (s1 c2),
  V[3,1]Vc[3,1]Vc[3,3] / (c1 c2),
  T32 c3 Vc[2,1]/(s1 c2 s3),
  -T32 Vc[3,2] / (c1 c2)}]];
RGE[phip] = 3yasq Im[Total[{
  s3 V[3,1]Vc[3,2] / c3,
  V[3,1]Vc[3,1]Vc[3,3] / (c1 c2),
  -s3 V[3,1]Vc[2,2]Vc[3,3] / (c1 c2 c3),
  -T32 Vc[2,2] / (s1 c2),
  -T32 Vc[3,2] / (c1 c2)}]];


Vmatrix = Module[{e = Exp[I delta]},{
  {c2 c3,                s3 c2,              s2/e },
  {-c1 s3 - s1 s2 c3 e,  c1 c3 - s1 s2 s3 e, s1 c2},
  { s1 s3 - c1 s2 c3 e, -s1 c3 - c1 s2 s3 e, c1 c2}}];
V[i:1|2|3, j:1|2|3] := Vmatrix[[i,j]];
Vc[i:1|2|3, j:1|2|3] := Conjugate[V[i,j]]
T32 = Exp[I phip]Vc[3,2]V[3,3](2m\[Nu][2]m\[Nu][3])/(m\[Nu][2]^2-m\[Nu][3]^2) + V[3,2]Vc[3,3](m\[Nu][2]^2+m\[Nu][3]^2)/(m\[Nu][2]^2-m\[Nu][3]^2);


(* ::Input:: *)
(*(* Verify V is unitary *)*)
(*Vmatrix.Conjugate[Transpose[Vmatrix]] // FullSimplify*)


RGE[phi]//FullSimplify
