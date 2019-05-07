(* ::Package:: *)

SetDirectory[Which[$FrontEnd=!=Null,NotebookDirectory[],$InputFileName=!="", DirectoryName[$InputFileName], True, "."]];
$ProjectRoot = RunProcess[{"git","rev-parse" , "--show-toplevel"}]["StandardOutput"] // StringTrim;

Needs["ResLept14041003`", "1404.1003.m"];
Needs["NuFIT`", "nufit40.m"];
msq = Get[FileNameJoin[{$ProjectRoot, "rge", "msq.dat"}]];

RRep["NH", z_, \[Zeta]_] := {{0, Cos[z], \[Zeta] Sin[z]}, {0, -Sin[z], \[Zeta] Cos[z]}};
RRep["IH", z_, \[Zeta]_] := {{Cos[z], \[Zeta] Sin[z], 0}, {-Sin[z], \[Zeta] Cos[z], 0}};
MajoranaPhase[\[Sigma]_] := DiagonalMatrix[{1, Exp[I \[Sigma]], 1}]; 

RemC[exp_] := exp //. Complex[a_Real, 0.`]:>a;


$Assumptions={\[CapitalDelta]M>0, lz\[Element]Reals};
vev = 246;

may6 = <|
  "z"->0.3 + 0.8I,
  "\[Zeta]"->1,
  "\[Delta]CP"->\[Pi]/2,
  "Sin2\[Theta]12"->0.31,
  "Sin2\[Theta]13"->0.02240,
  "Sin2\[Theta]23"->0.582,
  "\[CapitalDelta]msq21"->7.39*^-23,
  "\[CapitalDelta]msq3l"->2.525*^-21,
  "M1"->6.87433*^6,
  "M2"->6.87433*^6 + 277.251
|>;
M = {#M1, #M2}&[may6];
m = NuFIT`ToMasses[may6];
U = NuFIT`ToPMNS[may6];
R = RRep["NH", #z, #\[Zeta]]&[may6];
Y = (Sqrt[2]I/vev) DiagonalMatrix[Sqrt[M]].R.DiagonalMatrix[Sqrt[m]].ConjugateTranspose[U]//FullSimplify;
{msq[M[[1]]]/2, M[[1]]^2/(16\[Pi]^2) Tr[Y.ConjugateTranspose[Y]]}


res = TwoGenResLept[Transpose[Y], M];
res["\[Epsilon]mix(full)"]//Chop (*2.22 of 1404.1003 *)
res["\[Epsilon]mix"]//Chop  (*6.1 of 1608 or A.2 of 1404 *)
res["\[Epsilon]osc"]//Chop  (*6.3 of 1608 *)
res["\[Epsilon]"]//Chop     (* mix + osc *)


gstar = (6+3+3+2+1)*2*3 * (7/8) + 2*2 + (8 + 3 + 1)*2;
de = DiffEqs[res, gstar];

zl = Table[1.25Log[25de["Keff"][[l]]], {l,3}];
zc = M[[1]]/149;
(*2.33 of 1404*) (3/2)Sum[res["\[Epsilon]"][[l,a]]/(de["Keff"][[l]] 1.25Log[25de["Keff"][[l]]]), {l,3},{a,2}]
(*6.4  of 1608*) -(28/51)(1/27)(3/2)Sum[res["\[Epsilon]"][[l,a]]/(de["Keff"][[l]] Min[zl[[l]], zc]), {l,3},{a,2}]


