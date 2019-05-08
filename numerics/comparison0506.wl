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
CHOP[exp_] := exp // Chop[#, 1*^-30]& //N[#,10]&;


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
  "\[CapitalDelta]msq3l"->2.525*^-21 (* + 7.39*^-23 *),
  "M1"->6.87433*^6,
  "M2"->6.87433*^6 + 277.251
|>;
M = {#M1, #M2}&[may6] // SetPrecision[#, 50] &;
m = NuFIT`ToMasses[may6];
U = NuFIT`ToPMNS[may6];
R = RRep["NH", #z, #\[Zeta]]&[may6];
Y = (Sqrt[2]I/vev) DiagonalMatrix[Sqrt[M]].R.DiagonalMatrix[Sqrt[m]].ConjugateTranspose[U]//FullSimplify // SetPrecision[#, 50]&;
{msq[M[[1]]]/2, M[[1]]^2/(16\[Pi]^2) Tr[Y.ConjugateTranspose[Y]]} // CHOP


(* Verification of CI parameterization *)
Mnu = ArrayFlatten[{{0, Transpose[Y]vev/Sqrt[2]},{Y vev/Sqrt[2], {{M[[1]], 0}, {0, M[[2]]}}}}];
Eigenvalues[Mnu.ConjugateTranspose[Mnu]] //CHOP//Sqrt
Sort[%] / {1, m[[2]],m[[3]],M[[1]],M[[2]]}


(* Convention: h in 1404.1003 is equal to Y in 1611.03827. *)
res = TwoGenResLept[ConjugateTranspose[Y], M];
res["\[Epsilon]mix(full)"]//CHOP (*2.22 of 1404.1003 *)
res["\[Epsilon]mix"]//CHOP  (*6.1 of 1608 or A.2 of 1404 *)
res["\[Epsilon]osc"]//CHOP  (*6.3 of 1608 *)
res["\[Epsilon]"]//CHOP     (* mix + osc *)
%//Flatten//Total


gstar = (6+3+3+2+1)*2*3 * (7/8) + 2*2 + (8 + 3 + 1)*2;
de = DiffEqs[res, gstar];

zl = Table[1.25Log[25de["Keff"][[l]]], {l,3}];
zc = M[[1]]/149;
(*2.33 of 1404*) (3/2)Sum[res["\[Epsilon]"][[l,a]]/(de["Keff"][[l]] 1.25Log[25de["Keff"][[l]]]), {l,3},{a,2}]
(*6.4  of 1608*) -(28/51)(1/27)(3/2)Sum[res["\[Epsilon]"][[l,a]]/(de["Keff"][[l]] Min[zl[[l]], zc]), {l,3},{a,2}]


{lz0, lzSph} = {-5, Log[M[[1]] / 149]};
\[Eta]Init = {\[Eta]L[1][lz0]==1, \[Eta]L[2][lz0]==1};

\[Eta]DE = de["\[Eta]DE"][\[Eta],#,z] &/@ {1, 2} // Simplify
Chop[%] /. \[Eta][a_]->Function[z, \[Eta]L[a][Log[z]]] /. z->Exp[lz] // FullSimplify
\[Eta]sol = NDSolve[{%, \[Eta]Init}, {\[Eta]L[1], \[Eta]L[2]}, {lz, lz0, lzSph}];
{\[Eta]L[1][lzSph], \[Eta]L[2][lzSph], \[Eta]L[1][lzSph]+\[Eta]L[2][lzSph]} /. \[Eta]sol

LogPlot[{
    {\[Eta]L[1][lz], \[Eta]L[2][lz]}/.\[Eta]sol,
    de["\[Eta]SolApprox"][Exp[lz]]
 } //Evaluate, {lz, lz0, lzSph}, PlotLegends->{"\!\(\*SubscriptBox[\(\[Eta]\), \(1\)]\)", "\!\(\*SubscriptBox[\(\[Eta]\), \(2\)]\)", "\!\(\*SubsuperscriptBox[\(\[Eta]\), \(1\), \(approx\)]\)(2.25)", "\!\(\*SubsuperscriptBox[\(\[Eta]\), \(2\), \(approx\)]\)(2.25)"}]


{de["\[Delta]DE"][\[Eta],\[Delta],1,z], de["\[Delta]DE"][\[Eta],\[Delta],2,z], de["\[Delta]DE"][\[Eta],\[Delta],3,z]}//Simplify
\[Delta]DEtoSolve = % //. {
  \[Eta][a_]:>Function[z, \[Eta]L[a][Log[z]]],
  \[Delta][a_]:>Function[z, \[Delta]L[a][Log[z]]],
  z->Exp[lz]
} //. \[Eta]sol // Simplify;
\[Delta]sol = NDSolve[{\[Delta]DEtoSolve, (\[Delta]L[#][lz0]==1*^-18) &/@ {1,2,3}}, {\[Delta]L[1], \[Delta]L[2], \[Delta]L[3]}, {lz, lz0, lzSph}];
{\[Delta]L[1][lzSph], \[Delta]L[2][lzSph], \[Delta]L[3][lzSph], \[Delta]L[1][lzSph] + \[Delta]L[2][lzSph] + \[Delta]L[3][lzSph]} /. \[Delta]sol

LogPlot[{
    \[Delta]L[#][lz]&/@{1,2,3} /. \[Delta]sol,
    de["\[Delta]SolApprox"][Exp[lz]]
  } // Abs // Evaluate, {lz, lz0, lzSph},
  PlotLegends->{"\!\(\*SubscriptBox[\(\[Delta]\[Eta]\), \(1\)]\)", "\!\(\*SubscriptBox[\(\[Delta]\[Eta]\), \(2\)]\)", "\!\(\*SubscriptBox[\(\[Delta]\[Eta]\), \(3\)]\)", "\!\(\*SubsuperscriptBox[\(\[Delta]\[Eta]\), \(1\), \(approx\)]\)(2.33)", "\!\(\*SubsuperscriptBox[\(\[Delta]\[Eta]\), \(2\), \(approx\)]\)(2.33)", "\!\(\*SubsuperscriptBox[\(\[Delta]\[Eta]\), \(2\), \(approx\)]\)(2.33)"}]
