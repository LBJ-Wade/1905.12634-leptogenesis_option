#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Section:: *)
(*Casas\[Dash]Ibarra Parameterization [0103065] for two heavy neutrinos*)


VarSub[v_] := ToString[v] // StringReplace[#, {RegularExpression["\$\\d+"] -> "", RegularExpression["\\w+`"] -> "", " " -> ""}] &
SubsuperscriptBoxSub[z_,{},    {},    f_]:=MakeBoxes[z,f]
SubsuperscriptBoxSub[z_,a_List,{},    f_]:=SubscriptBox     @@(MakeBoxes[#,f] &/@ {z,Row[VarSub/@a]})
SubsuperscriptBoxSub[z_,{},    b_List,f_]:=SuperscriptBox   @@(MakeBoxes[#,f] &/@ {z,Row[VarSub/@b]})
SubsuperscriptBoxSub[z_,a_List,b_List,f_]:=SubsuperscriptBox@@(MakeBoxes[#,f] &/@ {z,Row[VarSub/@a],Row[VarSub/@b]})
MakeBoxes[obj:Cos[t12],f:StandardForm|TraditionalForm]:=InterpretationBox[#,obj]&@SubsuperscriptBoxSub["c",{"12"},{},f]
MakeBoxes[obj:Cos[t23],f:StandardForm|TraditionalForm]:=InterpretationBox[#,obj]&@SubsuperscriptBoxSub["c",{"23"},{},f]
MakeBoxes[obj:Cos[t13],f:StandardForm|TraditionalForm]:=InterpretationBox[#,obj]&@SubsuperscriptBoxSub["c",{"13"},{},f]
MakeBoxes[obj:Sin[t12],f:StandardForm|TraditionalForm]:=InterpretationBox[#,obj]&@SubsuperscriptBoxSub["s",{"12"},{},f]
MakeBoxes[obj:Sin[t23],f:StandardForm|TraditionalForm]:=InterpretationBox[#,obj]&@SubsuperscriptBoxSub["s",{"23"},{},f]
MakeBoxes[obj:Sin[t13],f:StandardForm|TraditionalForm]:=InterpretationBox[#,obj]&@SubsuperscriptBoxSub["s",{"13"},{},f]

$Assumptions = {t12>=0, t23>=0, t13>=0, \[Alpha]21>=0, \[Alpha]31>=0, \[Delta]>=0, v2>0, M[_]>0, m[_]>=0, GeV>0, eV>0};

GeV::usage = "Unit for Giga-electronvolt.";
eV::usage = "Unit for electronvolt.";

Unprotect[GeV, eV];
eV=10^-9 GeV;
N[GeV]=1;
Protect[GeV, eV];


(* PMNS matrix with notation in PDG2016(14.6) *)
V23 = {{1,0,0},{0,Cos[t23],Sin[t23]},{0,-Sin[t23],Cos[t23]}};
V13 = {{Cos[t13],0,Sin[t13]*Exp[-I*\[Delta]]},{0,1,0},{-Sin[t13]*Exp[I*\[Delta]],0,Cos[t13]}};
V12 = {{Cos[t12],Sin[t12],0},{-Sin[t12],Cos[t12],0},{0,0,1}};
V0  = V23.V13.V12;
V = V0.DiagonalMatrix[Exp[{0, I \[Alpha]21/2, I \[Alpha]31/2}]];

(* CI (5) and (10) *)
(* v2 ~ 174 GeV is the down-type higgs vev (and vSM/Sqrt[2]). *)
Dm = DiagonalMatrix[{m[1], m[2], m[3]}];
D\[Kappa] = Dm / v2^2;
DM = DiagonalMatrix[{M[1], M[2]}];

(* CI (13) *)
Y\[Nu]::usage = "Neutrino Yukawa matrix with (N,3)-shape.";
Y\[Nu] = Sqrt[DM].RR.Sqrt[D\[Kappa]].ConjugateTranspose[V];

(* parameterization; values based on PDG2018 *)
NH[z_] := NH[z, 1];
IH[z_] := IH[z, 1];
NH[z_, quadrant:1|2] := WithData[{
  RR -> {{0, Sqrt[1-z^2], z}, {0, -z, Sqrt[1-z^2]}},
  m[1]->0, m[2]->Sqrt[\[CapitalDelta]m21sq], m[3]->Sqrt[\[CapitalDelta]m32sqNH+\[CapitalDelta]m21sq], t23->t23NH}, quadrant]
IH[z_, quadrant:1|2] := WithData[{
  RR -> {{z, Sqrt[1-z^2], 0}, {Sqrt[1-z^2], -z, 0}},
  m[1]->Sqrt[-\[CapitalDelta]m32sqIH-\[CapitalDelta]m21sq], m[2]->Sqrt[-\[CapitalDelta]m32sqIH], m[3]->0, t23->t23IH}, quadrant]
WithData[rules_, quadrant:1|2] := Join[
  rules //. {
    \[CapitalDelta]m21sq->7.53*^-5 eV^2,
    \[CapitalDelta]m32sqIH->-2.56*^-3 eV^2,
    \[CapitalDelta]m32sqNH->2.51*^-3 eV^2,
    t23NH->{ArcSin[Sqrt[0.417]], \[Pi]-ArcSin[Sqrt[0.597]]}[[quadrant]],
    t23IH->{ArcSin[Sqrt[0.421]], \[Pi]-ArcSin[Sqrt[0.592]]}[[quadrant]]
  },{
    t12->ArcSin[Sqrt[0.307]],
    t13->ArcSin[Sqrt[2.12*^-2]],
    v2->246GeV,
    \[Delta]->1.5
  }
];


Y\[Nu]//.IH[z]//Simplify


masses = Join[NH[z], {\[Alpha]21->0, \[Alpha]31->0, M[1]->10^10 GeV, M[2]->1.01M[1]}];
Y\[Nu] //. masses // Simplify


massmatrix = ArrayFlatten[{{Table[0,3,3], Transpose[Y\[Nu]] v2}, {Y\[Nu] v2, DM}}//. masses // N]//Simplify
massmatrix /. z->0.3+5 I // Rationalize[#, 0] & // SingularValueDecomposition[#][[2]]& // N // Diagonal // Sort
{M[1], M[2], m[1], m[2], m[3]} //. masses // N // Sort


(* ::Subchapter:: *)
(*Resonant Leptogenesis*)


masses = Join[NH[RandomReal[] + RandomReal[] I], {\[Alpha]21->0, \[Alpha]31->0, M[1]->10^15 GeV, M[2]->1.1M[1]}];
Dag[h_] := ConjugateTranspose[Mat[h]]

(#[p__Integer] := Mat[#][[p]]) &/@ {h, hs, hdh, A, H, HC, Hs, H$, HdH, HCdHC};

Mat[h] = Transpose[Y\[Nu]] //. masses // Simplify
Mat[hs] = Conjugate[Mat[h]]

(* hdh == h^\[Dagger]h *)
Mat[hdh] = Dag[h].Mat[h]
Mat[A] = Conjugate[Mat[hdh]] / (16\[Pi])

Mat[H] = Module[{b=3-a}, Table[h[l,a]-I h[l,b]((M[a](M[a]A[a,b]+M[b]A[b,a]))/(M[a]^2-M[b]^2+2I M[a]^2A[b,b])), {l, 3}, {a, 2}]] //. masses
Mat[HC] = Module[{b=3-a}, Table[hs[l,a]-I hs[l,b]((M[a](M[a]A[b,a]+M[b]A[a,b]))/(M[a]^2-M[b]^2+2I M[a]^2A[b,b])), {l, 3}, {a,2}]] //. masses
Mat[Hs] = Conjugate[Mat[H]]
Mat[H$] = Conjugate[Mat[HC]]
Mat[HdH] = Dag[H].Mat[H]
Mat[HCdHC] = Dag[HC].Mat[HC]

eps[l_,a_]:=(Abs[H[l,a]]^2-Abs[HC[l,a]]^2)/Sum[Abs[H[k,a]]^2+Abs[HC[k,a]]^2,{k,3}]

freg[a:1|2]:=Module[{b=3-a},((M[a]^2-M[b]^2)M[a]\[CapitalGamma]0[b])/((M[a]^2-M[b]^2)^2+(M[a]\[CapitalGamma]0[b])^2)]
\[CapitalGamma]0[a:1|2]:=2M[a]A[a,a]
\[CapitalGamma][l_,a:1|2]:=M[a] 1/(16\[Pi]) H[l,a]Conjugate[H[l,a]]
\[CapitalGamma]C[l_,a:1|2]:=M[a] 1/(16\[Pi]) HC[l,a]Conjugate[HC[l,a]]

(*2.21*)
\[CapitalGamma]N[a:1|2]:=M[a]/(16\[Pi]) (HdH[a,a]+HCdHC[a,a])
\[Gamma]NL[a_]:=M[a]^3/(\[Pi]^2 z) BesselK[1,z]\[CapitalGamma]N[a]
\[Gamma]LL[k_,l_,a_,b_]:=(\[Gamma]NL[a]+\[Gamma]NL[b])/(1-2I (M[a]-M[b])/(\[CapitalGamma]N[a]+\[CapitalGamma]N[b]))*2(
Hs[l,a]H$[k,a]H[l,b]HC[k,b]+H$[l,a]Hs[k,a]HC[l,b]H[k,b]
)/(
HdH[a,a]+HCdHC[a,a]+HdH[b,b]+HCdHC[b,b]
)^2
\[Gamma]LLC[k_,l_,a_,b_]:=(\[Gamma]NL[a]+\[Gamma]NL[b])/(1-2I (M[a]-M[b])/(\[CapitalGamma]N[a]+\[CapitalGamma]N[b]))*2(
Hs[l,a]Hs[k,a]H[l,b]H[k,b]+H$[l,a]H$[k,a]HC[l,b]HC[k,b]
)/(
HdH[a,a]+HCdHC[a,a]+HdH[b,b]+HCdHC[b,b]
)^2
kappa[l_,a_,b_]:=(Sum[\[Gamma]LLC[l,k,a,b]+\[Gamma]LL[l,k,a,b],{k,3}]+\[Gamma]LLC[l,l,a,b]-\[Gamma]LL[l,l,a,b])/Sum[\[Gamma]NL[c]B[l,c],{c,2}]
B[l_,a_]:=(\[CapitalGamma][l,a]+\[CapitalGamma]C[l,a])/\[CapitalGamma]N[a]

kappa[3,1,1] // Simplify // Chop
%//. masses






