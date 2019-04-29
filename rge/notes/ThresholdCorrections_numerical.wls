#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Section:: *)
(*Threshold corrections for see-saw EFT matching: Numerical analysis*)


(* ::Text:: *)
(*We assume two heavy flavors which have similar Majorana masses (M2 ~ M1).*)


(* ::Subsubsection:: *)
(*Set up*)


(* Neutrino Yukawa Matrix *)
YY = Table[Y[i, j], {i, 3}, {j, 2}];
(* 5x5 mass matrix *)
MN = ArrayFlatten[{{Table[0, 3, 3], YY \[Phi]/Sqrt[2]}, {Transpose[YY]\[Phi]/Sqrt[2], {{M1, 0}, {0, M2}}}}];

Dag = ConjugateTranspose;

(* Numerical example of Yukawa Matrix and majorana mass *)
numbers = Flatten[{Table[Y[i, j]->RandomReal[] +  I RandomReal[], {i, 3}, {j, 2}], M1->2000, M2->2000}];


(* ::Subsubsection:: *)
(*Mass calculation*)


(* neutrino mass (descending) with functions of \[Phi] *)
masssq[\[Phi]in_] := Eigenvalues[Dag[MN].MN //. numbers //. \[Phi]->\[Phi]in] // Sort // Reverse


(* \[Phi]-fixed; Method 1: SVD *)
SingularValueDecomposition[MN //. numbers /. {\[Phi]->246}][[2]] // Diagonal
(* \[Phi]-fixed; Method 2: Eigenvalue *)
Eigenvalues[Dag[MN].MN //. numbers /. {\[Phi]->246}] // Sqrt
(* This line verifies that masssq gives accurate results at least for first two. *)
masssq[246] // Sqrt


(* ::Subsubsection:: *)
(*Definitions of Effective potential*)


(* Full Effective Potential *)
\[CapitalDelta]V[\[Phi]in_Real, Qin_] := Module[{mNsq = masssq[\[Phi]in]}, Sum[-2/(64\[Pi]^2) msq^2 (Log[msq / Qin^2]-3/2), {msq, {mNsq[[1]], mNsq[[2]]}}]]

(* My expression expanded with \[Phi] *)
\[CapitalDelta]m2a := (-M1^2/(4\[Pi]^2))(Norm[\[Omega]1]^2+Norm[\[Omega]2]^2) //. {\[Omega]1->YY[[All,1]], \[Omega]2->YY[[All,2]]} //. numbers;
\[CapitalDelta]\[Lambda]a := -(1/(16\[Pi]^2))Sum[(\[Omega]p.Conjugate[\[Omega]q])(\[Omega]p.Conjugate[\[Omega]q]), {\[Omega]p, {\[Omega]1, \[Omega]2}}, {\[Omega]q, {\[Omega]1, \[Omega]2}}]  //. {\[Omega]1->YY[[All,1]], \[Omega]2->YY[[All,2]]} //. numbers // Chop;

(* Brivio-Trott expression expanded with \[Phi] *)
\[CapitalDelta]m2b := \[CapitalDelta]m2a
\[CapitalDelta]\[Lambda]b := -(1/(16\[Pi]^2))Sum[(\[Omega]p.Conjugate[\[Omega]q])(\[Omega]q.Conjugate[\[Omega]p]), {\[Omega]p, {\[Omega]1, \[Omega]2}}, {\[Omega]q, {\[Omega]1, \[Omega]2}}]  //. {\[Omega]1->YY[[All,1]], \[Omega]2->YY[[All,2]]} //. numbers // Chop;

\[CapitalDelta]Va[\[Phi]1_Real, Q_] := Module[{\[Phi]0=246}, \[CapitalDelta]V[\[Phi]0, Q] + (-\[CapitalDelta]m2a/4) (\[Phi]1^2-\[Phi]0^2) + (\[CapitalDelta]\[Lambda]a / 4) (\[Phi]1^4-\[Phi]0^4)]
\[CapitalDelta]Vb[\[Phi]1_Real, Q_] := Module[{\[Phi]0=246}, \[CapitalDelta]V[\[Phi]0, Q] + (-\[CapitalDelta]m2b/4) (\[Phi]1^2-\[Phi]0^2) + (\[CapitalDelta]\[Lambda]b / 4) (\[Phi]1^4-\[Phi]0^4)]


(* ::Subsubsection:: *)
(*Comparison*)


Plot[{\[CapitalDelta]Va[\[Phi]in, 2000]/\[CapitalDelta]V[\[Phi]in, 2000], \[CapitalDelta]Vb[\[Phi]in, 2000]/\[CapitalDelta]V[\[Phi]in, 2000]}, {\[Phi]in, 230, 250},
  PlotLegends->{"\[CapitalDelta]Va/\[CapitalDelta]V", "\[CapitalDelta]Vb/\[CapitalDelta]V"},
  AxesLabel->{"\[Phi]", None}]


(* ::Text:: *)
(*You may evaluate the next block repeatedly to see more examples.*)
(*Always my result is more accurate.*)
(*Note the difference from unity comes from \[Phi]^6 terms etc.*)


numbers = Flatten[{Table[Y[i, j]->RandomReal[] +  I RandomReal[], {i, 3}, {j, 2}], M1->2000, M2->2000}];
{\[CapitalDelta]\[Lambda]a, \[CapitalDelta]\[Lambda]b}
Plot[{\[CapitalDelta]Va[\[Phi]in, 2000]/\[CapitalDelta]V[\[Phi]in, 2000], \[CapitalDelta]Vb[\[Phi]in, 2000]/\[CapitalDelta]V[\[Phi]in, 2000]}, {\[Phi]in, 230, 250},
  PlotLegends->{"\[CapitalDelta]Va/\[CapitalDelta]V", "\[CapitalDelta]Vb/\[CapitalDelta]V"},
  AxesLabel->{"\[Phi]", None}]
Plot[{\[CapitalDelta]V[\[Phi]in, 2000], \[CapitalDelta]Va[\[Phi]in, 2000], \[CapitalDelta]Vb[\[Phi]in, 2000]}, {\[Phi]in, 241, 251},
  PlotLegends->{"\[CapitalDelta]V", "\[CapitalDelta]Va", "\[CapitalDelta]Vb"},
  AxesLabel->{"\[Phi]", None}]


(* ::Text:: *)
(*You may want to enlarge the second plot...*)
