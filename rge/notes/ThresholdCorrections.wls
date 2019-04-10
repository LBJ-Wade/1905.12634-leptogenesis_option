#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Section:: *)
(*Threshold corrections for see-saw EFT matching*)


(* ::Text:: *)
(*We assume two heavy flavors which has similar Majorana masses (M2 ~ M1)*)


$Assumptions = {M2 > M1 > 0, \[Phi] > 0, v > 0};
definitions = {};
Dagger[matrix_List] := Transpose[Conjugate[matrix]]

V0 = -(msq/2) Hsq + \[Lambda] Hsq^2; (* notations following Brivio-Trott 1809.03450 *)
Hsq = (\[Phi]/Sqrt[2])^2;       (* \[Phi] follows the notation of DHM Appendix J *)
(* note that <\[Phi]> corresponds to Sqrt[2]v in DHM 0812.1594 Appendix J *)

MM = {{M1, 0}, {0, M2}};  (* majorana neutrino mass matrix *)
\[Omega]\[Omega] = Transpose[YY];       (* Neutrino Yukawa matrix in BT *)
YY = Table[Y[i, j], {i, 3}, {j, 2}];  (* Neurino Yukawa matrix in DHM/J *)
mD = v YY;                             (* DHM (J.2.20), where v = <\[Phi]>/Sqrt[2] *)

(* Then DHM (J.2.25) and (J.2.28) give the neutrino mass matrices before diagonalization *)
MnuL = -mD.Inverse[MM].Transpose[mD];
MnuH = MM + 1/2 (Inverse[MM].Dagger[mD].mD + Transpose[mD].Conjugate[mD].Inverse[MM]);
(* We define RR := Inverse[MM].Dagger[YY].YY *)
RR = Table[R[i, j], {i, 2}, {j, 2}];
AppendTo[definitions, R[i_, j_] :> (Inverse[MM].Dagger[YY].YY)[[i,j]]];
MnuH == MM + v^2(RR + Transpose[RR])/2 //. definitions //Simplify (* \[Rule] output "True" *)
MnuH = MM + v^2(RR + Transpose[RR])/2;

(* Diagonal elements of R are real. *)
# == Conjugate[#] &/@ Diagonal[RR/.definitions]//FullSimplify (* \[Rule] output "{True, True}" *)
AppendTo[$Assumptions, R[i_, i_] \[Element] Reals];


(* ::Text:: *)
(*With this set-up, the heavy neutrino masses-squared are given by:*)


mnuHsq = Simplify[Eigenvalues[MnuH.Dagger[MnuH]]] // Simplify;
(* expansion *)
mnuHsqApprox = Simplify[Series[mnuHsq, {v, 0, 4}]] //. {R[1,2]+R[2,1] :> rs, R[1,2]-R[2,1] :> rd} // Simplify
AppendTo[definitions, rs -> R[1,2]+R[2,1]];
AppendTo[definitions, rd -> R[1,2]-R[2,1]];


(* ::Text:: *)
(*So the effective potential difference \[CapitalDelta]V := V(full) - V(EFT) is given by, with \[Phi]-expansion,*)


\[CapitalDelta]V = Series[(-2/(64\[Pi]^2) (#^2 (Log[#/Q^2]-3/2)) &/@ Normal[mnuHsqApprox] // Total) /. v->\[Phi]/Sqrt[2], {\[Phi], 0, 4}]


(* ::Text:: *)
(*The threshold correction is given by, where the first (second) elements are couplings for EFT (full theory); they should be match at, with MSbar matching, Q=M.*)


Vlow = V0 /. {msq->msqEFT, \[Lambda]->\[Lambda]EFT};
Vhigh = V0 + Normal[\[CapitalDelta]V] /. {msq->msqFull, \[Lambda]->\[Lambda]Full};
Limit[-4*Coefficient[{Vlow, Vhigh}, \[Phi], 2], {M2->M1,Q->M1}]
Limit[+4*Coefficient[{Vlow, Vhigh}, \[Phi], 4], {M2->M1,Q->M1}]


(* ::Section:: *)
(*Validation of BT paper*)


(* ::Text:: *)
(*If we perform this matching at Q = M Exp[-3/4], we should reproduce BT results (3.11) and (3.12).*)


{\[Omega]1, \[Omega]2} = \[Omega]\[Omega];
Limit[-4*Coefficient[{Vlow, Vhigh}, \[Phi], 2],{M2->M1,Q->M1 Exp[-3/4]}] /. definitions /. M2->M1 // FullSimplify
%[[2]] == msqFull + (M1^2/(8\[Pi]^2))(Norm[\[Omega]1]^2 + Norm[\[Omega]2]^2)//FullSimplify

Limit[4*Coefficient[{Vlow, Vhigh}, \[Phi], 4],{M2->M1,Q->M1 Exp[-3/4]}]
%[[2]] == \[Lambda]Full - (7/(32\[Pi]^2))(Norm[\[Omega]1]^4 + Norm[\[Omega]2]^4 + 2(Re[\[Omega]1.Conjugate[\[Omega]2]])^2) //. definitions /. {M2->M1} // FullSimplify


(* ::Text:: *)
(*\[CapitalDelta]m^2 agrees, but \[CapitalDelta]\[Lambda] does not.*)
