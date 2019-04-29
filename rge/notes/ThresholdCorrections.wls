#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Section:: *)
(*Threshold corrections for see-saw EFT matching*)


(* ::Text:: *)
(*We assume two heavy flavors which have similar Majorana masses (M2 ~ M1)*)


$Assumptions = {M2 > M1 > 0, \[Phi] > 0, v > 0};
definitions = {};
Dagger[matrix_List] := Transpose[Conjugate[matrix]]

V0 = -(msq/2) Hsq + \[Lambda] Hsq^2; (* notations following Brivio-Trott 1809.03450 *)
Hsq = (\[Phi]/Sqrt[2])^2;         (* \[Phi] follows the notation of DHM Appendix J *)
(* note that <\[Phi]> corresponds to Sqrt[2]v in DHM 0812.1594 Appendix J *)

MM = {{M1, 0}, {0, M2}};  (* majorana neutrino mass matrix *)
\[Omega]\[Omega] = Transpose[YY];       (* Neutrino Yukawa matrix in BT *)
YY = Table[Y[i, j], {i, 3}, {j, 2}];  (* Neurino Yukawa matrix in DHM/J *)
mD = v YY;                             (* DHM (J.2.20), where v = <\[Phi]>/Sqrt[2] *)

Y  /: Conjugate[Y[a_,b_]]  := Yc[a,b]
Yc /: Conjugate[Yc[a_,b_]] := Y[a,b]
Y  /: Abs[Y[a_,b_]] := Sqrt[Y[a,b]Yc[a,b]]

(* We define RR := Dagger[YY].YY and r1-r4 as follows *)
AppendTo[definitions, R[i_, j_] :> (Dagger[YY].YY)[[i,j]]];
rRules = {
  Subscript[r, 1] -> (\[Phi]^2/2) R[1, 1],
  Subscript[r, 2] -> (\[Phi]^2/2) R[2, 2],
  Subscript[r, 3] -> (\[Phi]^2/2) Re[R[1,2]],
  Subscript[r, 4] -> (\[Phi]^2/2) Im[R[1,2]]
};
(* Effective potential calculated in NeutrinoMassMatrix.nb *)
\[CapitalDelta]Vraw = ((3+4 c) M1^4+4 (1+2 c) M1^2 Subscript[r, 1]+(-2+4 c) \!\(\*SubsuperscriptBox[\(r\), \(1\), \(2\)]\)+4 (1+2 c) M1^2 Subscript[r, 2]+(-2+4 c) \!\(\*SubsuperscriptBox[\(r\), \(2\), \(2\)]\)+4 (-1+2 c) \!\(\*SubsuperscriptBox[\(r\), \(3\), \(2\)]\)+4 \!\(\*SubsuperscriptBox[\(r\), \(4\), \(2\)]\)+8 c \!\(\*SubsuperscriptBox[\(r\), \(4\), \(2\)]\))/(32 \[Pi]^2);
\[CapitalDelta]V = \[CapitalDelta]Vraw /. rRules


(* Matching with MSbar scheme (c=0) *)
Vlow = V0 /. {msq->msqEFT, \[Lambda]->\[Lambda]EFT};
Vhigh = V0 + Normal[\[CapitalDelta]V] /. {msq->msqFull, \[Lambda]->\[Lambda]Full};
-4*Coefficient[{Vlow, Vhigh}, \[Phi], 2] /. c->0 // FullSimplify
+4*Coefficient[{Vlow, Vhigh}, \[Phi], 4] /. c->0 // FullSimplify


\[CapitalDelta]\[Lambda] = -(1/(16\[Pi]^2))(Norm[\[Omega]1]^4+Norm[\[Omega]2]^4- 2 Abs[\[Omega]1.Conjugate[\[Omega]2]]^2+4 Re[\[Omega]1.Conjugate[\[Omega]2]]^2);
\[CapitalDelta]msq = (-M1^2/(4\[Pi]^2))(Norm[\[Omega]1]^2+Norm[\[Omega]2]^2);
-4Coefficient[Vhigh, \[Phi], 2] == msqFull + \[CapitalDelta]msq //. c->0 //. definitions // Simplify
4Coefficient[Vhigh, \[Phi], 4] == \[Lambda]Full + \[CapitalDelta]\[Lambda] //. c->0 //. definitions // FullSimplify


(* ::Section:: *)
(*Validation of BT paper*)


(* ::Text:: *)
(*If we perform this matching at Q = M Exp[-3/4], we should reproduce BT results (3.11) and (3.12).*)


{\[Omega]1, \[Omega]2} = \[Omega]\[Omega];
-4*Coefficient[{Vlow, Vhigh}, \[Phi], 2] /. c->-3/4 // Simplify
%[[2]] == msqFull + (M1^2/(8\[Pi]^2))(Norm[\[Omega]1]^2 + Norm[\[Omega]2]^2) //. definitions // Simplify

4*Coefficient[{Vlow, Vhigh}, \[Phi], 4] /. c->-3/4 // ExpandAll
%[[2]] == \[Lambda]Full - (5/(32\[Pi]^2))(Norm[\[Omega]1]^4+Norm[\[Omega]2]^4+6Abs[\[Omega]1.Conjugate[\[Omega]2]]^2 + 4 Re[\[Omega]1.\[Omega]2]^2)
% //. definitions //. {
  Abs[a_]^2:>a Conjugate[a],
  Re[a_]:>(a+Conjugate[a])/2,
  Im[a_]:>(a-Conjugate[a])/(2I)
} // ExpandAll // FullSimplify


(* ::Text:: *)
(*\[CapitalDelta]m^2 does match, but \[CapitalDelta]\[Lambda] does not. In fact, \[CapitalDelta]\[Lambda] with this matching scale is given by:*)


\[CapitalDelta]\[Lambda]BTcorrected = -(1/(32\[Pi]^2))(5Norm[\[Omega]1]^4+5Norm[\[Omega]2]^4+2Abs[\[Omega]1.Conjugate[\[Omega]2]]^2 + 8Re[\[Omega]1.Conjugate[\[Omega]2]]^2);
\[CapitalDelta]msqBT = (M1^2/(8\[Pi]^2))(Norm[\[Omega]1]^2+Norm[\[Omega]2]^2);

-4Coefficient[Vhigh, \[Phi], 2] == msqFull + \[CapitalDelta]msqBT //. c->-3/4 //. definitions // Simplify
4Coefficient[Vhigh, \[Phi], 4] == \[Lambda]Full + \[CapitalDelta]\[Lambda]BTcorrected //. c->-3/4 //. definitions // FullSimplify


(* ::Text:: *)
(*Note also that (3.9) and (3.11) are not in agreement.*)
