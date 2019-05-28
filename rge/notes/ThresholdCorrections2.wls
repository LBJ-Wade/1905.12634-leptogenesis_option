#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Section:: *)
(*Threshold corrections for see-saw EFT matching / Another method*)


$Assumptions = {M2 > M1 > 0, \[Phi] > 0, v > 0, \[Epsilon] > 0, \[Phi]2sq>0};
definitions = <||>;
Dag = ConjugateTranspose;

(* all notations following Dev+ 1611.03827 *)
V0 = - msq/2 Hsq + \[Lambda] Hsq^2;  (* Higgs potential *)
Hsq = \[Phi]2sq;

MM = {{M1, 0}, {0, M2}};  (* majorana neutrino mass matrix *)
\[Omega]\[Omega] = Transpose[YY];       (* Neutrino Yukawa matrix in BT *)
YY = Table[Y[i, j], {i, 2}, {j, 3}];
mD = \[Phi]2 YY;

Y  /: Conjugate[Y[a_,b_]]  := Yc[a,b]
Yc /: Conjugate[Yc[a_,b_]] := Y[a,b]
Y  /: Abs[Y[a_,b_]] := Sqrt[Y[a,b]Yc[a,b]]

MN = ArrayFlatten[{{0, Transpose[mD]}, {mD, MM}}] (*Full mass matrix*)


MN\[Epsilon] = MN/.\[Phi]2->\[Epsilon] \[Phi]2;
eigen = Simplify[Eigenvalues[MN\[Epsilon].Dag[MN\[Epsilon]]]]//. Conjugate[\[Phi]2]->\[Phi]2sq/\[Phi]2;
equation = Cases[eigen,Root[__]][[All,1]] // DeleteDuplicates;
Length[equation]
feq = (equation[[1]])[x]/.\[Epsilon]->1;
Collect[feq, {x, \[Phi]2sq, M1, M2}, Simplify] //. Y[a_,1]Yc[b_,1]+Y[a_,2]Yc[b_,2]+Y[a_,3]Yc[b_,3] :> YYD[a,b]


feq2 = x^4+\[Phi]2sq^4 (R1)^2+x^3 (-M1^2-M2^2-2 \[Phi]2sq (YYD[1,1]+YYD[2,2]))+x^2 (M1^2 M2^2+\[Phi]2sq^2 (R2)+\[Phi]2sq 2(M2^2 YYD[1,1]+ M1^2 YYD[2,2]))+x (\[Phi]2sq^2 (M1 M2 (R3)-M2^2 YYD[1,1]^2-M1^2 YYD[2,2]^2)-2 \[Phi]2sq^3 (R4));
definitions[R1] = Y[1,1] (Y[2,2] (-Yc[1,2] Yc[2,1]+Yc[1,1] Yc[2,2])+Y[2,3] (-Yc[1,3] Yc[2,1]+Yc[1,1] Yc[2,3]))+Y[1,3] (Y[2,1] (Yc[1,3] Yc[2,1]-Yc[1,1] Yc[2,3])+Y[2,2] (Yc[1,3] Yc[2,2]-Yc[1,2] Yc[2,3]))+Y[1,2] (Y[2,1] (Yc[1,2] Yc[2,1]-Yc[1,1] Yc[2,2])+Y[2,3] (-Yc[1,3] Yc[2,2]+Yc[1,2] Yc[2,3]));
definitions[R2] = Y[1,1]^2 Yc[1,1]^2+Y[1,2]^2 Yc[1,2]^2+Y[1,3]^2 Yc[1,3]^2+4 Y[1,3] Y[2,1] Yc[1,3] Yc[2,1]+Y[2,1]^2 Yc[2,1]^2+4 Y[1,3] Y[2,2] Yc[1,3] Yc[2,2]+2 Y[2,1] Y[2,2] Yc[2,1] Yc[2,2]+Y[2,2]^2 Yc[2,2]^2-2 Y[1,3] Y[2,1] Yc[1,1] Yc[2,3]-2 Y[1,3] Y[2,2] Yc[1,2] Yc[2,3]+2 Y[1,3] Y[2,3] Yc[1,3] Yc[2,3]+2 Y[2,1] Y[2,3] Yc[2,1] Yc[2,3]+2 Y[2,2] Y[2,3] Yc[2,2] Yc[2,3]+Y[2,3]^2 Yc[2,3]^2+2 Y[1,1] (Y[1,2] Yc[1,1] Yc[1,2]+Y[1,3] Yc[1,1] Yc[1,3]+Y[2,1] Yc[1,1] Yc[2,1]-Y[2,2] Yc[1,2] Yc[2,1]-Y[2,3] Yc[1,3] Yc[2,1]+2 Y[2,2] Yc[1,1] Yc[2,2]+2 Y[2,3] Yc[1,1] Yc[2,3])+2 Y[1,2] (Y[1,3] Yc[1,2] Yc[1,3]+Y[2,2] Yc[1,2] Yc[2,2]-Y[2,3] Yc[1,3] Yc[2,2]+Y[2,1] (2 Yc[1,2] Yc[2,1]-Yc[1,1] Yc[2,2])+2 Y[2,3] Yc[1,2] Yc[2,3]);
definitions[R3] = -Y[2,1]^2 Yc[1,1]^2-Y[2,2]^2 Yc[1,2]^2-2 Y[2,2] Y[2,3] Yc[1,2] Yc[1,3]-Y[2,3]^2 Yc[1,3]^2-2 Y[2,1] Yc[1,1] (Y[2,2] Yc[1,2]+Y[2,3] Yc[1,3])-Y[1,1]^2 Yc[2,1]^2-2 Y[1,1] Y[1,2] Yc[2,1] Yc[2,2]-Y[1,2]^2 Yc[2,2]^2-2 Y[1,1] Y[1,3] Yc[2,1] Yc[2,3]-2 Y[1,2] Y[1,3] Yc[2,2] Yc[2,3]-Y[1,3]^2 Yc[2,3]^2;
definitions[R4] = Y[1,1]^2 Yc[1,1] (Y[2,2] (-Yc[1,2] Yc[2,1]+Yc[1,1] Yc[2,2])+Y[2,3] (-Yc[1,3] Yc[2,1]+Yc[1,1] Yc[2,3]))+(Y[1,3] (Y[2,1] (Yc[1,3] Yc[2,1]-Yc[1,1] Yc[2,3])+Y[2,2] (Yc[1,3] Yc[2,2]-Yc[1,2] Yc[2,3]))+Y[1,2] (Y[2,1] (Yc[1,2] Yc[2,1]-Yc[1,1] Yc[2,2])+Y[2,3] (-Yc[1,3] Yc[2,2]+Yc[1,2] Yc[2,3]))) (Y[1,2] Yc[1,2]+Y[1,3] Yc[1,3]+YYD[2,2])-Y[1,1] (Y[1,2] (Y[2,2] Yc[1,2] (Yc[1,2] Yc[2,1]-Yc[1,1] Yc[2,2])+Y[2,1] Yc[1,1] (-Yc[1,2] Yc[2,1]+Yc[1,1] Yc[2,2])+Y[2,3] (Yc[1,1] Yc[1,3] Yc[2,2]+Yc[1,2] (Yc[1,3] Yc[2,1]-2 Yc[1,1] Yc[2,3])))+Y[1,3] (Y[2,3] Yc[1,3] (Yc[1,3] Yc[2,1]-Yc[1,1] Yc[2,3])+Y[2,1] Yc[1,1] (-Yc[1,3] Yc[2,1]+Yc[1,1] Yc[2,3])+Y[2,2] (-2 Yc[1,1] Yc[1,3] Yc[2,2]+Yc[1,2] (Yc[1,3] Yc[2,1]+Yc[1,1] Yc[2,3])))+(Y[2,2] (Yc[1,2] Yc[2,1]-Yc[1,1] Yc[2,2])+Y[2,3] (Yc[1,3] Yc[2,1]-Yc[1,1] Yc[2,3])) YYD[2,2]);


feq - feq2 //. definitions //. {YYD[a_,b_]:>(YY.Dag[YY])[[a,b]]} // Simplify


Table[Root[(feq2/.x->#)&, i], {i, 4}];
masses = Series[%, {\[Phi]2sq, 0, 2}]


\[CapitalDelta]V = -2/(64\[Pi]^2) #^2 (Log[#/Q2]-3/2) &/@ masses[[3;;4]] // Total;
Coefficient[Normal[\[CapitalDelta]V],\[Phi]2sq,2] //. {Log[M1^2/Q2]->L1,Log[M2^2/Q2]->L2} // FullSimplify // Collect[#, {L1, L2}, FullSimplify]&


(* ::Subsection:: *)
(*Compare with old result*)


\[CapitalDelta]Vraw = ((3+4 c) M1^4+4 (1+2 c) M1^2 Subscript[r, 1]+(-2+4 c) \!\(\*SubsuperscriptBox[\(r\), \(1\), \(2\)]\)+4 (1+2 c) M1^2 Subscript[r, 2]+(-2+4 c) \!\(\*SubsuperscriptBox[\(r\), \(2\), \(2\)]\)+4 (-1+2 c) \!\(\*SubsuperscriptBox[\(r\), \(3\), \(2\)]\)+4 \!\(\*SubsuperscriptBox[\(r\), \(4\), \(2\)]\)+8 c \!\(\*SubsuperscriptBox[\(r\), \(4\), \(2\)]\))/(32 \[Pi]^2) //. {
  Subscript[r, 1] -> \[Epsilon] YYD[1, 1],
  Subscript[r, 2] -> \[Epsilon] YYD[2, 2],
  Subscript[r, 3] -> \[Epsilon] Re[YYD[1,2]],
  Subscript[r, 4] -> \[Epsilon] Im[YYD[1,2]],
  \[Epsilon]->\[Phi]2sq
};
FullSimplify[(\[CapitalDelta]Vraw/.c->0)==Normal[Series[\[CapitalDelta]V //. {Q2 -> Q^2, Q->M1}, {M2, M1, 0}]]] //. definitions //. {YYD[a_,b_]:>(YY.Dag[YY])[[a,b]]} // FullSimplify


(* ::Subsection:: *)
(*Analytic Expression*)


R3-YYD[1,1]^2-YYD[2,2]^2 == -Tr[YY.Dag[YY].Conjugate[YY].Transpose[YY]] //. definitions  //. YYD[a_,b_]:>(YY.Dag[YY])[[a,b]] // FullSimplify
R2-2 (YYD[1,1]+YYD[2,2])^2 == -Tr[YY.Dag[YY].YY.Dag[YY]] //. definitions  //. YYD[a_,b_]:>(YY.Dag[YY])[[a,b]] // Expand


(* ::Text:: *)
(*For M2-M1 series expansion*)


\[CapitalDelta]Vdegenerate = \[CapitalDelta]V //. {
  R3 -> -TrYDCT + YYD[1,1]^2 + YYD[2,2]^2,
  R2 -> -TrYDYD +2 (YYD[1,1]+YYD[2,2])^2,
  YYD[1,1]->H1,
  YYD[2,2]->H2
} // FullSimplify
Series[\[CapitalDelta]Vdegenerate,{M2, M1, 0}]//. Log[M1^2/Q2]->L1 //FullSimplify


(* ::Text:: *)
(*Without expansion:*)


\[CapitalDelta]Vwoe= \[CapitalDelta]V //. {
  Log[M1^2/Q2] -> L1, Log[M2^2/Q2] -> L2  (* do not use this replacement in M2-M1 series expansion! *)
} //. {
  R3 -> -TrYDCT + YYD[1,1]^2 + YYD[2,2]^2,
  R2 -> -TrYDYD +2 (YYD[1,1]+YYD[2,2])^2,
  YYD[1,1]->H1,
  YYD[2,2]->H2
} // FullSimplify;
Collect[-8\[Pi]^2 Coefficient[\[CapitalDelta]Vwoe, \[Phi]2sq, 1],{TrYDCT, TRYDYD}, FullSimplify]
Collect[-16\[Pi]^2  Coefficient[\[CapitalDelta]Vwoe, \[Phi]2sq, 2],{TrYDCT, TRYDYD, H1, H2}, FullSimplify]


(* ::Text:: *)
(*Comparison with Expressions in our paper*)


\[CapitalDelta]Vours = Simplify[\[CapitalDelta]Vwoe //. {
   L1 -> Log[M1sq/QMSbarsq],
   L2 -> Log[M2sq/QMSbarsq],
   QMSbarsq->Qsq Exp[-3/2]
}] //. { Log[M1sq/Qsq]->lnMsqQsq[1], Log[M2sq/Qsq]->lnMsqQsq[2] };
muexpression = Collect[-8\[Pi]^2 Coefficient[\[CapitalDelta]Vours, \[Phi]2sq, 1],{TrYDCT, TRYDYD}, FullSimplify]
lamexpression = Collect[+8\[Pi]^2 Coefficient[\[CapitalDelta]Vours, \[Phi]2sq, 2],{TrYDCT, TRYDYD, H1, H2}, FullSimplify]
(lamexpression == Sum[H[i]^2 L1[i,j]-1/2 TrYDCT L2[i,j]+1/2 TrYDYD L3[i,j] /. j->3-i,{i,1,2}]) //.{
   L1[i_,j_]:>(M[j]-M[i])^-1 M[j] Log[M[j]/M[i]] - 1,
   L2[i_,j_]:>(M[j]^2-M[i]^2)^-1 M[i]M[j] Log[M[j]/M[i]],
   L3[i_,j_]:>(M[j]^2-M[i]^2)^-1 M[i]^2 lnMsqQsq[i] - 1/4,
   M[1]->M1, M[2]->M2, H[1]->H1, H[2]->H2,
   Log[M2/M1]->(lnMsqQsq[2]-lnMsqQsq[1])/2,
   Log[M1/M2]->(lnMsqQsq[1]-lnMsqQsq[2])/2
   } // FullSimplify//Expand
lamexpression /. {lnMsqQsq[2]->Log[M2^2/Qsq], lnMsqQsq[1]->Log[M1^2/Qsq]} // Limit[#, M2->M1]& // Expand



