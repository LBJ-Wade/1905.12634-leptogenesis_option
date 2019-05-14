(* ::Package:: *)

(* ::Subsubsection:: *)
(*Packages*)


SetDirectory[Which[$FrontEnd=!=Null,NotebookDirectory[],$InputFileName=!="", DirectoryName[$InputFileName], True, "."]];
$ProjectRoot = RunProcess[{"git","rev-parse" , "--show-toplevel"}]["StandardOutput"] // StringTrim;

Needs["ResLept14041003`", "1404.1003.m"];
Needs["NuFIT`", "nufit40.m"];
Get["plot_tools.wl"];
msq = Get[FileNameJoin[{$ProjectRoot, "rge", "msq.dat"}]];

$Assumptions = {};

AddAssumptions[exp__] := (AppendTo[$Assumptions, #] &/@ {exp}; DeleteDuplicates[$Assumptions])
Validate[rules_, revertRules_] := rules //. Rule->Equal //. revertRules // FullSimplify

Unprotect[ConjugateTranspose];
ConjugateTranspose[ConjugateTranspose[a_]] := a;
ConjugateTranspose[Dot[a_, b__]] := ConjugateTranspose[Dot[b]].ConjugateTranspose[a];
ConjugateTranspose[Times[a___, x_, b___]] := Conjugate[x] ConjugateTranspose[Times[a,b]] /; Or[
  NumericQ[x],
  MatchQ[x, _rM | _rm | Subscript[M, _] | Subscript[m, _] | vev | vev^_Integer]];
Dagger = ConjugateTranspose;

Unprotect[Dot];
Dot[a___, Times[b___, x_, c___], d___] := x Dot[a, Times[b, c], d] /; Or[
  NumericQ[x],
  MatchQ[x, _rM | _rm | Subscript[M, _] | Subscript[m, _] | vev | vev^_Integer]];  


(* ::Subsubsection:: *)
(*Representations*)


RRep := Module[{z=w + I x}, Which[
  hierarchy==="NH", {{0, Cos[z], \[Zeta] Sin[z]}, {0, -Sin[z], \[Zeta] Cos[z]}},
  hierarchy==="IH", {{Cos[z], \[Zeta] Sin[z], 0}, {-Sin[z], \[Zeta] Cos[z], 0}}]];
  
URep = {
  { c[12] c[13],  s[12] c[13], s[13] Conjugate[e]},
  {-s[12] c[23] - c[12] s[23] s[13] e,  c[12] c[23] - s[12] s[23] s[13] e, s[23] c[13]},
  { s[12] s[23] - c[12] c[23] s[13] e, -c[12] s[23] - s[12] c[23] s[13] e, c[23] c[13]}
}.DiagonalMatrix[{1,Exp[I \[Sigma]],1}]//.{
  e->Exp[I \[Delta]], c[a_]:>Subscript[co,a], s[a_]:>Subscript[si,a]
};
U /: ConjugateTranspose[U].U := DiagonalMatrix[{1,1,1}];
U /: U.ConjugateTranspose[U] := DiagonalMatrix[{1,1,1}];

AddAssumptions[w\[Element]Reals, x\[Element]Reals, \[Zeta]\[Element]Reals, \[Zeta]^2==1, \[Sigma]\[Element]Reals, \[Delta]\[Element]Reals];
Do[AddAssumptions[-1<=Subscript[co, a]<=1, -1<=Subscript[si, a]<=1, Subscript[co, a]^2+Subscript[si, a]^2==1], {a, {12, 23, 13}}]


MM = DiagonalMatrix[{Subscript[M, 1], Subscript[M, 2]}];
rMM = DiagonalMatrix[{rM[1], rM[2]}];
mm = DiagonalMatrix[{Subscript[m, 1], Subscript[m, 2], Subscript[m, 3]}];
rmm = DiagonalMatrix[{rm[1], rm[2],rm[3]}];
rM /: rM[a_]^n_ /; n >= 2    := Subscript[M, a] rM[a]^(n-2)
rM /: rM[a_]^n_ /; n <= -2  := rM[a]^(n+2)/Subscript[M, a]
rm /: rm[a_]^n_ /; n >= 2    := Subscript[m, a] rm[a]^(n-2)
rm /: rm[a_]^n_ /; n <= -2  := rm[a]^(n+2)/Subscript[m, a]

y := I(Sqrt[2]/vev)rMM.R.rmm.Dagger[U];
AddAssumptions[vev>0, Subscript[m, 1]>=0, Subscript[m, 2]>0, Subscript[m, 3]>=0, Subscript[M, 2]>Subscript[M, 1]>0, rm[1]>=0, rm[2]>0, rm[3]>=0, rM[2]>rM[1]>0];


MRules1 = {Subscript[M, 2] -> (1+\[Delta]M)Subscript[M, 1]};
MRules2 = {Subscript[M, 2] -> Subscript[M, 1] Sqrt[1+\[Delta]M2]};
MRevert = {\[Delta]M -> (Subscript[M, 2]-Subscript[M, 1])/Subscript[M, 1], \[Delta]M2 -> (Subscript[M, 2]^2-Subscript[M, 1]^2)/Subscript[M, 1]^2};
AddAssumptions[\[Delta]M>0, \[Delta]M2>0];

mRules := Which[hierarchy==="NH", {Subscript[m, 2]|mL->mtot (1-\[Rho]m)/2, Subscript[m, 3]|mH->mtot (1+\[Rho]m)/2},
                hierarchy==="IH", {Subscript[m, 1]|mL->mtot (1-\[Rho]m)/2, Subscript[m, 2]|mH->mtot (1+\[Rho]m)/2}];
mRevert := {\[Rho]m->(mH-mL)/(mH+mL), mtot->mH+mL};
mValues := Which[hierarchy==="NH", {Subscript[m, 2]|mL->ToMasses[BestFit["NH"]][[2]], Subscript[m, 3]|mH->ToMasses[BestFit["NH"]][[3]]},
                 hierarchy==="IH", {Subscript[m, 1]|mL->ToMasses[BestFit["IH"]][[1]], Subscript[m, 2]|mH->ToMasses[BestFit["IH"]][[2]]}];
AddAssumptions[mtot>0, 0<\[Rho]m<1, 0<mL<mH];


(* ::Input:: *)
(*Validate[{MRules1,MRules2}, MRevert]*)
(*hierarchy="NH"; {Validate[mRules, mRevert],{\[Rho]m, mtot} //. mRevert //. mValues}*)
(*hierarchy="IH"; {Validate[mRules, mRevert],{\[Rho]m, mtot} //. mRevert //. mValues}*)


W[1, 1] := Cosh[2x] - \[Rho]m Cos[2w];
W[2, 2] := Cosh[2x] + \[Rho]m Cos[2w];
W[1, 2] := +I Sinh[2x] + \[Rho]m Sin[2w];
W[2, 1] := -I Sinh[2x] + \[Rho]m Sin[2w];
\[Mu][i:1|2] := (mtot Subscript[M, i])/(8\[Pi] vev^2);
\[CapitalGamma][i:1|2] := HoldForm[\[Mu]][i]HoldForm[W][i,i]Subscript[M, i];
yydag[i:1|2, j:1|2] := (mtot rM[i]rM[j])/vev^2 W[i,j]


(* ::Input:: *)
(*Do[*)
(*hierarchy=hierarchy$tmp;*)
(*\[Mu][1]//.mRevert//.{vev->246}//.mValues//Print;*)
(*FullSimplify[(y.Dagger[y]==Table[yydag[a,b],{a,2},{b,2}])//.R->RRep/.a_Conjugate:>FullSimplify[a]]//.mRules//ReleaseHold//Simplify//Print;*)
(*With[{Gamma=Function[i, (yydag[i,i])Subscript[M, i]/(8\[Pi])]},*)
(*{Gamma[1]==\[CapitalGamma][1],Gamma[2]==\[CapitalGamma][2]}]//ReleaseHold//Simplify//Print;*)
(*,{hierarchy$tmp,{"NH","IH"}}];*)


(* ::Section:: *)
(*Neutrino-option condition*)


RHS = (vev^2/2) (TrYY/M1);
fneutopt[M_] := 8\[Pi]^2 vev^2 \[Mu]sq[M Exp[-3/4]] /M^3;
fneutopt[M1] == RHS // Solve[#, \[Mu]sq[_]]&
RHS //. {M1->Subscript[M, 1], TrYY->yydag[1,1]+yydag[2,2]} //. MRules1  // Expand


Plot[{
    fneutopt[M1] /. {\[Mu]sq[q_]:>msq[q]/2, vev->246},
    Total[ToMasses[BestFit["NH"]]],
    Total[ToMasses[BestFit["IH"]]]} *10^9/. M1->M1*10^6 // Evaluate,  {M1, 6, 14},
    FrameLabel->{"\!\(\*SubscriptBox[\(M\), \(1\)]\) /PeV", None},
    PlotLegends->{"f(\!\(\*SubscriptBox[\(M\), \(1\)]\)) /eV", "\!\(\*SubscriptBox[\(m\), \(tot\)]\) /eV (IH)",  "\!\(\*SubscriptBox[\(m\), \(\(tot\)\(\\\ \)\)]\)/eV (NH)"}, PlotRange->{{6,14},{0.03,0.13}}
    ]


LogLogPlot[{
      (fneutopt[M1] /. {\[Mu]sq[q_]:>msq[q]/2, vev->246}) / Total[ToMasses[BestFit["NH"]]],
      (fneutopt[M1] /. {\[Mu]sq[q_]:>msq[q]/2, vev->246}) / Total[ToMasses[BestFit["IH"]]]
    },  {M1, 5*^5, 10^7},
    FrameLabel->{"\!\(\*SubscriptBox[\(M\), \(1\)]\) /GeV", None},
    PlotLegends->{"f(\!\(\*SubscriptBox[\(M\), \(1\)]\))/\[Sum] \!\(\*SubscriptBox[\(m\), \(i\)]\) (NH)", "f(\!\(\*SubscriptBox[\(M\), \(1\)]\))/\[Sum] \!\(\*SubscriptBox[\(m\), \(i\)]\) (IH)"},
    PlotRange->{{8*^5, 10^7}, {1, 1000}}
    ]
    


tm[i:1|2] := vev^2/2 yydag[i,i]/Subscript[M, i]
tm[1]
tm[2]


(fneutopt[M1] /. {\[Mu]sq[q_]:>msq[q]/2, vev->246}) / Total[ToMasses[BestFit["NH"]]]//.M1->2*^5
