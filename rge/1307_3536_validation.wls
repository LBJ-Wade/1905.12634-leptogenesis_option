#!wolframscript -script
(* ::Package:: *)

(* ::Section:: *)
(*Validation of 1307_3536.wl*)


(* ::Input::Initialization:: *)
SetDirectory[NotebookDirectory[]];
$ProjectRoot = RunProcess[{"git","rev-parse" , "--show-toplevel"}]["StandardOutput"] // StringTrim;
$Path = Append[$Path, FileNameJoin[{$ProjectRoot, "vendor", "LoopTools", "build"}]] // DeleteDuplicates;


(* ::Input::Initialization:: *)
Install["LoopTools"];
Needs["SMRGE`", "1307_3536.wl"];
<<plot_tools.wl


(* ::Subsection:: *)
(*Table 3*)


table2=<|
  "MW"->80.384,
  "MZ"->91.1876,
  "MH"->125.15,
  "MT"->173.34,
  "Gmu"->1/(v^2 Sqrt[2]) /. {v->246.21971},
  "as(MZ)"->0.1184,
  "MB"->4.2,
  "MTau"->1.777|>;

table3[order_]:={
  #[lam],
  Sqrt[#[ytsq]],
  Sqrt[#[g2sq]],
  Sqrt[#[g1sq]*3/5],
  Sqrt[#[msq]]
} &@ getCouplingsAtMt[table2, order]

table3[0]
table3[1]
table3[2]


table2aNeutOpt = <|
        "MW" -> 80.387,      "MZ" -> 91.1875,   "MH" -> 125.09,
        "MT" -> 173.2,       "MB" -> 4.18,      "MTau" -> 1.776,
        "as(MZ)" -> 0.1185,  "Gmu" -> 1.1663787*^-5 |>;
table2b[order_] := {
#[lam], Sqrt[#[msq]], Sqrt[#[g1sq]], Sqrt[#[g2sq]], Sqrt[#[g3sq]], Sqrt[#[ytsq]], Sqrt[#[ybsq]], Sqrt[#[yasq]]
    } &@ getCouplingsAtMt[table2aNeutOpt, order]


table2b[0]
table2b[1]
table2b[2]


(* ::Subsection:: *)
(*RGE validation against 1212.6829*)


<<"./refs/1212.6829/betayu.m";
betaCompare[rge12126829_] := Normal[Series[rge12126829, {h, 0, 4}]] * FourPi^2 /. {
  lam -> lam/FourPi^2
} //. {
  a1->g1sq/FourPi^2,
  a2->g2sq/FourPi^2,
  as->g3sq/FourPi^2,
  at->ytsq/FourPi^2,
  ab->ybsq/FourPi^2,
  atau->yasq/FourPi^2,
  h^6->0,
  h->1,
  NG->3
};
RGE[yasq] == betaCompare[batau]//Simplify
RGE[ybsq] == betaCompare[bab]//Simplify
RGE[ytsq, 2] == betaCompare[bat]//Simplify


(* ::Subsection:: *)
(*RGE plot*)


(* Returns {RGEs, params, {t, t_lowEnergy, t_highEnergy}} *)
getRGEs[couplingsAtMt_, mt_, mHigh_, order_] := Module[{
    t,
    t0 = Log[mt^2//N],
    t1 = Log[mHigh^2//N],
    couplingsToFunctions, boundary, diffeq
  },
  couplingsToFunctions = KeyValueMap[#1->#1[t]&, couplingsAtMt];
  boundary = KeyValueMap[#1[t0]==#2&, couplingsAtMt];
  diffeq = KeyValueMap[D[#1[t], t]==(RGE[#1, order] /. FourPi->4\[Pi] /. couplingsToFunctions)&, couplingsAtMt];
  Return[{Flatten[{boundary, diffeq}], couplings, {t, t0, t1}}]
]


rge1 = NDSolve@@getRGEs[getCouplingsAtMt[table2aNeutOpt, 0], table2aNeutOpt["MT"], 10^20, 1];
rge2 = NDSolve@@getRGEs[getCouplingsAtMt[table2aNeutOpt, 1], table2aNeutOpt["MT"], 10^20, 2];
rge3 = NDSolve@@getRGEs[getCouplingsAtMt[table2aNeutOpt, 2], table2aNeutOpt["MT"], 10^20, All];
LogLinearPlot[lam[Log[q^2]]/.{rge1, rge2, rge3}//Evaluate, {q, table2aNeutOpt["MT"], 10^19}, PlotRange->{{10^2, 10^20}, {-0.10, 0.15}}]


LogLinearPlot[{
  Sqrt[g1sq[Log[q^2]]],
  Sqrt[g2sq[Log[q^2]]],
  Sqrt[g3sq[Log[q^2]]],
  Sqrt[ytsq[Log[q^2]]],
  10Abs[lam[Log[q^2]]],
  Sqrt[msq[Log[q^2]]]/125
}/.rge3//Evaluate, {q, table2aNeutOpt["MT"], 10^19}, PlotLegends->{g1, g2, g3, yt, "10|\[Lambda]|", "mH/125GeV"}]
