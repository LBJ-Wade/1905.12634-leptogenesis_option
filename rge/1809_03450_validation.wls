#!wolframscript -script
(* ::Package:: *)

(* ::Section:: *)
(*Validation of 1809_03450 RGEs*)


(* ::Input::Initialization:: *)
SetDirectory[NotebookDirectory[]];
$ProjectRoot = RunProcess[{"git","rev-parse" , "--show-toplevel"}]["StandardOutput"] // StringTrim;
$Path = Append[$Path, FileNameJoin[{$ProjectRoot, "vendor", "LoopTools", "build"}]] // DeleteDuplicates;


(* ::Input::Initialization:: *)
Install["LoopTools"];
Needs["SMRGE`", "1307_3536.wl"];
<<plot_tools.wl

IntroduceTypo[introduce_] := Module[{typo = If[introduce===False, "NoTypo", "WithTypo"]},
  SMRGE`Private`WeakScaleThreshold["gY", 1] := SMRGE`Private`WeakScaleThreshold["gY", 1, typo]]  
getCouplingsWrongly[args__] := Module[{
  orig = getCouplingsAtMt[args]},
  orig[g3sq] = 4Pi {args}[[1]]["as(MZ)"];
  orig]


(* ::Subsection:: *)
(*Table 2*)


table2a=<|
  "Gmu"->1.1663787*^-5,
  "as(MZ)"->0.1185,
  "MZ"->91.1875,
  "MW"->80.387,
  "MH"->125.09,
  "MT"->173.2,
  "MB"->4.18,
  "MTau"->1.776|>;

table2b[p_]:={
  #[lam],
  Sqrt[#[msq]],
  Sqrt[#[g1sq]],
  Sqrt[#[g2sq]],
  Sqrt[#[g3sq]],
  Sqrt[#[ytsq]],
  Sqrt[#[ytsq]],
  Sqrt[#[ybsq]],
  Sqrt[#[yasq]]
} &@ p

IntroduceTypo[False]
table2b[getCouplingsAtMt[table2a, 0]]
table2b[getCouplingsAtMt[table2a, 1]]
table2b[getCouplingsAtMt[table2a, 2]]


(* Introduces the bugs found in 1809.03450. Yet, g2 values are slightly different. *)
IntroduceTypo[True]
table2b[getCouplingsWrongly[table2a, 0]]
table2b[getCouplingsWrongly[table2a, 1]]
table2b[getCouplingsWrongly[table2a, 2]]


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


(* All Possibilities *)
IntroduceTypo[False]
rgesA = (NDSolve@@#)[[1]]&/@Flatten[Table[getRGEs[getCouplingsAtMt[table2a, i], table2a["MT"], 10^20, j], {i, {0, 1, 2}}, {j, {1, 2, All}}],1];
rgesB = (NDSolve@@#)[[1]]&/@Flatten[Table[getRGEs[getCouplingsWrongly[table2a, i], table2a["MT"], 10^20, j], {i, {0, 1, 2}}, {j, {1, 2, All}}],1];
IntroduceTypo[True]
rgesC = (NDSolve@@#)[[1]]&/@Flatten[Table[getRGEs[getCouplingsAtMt[table2a, i], table2a["MT"], 10^20, j], {i, {0, 1, 2}}, {j, {1, 2, All}}],1];
rgesD = (NDSolve@@#)[[1]]&/@Flatten[Table[getRGEs[getCouplingsWrongly[table2a, i], table2a["MT"], 10^20, j], {i, {0, 1, 2}}, {j, {1, 2, All}}],1];
LogLinearPlot[lam[Log[q^2]]/.{rgesA, rgesB, rgesC, rgesD}//Evaluate, {q, table2a["MT"], 10^19}, PlotRange->{{10^5, 10^12}, {-0.05, 0.06}}]


rge1 = NDSolve@@getRGEs[getCouplingsAtMt[table2a, 0], table2a["MT"], 10^20, 1];
rge2 = NDSolve@@getRGEs[getCouplingsAtMt[table2a, 1], table2a["MT"], 10^20, 2];
rge3 = NDSolve@@getRGEs[getCouplingsAtMt[table2a, 2], table2a["MT"], 10^20, All];
LogLinearPlot[lam[Log[q^2]]/.{rge1, rge2, rge3}//Evaluate, {q, table2a["MT"], 10^19}, PlotRange->{{10^5, 10^12}, {-0.05, 0.06}}]
