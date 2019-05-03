#!wolframscript -script
(* ::Package:: *)

(* ::Section:: *)
(*Higgs potential running*)


(* ::Input::Initialization:: *)
SetDirectory[Which[$FrontEnd=!=Null,NotebookDirectory[],$InputFileName=!="", DirectoryName[$InputFileName], True, "."]];
$ProjectRoot = RunProcess[{"git","rev-parse" , "--show-toplevel"}]["StandardOutput"] // StringTrim;
$Path = Append[$Path, FileNameJoin[{$ProjectRoot, "vendor", "LoopTools", "build"}]] // DeleteDuplicates;


(* ::Input::Initialization:: *)
Install["LoopTools"];
Needs["SMRGE`", "1307_3536.wl"];
<<plot_tools.wl


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


EWInput =  <| "MW" -> 80.387,      "MZ" -> 91.1875,   "MH" -> 125.09,
              "MT" -> 173.2,       "MB" -> 4.18,      "MTau" -> 1.776, \
              "as(MZ)" -> 0.1185,  "Gmu" -> 1.1663787*^-5 |>;
SMThresholdCorrectionOrder = 1; (* regardless of this value, g3 is set by Eq. (60) of 1307.3536 *)
SMRGEOrder = 2;

EWatMT = getCouplingsAtMt[EWInput, SMThresholdCorrectionOrder]
rges = NDSolve @@ getRGEs[EWatMT, EWInput["MT"], 10^20, SMRGEOrder] // First


ExportTicks = Prepend[10^Range[5/2, 10, 1/4], EWInput["MT"]];
table = Table[
  {Log10[q], lam[Log[q^2]], Sqrt[msq[Log[q^2]]]} /. rges // SetPrecision[#,8]&,
  {q, ExportTicks}]
headers = If[Length[#]==Length[table[[1]]], StringRiffle[#, "\t"] <> "\n", Print["Length mismatch."]; Abort[]]&/@{
  {"# V = lam |H|^4 - (m^2/2) |H|^2;", "n_threshold=1;", "n_RGE=2"},
  {"# Log10[scale]", "lam[scale]", "m[scale]"}
};
Export["higgs_potential_running.tsv", headers <> ExportString[table, "TSV"], "Text"]  


Clear[msq]
f = msq /. rges;
msq = Function @@ {Q, f[2Log[Q]]};
Put[msq, "msq.dat"]
Clear[msq];


(* ::Subsection:: *)
(*RGE plot*)


rge1 = NDSolve@@getRGEs[getCouplingsAtMt[EWInput, 0], EWInput["MT"], 10^20, 1];
rge2 = NDSolve@@getRGEs[getCouplingsAtMt[EWInput, 1], EWInput["MT"], 10^20, 2];
rge3 = NDSolve@@getRGEs[getCouplingsAtMt[EWInput, 2], EWInput["MT"], 10^20, All];
LogLinearPlot[lam[Log[q^2]]/.{rge1, rge2, rge3}//Evaluate, {q, EWInput["MT"], 10^19}, PlotRange->{{10^2, 10^20}, {-0.10, 0.15}}]
