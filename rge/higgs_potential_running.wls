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
running = LogLinearPlot[
  {lam[Log[q^2]], Sqrt[msq[Log[q^2]]/2]/1000} /.{rge1, rge2, rge3} //Evaluate,
  {q, EWInput["MT"], 10^19},
    PlotRange->{{10^2, 10^18}, {-0.05, 0.15}},
  PlotStyle->{
    {color[1], Dotted},
    {color[3], Dotted},
    {color[1]},
    {color[3]},
    {color[1], Dashed},
    {color[3], Dashed}
  },
  FrameLabel->{MaTeX["\\bar Q\\,\\textrm{/\\,GeV}"], None},
  PlotLegends->Placed[
    LineLegend[
      {Dotted, Black, Dashed},
      {"1\[Hyphen]loop", "2\[Hyphen]loop", "3\[Hyphen]loop"},
      LegendFunction->(Framed[#, FrameStyle->None, FrameMargins->0, Background->Directive[White,Opacity[0]]]&),
      LabelStyle->{ FontFamily->"Times New Roman", FontSize->15}
    ], {0.81,0.45}],
  Epilog->{
    Inset[MaTeX["V=\\lambda |\\phi|^4 - \\mu^2|\\phi|^2"],{33,0.13}],
    Inset[MaTeX[ColorToTeX[color[3], "one"]<>"\\color{one}{\\mu(\\bar Q)\\,/\\mathrm{TeV}}"], {16, 0.112}],
    Inset[MaTeX[ColorToTeX[color[1], "one"]<>"\\color{one}{\\lambda(\\bar Q)}"], {18, 0.03}]
  }] //. {Charting`ScaledTicks[arg_]:>MyChartingScaledTicks[arg]}
outputPDF[running];


(* ::Subsection:: *)
(*For manuscript*)


KeyValueMap["$"<>#1<>"$ & " <> TextString[#2] <> "\\\\\n"&][EWInput] // StringJoin


getCouplingsAtMt[EWInput, 1] // {
  "g_1" -> (#[g1sq] // Sqrt),
  "g_2" -> (#[g2sq] // Sqrt),
  "g_3" -> (#[g3sq] // Sqrt),
  "\\lambda" -> #[lam],
  "\\mu" -> (#[msq] / 2 // Sqrt),
  "y_t" -> (#[ytsq] // Sqrt),
  "y_b" -> (#[ybsq] // Sqrt),
  "y_\\tau" -> (#[yasq] // Sqrt)
} &


NumToStr[a_] := NumberForm[a, 4] // TextString
{
  "Q" -> q,
  "g_1" -> (g1sq[Log[q^2]] // Sqrt),
  "g_2" -> (g2sq[Log[q^2]] // Sqrt),
  "g_3" -> (g3sq[Log[q^2]] // Sqrt),
  "\\lambda" -> lam[Log[q^2]],
  "\\mu" -> (msq[Log[q^2]] / 2 // Sqrt),
  "y_t" -> (ytsq[Log[q^2]] // Sqrt),
  "y_b" -> (ybsq[Log[q^2]] // Sqrt),
  "y_\\tau" -> (yasq[Log[q^2]] // Sqrt)
} //. rge2 //. q->{173.2, 10^7} //. Rule[a_, b_] :> StringJoin["$", a, "$ & ", NumToStr[b[[1]]], " & ", NumToStr[b[[2]]] , "\\\\\n"] // StringJoin 


NumToStr[a_] := NumberForm[a, 4] // TextString
{
  "Q" -> q,
  "g_1" -> (g1sq[Log[q^2]] // Sqrt),
  "g_2" -> (g2sq[Log[q^2]] // Sqrt),
  "g_3" -> (g3sq[Log[q^2]] // Sqrt),
  "\\lambda" -> lam[Log[q^2]],
  "\\mu" -> (msq[Log[q^2]] / 2 // Sqrt),
  "y_t" -> (ytsq[Log[q^2]] // Sqrt),
  "y_b" -> (ybsq[Log[q^2]] // Sqrt),
  "y_\\tau" -> (yasq[Log[q^2]] // Sqrt)
} //. rge3 //. q->{173.2, 10^7} //. Rule[a_, b_] :> StringJoin["$", a, "$ & ", NumToStr[b[[1]]], " & ", NumToStr[b[[2]]] , "\\\\\n"] // StringJoin 
