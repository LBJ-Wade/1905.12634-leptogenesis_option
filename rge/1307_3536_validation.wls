#!/usr/bin/env wolframscript
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



