(* ::Package:: *)

BeginPackage["NuFIT`"];
NuFITVersion = "4.0";


data::usage = "Raw data from NuFIT collaboration.";
MCMC::usage = "Generate an object for MCMC chain.";
BestFit::usage = "Give the best fit value.";
CI1::usage = "Give 1-sigma confidence interval.";
CI3::usage = "Give 3-sigma confidence interval.";
ToPMNS::usage = "Construct PMNS matrix from input values.";
ToMasses::usage = "Construct neutrino mass array from input values.";

NuFIT::InvalidKey = "Invalid key or hierarchy is specified.";


Begin["`Private`"]


ToAssociation[{bfp_, p1_, m1_, ci3min_, ci3max_}] := Association[
  "bfp"->bfp, "1\[Sigma]"->Interval[{bfp-Abs[m1], bfp+Abs[p1]}], "3\[Sigma]"->Interval[{ci3min, ci3max}]
]
keys = {"Sin2\[Theta]12", "Sin2\[Theta]23", "Sin2\[Theta]13", "\[Delta]CP", "\[CapitalDelta]msq21", "\[CapitalDelta]msq3l"};
data["NH", "Sin2\[Theta]12"] = {0.310, 0.013, -0.012, 0.275, 0.350} // ToAssociation;
data["IH", "Sin2\[Theta]12"] = {0.310, 0.013, -0.012, 0.275, 0.350} // ToAssociation;
data["NH", "Sin2\[Theta]23"] = {0.580, 0.017, -0.021, 0.418, 0.627} // ToAssociation;
data["IH", "Sin2\[Theta]23"] = {0.584, 0.016, -0.020, 0.423, 0.629} // ToAssociation;
data["NH", "Sin2\[Theta]13"] = {0.02241, 0.00065, -0.00065, 0.02045, 0.02439} // ToAssociation;
data["IH", "Sin2\[Theta]13"] = {0.02264, 0.00066, -0.00066, 0.02068, 0.02463} // ToAssociation;
data["NH", "\[Delta]CP"]    = {215, 40, -29, 125, 392} Degree // N // ToAssociation;
data["IH", "\[Delta]CP"]    = {284, 27, -29, 196, 360} Degree // N // ToAssociation;
data["NH", "\[CapitalDelta]msq21"] = {7.39, 0.21, -0.20, 6.79, 8.01} * 10^(-23) // ToAssociation;
data["IH", "\[CapitalDelta]msq21"] = {7.39, 0.21, -0.20, 6.79, 8.01} * 10^(-23) // ToAssociation;
data["NH", "\[CapitalDelta]msq3l"] = {+2.525, 0.033, -0.032, +2.427, +2.625} * 10^(-21) // ToAssociation;
data["IH", "\[CapitalDelta]msq3l"] = {-2.512, 0.034, -0.032, -2.611, -2.412} * 10^(-21) // ToAssociation;
data[__] := (Message[NuFIT::InvalidKey]; Abort[])


(* Prepare normal distribution *)
BestFit[hier_String] := Association @@ (#->data[hier, #]["bfp"]& /@ keys);
BestFit[hier_String, key_String] := data[hier, key]["bfp"];
CI1[hier_String, key_String] := data[hier, key]["1\[Sigma]"];
CI3[hier_String, key_String] := data[hier, key]["3\[Sigma]"];
AvgSD[hier_String, key_String] := Mean[Abs[MinMax[data[hier, key]["1\[Sigma]"]] - BestFit[hier, key]]];
Table[
  (ND[h, k] = NormalDistribution[BestFit[h, k], AvgSD[h, k]]),
  {h, {"NH", "IH"}}, {k, keys}
];
ToPMNS[data_Association] := Module[{
    s13=Sqrt[data["Sin2\[Theta]13"]],
    s23=Sqrt[data["Sin2\[Theta]23"]],
    s12=Sqrt[data["Sin2\[Theta]12"]],
    c13=Sqrt[1-data["Sin2\[Theta]13"]],
    c23=Sqrt[1-data["Sin2\[Theta]23"]],
    c12=Sqrt[1-data["Sin2\[Theta]12"]],
    e = Exp[I data["\[Delta]CP"]]},
  {{ c12 c13,  s12 c13, s13 Conjugate[e]},
   {-s12 c23 - c12 s23 s13 e,  c12 c23 - s12 s23 s13 e, s23 c13},
   { s12 s23 - c12 c23 s13 e, -c12 s23 - s12 c23 s13 e, c23 c13}
  }]
ToMasses[data_Association] := Module[{d21=data["\[CapitalDelta]msq21"], d3l=data["\[CapitalDelta]msq3l"]},
  If[d3l>0, {0, Sqrt[d21], Sqrt[d3l]}, {Sqrt[-d3l-d21], Sqrt[-d3l], 0}]];


MCMC[hier_String, fixed_Association] := MCMC[hier, fixed, 1/3];
MCMC[hier_String, fixed_Association, stepFactor_] := Module[{
  variables = Complement[keys, Keys[fixed]],
  MCMCObject
  },
  MCMCObject["values"] := Association @@ (Rule[#, MCMCObject[#]] &/@ keys);
  (* set fixed values and properties*)
  KeyValueMap[(MCMCObject[#1] = #2) &, fixed];
  MCMCObject["hierarchy"] = hier;
  MCMCObject["variables"] = variables;
  MCMCObject["stepFactor"] = stepFactor;
  MCMCObject["pdf"] = Function[assoc, Times @@ (PDF[ND[MCMCObject["hierarchy"], #]][assoc[#]]&/@ MCMCObject["variables"])];
  (* give random initial values *)
  (MCMCObject[#] = RandomVariate[ND[hier, #]]) &/@ variables;
  MCMCObject["currentP"] = MCMCObject["pdf"][MCMCObject["values"]];
  (* iterator *)
  MCMCObject["next"] := Module[{
    nextCandidate = Association @@ (#->RandomVariate[NormalDistribution[MCMCObject[#], MCMCObject["stepFactor"]*AvgSD[MCMCObject["hierarchy"], #]]] &/@ MCMCObject["variables"]),
    nextP
    },
    nextP = MCMCObject["pdf"][nextCandidate];
    If[nextP / MCMCObject["currentP"] > Random[],
      (MCMCObject[#] = nextCandidate[#]) &/@ variables;
      MCMCObject["currentP"] = nextP;
    ];
    MCMCObject["values"]];
  MCMCObject]


End[];


EndPackage[];
