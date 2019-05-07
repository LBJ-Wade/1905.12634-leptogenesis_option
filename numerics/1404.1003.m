(* ::Package:: *)

BeginPackage["ResLept14041003`"];


TwoGenResLept::usage = "Create objects from Yukawa(3,2) and Mass (2) matrices.";
DiffEqs::usage = "Construct differential equations from TwoGenResLept object and g*.";
B1dB2::usage = "Equals to K1[z] / K2[z]."

TwoGenResLept::InvalidDimension = "Invalid dimension of h or M.";


Begin["`Private`"];


Dag = ConjugateTranspose;
MPlanck = 1.2209102*^19; (*2.2*)
(* numerical expression for Bessel function ratio *)
B1dB2[z_] /; z<1*^-2 := 1/4 z (2+z^2 (EulerGamma+Log[z/2]));
B1dB2[z_] /; z>100 := 1+15/(8 z^2)-3/(2 z);
B1dB2[z_] /; NumericQ[z] := BesselK[1,z]/BesselK[2,z];
TwoGenResLept[hMat_, MMat_] := Module[{
     h, hc, M,
     AMat,  A,
     \[CapitalGamma]0Mat, \[CapitalGamma]0,
     obj,
     NL=3, NN=2
  },
  If[Dimensions[hMat] != {NL, NN} || Dimensions[MMat] != {NN},
    Message[TwoGenResLept::InvalidDimension];
    Abort[]];
  AMat = Conjugate[Dag[hMat].hMat]/(16\[Pi]);
  \[CapitalGamma]0Mat = Table[2 MMat[[a]] AMat[[a, a]], {a, NN}];
  obj["h"] = hMat;
  obj["hc"] = Conjugate[hMat];
  obj["M"] = MMat;
  obj["A"] = AMat;
  obj["\[CapitalGamma]0"] = \[CapitalGamma]0Mat;

  (* prepare function form to reduce brackets *)
  h  = Function[{i,j}, hMat[[i,j]]];
  M  = Function[{i},   MMat[[i]]];
  A  = Function[{i,j}, AMat[[i,j]]];
  \[CapitalGamma]0 = Function[{i},   \[CapitalGamma]0Mat[[i]]];
  hc = Function[{i,j}, Conjugate[hMat[[i,j]]]];

  obj["H"] = Table[Module[{b=3-a},
    h[l, a] - I h[l,b]M[a](M[a]A[a,b]+M[b]A[b,a])/(M[a]^2-M[b]^2+2I M[a]^2A[b,b])
  ], {l,3}, {a,2}]; (*A.1 or 3.31 of 0309342*)
  obj["Hc"] = Table[Module[{b=3-a},
    hc[l, a] - I hc[l,b]M[a](M[a]Conjugate[A[a,b]]+M[b]Conjugate[A[b,a]])/(M[a]^2-M[b]^2+2I M[a]^2Conjugate[A[b,b]])
  ], {l,3}, {a,2}]; (*cf. 3.32 of 0309342*)

  obj["\[CapitalGamma]"] = Table[(MMat[[a]]/(16\[Pi])) ((Dag[obj["H"]].obj["H"])[[a,a]]+(Dag[obj["Hc"]].obj["Hc"])[[a,a]]), {a,2}];

  obj["B"] = Module[{
      s = (Dag[obj["H"]].obj["H"]+Dag[obj["Hc"]].obj["Hc"])
    },
    Table[(Abs[obj["H"][[l,a]]]^2+Abs[obj["Hc"][[l,a]]]^2)/s[[a,a]], {l, 3}, {a, 2}]];
  
  obj["\[Epsilon]mix(full)"] = Table[Module[{H=obj["H"], Hc=obj["Hc"]},
    (Abs[H[[l,a]]]^2 - Abs[Hc[[l,a]]]^2) / ( (Dag[H].H)[[a,a]] + (Dag[Hc].Hc)[[a,a]] )
  ], {l,3}, {a,2}]; (*2.22*)

  obj["\[Epsilon]mix"] = Table[Module[{b=3-a, hdh=Dag[hMat].hMat}, 
    (Im[hc[l,a]h[l,b]hdh[[a,b]]] + (M[a]/M[b])Im[hc[l,a]h[l,b]hdh[[b,a]]]) / (hdh[[a,a]]hdh[[b,b]]) *
    (M[a]^2-M[b]^2)M[a]\[CapitalGamma]0[b] / ((M[a]^2-M[b]^2)^2 + (M[a]\[CapitalGamma]0[b])^2)
  ], {l,3}, {a,2}]; (*2.22\[Rule]A.2*)

  (* borrowed from 6.3 of 1611.03287; note that Y=Dagger[h] *)
  obj["\[Epsilon]osc"] = Table[Module[{b=3-a, hdh=Dag[hMat].hMat, YYD=Dag[hMat].hMat},
    (Im[hc[l,a]h[l,b]hdh[[a,b]]] + (M[a]/M[b])Im[hc[l,a]h[l,b]hdh[[b,a]]]) / (hdh[[a,a]]hdh[[b,b]]) *
    (M[a]^2-M[b]^2)M[a]\[CapitalGamma]0[b] / ((M[a]^2-M[b]^2)^2 + (M[a]\[CapitalGamma]0[a]+M[b]\[CapitalGamma]0[b])^2 Det[Re[YYD]]/(YYD[[a,a]]YYD[[b,b]]))
  ], {l,3}, {a,2}]; (*2.22\[Rule]A.2*)

  obj["\[Epsilon]"] = obj["\[Epsilon]mix"] + obj["\[Epsilon]osc"];
  
  (*kappa*)
  Module[{
      H=obj["H"],   Hs=Conjugate[obj["H"]],
      Hc=obj["Hc"], Hcs=Conjugate[obj["Hc"]],
      \[CapitalGamma]=obj["\[CapitalGamma]"],
      rs, rd, hsum
    },
    rs[p_, a_, b_] := H[[p,a]]Hs[[p,b]]+Hc[[p,a]]Hcs[[p,b]];
    rd[p_, a_, b_] := H[[p,a]]Hs[[p,b]]-Hc[[p,a]]Hcs[[p,b]];
    hsum[a_, b_] := (Dag[Hc].Hc + Dag[H].H)[[a,b]];
    obj["\[Kappa]"] = Table[Sum[
      ((M[a]^2 M[b]^2 (rd[p,b,a]^2+hsum[a,b] rs[p,b,a]) (M[a]^3 \[CapitalGamma][[a]]+M[b]^3 \[CapitalGamma][[b]]))/(8\[Pi] (M[1]^4 rs[p,1,1]+M[2]^4 rs[p,2,2]) (M[b] \[CapitalGamma][[a]]+M[a] \[CapitalGamma][[b]])^2))/(1-(2 I (M[a]-M[b]))/(\[CapitalGamma][[a]]+\[CapitalGamma][[b]])), {a,2}, {b,2}
    ],{p,3}];
  ];
  obj
]


DiffEqs[twoGenResLeptObj_, gstar_] := Module[{
    obj
  },
  obj["rl"] = twoGenResLeptObj;
  obj["gstar"] = gstar;
  obj["HN"] = Sqrt[4\[Pi]^3 gstar/45] Min[twoGenResLeptObj["M"]]^2 / MPlanck; (*2.2*)
  obj["Ka"] = Table[twoGenResLeptObj["\[CapitalGamma]"][[a]] / (obj["HN"] Zeta[3]), {a,2}];
  obj["Keff"] = Table[twoGenResLeptObj["\[Kappa]"][[l]] Sum[obj["Ka"][[a]] twoGenResLeptObj["B"][[l,a]], {a,2}], {l,3}];
  obj["\[Eta]DE"][\[Eta]_, a:1|2, z_] := Module[{},
    D[\[Eta][a][z],z]==B1dB2[z](1+(1-obj["Ka"][[a]] z)\[Eta][a][z])
  ];
  obj["\[Delta]DE"][\[Eta]_, \[Delta]_, l:1|2|3, z_] := Module[{rl=obj["rl"]},
    D[\[Delta][l][z],z]==z^3 BesselK[1,z] Sum[obj["Ka"][[a]](rl["\[Epsilon]"][[l,a]] \[Eta][a][z] -(2/3)rl["B"][[l,a]]rl["\[Kappa]"][[l]]\[Delta][l][z]), {a,2}]];
  obj
]  


End[];
EndPackage[];
