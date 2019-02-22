(* ::Package:: *)

(* Standard Model RGE at the two-loop level, amended by numerical three-loop results of gauge, top-yukawa, and Higgs couplings. *)
(* Prepared by Sho Iwamoto based on arXiv:1307.3536 (v4). *)


SMRGE::MissingLoopTools = "Install LoopTools before loading this module.";

(* Needs LoopTools before loading this file. *)
If[Not[MemberQ[$Packages, "LoopTools`"]], Message[SMRGE::MissingLoopTools]; Abort[]];


BeginPackage["SMRGE`"]


UsageFormat[str_] := Module[{func = {
  StringReplace[RegularExpression["arg\|(.+?)\|"] -> "\*StyleBox[\!\($1\), \"TI\"]"],
  StringReplace[RegularExpression["frac\|(.+?)\|(.+?)\|"] -> "\*FractionBox[\!\($1\), \!\($2\)]"],
  StringReplace[RegularExpression["\\*\\*(.+?)\\*\\*"] -> "\*StyleBox[$1, Bold]"],
  StringReplace[RegularExpression["\_(\\d+)"] -> "\_\!\(\*StyleBox[$1, \"TR\"]\)"],
  StringReplace[RegularExpression["\^(\\d+)"] -> "\^\!\(\*StyleBox[$1, \"TR\"]\)"]
}}, "\!\(" <> FixedPoint[Composition@@Reverse[func], str] <> "\)"]

couplings::usage = "The list of couplings for RGEs.";
observables::usage = "The list of observables for threshold corrections.";
FourPi::usage = "Numerical constant 4\[Pi].";

g1sq::usage = "Square of GUT-normalized U(1)-hypercharge coupling.";
g2sq::usage = "Square of SU(2)-weak coupling.";
g3sq::usage = "Square of SU(3)-color coupling.";
lam::usage = "Quartic coupling in the Higgs potential.";
msq::usage = "Quadratic coupling in the Higgs potential.";
ytsq::usage = "Square of 3rd-generation up-quark Yukawa coupling.";
ybsq::usage = "Square of 3rd-generation down-quark Yukawa coupling.";
yasq::usage = "Square of 3rd-generation charged-lepton Yukawa coupling.";

RGE::usage = "RGE[arg|p|] gives the RGE for the parameter arg|p| defined by frac|d arg|p||d log Q\^2|." // UsageFormat;

getCouplingsAtMt::usage = "hogehoge";
getCouplingsAtMt::InvalidOrder = "Invalid order `1` is specified.";


Begin["Private`"]


RGEList[g1sq] = {(41 g1sq^2)/(10 FourPi^2),(g1sq^2 ((199 g1sq)/50+(27 g2sq)/10+(44 g3sq)/5-(3 yasq)/2-ybsq/2-(17 ytsq)/10))/FourPi^4,(g1sq^2 (-((388613 g1sq^2)/24000)+(123 g1sq g2sq)/160+(789 g2sq^2)/64-(137 g1sq g3sq)/75-(3 g2sq g3sq)/5+(297 g3sq^2)/5+((27 g1sq)/50+(9 g2sq)/10-(9 lam)/5) lam+ytsq (-((2827 g1sq)/800)-(471 g2sq)/32-(29 g3sq)/5+(189 ytsq)/16)))/FourPi^6};
RGEList[g2sq] = {-((19 g2sq^2)/(6 FourPi^2)),(g2sq^2 ((9 g1sq)/10+(35 g2sq)/6+12 g3sq-yasq/2-(3 ybsq)/2-(3 ytsq)/2))/FourPi^4,(g2sq^2 (-((5597 g1sq^2)/1600)+(873 g1sq g2sq)/160+(324953 g2sq^2)/1728-(g1sq g3sq)/5+39 g2sq g3sq+81 g3sq^2+((3 g1sq)/10+(3 g2sq)/2-3 lam) lam+ytsq (-((593 g1sq)/160)-(729 g2sq)/32-7 g3sq+(147 ytsq)/16)))/FourPi^6};
RGEList[g3sq] = {-((7 g3sq^2)/FourPi^2),(g3sq^2 ((11 g1sq)/10+(9 g2sq)/2-26 g3sq-2 ybsq-2 ytsq))/FourPi^4,(g3sq^2 (-((523 g1sq^2)/120)-(3 g1sq g2sq)/40+(109 g2sq^2)/8+(77 g1sq g3sq)/15+21 g2sq g3sq+(65 g3sq^2)/2+ytsq (-((101 g1sq)/40)-(93 g2sq)/8-40 g3sq+15 ytsq)))/FourPi^6,-((2472.28` g3sq^5)/FourPi^8)};
RGEList[lam]  = {((27 g1sq^2)/400+(9 g1sq g2sq)/40+(9 g2sq^2)/16-yasq^2-3 ybsq^2-3 ytsq^2+lam (-((9 g1sq)/10)-(9 g2sq)/2+12 lam+2 yasq+6 ybsq+6 ytsq))/FourPi^2,1/FourPi^4 (-((3411 g1sq^3)/4000)-(1677 g1sq^2 g2sq)/800-(289 g1sq g2sq^2)/160+(305 g2sq^3)/32+((1887 g1sq^2)/400+(117 g1sq g2sq)/40-(73 g2sq^2)/16) lam+(-((9 g1sq^2)/8)+(33 g1sq g2sq)/20-(3 g2sq^2)/8) yasq+lam ((15 g1sq)/4+(15 g2sq)/4-yasq/2) yasq+yasq^2 (-((6 g1sq)/5)+5 yasq)+((9 g1sq^2)/40+(27 g1sq g2sq)/20-(9 g2sq^2)/8) ybsq+lam ((5 g1sq)/4+(45 g2sq)/4+40 g3sq-(3 ybsq)/2) ybsq+lam^2 ((54 g1sq)/5+54 g2sq-156 lam-24 yasq-72 ybsq-72 ytsq)+ybsq^2 ((2 g1sq)/5-16 g3sq+15 ybsq-3 ytsq)+(-((171 g1sq^2)/200)+(63 g1sq g2sq)/20-(9 g2sq^2)/8) ytsq+lam ((17 g1sq)/4+(45 g2sq)/4+40 g3sq-21 ybsq-(3 ytsq)/2) ytsq+ytsq^2 (-((4 g1sq)/5)-16 g3sq-3 ybsq+15 ytsq)),1/FourPi^6 (-1.508` g1sq^4-1.543` g1sq^3 g2sq+6.5` g1sq^2 g2sq^2-37.889` g1sq g2sq^3-114.091` g2sq^4+(0.663` g1sq^3+1.105` g1sq^2 g2sq+1.507` g1sq g2sq^2+7.536` g2sq^3) g3sq+g2sq^2 (79.638` g1sq+865.483` g2sq-57.144` g3sq) lam+g1sq^2 (28.168` g1sq+61.753` g2sq-8.381` g3sq) lam+(-185.532` g1sq^2-316.64` g1sq g2sq-790.28` g2sq^2) lam^2+g1sq^2 (11.117` g1sq+10.627` g2sq) ytsq+g2sq^2 (13.041` g1sq+62.5` g2sq) ytsq+(1.016` g1sq^2+11.386` g1sq g2sq+16.464` g2sq^2) g3sq ytsq+(-74.8599` g1sq^2+5.615` g1sq g2sq-319.664` g2sq^2+17.454` g1sq g3sq+15.1443` g2sq g3sq+356.968` g3sq^2) lam ytsq+(15.948` g1sq^2-70.356` g1sq g2sq+15.884` g2sq^2+17.57` g1sq g3sq+13.349` g2sq g3sq-50.201` g3sq^2) ytsq^2+lam (-21.015` g1sq-5.47` g2sq-662.866` g3sq-223.382` ytsq) ytsq^2+(33.93` g1sq+74.138` g2sq+250.494` g3sq-243.149` ytsq) ytsq^3+lam^3 (-77.49` g1sq-387.452` g2sq+6011.35` lam+873 ytsq)+lam^2 ytsq (-63.869` g1sq-359.539` g2sq+160.77` g3sq+1768.26` ytsq))};
RGEList[msq]  = {(msq (-((9 g1sq)/20)-(9 g2sq)/4+6 lam+yasq+3 ybsq+3 ytsq))/FourPi^2,(msq ((1671 g1sq^2)/800+(9 g1sq g2sq)/16-(145 g2sq^2)/32+((15 g1sq)/8+(15 g2sq)/8-(9 yasq)/4) yasq+((5 g1sq)/8+(45 g2sq)/8+20 g3sq-(27 ybsq)/4) ybsq+lam ((36 g1sq)/5+36 g2sq-30 lam-12 yasq-36 ybsq-36 ytsq)+((17 g1sq)/8+(45 g2sq)/8+20 g3sq-(21 ybsq)/2-(27 ytsq)/4) ytsq))/FourPi^4,1/FourPi^6 msq (g2sq^2 (9.931` g1sq+301.724` g2sq-28.572` g3sq)+g1sq^2 (8.378` g1sq+9.778` g2sq-4.191` g3sq)+(-65.8056` g1sq^2-37.8231` g1sq g2sq-64.5145` g2sq^2) lam+(-27.721` g1sq^2+11.47` g1sq g2sq-102.627` g2sq^2+8.727` g1sq g3sq+7.572` g2sq g3sq+178.484` g3sq^2) ytsq+lam^2 (-38.564` g1sq-192.822` g2sq+1026 lam+(297 ytsq)/2)+ytsq^2 (-7.50769` g1sq-3.82928` g2sq-209.24` g3sq+154.405` ytsq)+lam ytsq (-59.699` g1sq-318.591` g2sq+80.385` g3sq+347.394` ytsq))};
RGEList[ytsq] = {(ytsq (-((17 g1sq)/20)-(9 g2sq)/4-8 g3sq+yasq+(3 ybsq)/2+(9 ytsq)/2))/FourPi^2,(ytsq ((1187 g1sq^2)/600-(9 g1sq g2sq)/20-(23 g2sq^2)/4+(19 g1sq g3sq)/15+9 g2sq g3sq-108 g3sq^2+6 lam^2+((15 g1sq)/8+(15 g2sq)/8-(9 yasq)/4) yasq+((7 g1sq)/80+(99 g2sq)/16+4 g3sq+(5 yasq)/4-ybsq/4) ybsq+((393 g1sq)/80+(225 g2sq)/16+36 g3sq-12 lam-(9 yasq)/4-(11 ybsq)/4-12 ytsq) ytsq))/FourPi^4,1/FourPi^6 ytsq (16.099` g1sq^3-4.442` g1sq^2 g2sq-4.743` g1sq g2sq^2+169.829` g2sq^3-22.319` g1sq^2 g3sq-(321 g1sq g2sq g3sq)/20-21.072` g2sq^2 g3sq-15.096` g1sq g3sq^2+73.654` g2sq g3sq^2-619.35` g3sq^3+(-((1089 g1sq^2)/400)+(117 g1sq g2sq)/40-(171 g2sq^2)/16) lam+(9 g1sq+45 g2sq-36 lam) lam^2+(-24.422` g1sq^2+34.829` g1sq g2sq+16.99` g2sq^2+18.074` g1sq g3sq+48.37` g2sq g3sq+363.764` g3sq^2) ytsq+lam (-((127 g1sq)/10)-(135 g2sq)/2+16 g3sq+(15 lam)/4) ytsq+ytsq^2 (-((2437 g1sq)/80)-(1593 g2sq)/16-157 g3sq+198 lam+58.6028` ytsq))};
RGEList[ybsq] = {(ybsq (-(g1sq/4)-(9 g2sq)/4-8 g3sq+yasq+(9 ybsq)/2+(3 ytsq)/2))/FourPi^2,(ybsq (-((127 g1sq^2)/600)-(27 g1sq g2sq)/20-(23 g2sq^2)/4+(31 g1sq g3sq)/15+9 g2sq g3sq-108 g3sq^2+6 lam^2+((15 g1sq)/8+(15 g2sq)/8-(9 yasq)/4) yasq+((237 g1sq)/80+(225 g2sq)/16+36 g3sq-12 lam-(9 yasq)/4-12 ybsq) ybsq+((91 g1sq)/80+(99 g2sq)/16+4 g3sq+(5 yasq)/4-(11 ybsq)/4-ytsq/4) ytsq))/FourPi^4};
RGEList[yasq] = {(yasq (-((9 g1sq)/4)-(9 g2sq)/4+(5 yasq)/2+3 ybsq+3 ytsq))/FourPi^2,(yasq ((1371 g1sq^2)/200+(27 g1sq g2sq)/20-(23 g2sq^2)/4+6 lam^2+((537 g1sq)/80+(165 g2sq)/16-12 lam-3 yasq) yasq+((5 g1sq)/8+(45 g2sq)/8+20 g3sq-(27 yasq)/4-(27 ybsq)/4) ybsq+((17 g1sq)/8+(45 g2sq)/8+20 g3sq-(27 yasq)/4+(3 ybsq)/2-(27 ytsq)/4) ytsq))/FourPi^4};

couplings = {g1sq, g2sq, g3sq, lam, msq, ytsq, ybsq, yasq};
RGE[c_] := Total[RGEList[c]]


(* They use non-standard notation for A0 and B0. *)
(* A0 and B0 are the notation in LoopTools; TheirA0 and TheirB0 is the functions in 1307.3536. *)
(* In LoopTools, one should SetLambda[0] to match their definitions. *)
LoopTools`SetLambda[0];

ReTotal[exp__] := Total[Re/@List[exp]]
TheirA0[a_] := LoopTools`A0[a^2] /. {Mt^2->Mtsq, MW^2->MWsq, MZ^2->MZsq, Mh^2->Mhsq}
TheirB0[a_, b_, c_] := LoopTools`B0[a^2, b^2, c^2] /. {Mt^2->Mtsq, MW^2->MWsq, MZ^2->MZsq, Mh^2->Mhsq}
WeakScaleThreshold["lam", 0] := Gmu Mhsq / Sqrt[2];
WeakScaleThreshold["msq", 0] := Mhsq;
WeakScaleThreshold["yt", 0] := 2*Sqrt[Gmu Mtsq / Sqrt[2]];
WeakScaleThreshold["g2", 0] := 2 Sqrt[MWsq] Sqrt[Sqrt[2]Gmu];
WeakScaleThreshold["gY", 0] := 2 Sqrt[MZsq-MWsq] Sqrt[Sqrt[2]Gmu];

WeakScaleThreshold["yb", 0] := 2*Sqrt[Gmu Mbsq / Sqrt[2]];
WeakScaleThreshold["ya", 0] := 2*Sqrt[Gmu Masq / Sqrt[2]];

WeakScaleThreshold["lam", 1] := ReTotal[3*(Mhsq - 4*Mtsq)*Mtsq*TheirB0[Mh, Mt, Mt], 3*Mhsq*TheirA0[Mt], ((Mhsq^2 - 4*Mhsq*MZsq + 12*MZsq^2)*TheirB0[Mh, MZ, MZ])/ 4, (Mhsq*(7*MWsq - 4*MZsq)* TheirA0[MZ])/(2*(-MWsq + MZsq)), ((Mhsq^2 - 4*Mhsq*MWsq + 12*MWsq^2)* TheirB0[Mh, MW, MW])/ 2, (-3*Mhsq*MWsq* TheirA0[Mh])/(2*(Mhsq - MWsq)), (Mhsq*(-11 + (3*Mhsq)/(Mhsq - MWsq) - (3* MWsq)/(-MWsq + MZsq))*TheirA0[MW])/2, (9*Mhsq^2*TheirB0[Mh, Mh, Mh])/ 4, (Mhsq^2 + Mhsq*(-6*Mtsq + 2*MWsq + MZsq) - 8*(2*MWsq^2 + MZsq^2))/4]/(FourPi^2*V^4);
WeakScaleThreshold["msq", 1] := ReTotal[-3*Mhsq*TheirA0[Mh] + 24*Mtsq*TheirA0[Mt] - 2*(Mhsq + 6*MWsq)*TheirA0[MW] - (Mhsq + 6*MZsq)* TheirA0[MZ] + (9*Mhsq^2*TheirB0[Mh, Mh, Mh])/2 + 6*(Mhsq - 4*Mtsq)*Mtsq* TheirB0[Mh, Mt, Mt] + (Mhsq^2 - 4*Mhsq*MWsq + 12*MWsq^2)* TheirB0[Mh, MW, MW] + ((Mhsq^2 - 4*Mhsq*MZsq + 12*MZsq^2)*TheirB0[Mh, MZ, MZ])/ 2]/ (FourPi^2*V^2);
WeakScaleThreshold["yt", 1] := (g3sq*Sqrt[Mtsq]*(-8/3 - (8*TheirA0[Mt])/Mtsq))/(Sqrt[2]*FourPi^2*V) + (Sqrt[Mtsq]* ReTotal[(-Mhsq + 4*Mtsq)* TheirB0[Mt, Mh, Mt], ((-32*MWsq^2*MZsq + 40*MWsq*MZsq^2 - 17*MZsq^3 + Mtsq*(-64*MWsq^2 + 80*MWsq*MZsq - 7*MZsq^2))* TheirB0[Mt, Mt, MZ])/(9*Mtsq* MZsq), ((Mtsq^2 + Mtsq*MWsq - 2*MWsq^2)*TheirB0[Mt, 0, MW])/ Mtsq, (-10 + (3*Mhsq)/(Mhsq - MWsq) + (2*MWsq)/ Mtsq + (3*MWsq)/(MWsq - MZsq))* TheirA0[MW], (1 + (3*MWsq)/(-Mhsq + MWsq))* TheirA0[Mh], ((64*MWsq^2 + 36*Mtsq*MZsq - 56*MWsq*MZsq - 17*MZsq^2)*TheirA0[Mt])/(9*Mtsq* MZsq), (-3 + (3*MWsq)/(-MWsq + MZsq) + (32*MWsq^2 - 40*MWsq*MZsq + 17*MZsq^2)/(9*Mtsq*MZsq))*TheirA0[MZ], Mhsq/2 - 3*Mtsq - 9*MWsq, (7*MZsq)/ 18, (64*MWsq^2)/(9*MZsq)])/(Sqrt[2]*FourPi^2*V^3);
WeakScaleThreshold["g2", 1] := (2*Sqrt[MWsq]* ReTotal[((-2*Mhsq)/3 + Mhsq^2/(6*MWsq) + 2*MWsq)* TheirB0[MW, Mh, MW], (-Mtsq - Mtsq^2/MWsq + 2*MWsq)*TheirB0[MW, 0, Mt], ((-68*MWsq - (48*MWsq^2)/MZsq + 16*MZsq + MZsq^2/MWsq)* TheirB0[MW, MW, MZ])/6, ((-27 + Mhsq*(9/(Mhsq - MWsq) + MWsq^(-1)) + MWsq*(9/(MWsq - MZsq) + 48/MZsq) + MZsq/MWsq)*TheirA0[MW])/6, (2 - (Mhsq*(Mhsq + 8*MWsq))/(6*(Mhsq - MWsq)*MWsq))* TheirA0[Mh], (1 + Mtsq/MWsq)*TheirA0[Mt], ((-17 + (24*MWsq)/MZsq - MZsq/MWsq + (9*MWsq)/(-MWsq + MZsq))* TheirA0[MZ])/6, (-3*Mhsq + 18*Mtsq - 374*MWsq + (288*MWsq^2)/MZsq - 3*MZsq)/36])/ (FourPi^2*V^3);
WeakScaleThreshold["gY", 1] := (2*Sqrt[-MWsq + MZsq]* ReTotal[(88/ 9 - (124*MWsq)/(9*MZsq) + (Mhsq + 34*MWsq)/(6*(-MWsq + MZsq)))* TheirA0[MZ], ((Mhsq - 4*MWsq)*TheirA0[Mh])/(2*(Mhsq - MWsq)), (-7/9 + (64*MWsq)/(9*MZsq) - Mtsq/(-MWsq + MZsq))* TheirA0[Mt], ((Mhsq^2 + 2*MWsq*(MWsq - 15*MZsq) + 3*Mhsq*(2*MWsq + 7*MZsq))*TheirA0[MW])/ (6*(Mhsq - MWsq)*(MWsq - MZsq)), -(((Mtsq^2 + Mtsq*MWsq - 2*MWsq^2)* TheirB0[MW, 0, Mt])/(MWsq - MZsq)), -((Mhsq^2 - 4*Mhsq*MZsq + 12*MZsq^2)* TheirB0[MZ, Mh, MZ])/(6*(MWsq - MZsq)), ((Mhsq^2 - 4*Mhsq*MWsq + 12*MWsq^2)*TheirB0[MW, Mh, MW])/(6*(MWsq - MZsq)), ((-48*MWsq^3 - 68*MWsq^2*MZsq + 16*MWsq*MZsq^2 + MZsq^3)* TheirB0[MW, MW, MZ])/(6*(MWsq - MZsq)*MZsq), ((7*Mtsq - 23*MWsq - (64*Mtsq*MWsq)/MZsq + 17*MZsq - (9*(Mtsq - MWsq)*MWsq)/(-MWsq + MZsq))* TheirB0[MZ, Mt, Mt])/9, ((-48*MWsq^3 - 68*MWsq^2*MZsq + 16*MWsq*MZsq^2 + MZsq^3)* TheirB0[MZ, MW, MW])/(6*MZsq*(-MWsq + MZsq)), (-3*Mhsq - 242*MWsq + Mtsq*(82 - (256*MWsq)/MZsq) + (576*MWsq^2)/MZsq + 257*MZsq + (36*MWsq)/(-MWsq + MZsq))/36])/(FourPi^2*V^3);

(* These two-loop numerical values are valid for SetMudim[Mtsq]. (not Mt, but Mtsq) *)
WeakScaleThreshold["lam", 2] := Plus[
  (g3sq/FourPi^4)*(-23.88 + 0.12*(Sqrt[Mhsq]-125) - 0.64*(Sqrt[Mtsq]-173)),
  (1/FourPi^4)   *( -9.45 - 0.12*(Sqrt[Mhsq]-125) - 0.21*(Sqrt[Mtsq]-173))];
WeakScaleThreshold["msq", 2] := Plus[
  (g3sq Mhsq / FourPi^4)*(-140.50 + 2.91*(Sqrt[Mhsq]-125) - 3.72*(Sqrt[Mtsq]-173)),
  (Mhsq / FourPi^4)     *(-149.47 + 2.54*(Sqrt[Mhsq]-125) - 4.69*(Sqrt[Mtsq]-173))];
WeakScaleThreshold["yt", 2] := Plus[
  (1      / FourPi^4)*( 6.48 - 0.01*(Sqrt[Mhsq]-125) + 0.18*(Sqrt[Mtsq]-173)),
  (g3sq   / FourPi^4)*(-7.53 + 0.09*(Sqrt[Mhsq]-125) - 0.23*(Sqrt[Mtsq]-173)),
  (g3sq^2 / FourPi^4)*(-145.08 - 0.84*(Sqrt[Mtsq]-173))];
WeakScaleThreshold["g2", 2] := Plus[
  (1      / FourPi^4)*( 2.25 + 0.01*(Sqrt[Mhsq]-125) + 0.01*(Sqrt[Mtsq]-173)),
  (g3sq   / FourPi^4)*( 3.00 + 0.01(Sqrt[Mtsq]-173))];
WeakScaleThreshold["gY", 2] := Plus[
  (1      / FourPi^4)*(-7.55  - 0.01*(Sqrt[Mhsq]-125) - 0.11*(Sqrt[Mtsq]-173)),
  (g3sq   / FourPi^4)*(-14.66 - 0.14(Sqrt[Mtsq]-173))];

WeakScaleThreshold["g3", {0,1,2}] := 1.1666 + (0.00314/0.0007)(asMZ-0.1184) - 0.00046(Sqrt[Mtsq]-173.34);

getCouplingsAtMt[observableValues_, order_: 2] := Module[{
   repRule = {
     Mtsq -> observableValues["MT"]^2,
     Mhsq -> observableValues["MH"]^2,
     MZsq -> observableValues["MZ"]^2,
     MWsq -> observableValues["MW"]^2,
     g3sq -> WeakScaleThreshold["g3", {0,1,2}]^2 /. asMZ->observableValues["as(MZ)"],
     Mbsq -> observableValues["MB"]^2,
     Masq -> observableValues["MTau"]^2,
     Gmu -> observableValues["Gmu"],
     V -> 1/Sqrt[Sqrt[2]observableValues["Gmu"]],
     FourPi -> 4\[Pi]
   },
   sub = Switch[order,
     0, Function[{name, rep}, WeakScaleThreshold[name,0] //. rep],
     1, Function[{name, rep}, WeakScaleThreshold[name,0] + WeakScaleThreshold[name,1] //. rep],
     2, Function[{name, rep}, WeakScaleThreshold[name,0] + WeakScaleThreshold[name,1] + WeakScaleThreshold[name, 2] //. rep],
     _, Message[getCouplingsAtMt::InvalidOrder, order]; Abort[]]
},
  LoopTools`SetLambda[0];
  LoopTools`SetMudim[Mtsq/.repRule]; (* Not mt but mtsq *)
  <|
    g1sq -> (5/3)*sub["gY", repRule]^2,
    g2sq -> sub["g2", repRule]^2,
    g3sq -> (g3sq//.repRule),
    lam -> sub["lam", repRule],
    msq -> sub["msq", repRule],
    ytsq -> sub["yt", repRule]^2,
    ybsq -> (WeakScaleThreshold["yb", 0]^2 //. repRule),
    yasq -> (WeakScaleThreshold["ya", 0]^2 //. repRule)
  |>]


End[];
EndPackage[];


(* ::Input:: *)
(*initial = {0.01, 0.4, 3, 0.1, 80, 1, 0.04, 0.01};*)
(*sol = NDSolve[{{D[#[t], t] == Total[RGEList[#] /. CouplingsReplace[t]] & /@ couplings}, MapThread[#1[Log10[81.0]] == #2 &, {couplings, initial}]} /. FourPi -> 4 \[Pi], couplings, {t, Log10[81.0], 17} ];*)
(*LogPlot[Evaluate[#[x] & /@ sol[[1, All, 2]]], {x, 2, 17}]*)
