(* ::Package:: *)
(* Time-Stamp: <2019-05-15 18:17:47> *)

thisFile := If[$FrontEnd === Null, $Input, NotebookFileName[]];
SetAttributes[outputPDF, HoldFirst];
outputPDF[obj_] := outputPDF[TextString[HoldForm[obj]], obj]
(* magnify for a magical spell... *)
outputPDF[title_, obj_] := Export[FileBaseName[thisFile] <> "_" <> Evaluate[title] <> ".pdf", Magnify[obj, 1]];


(* Color Scheme good for color-blind and monochromatic; cf. https://github.com/misho104/scicolpick ; colordistance 29.67 *)
colors = RGBColor /@ {"#001b95", "#6e501f", "#d2454f", "#639bf3", "#00e47b"};
color[i_Integer] /; 1<=i<=9 := colors[[i]];
color[i_Integer] /; i>9 := (Print["Color undefined"]; Abort[];)
color[0] := RGBColor["#000000"];

ColorToTeX[RGBColor[c__], name_String] := StringJoin["\\definecolor{", name, "}{rgb}", TextString[{c}]]


(* Mathematica Default markers *)
markers = {{"\[FilledCircle]",8.96`},{"\[FilledSquare]",8.96`},{"\[FilledDiamond]",10.88`},{"\[FilledUpTriangle]",10.24`},{"\[FilledDownTriangle]",10.24`},{"\[EmptyCircle]",10.24`},{"\[EmptySquare]",10.24`},{"\[EmptyDiamond]",10.24`},{"\[EmptyUpTriangle]",11.136`},{"\[EmptyDownTriangle]",11.136`}};
marker[i_Integer] := markers[[Mod[i-1, Length[markers]]+1]];


Themes`AddThemeRules["MishoStyle", Join[
  Charting`ResolvePlotTheme["Detailed", "ListLogPlot"] /.  {
    List[Directive[__],___] -> Thread@Directive[colors, AbsoluteThickness[1.6`]]
  },
  {LabelStyle->Black,
    BaseStyle -> {Black, FontFamily->"Times New Roman", FontSize->15}
  }]];
SetOptions[#,
  PlotTheme -> "MishoStyle",
  ImageSize -> {Automatic, 250}
] &/@ {Plot, LogPlot, LogLogPlot, LogLinearPlot, ListPlot, ListLogPlot, ListLogLogPlot, ListLogLinearPlot};

<<MaTeX`
SetOptions[MaTeX, FontSize -> 16, "Preamble"->{"\\usepackage{newtxtext,newtxmath,color}"}, ContentPadding->False];

TeXParamAligned[params_List] := "\\begin{aligned}" <> StringRiffle[#[[1]] <> "&=" <> If[Head[#[[2]]]===String,#[[2]],MyTextString[#[[2]]]] &/@ params, "\\\\"] <> "\\end{aligned}"
TeXParamRow[params_List] := StringRiffle[#[[1]] <> "=" <> If[Head[#[[2]]]===String,#[[2]],MyTextString[#[[2]]]] &/@ params, ",\\ \\ "]

Options[MyChartingScaledTicks] := {"RawRange" -> {0.1, 10}, "Separator" -> "\[CenterDot]"};
MyChartingScaledTicks[arg_, OptionsPattern[]] := MapAt[Module[{form},
    form[n_] := If[Between[n, OptionValue["RawRange"]], TextString[n],
      Module[{e = Floor[Log10[n]], r},
        r = TextString[n/10^e] // StringReplace["." ~~ EndOfString -> ""];
        If[r == "1", Superscript[10, e], Row[{r, Superscript[10, e]}, OptionValue["Separator"]]]]];
    Replace[#, {
      n_?NumericQ | NumberForm[n_, _] :> form[n],
      Superscript[a_, b_] :> form[a^b],
      Row[{a_, Superscript[b_, c_]}, _] :> form[a*b^c]
    }]] &, Charting`ScaledTicks[arg][##], {All, 2}] &
