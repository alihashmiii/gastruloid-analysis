(* ::Package:: *)

(* ::Section:: *)
(*Begin Package*)


BeginPackage["PolElong`", {"principleDir`"}]


contrastElong::usage = "contElong[dir_, maskStr_,sigStr_,OptionsPattern[]][serial_] takes in a directory address as well as string-prefixes for the mask folder and
the signal.tif stack and the aggregate ID to output aggregate property:
ID->{perimeter,area,elongation,MinMax[diameter], MinMax[radius],MinMax[ contrast[p][[All,1]] ], MinMax[ contrast[p][[All,2]] ]}";


plotFn::usage = "plotFn[res_] uses the result of contElong to generate plot(s) for elongation/contrast vs time";


contrastTheta::usage = "aux Fn to calculate the contrast by sweeping pseudo-axis of symmetry";


Begin["`Private`"];


(* ::Section:: *)
(*Main Functions*)


(* ::Subsection:: *)
(*contrast \[Theta]*)


(* ::Subsubsection:: *)
(*Bresenham line algorithm*)


bresenhamLine[p0_,p1_]:=Module[{dx,dy,sx,sy,err,newp},{dx,dy}=Abs[p1-p0];
{sx,sy}=Sign[p1-p0];
err=dx-dy;
newp[{x_,y_}]:=With[{e2=2 err},{If[e2>-dy,err-=dy;x+sx,x],If[e2<dx,err+=dx;y+sy,y]}];
NestWhileList[newp,p0,#=!=p1&,1]
];


(* ::Subsubsection:: *)
(*contrast \[Theta]*)


(* Options[contrastTheta]={"imageSize"->{1024,1024}};
contrastTheta[boundarycoord_,image_,surfacepts_,nF_,OptionsPattern[]]:=With[{nearest=nF,
img=ConstantImage[0,OptionValue["imageSize"]]},
Block[{ptb1,ptb2,pos1,pos2,fmask,smask,bresen,arc1,arc2,r1,r2,contrast},
({ptb1,ptb2}=First@*nearest/@#;
bresen = bresenhamLine[ptb1,ptb2];
{pos1,pos2}=First/@ Position[boundarycoord,ptb1|ptb2][[ ;; 2]];
arc1 = boundarycoord[[pos1 ;; pos2]];
arc2 = Join[boundarycoord[[pos2 ;; -2]], boundarycoord[[;; pos1]]];
{fmask,smask}=FillingTransform@ReplacePixelValue[img,Join[#,bresen]-> 1]&/@{arc1,arc2};
{r1,r2}=ImageMeasurements[image*#,"TotalIntensity",Masking->#]&/@{fmask,smask};
contrast={(Abs[Subtract@@#]/Total[#]),Abs[Subtract@@#]}&@{r1,r2})&/@surfacepts]
]; *)


Options[contrastTheta]={"imageSize"->{1024,1024}};
contrastTheta[boundarycoord_,image_,surfacepts_,nF_,OptionsPattern[]]:=With[{nearest=nF, imagesize1 = First@OptionValue["imageSize"],
imagesize2=Last@OptionValue["imageSize"]},
Block[{ptb1,ptb2,pos1,pos2,fmask,smask,r1,r2,contrast},
({ptb1,ptb2}=First@*nearest/@#;
{pos1,pos2} = Position[boundarycoord,ptb1|ptb2][[;;2]];
fmask = PadRight[#,Length@#+1,#]&@Part[boundarycoord, First[pos1] ;; First[pos2]]//Graphics[{White,Line@#},
PlotRange->{{0,imagesize1},{0,imagesize2}}]&//Rasterize[#,"Image",Background->Black,
ImageSize->{{imagesize1},{imagesize2}}]&//FillingTransform@*Binarize;
(* smask=Join[boundarycoord[[pos2 ;; -2]], boundarycoord[[ ;; pos1]]]//Graphics[{White,Polygon@#},
PlotRange->{{0,imagesize1},{0,imagesize2}}]&//Rasterize[#,"Image",Background->Black,
ImageSize->{{imagesize1},{imagesize2}}]&//FillingTransform@*Binarize; *)
smask=Join[boundarycoord[[;;First@pos1+1]],boundarycoord[[First@pos2+1;;]]]//Graphics[{White,Line@#},
PlotRange->{{0,imagesize1},{0,imagesize2}}]&//Rasterize[#,"Image",Background->Black,
ImageSize->{{imagesize1},{imagesize2}}]&//FillingTransform@*Binarize;
{r1,r2}=ImageMeasurements[image*#,"TotalIntensity",Masking->#]&/@{fmask,smask};
contrast = {(Abs[Subtract@@#]/Total[#]),Abs[Subtract@@#]}&@{r1,r2})&/@surfacepts
]
];


(* ::Subsection:: *)
(*Mains*)


Options[contrastElong]={"div"-> 5, "printBisect" -> False};
contrastElong[dir_, maskStr_,sigStr_,OptionsPattern[]][serial_]:=With[{div=OptionValue["div"], 
printBisect=OptionValue["printBisect"]},
Block[{files,filesind,strserial = ToString[serial]},
files=Import[dir<>strserial<>"\\"<> maskStr <>strserial];
filesind=Sort@Flatten@StringCases[files,(x:DigitCharacter..):>FromDigits[x]];

Block[{mask,signal,ind=#,p,contrast,centroid,pt1,pt2,boundary,tour,boundaryc,nF,\[ScriptCapitalR]1,\[ScriptCapitalR]2,intersections,rotationTransform,pts,
\[Theta],r,infinitelines,contour,surfacepts,pos,surfacepts2,diameter,perimeter,area,elongation},
Print[ind];
mask = Import[dir<>strserial<>"\\"<> maskStr<>strserial<>"\\Mask"<>ToString[ind]<>First@global$ext];
signal = ImageAdjust@Import[dir<>strserial<>"\\"<>sigStr<>"\\"<>ToString[ind]<>Last@global$ext];
{centroid,pt1,pt2} = bisectMask[mask,printBisect];
boundary = PixelValuePositions[MorphologicalPerimeter@mask,1];
tour = Last@FindShortestTour[boundary];
boundaryc = Part[boundary,tour];
nF = Nearest[boundaryc]; (* Nearest Function *)
\[ScriptCapitalR]1 = BoundaryMeshRegion[boundary,Line@tour]; (* boundary *)
\[ScriptCapitalR]2 = InfiniteLine[{centroid,pt1}]; (* pseudo axis of symmetry *)
intersections = Cases[RegionIntersection[\[ScriptCapitalR]2,#]&/@MeshPrimitives[\[ScriptCapitalR]1,1],_Point]/.Point[{x_}]:> x;
rotationTransform = RotationTransform[div Degree,centroid];
pts = Rest@NestList[rotationTransform,pt1,360/(div*2)];
infinitelines = InfiniteLine[{centroid,#}]&/@pts;
contour = Polygon@boundaryc;
(*surfacepts=(x\[Function]RegionIntersection[x,contour])/@infinitelines/.Point[x_]:>x; *)
surfacepts = Graphics`Mesh`FindIntersections[Graphics[{FaceForm[],EdgeForm[Blue],contour,#}]]&/@infinitelines;
pos = Position[surfacepts,_?(Length@#>2&),{1}];
surfacepts2 = MapAt[Sequence@@TakeLargestBy[Permutations[#,{2}], EuclideanDistance[#[[1]],#[[2]]]&,1]&,surfacepts,pos];
\[Theta] = Array[# &, 360/div,{0, 2\[Pi]}];
p/:r[p]:= Join[#[[;;,1]],#[[;;,2]]]&@Map[EuclideanDistance[centroid,#]&,surfacepts2,{2}];
p/:contrast[p] = contrastTheta[boundaryc,signal,surfacepts2,nF];
diameter = Plus@@Partition[r[p],Length[r[p]]/2];
{perimeter,area,elongation} = Values@@ComponentMeasurements[mask,{"PerimeterLength","Area","Elongation"}];
ind -> {perimeter,area,elongation,MinMax@diameter, MinMax@r[p],MinMax[ contrast[p][[All,1]] ], MinMax[ contrast[p][[All,2]] ]}
]&/@filesind
]
];


(* ::Subsection:: *)
(*plot contrast/elongation for aggregate*)


Options[plotFn]={"contrast" -> "normalized","join"-> True,"filter"-> 5,
"filling"-> {1->{{2},{LightGreen,LightBlue}}},"imgsize"-> 200,"rescale"-> True};
plotFn[res_,harmonics_,OptionsPattern[]]:=With[{secondharmonic=harmonics,size=OptionValue["imgsize"],
rad=OptionValue["filter"],j=OptionValue["join"],fill=OptionValue["filling"]},
 Module[{contrast, radratio, g, ind1, len = Replace[res,{p:{HoldPattern[_-> _]..}:>Length@{p},x_ :> Length[x]}],
  results = res},
If[len == 1, results = {res}];
ind1 = Switch[OptionValue["contrast"],"normalized", -2, _ , -1];
Function[x,(contrast = Max /@ ((Values /@ results[[x]])[[All, ind1]]);
If[OptionValue["rescale"],
Show[
ListPlot[{Thread[{If[#==1, 72, N[#/6+72,4]]&/@Keys@results[[x]], #}]&[GaussianFilter[Rescale@contrast,rad]],
Thread[{Part[secondharmonic[[x]],All,1],GaussianFilter[Rescale@secondharmonic[[x]][[All,2]],rad]}]} ,
 Joined ->j, PlotStyle -> {{Thick,RGBColor[0.07924700513978256, 0.07394564469935314, 0.4499539624003866]},
 {Thick,RGBColor[0.10823426242933962`, 0.6204150625738422, 0.4678154145759483]}},Filling->fill,PlotRange->{{72,92},{0.,1.}},ImageSize->size,
 PlotLegends-> {"polarization","tip elongation"},AxesLabel->{"hrs pp",""},AxesStyle->Directive[Black,Bold,12]],
Graphics[{Opacity[0.5],Dashed,Thick,Black,Line[{{82.5,0.0},{82.5,1.0}}]}]
],
{ListPlot[Thread[{If[#==1, 72, N[#/6+72,4]]&/@Keys@results[[x]], #}]&[GaussianFilter[contrast,rad]],Joined->True,
PlotStyle->{Thick,RGBColor[0.07924700513978256, 0.07394564469935314, 0.4499539624003866]},AxesStyle->Directive[Black,Bold,12],AxesLabel->{"hrs pp","polarization"},Filling->Axis,
PlotRange->{{72,92},{0,Automatic}}],
ListPlot[Thread[{Part[secondharmonic[[x]],All,1],GaussianFilter[secondharmonic[[x]][[All,2]],rad]}],Joined->True,
PlotStyle->{Thick,RGBColor[0.10823426242933962`, 0.6204150625738422, 0.4678154145759483]},AxesStyle->Directive[Black,Bold,12],AxesLabel->{"hrs pp","elongation"},Filling->Axis,
PlotRange->{{72,92},{0,Automatic}}]}])]/@Range[Length@res]
 ]
]


global$ext = {".tif",".tif"};


(* ::Section:: *)
(*End Package*)


End[];


EndPackage[]
