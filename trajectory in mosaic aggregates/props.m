(* ::Package:: *)

BeginPackage["props`", {"principleDir`"}]


elongation::usage = "elongation[mask_,OptionsPattern[]] computes the maximal elongation of the aggregate mask";
MSD::usage = "MSD[ptslist_,ntime_:2] -> the mean square displacement of the particle trajectory. ntime is the number of times which the trajectory
length is reduced for calculating MSD.";
centroidsIntensity::usage = "measures the coordinates of the intensity weighted centroid";
coloredTraj::usage = "";
trajVid::usage = "";
histogramDist::usage = "";
extractSignal::usage = "";
ExtractSignal::usage = "";
centroidsIntensityM2::usage = "";


Begin["`Private`"]


Options[elongation]={"div"-> 5, "printBisect" -> False};
elongation[mask_,OptionsPattern[]]:=With[{div=OptionValue["div"], printBisect=OptionValue["printBisect"]},
Block[{p,centroid,pt1,pt2,boundary,tour,boundaryc,nF,\[ScriptCapitalR]1,\[ScriptCapitalR]2,intersections,rotationTransform,pts,
\[Theta],r,infinitelines,contour,surfacepts,pos,surfacepts2,surfacepts3,diameter,perimeter,area,elongation},
{centroid,pt1,pt2}=bisectMask[mask,printBisect];
boundary=PixelValuePositions[MorphologicalPerimeter@mask,1];
tour=Last@FindShortestTour[boundary];
boundaryc=Part[boundary,tour];
nF=Nearest[boundaryc]; (* Nearest Function *)
\[ScriptCapitalR]1=BoundaryMeshRegion[boundary,Line@tour]; (* boundary *)
\[ScriptCapitalR]2=InfiniteLine[{centroid,pt1}]; (* pseudo axis of symmetry *)
intersections=Cases[RegionIntersection[\[ScriptCapitalR]2,#]&/@MeshPrimitives[\[ScriptCapitalR]1,1],_Point]/.Point[{x_}]:> x;
rotationTransform=RotationTransform[div Degree,centroid];
pts=Rest@NestList[rotationTransform,pt1,360/(div*2)];
infinitelines=InfiniteLine[{centroid,#}]&/@pts;
contour= Polygon@boundaryc;
(*surfacepts=(x\[Function]RegionIntersection[x,contour])/@infinitelines/.Point[x_]:>x; *)
surfacepts = Graphics`Mesh`FindIntersections[Graphics[{FaceForm[],EdgeForm[Blue],contour,#}]]&/@infinitelines;
pos = Position[surfacepts,_?(Length@#>2&),{1}];
surfacepts2=MapAt[Sequence@@TakeLargestBy[Permutations[#,{2}],EuclideanDistance[#[[1]],#[[2]]]&,1]&,surfacepts,pos];
\[Theta] = Array[# &, 360/div,{0, 2\[Pi]}];
p/:r[p]:= Join[#[[;;,1]],#[[;;,2]]]&@Map[EuclideanDistance[centroid,#]&,surfacepts2,{2}];
diameter=Plus@@Partition[r[p],Length[r[p]]/2];
{perimeter,area,elongation}=Values@@ComponentMeasurements[mask,{"PerimeterLength","Area","Elongation"}];
{perimeter,area,elongation,MinMax@diameter, MinMax@r[p]}]
];


MSD[ptslist_, ntime_:2]:= Block[{},
Map[x \[Function] With[{timesteps = Round[Length@ptslist[[x]]/ntime]},
Module[{threadedlist},
threadedlist = Table[Thread[{ptslist[[x]],
RotateLeft[ptslist[[x]],i]}][[;;-(i+1)]],{i,timesteps}];
Mean/@Map[Function[x,Power[#,2]&@
(Norm@(Subtract@@Reverse@x))],threadedlist,{2}]
]
], Range@Length@ptslist]
];


MSD[pos_,t_:50,"M2"]:=With[{timesteps=t},
Module[{transposedlist,length},
transposedlist=Table[
Thread[{pos,RotateLeft[pos,i]}][[;;-(i+1)]],{i,timesteps}];
length=(Dimensions/@transposedlist);
Mean/@Table[
Module[{listreshape,lr},
listreshape=ArrayReshape[transposedlist[[i]],length[[i]]];
lr=Map[Function[x,Power[#,2]&@(EuclideanDistance@@Reverse@x)],listreshape[[;;All]]]
],{i,timesteps}]
]
];


extractSignal[images:{__Image}]:=Block[{threshim,img=#,immean,im,imstd,imbin},
threshim=Threshold[img,{"Soft","Cluster"}];
immean=Mean@threshim; imstd=StandardDeviation@img;
im = MedianFilter[ImageSubtract[threshim,immean],2];
imbin=MorphologicalBinarize@im;
DeleteBorderComponents[imbin]*im
]&/@images;


ExtractSignal[images:{__Image}]:=
Block[{threshim,immean,im,imstd,imbin,masks,signal,comp,mediansize},
{masks,signal}=Transpose[(threshim=Threshold[#,{"Soft","Cluster"}];
immean=Mean@threshim; imstd=StandardDeviation@threshim;
im = MedianFilter[ImageSubtract[threshim,immean+0.5imstd],2];
imbin=MorphologicalBinarize@im;
{DeleteBorderComponents[imbin],im})&/@images];
comp=MorphologicalComponents[#]&/@masks;
mediansize=Median@Flatten[Values@ComponentMeasurements[#,"Count"]&/@comp,1];
MapThread[#1*#2&,{DeleteSmallComponents[#,mediansize]&/@masks,signal}]
];


centroidsIntensityM2[images_,masks_,NLargestValues_]:= MapThread[
ImageMeasurements[Threshold[(#1*#2),{"LargestValues",NLargestValues}], "IntensityCentroid", Masking->#2]&,
{images,masks}];


centroidsIntensity[images:{_Image..}]:=Block[{img=#},
ImageMeasurements[img,"IntensityCentroid"]
]&/@images;


coloredTraj[pts_,ptsSize_]:= Graphics[{PointSize[ptsSize],
Point[#, VertexColors -> ColorData["Rainbow"]/@Rescale@Range@Length@#]&@pts,Line@pts}];


trajVid[images_,centroids_]:= Module[{comb},
comb = MapThread[
Show[#1,Graphics[{Green,PointSize[0.025],Point[#2]}]]&,{images,centroids}];
ListAnimate@comb
];


histogramDist[disp_,dist_,range_:1]:=Show[
Histogram[disp,Automatic,"PDF",PlotRange->{{0,range},{0,Automatic}}]/.RGBColor[__]-> Opacity[0.15,Blue],
Plot[PDF[dist,x],{x,0,range},PlotStyle->Opacity[0.2,Purple],Filling-> Bottom,PlotRange->All],
AxesLabel->{Style["velocity (um/s)",{Black,Bold,16}],Style["PDF",{Black,Bold,16}]},
AxesStyle->Directive[Black, 16],ImageSize-> Large
];


End[];
EndPackage[];
