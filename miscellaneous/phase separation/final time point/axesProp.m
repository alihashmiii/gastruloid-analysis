(* ::Package:: *)

BeginPackage["axes`",{"principleDir`"}]


bresenhamLine::usage = "brensenhamLine[pt1,pt2] uses two pixel points to draw a line between them";
elongationVector::usage = "elongationVector[mask_,div_:5] yields opposing points on the aggregate boundary with maximal elongation.
A default angular sweep of 5 degree is used";
bresenhamVector::usage = "gives the bresenham pixel points along which signal is read";


Begin["`Private`"]


bresenhamLine[p0_,p1_]:=Module[{dx,dy,sx,sy,err,newp},
{dx,dy}=Abs[p1-p0];
{sx,sy}=Sign[p1-p0];
err=dx-dy;
newp[{x_,y_}]:=With[{e2=2 err},
{If[e2>-dy,err-=dy;x+sx,x], If[e2<dx,err+=dx;y+sy,y]}
];
NestWhileList[newp,p0,#=!=p1&,1]
];


elongationVector[mask_,div_:5]:=Module[{centroid,pt1,pt2,boundary,tour,boundaryc,nF,\[ScriptCapitalR]1,\[ScriptCapitalR]2,
intersections,rotationTransform,pts,\[Theta],r,infinitelines,contour,surfacepts,pos,surfacepts2,
diameter,p,maxD},
{centroid,pt1,pt2}=bisectMask[mask];
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
surfacepts = Graphics`Mesh`FindIntersections[Graphics[{FaceForm[],EdgeForm[Blue],contour,#}]]&/@infinitelines;
pos = Position[surfacepts,_?(Length@#>2&),{1}];
surfacepts2=MapAt[Sequence@@TakeLargestBy[Permutations[#,{2}],EuclideanDistance[#[[1]],#[[2]]]&,1]&,surfacepts,pos];
MaximalBy[surfacepts2,Apply[EuclideanDistance]]
]


bresenhamVector[mask_]:=Module[{pts,imgPeripts,nf,ptsB,perms,extrema,pointsL,subdiv},
pts=bisectMask[mask];
imgPeripts=First@@Values@ComponentMeasurements[mask,"PerimeterPositions"];
nf=Nearest[imgPeripts];
ptsB={Round@*First@pts}~Join~Flatten[Round@*nf/@Rest[pts],1];
perms=Permutations[ptsB,{2}];
extrema = Extract[perms,FirstPosition[#,Max@#]&[N@EuclideanDistance[Sequence@@#]&/@perms]];
pointsL=bresenhamLine[First@extrema,Last@extrema];
subdiv=N@Subdivide[Length@pointsL];
{subdiv,pointsL}
];


End[];
EndPackage[];
