(* ::Package:: *)

BeginPackage["masktoRegion`"]


findSym::usage = "findSym[mask_,option_:False] takes a binarized mask of the aggregate and yields center of mass and the eigenvectors";
bisectMask::usage = "bisectMask[mask_,option_:False] "
rotateMask::usage = "rotateMask[mask_?ImageQ] takes the mask and rotates the mask for alignment"
extractRadius::usage = "extractRadius[mask_?ImageQ] rotates the 2D image to generate and 3D mesh. Computes the equivalent radius"
equivalentRadii::usage = "equivalentRadii[directory_String] takes the name of the directory containing the masks and outputs the equivalent radius 
along with the respective timestamp";


Begin["`Private`"]


findSym[mask_,option_:False]:=Block[{imd,nx,ny,points,weights,mass,centerOfMass,pointsCom,diag,inertia,v1,v2},
imd=Transpose@ImageData[mask];(* here we transpose the image data *)
{nx,ny}=Dimensions[imd];(* dimensions of the image *)
points=Tuples[{Range[nx],Reverse@Range[ny]}]; (* x-y coordinates on the image *)
weights=Flatten[imd]; (* vector of imagedata *)

centerOfMass=ImageMeasurements[mask,"IntensityCentroid"];(* here we calculate the x-y coordinate of the center of mass *)

pointsCom=(Subtract[#,centerOfMass]&/@points) weights; (* we subtract from each point the center of mass. This gives us a vector of n number of points the {x,y}
 coordinates of which are displacements from the centroids. We multiply this resulting matrix with the intensity vector. This gives us a vector consisting of two
 coordinates which are x y displacments from centroid * intensity for each point in the grid *)

diag=Total[(#.#)&/@pointsCom]; (* sum of the {x,y}.{x,y} from above. gives a single values *)
inertia=Table[KroneckerDelta[i,j] diag-Total[(#[[i]] #[[j]])&/@pointsCom],{i,2},{j,2}]; (* inertia *)

{v1,v2}=Eigenvectors[inertia]; (* extracting eigenvectors from inertia *)
(* Print["center of Mass: ", centerOfMass,"\n","pt1 : ",v1 + centerOfMass,"\n","pt2 : ",v2 + centerOfMass]; *)

(* display *)
If[option,
Print@Show[mask,Graphics[{{Thick,Dashed,XYZColor[0,0,1,0.6],InfiniteLine[{v2+centerOfMass,centerOfMass}]},
{Thick,Dashed,XYZColor[1,0,0,0.6],InfiniteLine[{v1+centerOfMass,centerOfMass}]},{Darker@XYZColor[0,1,0,0.5],PointSize[Large],Point[centerOfMass]},
Red,Point[v1+centerOfMass],Blue,Point[v2+centerOfMass]}]]
];
{centerOfMass,v1,v2}
];


bisectMask[mask_,option_:True]:= Block[{centerofMass,v1,v2,pixelordered,\[ScriptCapitalR]1,\[ScriptCapitalR]2,temp={},pt1,pt2,
list,trimList,npixelordered,nf, ptslocation,firstPreMask,secondPreMask},
{centerofMass,v1,v2}=findSym[mask];
\[ScriptCapitalR]1=ImageMesh[mask];
pixelordered=MeshPrimitives[\[ScriptCapitalR]1,0]/.Point[x_]:> x;
\[ScriptCapitalR]2=InfiniteLine[{centerofMass,centerofMass+v2}];
list=(RegionIntersection[\[ScriptCapitalR]2,#]&/@MeshPrimitives[\[ScriptCapitalR]1,1]//Cases[#,Point[arg_]:> arg]&);

trimList[ls_]:= Module[{p=ls,pairs,firstelem = First@ls,dist,maxdist},
 pairs={firstelem,#}&/@p;
 dist=EuclideanDistance[##]&@@@pairs;
 maxdist=Max@dist;
 pairs[[First@@Position[dist,maxdist]]]
]/;Length@ls>2;

{pt1,pt2}=trimList[list]/.trimList[arg__]:> arg;

If[option,
temp = Show[HighlightMesh[\[ScriptCapitalR]1,{Style[1,Directive[XYZColor[0,0,1,0.9]]],Style[2,Directive[XYZColor[0,0,1,0.18]]]}],
Graphics[{Red,PointSize[Large],Point@pt1,PointSize[Large],Darker@Green,Point@pt2,Thick,Dashed,XYZColor[1,0,0,0.5],\[ScriptCapitalR]2,Blue,Point@centerofMass}]]
];

npixelordered = N@pixelordered;

nf=Nearest[npixelordered];
ptslocation=Flatten@Position[npixelordered,Flatten@nf[pt1]|Flatten@nf[pt2]];
Join[{temp},{centerofMass},pt1,pt2]
];


rotateMask[mask_?ImageQ]:= Module[{temp,centroid,pts1,pts2,grad,perimeterpts,norm,\[Theta],
rotateFn,canvas,rotatedpts,dim = ImageDimensions@mask,translationFn,transpts},
{temp,centroid,pts1,pts2}=bisectMask[mask];
perimeterpts= PixelValuePositions[MorphologicalPerimeter@mask,1];
norm = Normalize[Subtract[pts2,pts1]];
If[First@pts2 == First@pts1, Abort[]];
grad=(pts2[[2]]-pts1[[2]])/(pts2[[1]]-pts1[[1]]);
\[Theta] = VectorAngle[norm,{0,1}](180/\[Pi]);
\[Theta] = Which[grad>=0, If[\[Theta]>90,-\[Theta],\[Theta]], 
 grad<0, If[\[Theta]<90,-\[Theta],\[Theta]]];
rotateFn=RotationTransform[\[Theta] Degree,norm];
canvas=ConstantImage[0,dim];
rotatedpts=rotateFn[perimeterpts~Append~centroid];
translationFn = TranslationTransform[dim/2-Last@rotatedpts];
transpts= Abs[translationFn@Most@rotatedpts];
ReplacePixelValue[canvas,transpts-> 1]//FillingTransform
];


extractRadius[mask_?ImageQ]:=Block[{mesh,centroid,mesh2Dpts,rotateTransform,mesh3Dpts,hullmesh,r},
mesh = ImageMesh[mask];
centroid = Append[RegionCentroid[mesh],0];
mesh2Dpts = MeshPrimitives[mesh,0]/.Point[x_]:> x;
rotateTransform = RotationTransform[3.Degree,{0,1,0},centroid];
mesh3Dpts = NestList[rotateTransform,ArrayFlatten[{{mesh2Dpts,0}}],59];
hullmesh = ConvexHullMesh[Flatten[mesh3Dpts,1]];
Print[mask,BoundaryMeshRegion[hullmesh,ImageSize->Tiny,Boxed->True]];
Values@(First@@Solve[4/3Pi r^3 ==Volume[hullmesh],r\[Element]Reals]) (* equivalent radius *)
];


Fn[arg_]:=With[{pixelsize=0.633},
MapThread[{#1, First@*Values@ComponentMeasurements[MorphologicalComponents@#2,"EquivalentDiskRadius"]*pixelsize}&,arg]
];


equivalentRadii[directory_String]:=With[{pixelsize=0.633}, 
Module[{files,filteredfiles,timestamps,masks,deformedmaskpos,deformedmasks,rotatedmasks,
maskseries1,maskseries2,timestamps1,timestamps2,aggregate1R,aggregate2R,repmasks},
files = Import[directory];
filteredfiles = Flatten@StringCases[files,__~~"Mask.tif"];
timestamps = Flatten[StringCases[filteredfiles,x:__ ~~" hours"~~___ :> ToExpression@x]];
masks = Import[directory<>#]&/@(filteredfiles);
deformedmaskpos = Position[timestamps,_?(#>72&)];
deformedmasks = Extract[masks,deformedmaskpos];
Which[
deformedmaskpos!={}, (rotatedmasks=rotateMask/@deformedmasks;
repmasks=ReplacePart[masks,Thread[deformedmaskpos-> rotatedmasks]];
{maskseries1,timestamps1}=Take[#,Min@deformedmaskpos-1]&/@{repmasks,timestamps};
{maskseries2,timestamps2} = Drop[#,Min@deformedmaskpos-1]&/@{repmasks,timestamps};
aggregate1R=Fn[{timestamps1,maskseries1}];
aggregate2R=Thread[{timestamps2,Composition[Times[#,pixelsize]&,extractRadius]/@rotatedmasks}];
aggregate1R~Join~aggregate2R),
True,
Fn[{timestamps,masks}]]
 ]
];


End[];
EndPackage[];
