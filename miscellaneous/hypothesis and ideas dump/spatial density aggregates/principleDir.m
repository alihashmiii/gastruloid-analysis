(* ::Package:: *)

BeginPackage["principleDir`"]


maskToRegion::usage = "given a mask, maskToRegion[mask_] outputs an ordered list of pixels of the mask boundary 
as well as converts the mask to a region.";
findSym::usage = "findSym[mask_,option_:False] takes in a mask and outputs the principle directions for the mask as
{centroid, boundarypt1, boundarypt2}. Graphic displayed with option set to True."
bisectMask::usage = "bisectMask[mask_,option_:False] splits the aggregate mask into two halves defined by the principle directions. 
With option set to True the graphic is displayed."


Begin["`Private`"]


findSym[mask_,option_:False]:=Block[{imd,nx,ny,points,weights,mass,centerOfMass,pointsCom,diag,inertia,v1,v2},
imd = Transpose@ImageData[mask];(* here we transpose the image data *)
{nx,ny} = Dimensions[imd];(* dimensions of the image *)
points = Tuples[{Range[nx],Reverse@Range[ny]}]; (* x-y coordinates on the image *)
weights = Flatten[imd]; (* vector of imagedata *)
centerOfMass=ImageMeasurements[mask,"IntensityCentroid"];(* here we calculate the x-y coordinate of the center of mass *)

pointsCom=(Subtract[#,centerOfMass]&/@points) weights; (* we subtract from each point the center of mass. This gives us
 a vector of n number of points the {x,y} coordinates of which are displacements from the centroids. We multiply this resulting
 matrix with the intensity vector. This gives us a vector consisting of two coordinates which are x,y displacments from centroid 
 for each point in the grid *)

diag=Total[(#.#)&/@pointsCom]; (* sum of the {x,y}.{x,y} from above. gives a single values *)
inertia=Table[KroneckerDelta[i,j] diag-Total[(#[[i]] #[[j]])&/@pointsCom],{i,2},{j,2}]; (* inertia *)

{v1,v2}=Eigenvectors[inertia]; (* extracting eigenvectors from inertia *)
(* Print["center of Mass: ", centerOfMass,"\n","pt1 : ",v1 + centerOfMass,"\n","pt2 : ",v2 + centerOfMass]; *)

If[option, Print@Show[mask, Graphics[{{Thick,Dashed,XYZColor[0,0,1,0.6],InfiniteLine[{v2+centerOfMass,centerOfMass}]},
{Thick,Dashed,XYZColor[1,0,0,0.6],InfiniteLine[{v1+centerOfMass,centerOfMass}]},{Darker@XYZColor[0,1,0,0.5],PointSize[Large],
Point[centerOfMass]},Red,Point[v1+centerOfMass],Blue,Point[v2+centerOfMass]}]]
];
{centerOfMass,v1,v2}
];


maskToRegion[mask_]:= Block[{pixelvaluepos,tour,pixelposordered},
pixelvaluepos= PixelValuePositions[MorphologicalPerimeter@mask,1];
tour=Last@FindShortestTour[pixelvaluepos];
pixelposordered = Part[pixelvaluepos,tour];
{pixelposordered,BoundaryMeshRegion[N@pixelvaluepos,Line@tour]}
];


bisectMask[mask_,option_:False]:= Block[{centerofMass,v1,v2,pixelordered,\[ScriptCapitalR]1,\[ScriptCapitalR]2,temp,pt1,pt2,
list,trimList,npixelordered,nf,ptslocation,firstPreMask,secondPreMask,dim},
dim = First@ImageDimensions[mask];
{centerofMass,v1,v2} = findSym[mask];
{pixelordered,\[ScriptCapitalR]1} = maskToRegion[mask];
\[ScriptCapitalR]2=InfiniteLine[{centerofMass,centerofMass+v2}];
list=(RegionIntersection[\[ScriptCapitalR]2,#]&/@MeshPrimitives[\[ScriptCapitalR]1,1]//Cases[#,Point[arg_]:> arg]&);

trimList[ls_]:= Module[{p=ls,pairs,firstelem = First@ls,dist,maxdist},
pairs={firstelem,#}&/@p;
dist=EuclideanDistance[##]&@@@pairs;
maxdist=Max@dist;
pairs[[First@@Position[dist,maxdist]]]
]/;Length@ls>2;

{pt1,pt2}=trimList[list]/.trimList[arg__]:> arg;

If[option,temp=Show[HighlightMesh[\[ScriptCapitalR]1,{Style[1,Directive[XYZColor[0,0,1,0.9]]],Style[2,Directive[XYZColor[0,0,1,0.18]]]}],
Graphics[{Red,PointSize[Large],Point@pt1,PointSize[Large],Darker@Green,
Point@pt2,Thick,Dashed,XYZColor[1,0,0,0.5],\[ScriptCapitalR]2,Blue,Point@centerofMass}]
]
];

npixelordered = N@pixelordered;
nf=Nearest[npixelordered];
ptslocation=Flatten@Position[npixelordered,Flatten@nf[pt1]|Flatten@nf[pt2]];

firstPreMask=PadRight[#,Length@#+1,#]&@Part[pixelordered,First@ptslocation;;Last@ptslocation]//Graphics[{White,Line@#},
PlotRange->{{0,dim},{0,dim}}]&//Rasterize[#,"Image",Background->Black,ImageSize->{{dim},{dim}}]&//FillingTransform@*Binarize;
secondPreMask=Join[pixelordered[[;;First@ptslocation+1]],pixelordered[[Last@ptslocation+1;;]]]//Graphics[{White,Line@#},
PlotRange->{{0,dim},{0,dim}}]&//Rasterize[#,"Image",Background->Black,ImageSize->{{dim},{dim}}]&//FillingTransform@*Binarize;

Switch[option,True,Join[{temp},{centerofMass},pt1,pt2],
_, Join[{centerofMass},pt1,pt2]]
];


End[];
EndPackage[];
