(* ::Package:: *)

(* ::Chapter:: *)
(*ARAPlib.m (Mathematica Modules 2016/12/20)*)
(*Kyushu University 2016*)


(* ::Section::Closed:: *)
(*Basic Functions*)


MakePolygon[V_,tindex_]:=Module[{n,cl},
n=Length[tindex];
cl=ColorData[35,"ColorList"];
MapIndexed[{ColorData["TemperatureMap"][If[n ==1,1,(#2[[1]]-1)/(n-1)]],
Polygon[{V[[#1[[1]]]],V[[#1[[2]]]],V[[#1[[3]]]]}]}&,tindex]];
(*HorizontalSnake[n_]:=Flatten[Table[{{i,0},{i,1}},{i,0,n-1}],1];
VerticalSnake[n_]:=Flatten[Table[{{1,i},{0,i}},{i,0,n-1}],1];
SnakeTriangle[n_]:=Flatten[Table[{{2i+1,2i+2,2(i+1)+1},{2(i+1)+1,2(i+1),2(i+2)}},
{i,0,n-2}],1];*)
AnimationRange[conf_]:={Min[#],Max[#]}&/@Transpose[conf[[1]]~Join~conf[[2]]];
(*l:list*)
AnimationRange2[range1_,range2_]:=AnimationRange[{Transpose[range1]~Join~Transpose[range2],{},{}}];


(*\:6975\:5ea7\:6a19\:5909\:63db\:3000a=(x,y)*)
Polar[a_]:=Module[{x,y,r,s},
x=a[[1]];y=a[[2]];
r=Sqrt[x*x+y*y];
s=If[x==0,If[y>=0,\[Pi]/2,-( \[Pi]/2)],ArcTan[x,y]];{r,s}];
(* \:30c7\:30ab\:30eb\:30c8\:5ea7\:6a19\:5909\:63db a=(r,s)*)
Cartesian[a_]:=Module[{x,y,r,s},
r=a[[1]];s=a[[2]];
x=r Cos[s];y=r Sin[s];{x,y}];
(* \:7dda\:5f62\:88dc\:9593 \:59cb\:70b9 a \:7d42\:70b9 b \:30d1\:30e9\:30e1\:30fc\:30bf t=0->1 *)
LinearInterpolate[a_,b_,t_]:=a+t(b-a);
LinearInterpolation[p_,q_,t_]:=Map[LinearInterpolate[#[[1]],#[[2]],#[[3]]]&,
Transpose[{p,q,Table[t,{Length[p]}]}]];
(* \:6975\:5ea7\:6a19\:88dc\:9593 \:59cb\:70b9 a \:7d42\:70b9 b \:30d1\:30e9\:30e1\:30fc\:30bf t=0->1 *)
PolarInterpolate[a_,b_,t_]:=Module[{pa,pb},
	pa=Polar[a];pb=Polar[b];
	If[pa[[2]]-pb[[2]]>\[Pi],pb[[2]]=pb[[2]]+2\[Pi],
    If[pb[[2]]-pa[[2]]>\[Pi],pb[[2]]=pb[[2]]-2\[Pi]]];
	Cartesian[LinearInterpolate[pa,pb,t]]];
PolarInterpolation[p_,q_,t_]:=Map[PolarInterpolate[#[[1]],#[[2]],#[[3]]]&,
  Transpose[{p,q,Table[t,{Length[p]}]}]];
(* \:70b9\:5217 pl \:3092 p \:3060\:3051\:5e73\:884c\:79fb\:52d5\:3059\:308b *)
pt[p_,pl_]:=Map[p+#&,pl];
(*LinearInterpolateSnake2D[n_,t_]:=MakePolygon[
LinearInterpolation[pt[{1,1},HorizontalSnake[n]],
pt[{1,1},VerticalSnake[n]],t],SnakeTriangle[n]];
PolarInterpolateSnake2D[n_,t_]:=MakePolygon[
PolarInterpolation[pt[{1,1},HorizontalSnake[n]],
pt[{1,1},VerticalSnake[n]],t],SnakeTriangle[n]]*)


(**)
PolarDecomposition[m_]:=Module[{u,w,v},
{u,w,v}=SingularValueDecomposition[m]//N;
{u.Transpose[v],v.w.Transpose[v]}];
PolarDecompositionPlus[m_]:=Module[{p,q,s,pt,qt},
{p,q}=PolarDecomposition[m];
pt=Transpose[p];qt=Transpose[q];
s=ArcTan[pt[[1]][[1]],pt[[1]][[2]]];
(* RotationMatrix[s]={{Cos[s],-Sin[s]},{Sin[s],Cos[s]}} *)
{RotationMatrix[s],Transpose[qt]} 
];
RotateAngle[m_]:=ArcTan[m[[1]][[1]],m[[2]][[1]]];
(* PolarDecomposition\:306e\:56de\:8ee2\:884c\:5217\:3092\:89d2\:5ea6\:3067\:8fd4\:3059\:95a2\:6570 *)


(**)
CogTrans[pl_]:=Map[-Cog[pl]+#&,pl];
PolygonToTriangles[l_]:=Map[CogTrans,Map[#[[2]][[1]]&,l]];
(*\:5e73\:884c\:79fb\:52d5*)
TrianglesToTriangles[offset_,cog_]:=
Map[pt[#[[1]],#[[2]]]&,Transpose[{cog,offset}]];
TrianglesToPolygon[tl_]:=Module[{n,cl},
n=Length[tl];
cl=ColorData[35,"ColorList"];
MapIndexed[{ColorData["TemperatureMap"][If[n>1,(#2[[1]]-1)/(n-1),1]],
Polygon[#]}&,tl]];


(* valiable names list with size n *)
ValNames[n_]:=Transpose[{Table[ToExpression["v"~StringJoin~ToString[i]~StringJoin~"x"],{i,1,n}],
Table[ToExpression["v"~StringJoin~ToString[i]~StringJoin~"y"],{i,1,n}]}];
ValNames[s_,n_]:=Transpose[{Table[ToExpression[s~StringJoin~ToString[i]~StringJoin~"x"],{i,1,n}],
Table[ToExpression[s~StringJoin~ToString[i]~StringJoin~"y"],{i,1,n}]}];
Const1[n_]:=Function[{conf,t},Total[#^2&/@((1-t)conf[[1]][[n]]+t conf[[2]][[n]]
-ValNames[Length[conf[[1]]]][[n]])]];
(**)
NormF[m_]:=Total[Map[#^2&,Flatten[m]]];


(* ::Section::Closed:: *)
(*QuadraticFormMatrix + LinearFormVector*)


(*poly:\:591a\:9805\:5f0f\:3001vl:\:5909\:6570\:5217*)
QuadraticFormVariableMatrix[vl_]:=Outer[#1 #2&,vl,vl,1]
Div2if[n_,l_]:=Table[If[n==i,l[[i]],l[[i]]/2],{i,1,Length[l]}];
Div2Matrix[m_]:=MapIndexed[Div2if[First[#2],#1]&,m];
(* QuadraticFormMatrix *)
QuadraticFormMatrix[poly_,vl_]:=
Div2Matrix[Partition[Map[Coefficient[poly,#]&,
Flatten[QuadraticFormVariableMatrix[vl],1]],Length[vl]]];
(* LinearFormVector *)
LinearFormVector[poly_,vl_]:=Module[{allzero},
allzero=Map[#->0&,vl];
Map[Coefficient[poly,#]/.allzero&,vl]]


(* ::Section::Closed:: *)
(*FindMatrix etc*)


(*\:4e09\:89d2\:5f62\:306e\:9802\:70b9\:306e\:5ea7\:6a19\:3092\:6c42\:3081\:308b\:95a2\:6570 VtoTriangle[V,tindex]*)
VtoTriangle[V_,tindex_]:= Part[V,tindex];
VtoTriangles[V_,tindex_]:=Map[VtoTriangle[V,#]&,tindex];
(*P\:306e\:91cd\:5fc3\:3092\:6c42\:3081\:308b\:95a2\:6570*)
Cog[P_]:= Mean[P];
(*\:5ea7\:6a19P\:3092l=(x,y)\:3060\:3051\:5e73\:884c\:79fb\:52d5\:3055\:305b\:308b\:95a2\:6570*)
Trans[P_,l_] :=  {#[[1]]-l[[1]],#[[2]]-l[[2]]}&/@P;

(*\:91cd\:5fc3\:304c\:539f\:70b9\:3067\:3042\:308b\:5834\:5408 2\:901a\:308a\:306e\:884c\:5217\:3092\:6c42\:3081\:308b\:95a2\:6570*)
FindMatrix[Tri1_,Tri2_]:=(Join[Transpose[Tri2],{{1,1,1}}].Inverse[Join[Transpose[Tri1],{{1,1,1}}]])[[1;;2,1;;2]];
FindMatrix1[Tri1_,Tri2_]:=Transpose[Take[Tri2,2]].Inverse[Transpose[Take[Tri1,2]]];
FindMatrices[V1_,V2_,tindex_]:=FindMatrix1[
Trans[VtoTriangle[V1,#],Cog[VtoTriangle[V1,#]]],
Trans[VtoTriangle[V2,#],Cog[VtoTriangle[V2,#]]]]&/@tindex

(*\:91cd\:5fc3\:304c\:539f\:70b9\:3067\:306a\:3044\:5834\:5408 3x3\:884c\:5217*)
FindAffineMatrix[Tri1_,Tri2_]:= MatrixForm[
Join[Transpose[Tri2],{{1,1,1}}].Inverse[Join[Transpose[Tri1],{{1,1,1}}]]];
FindAffineMatrices[V1_,V2_,tindex_]:=FindAffineMatrix[VtoTriangle[V1,#],VtoTriangle[V2,#]]&/@ tindex


(* ::Section::Closed:: *)
(*EmbedMatrix,EmbedVector*)


(**)
(*En[P1_,P2_,A_]:= NormF[A-FindMatrix1[Trans[P1,Cog[P1]],Trans[P2,Cog[P2]]]];
En2[P1_,P2_,A_]:=Module[{B},
B=FindMatrix[Trans[P1,Cog[P1]],Trans[P2,Cog[P2]]];
NormF[B]-
1/NormF[Transpose[A]]
(NormF[B.Transpose[A]]+2 Det[B.Transpose[A]])];*)

(*F1[P1_,P2_,A_]:= QuadraticFormMatrix[En[P1,P2,A],Transpose[P2][[1]]];
F2[P1_,P2_,A_]:= QuadraticFormMatrix[En2[P1,P2,A],Flatten[P2]];*)

(*QuadraticFormMatrix\:306e\:7d50\:679c\:3092\:7528\:3044\:3066\:76f4\:63a5\:4f5c\:6210\:3057\:305f\:3082\:306e*)
F1a[{{a1x_,a1y_},{b1x_,b1y_},{c1x_,c1y_}},{{m11_,m12_},{m21_,m22_}}]:=
Module[{q,p12,p13,p23},
q=(a1y b1x-a1x b1y-a1y c1x+b1y c1x+a1x c1y-b1x c1y)^2;
p11=b1x^2+b1y^2-2 b1x c1x+c1x^2-2 b1y c1y+c1y^2;
p12=-(a1x b1x+a1y b1y-a1x c1x-b1x c1x+c1x^2-a1y c1y-b1y c1y+c1y^2);
p13=-b1x^2+a1x (b1x-c1x)+b1x c1x+(a1y-b1y) (b1y-c1y);
p22=a1x^2+a1y^2-2 a1x c1x+c1x^2-2 a1y c1y+c1y^2;
p23=-(a1x^2+a1y^2+b1x c1x-a1x (b1x+c1x)+b1y c1y-a1y (b1y+c1y));
p33=a1x^2+a1y^2-2 a1x b1x+b1x^2-2 a1y b1y+b1y^2;
{{p11/q,p12/q,p13/q},{p12/q,p22/q,p23/q},{p13/q,p23/q,p33/q}}];
(**)
F2a[{{a1x_,a1y_},{b1x_,b1y_},{c1x_,c1y_}},{{m11_,m12_},{m21_,m22_}}]:=
Module[{fm,cm,m1222,m1121,X,Y0,Y1,Y2,Y3,Y4,Y5,Y6},
fm=(m11^2+m12^2+m21^2+m22^2);
cm=(m11 m12+m21 m22);
m1121=(m11^2+m21^2);
m1222=(m12^2+m22^2);
X=(a1y (b1x-c1x)+b1y c1x-b1x c1y+a1x (-b1y+c1y));
Y0=(m12 m21-m11 m22);
Y1=-b1x^2 m1121+a1x (b1y cm-c1y cm+(b1x-c1x) m1121)+b1x ((a1y-2 b1y+c1y) cm+c1x m1121)
+(a1y-b1y) (-c1x cm+(b1y-c1y) m1222);
Y2=(-a1y b1x cm-a1y c1x cm+b1y c1x cm+b1x c1y cm+a1x^2 m1121+b1x c1x m1121
-a1x ((-2 a1y+b1y+c1y) cm+(b1x+c1x) m1121)+(a1y-b1y) (a1y-c1y) m1222);
Y3=((a1x-b1x) (2 a1y cm-2 b1y cm+(a1x-b1x) m1121)+(a1y-b1y)^2 m1222);
Y4=(-a1y b1x cm-a1x b1y cm+a1x c1y cm+b1x c1y cm-c1x^2 m11^2-a1x b1x m1121
+a1x c1x m1121+b1x c1x m1121+a1y c1x m11 m12+b1y c1x m11 m12
-2 c1x c1y m11 m12-a1y b1y m12^2+a1y c1y m12^2+b1y c1y m12^2-c1y^2 m12^2-c1x^2 m21^2
+c1x (a1y+b1y-2 c1y) m21 m22+(a1y-c1y) (-b1y+c1y) m22^2);
Y5=((b1x-c1x) (2 b1y cm-2 c1y cm+(b1x-c1x) m1121)+(b1y-c1y)^2 m1222);
Y6=((a1x-c1x) (2 a1y cm-2 c1y cm+(a1x-c1x) m1121)+(a1y-c1y)^2 m1222);
{{Y5/(fm X^2),0,Y4/(fm X^2),-(Y0/(fm X)),Y1/(fm X^2),Y0/(fm X)},{0,Y5/(fm X^2),Y0/(fm X),Y4/(fm X^2),-(Y0/(fm X)),Y1/(fm X^2)},
{Y4/(fm X^2),Y0/(fm X),Y6/(fm X^2),0,-(Y2/(fm X^2)),-(Y0/(fm X))},{-(Y0/(fm X)),Y4/(fm X^2),0,Y6/(fm X^2),Y0/(fm X),-(Y2/(fm X^2))},
{Y1/(fm X^2),-(Y0/(fm X)),-(Y2/(fm X^2)),Y0/(fm X),Y3/(fm X^2),0},{Y0/(fm X),Y1/(fm X^2),-(Y0/(fm X)),-(Y2/(fm X^2)),0,Y3/(fm X^2)}}];

(*EmbedMatrix*)
EmbedMatrix[n_,i_,j_,M_]:=Table[Switch[k,
i,Switch[l,i,M[[1,1]],j,M[[1,2]],_,0],
j,Switch[l,i,M[[2,1]],j,M[[2,2]],_,0],_,0],
{k,1,n},{l,1,n}];
EmbedMatrix[n_,i_,j_,k_,M_]:=Table[Switch[m,
i,Switch[l,i,M[[1,1]],j,M[[1,2]],k,M[[1,3]],_,0],
j,Switch[l,i,M[[2,1]],j,M[[2,2]],k,M[[2,3]],_,0],
k,Switch[l,i,M[[3,1]],j,M[[3,2]],k,M[[3,3]],_,0],_,0],
{m,1,n},{l,1,n}];

(*6\[Times]6\[Rule] 8\[Times]8*)
EmbedMatrix2[n_,i_,j_,k_,M_]:=Table[Switch[m,
2*i-1,Switch[l,2*i-1,M[[1,1]],2*i,M[[1,2]],
2*j-1,M[[1,3]],2*j,M[[1,4]],
2*k-1,M[[1,5]],2*k,M[[1,6]],_,0],
2*i,Switch[l,2*i-1,M[[2,1]],
2*i,M[[2,2]],2*j-1,M[[2,3]],2*j,M[[2,4]],
2*k-1,M[[2,5]],2*k,M[[2,6]],_,0],
2*j-1,Switch[l,2*i-1,M[[3,1]],2*i,M[[3,2]],
2*j-1,M[[3,3]],2*j,M[[3,4]],
2*k-1,M[[3,5]],2*k,M[[3,6]],_,0],
2*j,Switch[l,2*i-1,M[[4,1]],2*i,M[[4,2]],
2*j-1,M[[4,3]],2*j,M[[4,4]],
2*k-1,M[[4,5]],2*k,M[[4,6]],_,0],
2*k-1,Switch[l,2*i-1,M[[5,1]],
2*i,M[[5,2]],2*j-1,M[[5,3]],2*j,M[[5,4]],
2*k-1,M[[5,5]],2*k,M[[5,6]],_,0],
2*k,Switch[l,2*i-1,M[[6,1]],2*i,M[[6,2]],
2*j-1,M[[6,3]],2*j,M[[6,4]],2*k-1,M[[6,5]],
2*k,M[[6,6]],_,0],_,0],{m,1,2n},{l,1,2n}];


(*LineatFormVector\:306e\:7d50\:679c\:3092\:7528\:3044\:3066\:76f4\:63a5\:4f5c\:6210\:3057\:305f\:3082\:306e*)
F1v[{{a1x_,a1y_},{b1x_,b1y_},{c1x_,c1y_}},
{{m11_,m12_},{m21_,m22_}}]:=
Module[{q,p1,p2,p3,p4,p5,p6},
q=a1y b1x-a1x b1y-a1y c1x+b1y c1x+a1x c1y-b1x c1y;
p1=2 b1y m11-2 c1y m11-2 b1x m12+2 c1x m12;
p2=2 b1y m21-2 c1y m21-2 b1x m22+2 c1x m22;
p3=2 a1y m11-2 c1y m11-2 a1x m12+2 c1x m12;
p4=2 a1y m21-2 c1y m21-2 a1x m22+2 c1x m22;
p5=-2 a1y m11+2 b1y m11+2 a1x m12-2 b1x m12;
p6=-2 a1y m21+2 b1y m21+2 a1x m22-2 b1x m22;
{p1/q,p2/q,p3/-q,p4/-q,p5/-q,p6/-q}];
(*F2v[P1_,P2_,A_]:=LinearFormVector[En2[P1,P2,A],Flatten[P2]];*)


(*EmbedVector*)
EmbedVector[n_,i_,j_,k_,V_]:=
Table[Switch[m,2*i-1,V[[1]],2*i,V[[2]],
2*j-1,V[[3]],2*j,V[[4]],
2*k-1,V[[5]],2*k,V[[6]],_,0],
{m,1,2n}];


(*QuadraticMatrix\:3092\:57cb\:3081\:8fbc\:3093\:3060\:884c\:5217\:306e\:548c Fn*)
(*Fn[P_,V_,T_,A_]:=Total[Table[
EmbedMatrix[Length[P],T[[i,1]],T[[i,2]],T[[i,3]],
F1[VtoTriangle[P,T[[i]]],VtoTriangle[V,T[[i]]],A]],
{i,1,Length[T]}]];
(*F1a\:3092\:7528\:3044\:305f\:3082\:306e Fn2*)
Fn2[P_,V_,T_,A_]:=Total[Table[
EmbedMatrix[Length[P],T[[i,1]],T[[i,2]],T[[i,3]],
F1a[VtoTriangle[P,T[[i]]],A]],
{i,1,Length[T]}]];*)

(*LinearFormVector\:3092\:57cb\:3081\:8fbc\:3093\:3060\:884c\:5217\:306e\:548c Ln*)
(*Ln[P_,V_,T_,A_]:=Total[Table[
EmbedVector[Length[P],T[[i,1]],T[[i,2]],T[[i,3]],
F1v[VtoTriangle[P,T[[i]]],A]]
,{i,1,Length[T]}]];*)


(* ::Section::Closed:: *)
(*Local Interpolation*)


RotateAngle[m_,flag_]:=Module[{\[Theta]},
\[Theta]=ArcTan[m[[1,1]],m[[2,1]]];
Which[flag == 0, \[Theta],
flag > 0, If[\[Theta]>=0, \[Theta], \[Theta]+2\[Pi]],
flag < 0, If[\[Theta]<=0, \[Theta], \[Theta]-2\[Pi]]]];


NewFindMatrices[conf_]:=FindMatrices[conf[[1]],conf[[2]],conf[[3]]];
(*m:\:884c\:5217*)
(**)
LocalLinear[m_]:=Function[{t},(1-t) IdentityMatrix[2] + t m];
LocalAlexa[m_]:=Function[{t},
Module[{p,q,\[Theta]},
{p,q}=PolarDecompositionPlus[m];
\[Theta]=RotateAngle[p,0];
RotationMatrix[t \[Theta]].((1-t) IdentityMatrix[2] + t q)]];
LocalLogExp[m_]:=Function[{t},
Module[{p,q,s,\[Theta],evec,eval},
{p,q}=PolarDecompositionPlus[m];
{eval,evec}=Eigensystem[q];
\[Theta]=RotateAngle[p,0];
RotationMatrix[t \[Theta]].Transpose[evec].DiagonalMatrix[{eval[[1]]^t,eval[[2]]^t}].evec
]];
(*Flag\:306f\:4fdd\:7559?*)

LocalInterpolations[local_,conf_]:=Function[{t},
(local[#][t]&)/@ NewFindMatrices[conf]];


(* ::Section::Closed:: *)
(*Constraint Functions*)


Const1[n_]:=Function[{conf,t},Total[#^2&/@((1-t)conf[[1]][[n]]+t conf[[2]][[n]]
-ValNames[Length[conf[[1]]]][[n]])]];
ConstM[conf_,t_]:= Total[#^2&/@((1-t)Mean[conf[[1]]]+t Mean[conf[[2]]]
-Mean[ValNames[Length[conf[[1]]]]])];
Constfix[k_,l_]:=Function[{conf,t},Total[#^2&/@(ValNames[Length[conf[[1]]]][[k]]
-ValNames[Length[conf[[1]]]][[l]]-(conf[[1]][[k]] -conf[[1]][[l]]))]];
(*Constfix2[k_,l_]:=Function[{conf,t},Total[#^2&/@(ValNames[Length[conf[[1]]]][[k]]
-ValNames[Length[conf[[1]]]][[l]]-Ekl[k,l][conf,t])]];*)
Constfix2[k_,l_]:=Function[{conf,t},Module[{st,n},
st = conf[[1]];
n = Length[st];
Total[#^2&/@(ValNames[n][[k]]
-ValNames[n][[l]]-RotationMatrix[2\[Pi] t].(st[[k]]-st[[l]]))]]];


(*Const\:306evn^2\:306e\:4fc2\:6570,m:\:9802\:70b9\:9078\:629e*)
ConstMatrix[m_,st_]:=Table[Switch[i,
2m-1,Switch[j,2m-1,1,_,0],
2m,Switch[j,2m,1,_,0],_,0],
{i,1,2Length[st]},{j,1,2Length[st]}]
(*Const\:306evn\:306e\:4fc2\:6570,m:\:9802\:70b9\:9078\:629e*)
(*ConstVector[m_,st_,en_,t_]:=Table[Switch[i,
2m-1,-2((1-t)st[[m,1]]+t en[[m,1]]),
2m,-2((1-t)st[[m,2]]+t en[[m,2]]),_,0],{i,1,2Length[st]}];*)

(*\:91cd\:5fc3\:3092\:9078\:629e*)
ConstMatrixM[conf_,t_]:=DoubleMatrix[Table[1,{i,1, Length[conf[[1]]]},{j,1 , Length[conf[[1]]]}]/(Length[conf[[1]]]^2)];
ConstVectorM[conf_,t_]:=Module[{st,en,p,x,y},
st = conf[[1]];
en = conf[[2]];
p=(1/Length[st]) Map[-2((1-t)#[[1]]+t #[[2]]){Table[1,{i,1,Length[st]}]}&,Transpose[{Mean[st],Mean[en]}]];
x=Flatten[Transpose[{Flatten[p[[1]]],Table[0,{Length[st]}]}]];
y=Flatten[Transpose[{Table[0,{Length[st]}],Flatten[p[[2]]]}]];
x+y];
(*ConstMatrixM[st_,en_,tri_,t_]:=DoubleMatrix[Table[1,{i,1, Length[st]},{j,1 , Length[st]}]/(Length[st]^2)];
ConstVectorM[st_,en_,tri_,t_]:=Module[{p,x,y},
p=(1/Length[st]^2) Map[-2((1-t)#[[1]]+t #[[2]]){Table[1,{i,1,Length[st]}]}&,Transpose[{Mean[st],Mean[en]}]];
x=Flatten[Transpose[{Flatten[p[[1]]],Table[0,{Length[st]}]}]];
y=Flatten[Transpose[{Table[0,{Length[st]}],Flatten[p[[2]]]}]];
x+y];*)

(*Constfix*)
ConstfixMatrix[k_,l_]:= Function[{conf,t},DoubleMatrix[EmbedMatrix[Length[conf[[1]]],k,l,{{1,-1},{-1,1}}]]];
ConstfixVector[k_,l_]:=Function[{conf,t},Module[{st},
st = conf[[1]];
Table[Switch[i,2k-1,-2(st[[k,1]]-st[[l,1]]),2k,-2(st[[k,2]]-st[[l,2]]),
2l-1,2(st[[k,1]]-st[[l,1]]),2l,2(st[[k,2]]-st[[l,2]]),_,0],{i,1,2Length[st]}]]];
(*Constfix2:\:56de\:8ee2\:884c\:5217*)
(*Constfix2Matrix = ConstfixMatrix*)
Constfix2Vector[k_,l_]:=Function[{conf,t},Module[{st},
st = conf[[1]];
EmbedVector2[Length[st],k,l,
{-2 st[[k,1]]Cos[2 \[Pi] t]+2 st[[l,1]]Cos[2 \[Pi] t]+2 st[[k,2]]Sin[2 \[Pi] t]-2 st[[l,2]]Sin[2 \[Pi] t],
-2 st[[k,2]]Cos[2 \[Pi] t]+2 st[[l,2]]Cos[2 \[Pi] t]-2 st[[k,1]]Sin[2 \[Pi] t]+2 st[[l,1]]Sin[2 \[Pi] t],
2 st[[k,1]]Cos[2 \[Pi] t]-2 st[[l,1]]Cos[2 \[Pi] t]-2 st[[k,2]]Sin[2 \[Pi] t]+2 st[[l,2]]Sin[2 \[Pi] t],
2st[[k,2]]Cos[2 \[Pi] t]-2 st[[l,2]]Cos[2 \[Pi] t]+2 st[[k,1]]Sin[2 \[Pi] t]-2 st[[l,1]]Sin[2 \[Pi] t]}]]];

(*Ekl[k_,l_]:=Function[{conf,t},RotationMatrix[2\[Pi] t].(conf[[1,k]]-conf[[1,l]])];
ConstMatrix2[n_,k_,l_]:=Table[Switch[i,
2k-1,Switch[j,2l-1,1,_,0],
2k,Switch[j,2l,1,_,0],_,0],{i,1,2n},{j,1,2n}];
ConstfixMatrix[n_,k_,l_,st_]:=Module[{C1,C2},
C1 = Table[Switch[i,2k-1,Switch[j,2l-1,1,_,0],
2k,Switch[j,2l,1,_,0],_,0],{i,1,2Length[st]},{j,1,2Length[st]}];
C2 = Table[Switch[i,2l-1,Switch[j,2k-1,1,_,0],
2l,Switch[j,2k,1,_,0],_,0],{i,1,2Length[st]},{j,1,2Length[st]}];
n ConstMatrix[k,st]+n ConstMatrix[l,st]
-n C1-n C2];
ConstfixVector[n_,k_,l_,st_]:=Table[Switch[i,
2k-1,-2(st[[k,1]]-st[[l,1]]),2k,-2(st[[k,2]]-st[[l,2]]),
2l-1,2(st[[k,1]]-st[[l,1]]),2l,2(st[[k,2]]-st[[l,2]]),_,0],{i,1,2Length[st]}];
Constfix2Vector[n_,k_,l_,st_,t_]:=Table[Switch[i,
2k-1,-2 Ekl[st,k,l,t][[1]],2k,-2Ekl[st,k,l,t][[2]],
2 l-1,2Ekl[st,k,l,t][[1]],2 l ,2Ekl[st,k,l,t][[2]],_,0],{i,1,2Length[st]}];*)


(*ARAP\:306eConst*)
ConstPair[m_]:=Function[{conf,t},{
Table[Switch[i,
2m-1,Switch[j,2m-1,1,_,0],
2m,Switch[j,2m,1,_,0],_,0],{i,1,2Length[conf[[1]]]},{j,1,2Length[conf[[1]]]}],
Table[Switch[i,
2m-1,-2((1-t)conf[[1,m,1]]+t conf[[2,m,1]]),
2m,-2((1-t)conf[[1,m,2]]+t conf[[2,m,2]]),_,0],{i,1,2Length[conf[[1]]]}]
}]
ConstPair[m_,n_]:=Function[{conf,t},ConstPair[m][conf,t]+ConstPair[n][conf,t]];
ConstPairM:=Function[{conf,t},{ConstMatrixM[conf,t],ConstVectorM[conf,t]}];
ConstPairfix[k_,l_]:=Function[{conf,t},
{ConstfixMatrix[k,l][conf,t],ConstfixVector[k,l][conf,t]}];
ConstPairfix2[k_,l_]:=Function[{conf,t},
{ConstfixMatrix[k,l][conf,t],Constfix2Vector[k,l][conf,t]}];
(*ConstM,Constfix\:306f\:6e96\:5099\:4e2d*)


(**)
DoubleMatrix[m_]:=Module[{X0,Y0},
X0[l_]:=Flatten[Transpose[{l,Table[0,{Length[l]}]}]];
Y0[l_]:=Flatten[Transpose[{Table[0,{Length[l]}],l}]];
Flatten[Map[{X0[#],Y0[#]}&,m],1]];


(* ::Section::Closed:: *)
(*Global Interpolation*)
(*(Embed\:3092\:7528\:3044\:305f\:3082\:306e\:3001Alexa)*)


(*energy*)
(*Alexa\:306eMatrix,Vector*)
QuadraticFormAlexa[local_,conf_]:=Function[{t},{
DoubleMatrix[
Total[Map[EmbedMatrix[Length[conf[[1]]],#[[2,1]],#[[2,2]],#[[2,3]],#[[1]]]&,
Transpose[
{Map[F1a[#[[1]],#[[2]]]&,
Transpose[{VtoTriangles[conf[[1]],conf[[3]]],LocalInterpolations[local,conf][t]}]],
conf[[3]]}]]]],
Flatten[
Total[Map[EmbedVector[Length[conf[[1]]],#[[2,1]],#[[2,2]],#[[2,3]],#[[1]]]&,
Transpose[
{Map[F1v[#[[1]],#[[2]]]&,
Transpose[{VtoTriangles[conf[[1]],conf[[3]]],LocalInterpolations[local,conf][t]}]],
conf[[3]]}]]]]
}];
(*Sim\:306eMatrix,Vector*)
QuadraticFormSim[local_,conf_]:=Function[{t},{
Total[Map[
EmbedMatrix2[Length[conf[[1]]],#[[2,1]],#[[2,2]],#[[2,3]],#[[1]]]&,
Transpose[{Map[F2a[#[[1]],#[[2]]]&,
Transpose[{VtoTriangles[conf[[1]],conf[[3]]],LocalInterpolations[local,conf][t]}]],
conf[[3]]}]]],
Table[0,{i,1,2Length[conf[[1]]]}]}];
(**)
ARAP[local_,energy_,const_,conf_]:=Function[{t},
Module[{G,h},
G=energy[local,conf][t][[1]]
+const[conf,t][[1]];
h=energy[local,conf][t][[2]]
+const[conf,t][[2]];
{G,h}
]];


(* ::Section::Closed:: *)
(*Animation*)


(* *)
ShowStatus[st_,en_,tri_,plotrange_]:=
Show[Graphics[{MakePolygon[st,tri],MakePolygon[en,tri]}],
PlotRange->plotrange,Frame->True]
(* *)
DrawAnimation[local_,energy_,const_,conf_]:=
DynamicModule[{GU},
Animate[
GU=ARAP[local,energy,const,conf][x];
Show[Graphics[MakePolygon[
Partition[-0.5 Inverse[GU[[1]]].GU[[2]],2],conf[[3]]],
PlotRange -> AnimationRange[conf]]],{x,0,1}, 
AnimationRunning->False,AnimationRepetitions->1]];


(* Export["filename.mov",ListAnimation[k,local,energy,const,conf]] *)
(*ListAnimation[k_,local_,energy_,const_,conf_]:=
Module[{GU},
Table[
GU=ARAP[local,energy,const,conf][i/(k-1)];
Show[Graphics[MakePolygon[
Partition[-0.5 Inverse[GU[[1]]].GU[[2]],2],conf[[3]]],
PlotRange -> AnimationRange[conf]
]],{i,0,k-1}]];*)
ListAnimation[k_,local_,energy_,const_,conf_]:=
Module[{GU,V,range,table},
range = {conf[[1,1]],conf[[1,1]]};
table=Table[
GU=ARAP[local,energy,const,conf][i/(k-1)];
V=Partition[-0.5 Inverse[GU[[1]]].GU[[2]],2];
range=AnimationRange2[range,AnimationRange[{V,{},{}}]];
MakePolygon[V,conf[[3]]],{i,0,k-1}];
Map[Graphics[#,PlotRange->range]&,table]];
