clear all; close all; clc;
cd('C:\Users\aliha\Desktop\analysis notebooks\trajectory in mosaic aggregates\');
load("C:\Users\aliha\Desktop\analysis notebooks\trajectory in mosaic aggregates\motion.mat");

dT = 6; blur = 0; resbin = 2*0.633;

st1 = struct2cell(DSigma2_Est(Expression1*resbin,blur,dT,'DSigma2_MLE',[]));
st1d = st1{1};
st1v = norm(st1{3});

st2 = struct2cell(DSigma2_Est(Expression2*resbin,blur,dT,'DSigma2_MLE',[]));
st2d = st2{1};
st2v = norm(st2{3});

st3 = struct2cell(DSigma2_Est(Expression3*resbin,blur,dT,'DSigma2_MLE',[]));
st3d = st3{1};
st3v = norm(st3{3});

st4 = struct2cell(DSigma2_Est(Expression4*resbin,blur,dT,'DSigma2_MLE',[]));
st4d = st4{1};
st4v = norm(st4{3});

st5 = struct2cell(DSigma2_Est(Expression5*resbin,blur,dT,'DSigma2_MLE',[]));
st5d = st5{1};
st5v = norm(st5{3});

st6 = struct2cell(DSigma2_Est(Expression6*resbin,blur,dT,'DSigma2_MLE',[]));
st6d = st6{1};
st6v = norm(st6{3});

st7 = struct2cell(DSigma2_Est(Expression7*resbin,blur,dT,'DSigma2_MLE',[]));
st7d = st7{1};
st7v = norm(st7{3});

avg_V = mean([st1v,st2v,st3v,st4v,st5v,st6v,st7v])
avg_D = mean([st1d,st2d,st3d,st4d,st5d,st6d,st7d])

