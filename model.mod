set SATELITY;
param r := 6378137;
param x{SATELITY} >0;
param y{SATELITY} >0;
param z{SATELITY} >0;
param w := 0.00001;

param vst{SATELITY} >0;

var X >= -r, <=r;
var Y >= -r, <=r;
var Z >= -r, <=r;

minimize f:
    sum{i in SATELITY} (w*(X-x[i])^2+w*(Y-y[i])^2+w*(Z-z[i])^2-w*(vst[i])^2)^2
;
