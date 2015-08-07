function X = ihaar2d(H,nlevel)

% Author: Jun Li, more info@ http://goldensectiontransform.org/
%
% function X = ihaar2d(H,nlevel)
% 3-level inverse 2d haar wavelet transform of 8*8 image matrix

% nlevel = 3; % 3-level inverse transform for each 8*8 image block

[x,y] = size(H);
xx = x/(2^nlevel);
yy = y/(2^nlevel);

X=H;
for i=1:nlevel
   
   for j=1:yy*2
      ss = X(1:xx,j);
      dd = X(xx+1:xx*2,j);
      X(1:xx*2,j) = ihaar1d(ss',dd')'; % inverse column transform
   end
   
   for k=1:xx*2 
      ss = X(k,1:yy);
      dd = X(k,yy+1:yy*2);
      X(k,1:yy*2) = ihaar1d(ss,dd); % inverse row transform
   end
   
   xx = xx*2;
   yy = yy*2;
   
end


%% inverse 1d haar lifting scheme

function S = ihaar1d(ss,dd)

N = length(ss) + length(dd);

d1 = sqrt(2)*dd;
s1 = ss/sqrt(2);
d0 = s1 - 1/2*d1;
s0 = d1 + d0;

S(2:2:N) = d0;
S(1:2:N-1) = s0;
