function H = haar2d(X,nlevel)

% Author: Jun Li, more info@ http://goldensectiontransform.org/
%
% function H = haar2d(X,nlevel)
% 3-level 2d haar wavelet transform of 8*8 image matrix

% nlevel = 3; % 3-level transform for each 8*8 image block

[xx,yy] = size(X);

H=X;
for i=1:nlevel

   for j=1:xx
      [ss,dd] = haar1d(H(j,1:yy)); % row transform
      H(j,1:yy) = [ss,dd];
   end
   
   for k=1:yy
      [ss,dd] = haar1d(H(1:xx,k)'); % column transform
      H(1:xx,k) = [ss,dd]';
   end
   
   xx = xx/2;
   yy = yy/2;

end


%% 1d haar lifting scheme

function [ss,dd] = haar1d(S)

N = length(S);

s0 = S(1:2:N-1);  % S(1),S(3),S(5),S(7)...
d0 = S(2:2:N);    % S(2),S(4),S(6),S(8)...

d1 = s0 - d0;
s1 = d0 + 1/2*d1;
ss = sqrt(2)*s1;
dd = d1/sqrt(2);
