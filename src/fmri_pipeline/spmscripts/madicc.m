function rmad = madicc(x,y)
% Median Absolute Deviation Intraclass Correlation Coefficient
%
% Impliments Median Absolute Deviation Correlation Coefficient, (as
% described in Shevlyakov & Smirnov (2011)), modified to be the
% intraclass version of correlation.  The (non-intrasclass) 
% estimate is
%     r = ( mad(Sp)^2 - mad(Sm)^2 ) / ( mad(Sp)^2 + mad(Sm)^2 )
% where
%     Sp = (x-m(x))/mad(x) + (y-m(y))/mad(y);
%     Sm = (x-m(x))/mad(x) - (y-m(y))/mad(y);
% and m() is median and mad() is the median absolute deviation,
%     mad(x) = m(abs(x-m(x)))
%
% For intraclass correlation we assume mad(x)=mad(y) and so the divisors
% cancel; further, we find a common estimate of the median mm=m([x,y])
% can compute Sp & Sm as:
%     Sp = (x-mm) + (y-mm);
%     Sm = (x-mm) - (y-mm);
%
%
% REFERENCES
%
% Shevlyakov, G., & Smirnov, P. (2011). Robust estimation of a
% correlation coefficient: an attempt of survey. Australian & New
% Zealand Journal of Statistics, 40(1), 147156.
% 
% Kharin, Y. S., & Voloshko, V. A. (2011). Robust estimation of AR 
% coefficients under simultaneously influencing outliers and missing 
% values. Journal of Statistical Planning and Inference, 141(9),
% 3276288.
%
% 2014-07-08
% Thomas Nichols http://warwick.ac.uk/tenichols

I=find(all(~isnan([x(:) y(:)]),2));
if isempty(I)
  rmad=NaN;
else
  mx    = median(x(I));
  my    = median(y(I));
  Sp    = (x(I)-mx) + (y(I)-my);
  Sm    = (x(I)-mx) - (y(I)-my);
  madSp = median(abs(Sp-median(Sp)));
  madSm = median(abs(Sm-median(Sm)));
  if madSp==0 && madSm==0
    rmad = NaN;
  else
    rmad = (madSp^2 - madSm^2)/(madSp^2 + madSm^2);
  end
end
