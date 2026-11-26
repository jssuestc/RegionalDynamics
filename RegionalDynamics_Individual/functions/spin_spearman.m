function [rho,p,nullR] = spin_spearman(RAD,iRTD,SpinIndices)
% One-sided map rotation (iRTD only), two-sided statistical test.
% One-to-one permutations preserve ranks; vectorized rank correlation.
a=tiedrank(RAD(:)); b=tiedrank(iRTD(:));
a=a-mean(a); b=b-mean(b); den=norm(a)*norm(b);
assert(isfinite(den) && den>0,'Spatial correlation undefined for a constant map.');
rho=max(-1,min(1,a'*b/den));
nullR=max(-1,min(1,(a'*b(SpinIndices))/den));
p=(1+sum(abs(nullR)>=abs(rho)-1e-12))/(numel(nullR)+1);
end
