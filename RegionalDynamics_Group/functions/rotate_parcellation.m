function perm_id = rotate_parcellation(coord_l, coord_r, nrot)
%
% Usage:
%   perm_id = rotate_parcellation(coord_l, coord_r, nrot)
%
% Description:
%   Rotate parcellation (authors: @frantisekvasa, @saratheriver)
%
% Inputs:
%   coord_l (double array) - Coordinates of left hemisphere regions on the sphere, size = [m x 3]
%   coord_r (double array) - Coordinates of right hemisphere regions on the sphere, size = [m x 3]
%   nrot (int) - Number of rotations.
%
% Outputs:
%   perm_id (double array) - Array of permutations, size = [m x nrot]
%
% References:
%   Alexander-Bloch A, Shou H, Liu S, Satterthwaite TD, Glahn DC, Shinohara RT, 
%       Vandekar SN and Raznahan A (2018). On testing for spatial correspondence 
%       between maps of human brain structure and function. NeuroImage, 178:540-51.
%   Vása F, Seidlitz J, Romero-Garcia R, Whitaker KJ, Rosenthal G, Vértes PE, 
%       Shinn M, Alexander-Bloch A, Fonagy P, Dolan RJ, Goodyer IM, the NSPN 
%       consortium, Sporns O, Bullmore ET (2017). Adolescent tuning of association 
%       cortex in human structural brain networks. Cerebral Cortex, 28(1):281?294.
%
% Sara Lariviere  |  saratheriver@gmail.com

% Adaptation: vectorized distances and deterministic first-index tie breaking.
% The original Haar rotations, mirrored right hemisphere, greedy maximum
% of nearest distances, and exclusion of identity rotations are retained.
if size(coord_l,2)~=3 && size(coord_l,1)==3, coord_l=coord_l'; end
if size(coord_r,2)~=3 && size(coord_r,1)==3, coord_r=coord_r'; end
assert(size(coord_l,2)==3 && size(coord_r,2)==3,'Expected N x 3 sphere centroids.');
assert(all(isfinite(coord_l(:))) && all(isfinite(coord_r(:))),'Invalid centroids.');
nL=size(coord_l,1); nR=size(coord_r,1); n=nL+nR;
perm_id=zeros(n,nrot); r=0; I1=diag([-1 1 1]);
while r<nrot
    [TL,R]=qr(randn(3)); TL=TL*diag(sign(diag(R)));
    if det(TL)<0, TL(:,1)=-TL(:,1); end
    TR=I1*TL*I1;
    pL=greedy_match(coord_l,coord_l*TL);
    pR=greedy_match(coord_r,coord_r*TR);
    p=[pL;nL+pR];
    if all(p==(1:n)'), continue; end
    r=r+1; perm_id(:,r)=p;

end
end

function p=greedy_match(A,B)
% Squared Euclidean distances preserve the ordering of Euclidean distances.
D=max(0,sum(A.^2,2)+sum(B.^2,2)'-2*(A*B'));
n=size(A,1); p=zeros(n,1); active=true(n,1);
for k=1:n
    [nearest,j]=min(D,[],2);
    nearest(~active)=-Inf;
    [~,i]=max(nearest); col=j(i);
    p(i)=col; active(i)=false; D(:,col)=Inf;
end
end
