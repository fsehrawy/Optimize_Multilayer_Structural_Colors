function [R,T]=TMM_Calculate_R_T(layers,t,lambda,stepsize,Rgoal)

% Load in index of refraction for each material
n = zeros(size(layers,2),size(lambda,2));
for index = 1:size(layers,2)
    n(index,:) = Materials(layers{index},lambda);
end

% Calculate transfer matrices, and field at each wavelength and position
t(1)=0;
t_cumsum=cumsum(t);
x_pos =(stepsize/2):stepsize:sum(t); %positions to evaluate field
%x_mat specifies what layer number the corresponding point in x_pos is in:
x_mat = sum(repmat(x_pos,length(t),1)>repmat(t_cumsum',1,length(x_pos)),1)+1; 
R=lambda*0;
Enorm=zeros(length(x_pos),length(lambda));

for l = 1:length(lambda)
    % Calculate transfer matrices for incoherent reflection and transmission at the first interface
    S=I_mat(n(1,l),n(2,l));
    for matindex=2:(length(t)-1)
        S=S*L_mat(n(matindex,l),t(matindex),lambda(l))*I_mat(n(matindex,l),n(matindex+1,l));
    end
    R(l)=abs(S(2,1)/S(1,1))^2; %JAP Vol 86 p.487 Eq 9 Power Reflection from layers other than substrate
    T(l)=abs(1/S(1,1))^2; %JAP Vol 86 p.487 Eq 10 Power Transmission from layers other than substrate
end

