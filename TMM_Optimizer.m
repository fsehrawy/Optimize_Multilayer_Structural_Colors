function Rsum=TMM_Optimizer(layers,L,lambda,stepsize,Rgoal)

t = [100 L 100];  % thickness of each corresponding layer in nm (thicknesses of the first and last layers are irrelevant)

Nl=length(lambda); %initialise
R=zeros(Nl,1);
    
[R,T]=TMM_Calculate_R_T(layers,t,lambda,stepsize,Rgoal);

plot(lambda,R,lambda,Rgoal,'linewidth',1.5)
hold on
plot(lambda,T,'--','linewidth',1.5)
hold off
axis([lambda(1) lambda(end) 0 1])
title(['L: ' num2str(round(L))])
xlabel('Wavelength, nm')
ylabel('Reflectance, Transmittance')
pause(0.001)

spectra_corr=corrcoef(R,Rgoal);
[M_R,I_R] = max(R(:));
[M_Rgoal,I_Rgoal] = max(Rgoal(:));

Rsum=mean(sum((R-Rgoal').^2))*(1-spectra_corr(1,2));

end