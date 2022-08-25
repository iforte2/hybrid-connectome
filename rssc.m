%rs_SC = optimal hybrid connectome
%comp = all grid values
%SC_opt =  similarity with structure for optimal hybrid connectome
%FC_opt = maximum correlation between FC reconstructed and observed FC
%Optimal is the index list where the highest comp occurs
%beta_sim is the simulated temp/beta values
%optimum is the function that computes the optimzation values

%for SC with fiber counts, scale by the sum total of all fibers
%(sum(sum(SC));

function [rs_SC,comp,SC_opt,FC_opt,beta_sim,optimal] = rssc(SC,series,lambd_beta)
tic
[comp,SC_opt,FC_opt,beta_sim,optimal] = optimum(SC,series,lambd_beta);
toc

%create optimal rs-SC
for j = 1:size(SC,3)
        tic
k = optimal(j);% k gets the index
lam = lambd_beta(k,1);
bet = lambd_beta(k,2);

rs_SC(:,:,j) = FSE(SC(:,:,j),series(:,:,j),lam,bet);
toc

end
