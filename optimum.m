
%SC and FC are an n x n x subject count
%series is m x n x subject count, where m is the number of time points.
%Each subjects should have the same amount
%236x80 time series takes about 215 seconds to run with 144 combinations of
%parameter values

function [comp,SC_opt,FC_opt,beta_sim,optimal] = optimum(SC,series,lambd_beta,ts_sz)
%ts_sz = [number of time_points, number of ROIs) per subject; 
SC_comp = [];%similarity with SC
FC_comp = [];%FC reconstruction quality
maxi = [];%max temperature

sims = 2000;

%run every subject
      for z = 1:size(SC,3)

          %initiate the observed data (SC, FC, timeseries)
        SC1 = SC(:,:,z);
        SC1 = SC1 - diag(diag(SC1));%remove diagonal
        %scale the SC if using deterministic to the sum total of all fibers. Not necessary for probabilistic
        rowcol = ts_sz(:,:,z);
        row = rowcol(1);
        col = rowcol(2);
        ts = series(1:row,1:col,z);
        time_series1 = ts;
        FC = corrcoef(time_series1);
        FC1 = FC-diag(diag(FC));
        %for each subject, build grid search
for i = 1:length(lambd_beta)
    l = lambd_beta(i,1);
    k = lambd_beta(i,2);

    rssc = FSE(SC1,time_series1,l,k);
    SC_comp(i,1) = corr2(SC1, (sign(SC1).*(abs(rssc))));%similarity computation for Structure and magnitude of rs-sc
    
    [max_index,FC_recon] = ising(rssc,FC1,sims);%fc_recon is the reconstructed functional connectome, max index is the beta value where the reconstruction highest

    maxi(i,1) = max_index*0.2; %multiplied by 0.2 beause the beta goes from 0.2 to 3.0 with 0.2 increments and max index is a scalar, so we convert it.
    FC_comp(i,1) = corr2(FC1,FC_recon);%reconstruction quality

end
comp(:,z) = SC_comp+FC_comp;%combining max(fc) + Sm, giving equal weight to both

beta_sim(:,z) = maxi;%sanity check the simulated optimal temperature, to make sure that simulated temp is not hitting 3.0.
% if simualted temp is hitting 3.0, then most likely the structural
% connectivity needs to be scaled/normalized by the sum total.

optimal(z,1) = find(comp(:,z)==max(comp(:,z)));%find index of optimal params

SC_opt(z) = SC_comp(optimal(z,1));%correlation with SC at optimal params

FC_opt(z) = FC_comp(optimal(z,1));%reconstruction quality with optimal params

      end
     
end