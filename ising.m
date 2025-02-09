
%{ 
****NOTES****
[max_index,T,FC_recon,inter_corr] = ising_model(J,FC);

1) 
inter corr is the cost function comparing observed FC and reconstructed FC
%for each temperature

2) 
max_index is the x-axis location where the cost function is at its peak

3)
FC_recon is the reconstructed functional connectome that gives the peak of
inter_corr
%}


function [max_index,FC_recon,inter_corr,S_recon] = ising(J,FC,sims)
%tic
%% Set up the parameters and run the simulaiton 
% Set matrix of temperatures -- might have to break this up as you approach array size limits
FC(isnan(FC))=0;
FC(FC == diag(FC)) =0; % make sure diag is zero
N=size(FC,1); % number of regions

normJ = J;
normJ(isnan(normJ))=0;
temp = 0.4:0.2:1.6;
%temp = 1:1:10;
%the key is to capture where the maximum happens without increasing run time too much

time_points = sims;
%time_points = 2000; % number of samples for MCMC. Processing time increases the larger this is
%tested 200->300,000. After 50000, the only thing that increases is
%processing time, not the res1lt

rand_coord = randperm(N); % Random assignment of regions 1 through N
count = 0;
temp_length = length(temp);
S = zeros (N,time_points,temp_length); % complete spin matrix [3D array]
%%% 
   
for t = temp
    %% +1,-1
    spin_vec = (((rand (N,1) > .5)*2 -1)); % initialize spin vector --> an Nx1 array whose elements are +1 or -1                                       
    %spin_vec = randi(-5,N,5);
    count = count + 1;
    index = 0;
    
    for j = 1: time_points
        % the interior loop serves to execute the recursive formula : 
        % index = index + 1 a total of N times; 
        % --> at each step, index takes on a new value, and so, too, does
        % 'flip' ==> testing a new Hamiltonian at each step 
        %        ==> potentially getting a new spin vector 
        for i = 1:N
            index = index + 1;
            if index > N
                index = index - N;
                rand_coord = randperm(N);   % randperm() generates a matrix consisting of a random permutation of N elements (1:N)
         
            end
            flip = rand_coord(index);       % use 'flip' to pick out a random spin (among the N)
                % rand_coord defined OUTSIDE the loop --> ensures that the
                % same value for 'flip' will NOT be chosen more than once
         
                
            dE = 0;    % initialize energy as 0
                % hamiltonian --> w/ NO external (B) field term
               for k = 1:N

                if (k~=flip)
                    %dE = dE + (normJ(flip,k)*spin_vec(k));
                     dE = dE + (normJ(flip,k)*spin_vec(k));
                end
               end

           dE = 2*(dE*spin_vec(flip)); % completion of hamiltonian --> this is the TOTAL energy cost of the spin configuration
    
            if dE <= 0 %energy has to be non-negative
                spin_vec(flip) = -spin_vec(flip);
            elseif (rand <= exp (-dE/t)) % the probability of the configuration has to be non-random
                spin_vec(flip) = -spin_vec(flip);

            end
        end

        S(:,j,count) = spin_vec;    % store the (new) spin vector as an element of the 3D array, S, after each time point, j
    end

%excluded 05/15
    spin_prod_allT(:,:,count) = (S(:,:,count)*S(:,:,count)')./time_points; % compute the paired spin interaction
    %, which is just another version of  correlation instead of pearson
    corr_recon = spin_prod_allT(:,:,count);
    corr_recon(corr_recon == diag(corr_recon)) =0; % zero out the diagonal
    
    
    inter_corr(count) = corr2(corr_recon,FC); %compute correlation between the two correlations PSI and Pearson of observed

end %%% End Temperature loop
max_index = find(inter_corr == max(inter_corr(:)));
FC_recon = spin_prod_allT(:,:,max_index);
S_recon = S(:,:,max_index);
%FC_recon = FC_recon-diag(diag(FC_recon));
FC_recon(FC_recon == diag(FC_recon)) =0;

%toc
end
