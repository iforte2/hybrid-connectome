function [J] = pseudo(binarizedData,SC,lambda,B)

SC(SC == diag(SC)) =0; % make sure diagonal is zero for consistency
[nodeNumber,dataLength] = size(binarizedData);
iterationMax = 50000;
dt = .1;
%lambda = 1;

dataCorrelation = (binarizedData*binarizedData')/dataLength;
J = zeros(nodeNumber);

for t=1:iterationMax
 %gradient ascent   
    LdJ = -dataCorrelation + .5 * binarizedData * (tanh(B.*J * binarizedData))'/dataLength ...
  + .5 * (binarizedData * (tanh(B.*J * binarizedData))')'/dataLength; %based on Callen's identity
    LdJ = LdJ - diag(diag(LdJ));%eliminate diagonal entities (J_{ii}=0)
    %}
 
%}
%update
J = J - (dt.*LdJ)- ((dt.*lambda)*(J-(sign(J).*SC)));
 
if(t>1 && corr2(J,J-dt*LdJ-(dt.*lambda.*(J-(sign(J).*SC))))== 1)

        break
end

end


