clear phi X t sigmaINV
M2 = 5;
lambda2=0.4;
load synthetic_data.mat;

[row,col]=size(X);
N2=round(row*0.8);
E = N2;
D = col;
for i = (1:N2)
    trainInd2(i,1) = i;
end

NX = round(row*0.9) - N2;
i = N2 + 1;
for j = (1:NX)
    validInd2(j,1) = i;
    i=i+1;
end
target = t(1:N2);
Sigma2 = 1/3 * diag(var(X) + 10 ^ -2);
sigmaINV = inv(Sigma2);
Sigma2 = repmat(Sigma2,1,1,M2);

p = randperm(N2,M2);
mu2 = X(p,:);

phi(1:N2,1) = 1;
for i = (1:N2)
   for j = (2:M2)
      val = X(i,:) - mu2(j,:);
      phi(i,j) = exp(-1/2 * val * sigmaINV * transpose(val));
   end
end

%%%%%%%%Calculate w
transphi = transpose(phi);
I = eye(M2,M2);
inverseY = inv(( lambda2 * I ) + ( transphi * phi ));
w2 = inverseY * transphi * target;
transW2 = transpose(w2);
%%%%%%%

%%%%%%%Calculate Error trainPer
ED_W=0;
for i = (1:N2)
    multi = ( transW2 * transphi(:,i) );
    ED_W = ED_W + ( target(i,1) - multi )^2;
end
ED_W = 1/2 * ED_W;
EW_W = 1/2 * sum(power(w2,2));
E_W = ED_W + lambda2 * EW_W;
trainPer2 = sqrt(2 * E_W / N2)
%%%%%%%

save x.mat w2 lambda2 sigmaINV M2 mu2 D X t row col
validPer2 = syntheticValid

%%%%%%%%Stochastic

%%%%%%%%%%%%%%%%%Stochastic
clear target mu phi
p = randperm(E,M2);
target = t(1:E);
mu = X(p,:);
clear phi;
phi(1:E,1) = 1;
for i = (1:E)
   for j = (2:M2)
      val = X(i,:) - mu(j,:);
      phi(i,j) = exp(-1/2 * val * sigmaINV * transpose(val));
   end
end

%%%%%%%%Calculate w02
I = eye(M2,M2);
inverseY = inv( ( lambda2 * I ) +phi' * phi );
w02 = inverseY * phi' * target;

%%%%%%%

clear dw2 eta2 x Err
wx=w02;
eta2(1,1)=1.0000008765;
Err = 0;
for i = (1:E)
Err2 = t(i,1) - transpose(wx) * phi(i,:)';
if(i>1)
euclidean = norm(dw2(:,i-1));
    eta2(1,i)=eta2(1,i-1);
if( Err > 0)
    eta1(1,i)=eta1(1,i-1)/1.145;
else
    eta1(1,i)=eta1(1,i-1)*1.00001;
end

end
if i > 1 && norm(wx-w2) < 0.01
    eta2(1,i)=0;
end
Err = Err2;
DeltaE = -Err2 * phi(i,:)' + (1/E)*lambda2 * wx;

DeltaW = -eta2(1,i) * DeltaE;
dw2(:,i) = DeltaW;
wx = wx + DeltaW;
end
euclideanInit = norm(w02 - w2)
euclideanLast = norm(wx - w2)
euclideanInit-euclideanLast;
ED_W=0;EW_W=0;E_W=0;
for i = (1:E)
    
    multi = ( wx' * phi(i,:)' );
    ED_W = ED_W + ( target(i,1) - multi )^2;
    
end

ED_W = 1/2 * ED_W;
EW_W = 1/2 * sum(power(wx,2));
E_W = ED_W + lambda2 * EW_W;
RMSPer2 = sqrt(2 * E_W / E);

mu2=mu2';
%clear transW2 transWX RMSPer2 RMSPer2 mu Err E sigmaINV transW2 X wx val row col I inverseY DeltaE DeltaW DeltaEd DeltaEw E_W ED_W EW_W multi transphi phi X D NX i j p target t