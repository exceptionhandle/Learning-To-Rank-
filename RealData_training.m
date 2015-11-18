clear %phi X t sigmaINV
M1 = 2;
lambda1=0.001;

load readData.mat;

[row,col]=size(X);
N1=round(row*0.8);
E = row;%N1;
D = col;
for i = (1:N1)
    trainInd1(i,1) = i;
end

NX = round(row*0.9) - N1;
i = N1 + 1;
for j = (1:NX)
    validInd1(j,1) = i;
    i=i+1;
end
target = t(1:N1);
Sigma1 = 1/10 * diag(var(X) + 10 ^ -2);
sigmaINV = inv(Sigma1);
Sigma1 = repmat(Sigma1,1,1,M1);

p = randperm(N1,M1);
mu1 = X(p,:);

phi(1:N1,1) = 1;
for i = (1:N1)
   for j = (2:M1)
      val = X(i,:) - mu1(j,:);
      phi(i,j) = exp(-1/2 * val * sigmaINV * transpose(val));
   end
end

%%%%%%%%Calculate w
transphi = transpose(phi);
I = eye(M1,M1);
inverseY = inv(( lambda1 * I ) + ( transphi * phi ));
w1 = inverseY * transphi * target;
transW1 = transpose(w1);
%%%%%%%

%%%%%%%Calculate Error trainPer
ED_W=0;
for i = (1:N1)
    multi = ( transW1 * transphi(:,i) );
    ED_W = ED_W + ( target(i,1) - multi )^2;
end
ED_W = 1/2 * ED_W;
EW_W = 1/2 * sum(power(w1,2));
E_W = ED_W + lambda1 * EW_W;
trainPer1 = sqrt(2 * E_W / N1)
%%%%%%%

save x.mat w1 lambda1 sigmaINV M1 mu1 D X t row col
validPer1 = Validation

%%%%%%%%Stochastic

%%%%%%%%%%%%%%%%%Stochastic
p = randperm(E,M1);
target = t(1:E);
mu = X(p,:);
clear phi;
phi(1:E,1) = 1;
for i = (1:E)
   for j = (2:M1)
      val = X(i,:) - mu(j,:);
      phi(i,j) = exp(-1/2 * val * sigmaINV * transpose(val));
   end
end

%%%%%%%%Calculate w01
w01(1:M1,1) = 0;

%%%%%%%

wx=w01;
eta1(1,1)=1.00008;
Err = 0;
for i = (1:E)
Err2 = t(i,1) - (wx)' * phi(i,:)';
if(i>1)
euclidean = norm(dw1(:,i-1));
if( Err > 0)
    eta1(1,i)=eta1(1,i-1)/1.145;
else
    eta1(1,i)=eta1(1,i-1)*1.00001;
end
end
if i > 1 && norm(wx-w1) < 0.001
    eta1(1,i)=0;
end

Err = Err2;
DeltaE = -Err2 * phi(i,:)' + (1/E)*lambda1 * wx;

DeltaW = -eta1(1,i) * DeltaE;
dw1(:,i) = DeltaW;
wx = wx + DeltaW;
end
euclideanInit = norm(w01 - w1)
euclideanLast = norm(wx - w1)


ED_W=0;E_W=0;EW_W=0;
for i = (1:E)
    
    multi = ( (wx)' * phi(i,:)' );
    
    ED_W = ED_W + ( target(i,1) - multi )^2;
    
end

ED_W = 1/2 * ED_W;
EW_W = 1/2 * sum(power(wx,2));
E_W = ED_W + lambda1 * EW_W;
RMSPer1 = sqrt(2 * E_W / E);

mu1=transpose(mu1);




%clear transW1 transWX RMSPer1 RMSPer1 mu Err E sigmaINV transW1 X wx val row col I inverseY DeltaE DeltaW DeltaEd DeltaEw E_W ED_W EW_W multi transphi phi X D NX i j p target t