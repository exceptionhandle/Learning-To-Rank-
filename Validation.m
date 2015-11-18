
function RMS = Validation()
RMS = evaluate();
end

function RMS = evaluate(X,strt_N,end_N)

load x.mat

strt_N = round(row*0.8)+1;
end_N = round(row*0.9);
X = X(strt_N:end_N,:);
target = t(strt_N:end_N,1);

phi(1:end_N-strt_N,1)=1;
for i = (1:end_N - strt_N)
   for j = (2:M1)
      val = X(i,1:D) - mu1(j,1:D);
      phi(i,j) = exp(-1/2 * val * sigmaINV * transpose(val));
   end
end
transphi = transpose(phi);

ED_W=0;E_W=0;EW_W=0;
for i = (1:end_N - strt_N)
    ED_W = ED_W + (target(i,1) - transpose(w1) * transphi(:,i)) ^ 2;
end
ED_W = 1/2 * ED_W;
EW_W = 1/2 * sum(power(w1,2));
E_W = ED_W + lambda1 * EW_W;
RMS = sqrt(2 * E_W / (end_N - strt_N));

end
