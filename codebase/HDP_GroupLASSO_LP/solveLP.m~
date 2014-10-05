Y = load('Y');
Z = load('Z');
d_assign = load('dataset_assign');

[N,K] = size(Y);
D = max(d_assign);
n = N*K + D*K + K;

rho = 1;
lambda= 10;

H = sparse(n,n);
for i = 1:N*K
	H(i,i) = rho;
end

C = Y - rho*Z;
c = C(:);
c = [c ; lambda*ones(K,1)];

A = sparse(N*K, n);
b = zeros(N*K,1);
count = 1;
for j = 1:K
	for i = 1:N
		A(count, (j-1)*N+i) = 1;
		A(count, N*K + j ) = -1;
		b(count) = 0.0;
		count = count + 1;
	end
end

lb = zeros(n,1);
options = optimset('Display','iter','Algorithm','interior-point-convex','TolFun',1e-13,'MaxIter',200);
%options = optimoptions(@quadprog,'Algorithm','interior-point-convex');
[x,fval] = quadprog(H,c,A,b,[],[],lb,[],[],options);

W = reshape(x(1:N*K),N,K);

Wt = W';
fp = fopen('W_out','w');
for i = 1:N
	fprintf(fp,'%g ',Wt(:,i));
	fprintf(fp,'\n');
end
fclose(fp);

fp = fopen('xi_out','w');
fprintf(fp,'%g\n',x(N*K+1:end));
fclose(fp);

fp = fopen('fval','w');

%fval = fval - sum(sum(Y.*Z)) + rho*sum(sum(Z.*Z))/2;
%fprintf(fp,'fval=%g\n',fval);

xi_sum = 0.0;
for k= 1:K
	xi_sum = xi_sum + x(N*K+k);
end

term1 = lambda*xi_sum;
term2 = sum(sum(Y.*(W-Z)));
term3 = rho*sum(sum((W-Z).*(W-Z)))/2;
fval = term1 + term2 + term3 ;

fprintf(fp,'lambda*xi_sum = %g\n', term1);
fprintf(fp,'Y(W-Z)=%g\n', term2);
fprintf(fp,'rho/2*|W-Z|^2=%g\n', term3 );
fprintf(fp,'fval=%g\n',fval);

fclose(fp);
