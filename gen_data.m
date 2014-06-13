% K: number of clusters
% N: number of data
% D: number of dimension
% p_1~p_{K+1}: mix prob. of clusters (K+1 is noise)
% SIGMA: std among clusters
% sigma1~sigmaK: std within clusters
% mu1~muK: means of clusters

K = 4;
D = 2;
N = 10000;
p = [0.225 0.225 0.225 0.225    0.1];
SIGMA = 4;
sigma = ones(K,1)*1;
mu = randn(K,D)*SIGMA;
R = 8;

z = zeros(N,1);
X = zeros(N,D);
for i = 1:N
	r = rand;
	psum = p(1);
	k = 1;
	while r > psum 
		k = k + 1;
		psum = psum + p(k);
	end

	if k <= K
		z(i) = k;
		X(i,:) = randn(1,D)*sigma(k) + mu(k,:);
	else % it's noise
		z(i) = 0;
		X(i,:) = (rand(1,D)-0.5)*2*R ;
	end
end

f = fopen('simdata','w');
for i=1:N
	fprintf(f, '%d', z(i));
	for j=1:D
		fprintf(f, ' %d:%g', j, X(i,j));
	end
	fprintf(f, '\n');
end
fclose(f);
