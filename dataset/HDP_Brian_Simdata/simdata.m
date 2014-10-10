K = 15;
D = 2;
sigma = 0.1;
M = 50; %number of data set
Km = 5; %number of clusters per data set
Nmk = 5; % number of points for each topic in each data set 

MU = rand(D,K);

label = zeros(1,Nmk*Km*M);
dataset = zeros(1,Nmk*Km*M);
data = zeros(D,Nmk*Km*M);

i = 1;
for m = 1:M
	for k = 1:Km
		z = ceil(rand*K);
		for n = 1:Nmk
			label(:,i) = z;
			dataset(:,i) = m;
			data(:,i) = MU(:,z) + sigma*randn(D,1);
			i = i+1;
		end
	end
end

fp  = fopen('ground_truth','w');
fprintf(fp,'%d\n',label);
fclose(fp);

fp = fopen('data','w');
for i =1 : size(data,2)
	fprintf(fp,'%d %g %g\n',dataset(1,i), data(1,i), data(2,i));
end
fclose(fp);
