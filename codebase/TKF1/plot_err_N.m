s_g = 0.2;
lambda_g = 0.1818;
mu_g = 0.2;

A = load('err_N_noalign.txt');
N = A(:,1);

s = A(:,2);
lambda = A(:,3);
mu = A(:,4);

err_1 = abs(s-s_g);
err_2 = abs(lambda-lambda_g);
err_3 = abs(mu-mu_g);

logN = log10(N);

hold on;
legend('-DynamicLegend','Location','NorthEast');
plot(logN,err_1,'r-','DisplayName','SubRate');
legend('-DynamicLegend','Location','NorthEast');
plot(logN,err_2,'b-','DisplayName','InsRate');
legend('-DynamicLegend','Location','NorthEast');
plot(logN,err_3,'k-','DisplayName','DelRate');
hold off;

xlabel('log(N)');
ylabel('error');

saveas(gcf,'err_N.eps','epsc');
