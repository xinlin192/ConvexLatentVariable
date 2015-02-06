xlambdas = [50 55 60 65 ];
Xpoints = [232 229 218 211];

lambdas = [50 55   60 65   70 155 160 300];
redline = [232 229 218 211 2 2    1   1];

size(xlambdas)
size(Xpoints)
MYPLOT = @semilogx

XXplot = MYPLOT(xlambdas, Xpoints, 'rx', 'markersize', 15)
hold all
redplot = MYPLOT(lambdas, redline, 'r--')
ISplot = MYPLOT([70], [2], 'ks', [155], [2], 'ks', [160], [1], 'ks', [300],[1], 'ks', 'markersize', 15)
IRplot = MYPLOT([70 155], [2 2], 'k', [160 300], [1 1], 'k', 'linewidth', 3)

legend([XXplot(1) ISplot(1) IRplot(1)], 'Fractional Sample', 'Integer Sample', 'Integer Range') 


axis([40, 350, -5, 250])
xlabel('lambda', 'Fontsize', 18)
ylabel('K', 'Fontsize', 18)
title('DNA', 'Fontsize', 18)
set(gca, 'Fontsize', 18)
