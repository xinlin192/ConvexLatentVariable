xlambdas = [0.5 0.8        1.2      2.5     3     4     5       ];
Xpoints = [13.0006 10.3995 7.99974  5.99914 5.006 4.006 4.00138 ];

lambdas = [0.5 0.8 1 1.2 1.5 2 2.5 3 4 5 6 20 25 100 120 500];
redline = [13.0006 10.3995 10 7.99974 7 7 5.99914 5.006 4.006 4.00138 3 3 2 2 1 1];

size(xlambdas)
size(Xpoints)

XXplot = semilogx(xlambdas, Xpoints, 'rx', 'markersize', 15)
hold all
redplot = semilogx(lambdas, redline, 'r--')
ISplot = semilogx([1.5], [7], 'ks', [1], [10], 'ks', [2], [7], 'ks', [6 20], [3 3], 'ks', [25 100], [2 2], 'ks',[120, 500], [1 1],'ks', 'markersize', 15)
IRplot = semilogx([1.5 2], [7 7], 'k', [6 20], [3 3], 'k', [25 100], [2 2], 'k', [120, 500], [1 1], 'k', 'linewidth', 3)

legend([XXplot(1) ISplot(1) IRplot(1)], 'Fractional Sample', 'Integer Sample', 'Integer Range') 


axis([0.8, 600, 0, 12])
xlabel('lambda', 'Fontsize', 18)
ylabel('K', 'Fontsize', 18)
title('Iris', 'Fontsize', 18)
set(gca, 'Fontsize', 18)
