%xlambdas = [1 1.2 1.5 1.8 2 2.5 3 4 5 6 7 8 9 10 12 15 18 20 25 60 70 150 180 500];
%xplot = [54 47 35.333 30.5 28 22.9 18.7 14.3 11.4 10.25 8.5 7.4 6.9 6.0 5.5 4.3 4 4 3 3 2 2 1 1];

<<<<<<< HEAD
xlambdas = [1.5 1.8 2.5 3 4         5 6 7 8 9              10 12 15 ];
xplot = [35.333 30.5 22.9 18.7 14.3 11.4 10.25 8.5 7.4 6.9 6.0 5.5 4.3 ];

lambdas = [1 1.2 1.5 1.8 2 2.5 3 4 5 6 7 8 9 10 12 15  18 20 25 60 70  150 180 500];
redline = [54 47 35.333 30.5 28 22.9 18.7 14.3 11.4 10.25 8.5 7.4 6.9 6.0 5.5 4.3 4 4 3 3 2 2 1 1];

size(lambdas)
size(redline)
% redplot
XXplot = semilogx(xlambdas, xplot, 'rx', 'markersize', 15); 
hold all 
RLplot = semilogx(lambdas, redline, 'r--');
ISplot = semilogx([1 1.2 2 18 20 25 60 70 150 180 500], [54 47 28 4 4 3 3 2 2 1 1], 'ks', 'markersize', 15);
IRplot = semilogx([18 20], [4 4], 'k', [25 60], [3 3], 'k', [70 150], [2 2], 'k', [180, 500], [1 1], 'k', 'linewidth', 3);
hold off

legend([XXplot(1) ISplot(1) IRplot(1)], 'Fractional Sample', 'Integer Sample', 'Integer Range') 
%, , , IRplot, 'Integer Range'
% ,  , 
=======
lambdas = [1 1.2 1.5 1.8 2 2.5 3 4 5 6 7 8 9 10 12 15 18 20 25 60 70 150 180 500];
redplot = [54 47 35.333 30.5 28 22.9 18.7 14.3 11.4 10.25 8.5 7.4 6.9 6.0 5.5 4.3 4 4 3 3 2 2 1 1];
semilogx(lambdas, redplot, 'rx--', [18], [4], 'ks', 'markersize', 18)
hold on
legend('fractions', 'integers')
semilogx([1 1.2 2 18 20 25 60 70 150 180 500], [54 47 28 4 4 3 3 2 2 1 1], 'ks', 'markersize', 18)
semilogx([18 20], [4 4], 'k', [25 60], [3 3], 'k', [70 150], [2 2], 'k', [180, 500], [1 1], 'k', 'linewidth', 3)
>>>>>>> 5f3e434b71df31984b23b366d65898a50b666bc6

axis([0.8, 600, 0, 55])
xlabel('lambda', 'Fontsize', 18)
ylabel('K', 'Fontsize', 18)
title('Wine', 'Fontsize', 18)
set(gca, 'Fontsize', 25)
