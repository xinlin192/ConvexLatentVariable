
lambdas = [1 1.2 1.5 1.8 2 2.5 3 4 5 6 7 8 9 10 12 15 18 20 25 60 70 150 180 500];
redplot = [54 47 35.333 30.5 28 22.9 18.7 14.3 11.4 10.25 8.5 7.4 6.9 6.0 5.5 4.3 4 4 3 3 2 2 1 1];
semilogx(lambdas, redplot, 'rx--', [18], [4], 'ks', 'markersize', 18)
hold on
legend('fractions', 'integers')
semilogx([1 1.2 2 18 20 25 60 70 150 180 500], [54 47 28 4 4 3 3 2 2 1 1], 'ks', 'markersize', 18)
semilogx([18 20], [4 4], 'k', [25 60], [3 3], 'k', [70 150], [2 2], 'k', [180, 500], [1 1], 'k', 'linewidth', 3)

axis([0.8, 600, 0, 55])
xlabel('lambda', 'Fontsize', 18)
ylabel('K', 'Fontsize', 18)
title('Wine', 'Fontsize', 18)
set(gca, 'Fontsize', 25)
