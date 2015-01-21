
lambdas = [1 1.2 1.5 2 2.5 3 4 5 6 20 25 100 120 500];
redplot = [10 7.99974 7 7 5.99914 5.006 4.006 4.00138 3 3 2 2 1 1];
semilogx(lambdas, redplot, 'rx--', [1.5], [7], 'ks', 'markersize', 18)
hold on
legend('fractions', 'integers')
semilogx([1], [10], 'ks', [2], [7], 'ks', [6 20], [3 3], 'ks', [25 100], [2 2], 'ks',[120, 500], [1 1],'ks', 'markersize', 18)
semilogx([1.5 2], [7 7], 'k', [6 20], [3 3], 'k', [25 100], [2 2], 'k', [120, 500], [1 1], 'k', 'linewidth', 3)

axis([0.8, 600, 0, 11])
xlabel('lambda', 'Fontsize', 18)
ylabel('K', 'Fontsize', 18)
title('Iris', 'Fontsize', 18)
set(gca, 'Fontsize', 25)
