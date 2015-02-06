<<<<<<< HEAD
xlambdas = [1.2 1.5 1.8 2 2.5         6     25 ];
redplot = [23.1 22.9 15.75 14.34 13.9 8.5   3.4 ];

lambdas = [1 1.2 1.5 1.8 2 2.5 3 4 5 6 7 8 10 12 15 20 25 30 200 250 500];
redline = [28 23.1 22.9 15.75 14.34 13.9 13 11 9 8.5 8 6 6 5 4 4 3.4 2 2 1 1];

XXplot = semilogx(xlambdas, redplot, 'rx', 'markersize', 15)
hold on
semilogx(lambdas, redline, 'r--')
ISplot = semilogx([1 3 4 5 7 8 10 12 15 20 30 200 250 500], [28 13 11 9 8 6 6 5 4 4 2 2 1 1], 'ks', 'markersize', 15)
IRplot = semilogx([8 10], [6 6], 'k', [15 20], [4 4], 'k', [30 200], [2 2], 'k', [250, 500], [1 1], 'k', 'linewidth', 3)

legend([XXplot(1) ISplot(1) IRplot(1)], 'Fractional Sample', 'Integer Sample', 'Integer Range') 
=======
lambdas = [1 1.2 1.5 1.8 2 2.5 3 4 5 6 7 8 10 12 15 20 25 30 200 250 500];
redplot = [28 23.1 22.9 15.75 14.34 13.9 13 11 9 8.5 8 6 6 5 4 4 3.4 2 2 1 1];
semilogx(lambdas, redplot, 'rx--', [20], [4], 'ks', 'markersize', 18)
hold on
legend('fractions', 'integers')
semilogx([1 3 4 5 7 8 10 12 15 20 30 200 250 500], [28 13 11 9 8 6 6 5 4 4 2 2 1 1], 'ks', 'markersize', 18)
semilogx([8 10], [6 6], 'k', [15 20], [4 4], 'k', [30 200], [2 2], 'k', [250, 500], [1 1], 'k', 'linewidth', 3)
>>>>>>> 5f3e434b71df31984b23b366d65898a50b666bc6

axis([0.8, 600, 0, 30])
xlabel('lambda', 'Fontsize', 18)
ylabel('K', 'Fontsize', 18)
title('Glass', 'Fontsize', 18)
<<<<<<< HEAD
set(gca, 'Fontsize', 18)
=======
set(gca, 'Fontsize', 25)
>>>>>>> 5f3e434b71df31984b23b366d65898a50b666bc6
