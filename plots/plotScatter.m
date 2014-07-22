function plotScatter(fname)

color_arr = { 'm' 'k' 'r' 'g' 'b' 'c' };

A = load(fname);
[N,tmp] = size(A);

labels = A(:,1); 
pos = A(:,2:3);

hold on;
for i=1:N
	plot( pos(i,1), pos(i,2), ['x' color_arr{labels(i)+1}] );
end
%axis([0,1,0,1]);

saveas(gcf,[fname '.pdf'],'pdf');
exit(0);
