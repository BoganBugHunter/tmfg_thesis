%find the minimum area line


possibleDistances = zeros(size(x,1),size(x,1));
for i = 1:countAlong
    point_on_bottom_x = x(i,1);
    point_on_bottom_y = y(i,1);
    for j = 1:countAlong
        point_on_top_x = x(j,size(x,1));
        point_on_top_y = y(j,size(y,1));
        possibleDistances(i,j) = sqrt( (point_on_bottom_x-point_on_top_x)^2 +...
            (point_on_bottom_y-point_on_top_y)^2 );
    end
end

[mi,index] = min(possibleDistances);

[mi2,index2] = min(mi);

minimumDistance = possibleDistances(index(index2),index2);



line_x = [ x(index(index2),1) ; x(index(index2),size(x,1)) ];
line_y = [ y(index(index2),1) ; y(index(index2),size(y,1))  ];

figure(20)
hold on
axis equal
plot(x,y);
plot(line_x,line_y,'k')



