clc

[x1,y1]=ncgeodetictoxy(80,0,1);
[x2,y2]=ncgeodetictoxy(60,45,1);
[x3,y3]=ncgeodetictoxy(60,45,-1);

a = sprintf('x1 = %f, y1 = %f\n',x1,y1);
b = sprintf('x2 = %f, y2 = %f\n',x2,y2);
c = sprintf('x3 = %f, y3 = %f\n\n',x3,y3);

disp(a)
disp(b)
disp(c)

[lat1,lon1]=ncxytogeodetic(x1,...
    y1,1);
[lat2,lon2]=ncxytogeodetic(x2,y2,1);
[lat3,lon3]=ncxytogeodetic(x3,y3,-1);

a = sprintf('lat1 = %f, lon1 = %f\n',lat1,lon1);
b = sprintf('lat2 = %f, lon2 = %f\n',lat2,lon2);
c = sprintf('lat3 = %f, lon3 = %f\n\n',lat3,lon3);

disp(a)
disp(b)
disp(c)
