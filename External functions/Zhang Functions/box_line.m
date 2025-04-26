function H = box_line(x,y,e,color)
x = x(:); X = [x;flipud(x); x(1)];
y = y(:);
e = e(:);
Y = [y-e; flipud(y+e); y(1)-e(1)];

H = plot(X,Y,'color',color);
% hold on;
% plot(x,y,'color',color);
end