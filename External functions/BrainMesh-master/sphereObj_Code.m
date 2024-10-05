[x,y,z] = sphere;
x = (reshape(x,1,441)*80) + 675;
y = (reshape(y,1,441)*80) + 225;
z = (reshape(z,1,441)*80) + 900;
sphereObj = [x;y;z]';