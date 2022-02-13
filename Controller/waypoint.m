resolution = 4;


r = ones(360/resolution + 1, 1);
theta = (0:resolution:360)';
radian_divider = r*(2*pi)/360;

theta_radian = radian_divider.*theta;
[x, y] = pol2cart(theta_radian, r);
point_cart = [x, y];
point_cart2 = [point_cart; point_cart];