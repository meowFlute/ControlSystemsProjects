clear;
clc;

t = 90:-10:-270;
x = cosd(t);
y = sind(t);

%{
for i=1:numel(t);
    theta = pi/2 - atan2(y(i),x(i));
    if theta < 0, theta = theta + 2*pi; end
    disp(['(',num2str(mod(t(i),360)),') => [',num2str(theta*180/pi),']']);
end
%}

for i=1:numel(t);
    %atan2(y(i),x(i))*180/pi
    theta = mod(pi/2 + 2*pi - atan2(y(i),x(i)),2*pi);
    disp(['(',num2str(mod(t(i),360)),') => [',num2str(theta*180/pi),']']);
end