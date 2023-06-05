function x = projhalf(y, a, b)
% project an n-dim vector y to half space
% Dn = { x : a'x <= b}

% Xiaoqian Liu
% May 30, 2023.

if a'*y <= b
    x = y;
else
    x = y - ( (a'*y - b)/norm(a)^2 )*a;
end

end