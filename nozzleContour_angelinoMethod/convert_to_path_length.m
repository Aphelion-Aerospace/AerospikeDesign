function [x,s] = convert_to_path_length(x,y)
% given discrete points, (x_i,y_i), convert them to cumulative path length
% coordinates

y_dummy = y(2:end) - y(1:(end-1));
x_dummy = x(2:end) - x(1:(end-1));
s_dummy = sqrt(y_dummy.^2 + x_dummy.^2);

s_dummy = [0, s_dummy];

s = zeros(size(s_dummy));

for i = 2:length(s_dummy)
    s(i) = s_dummy(i) + s(i-1);
end

end