function [x,y] = line_int(m1,x1,y1,m2,x2,y2)
% calcs. intercept of two lines
  if (m1 == inf) | (m1 == -inf) | isnan(m1)
    x = x1;
  else
    x = (y2 - y1 - m2.*x2 + m1.*x1)./(m1-m2);
  end
  y = m2.*(x-x2) + y2;
end
