function max8 = max_8(in)
%UNTITLED4 此处提供此函数的摘要
%   此处提供详细说明
max_temp=in(1);
for max_i=1:8
    if in(max_i) > max_temp
			max_temp=in(max_i);	
    end
end

max8 = max_temp;

end
