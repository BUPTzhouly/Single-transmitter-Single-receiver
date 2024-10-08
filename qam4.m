%该程序用来模拟信号4QAM调制，属于《链路级仿真软件设计》程序二的调制解调模块中QAM调制部分
%作者：赵亚利  
%编程日期：2005－3－7


function [y]=qam4(x)
% y= qam4(x)
% x是1*(2(Lp+8)) 的向量，其中向量元素为[1] 或[0]，表示表示经过交织的结果数据
% y为x通过4qam调制后的输出信号,是1x(Lp+8) 的向量

% 建立符号映射关系
%(00->sqrt(2)/2+sqrt(2)/2*j;01->-sqrt(2)/2+sqrt(2)/2*j;11->-sqrt(2)/2-sqrt(2)/2*j;10->sqrt(2)/2-sqrt(2)/2*j)
mapping=[sqrt(2)/2+sqrt(2)/2*j,sqrt(2)/2-sqrt(2)/2*j,-sqrt(2)/2+sqrt(2)/2*j,-sqrt(2)/2-sqrt(2)/2*j]; 
% 取得输入二进制序列长度
len=length(x);
% 进行符号映射
y=zeros(1,len/2);
for I=1:len/2
    temp=x(2*I-1)+x(2*I)*2;
    y(I)=mapping(temp+1);
end
