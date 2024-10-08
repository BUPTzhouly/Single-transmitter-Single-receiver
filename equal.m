%该程序用来完成对输入信号进行均衡，属于《链路级仿真软件设计》程序二的均衡模块
%作者：赵亚利  
%编程日期：2005－3－7


function [y]=equal(x,h,U,S)
%[y]=equal(x,h)
%S    =2;               % number of receiving antennas
%U    = 2;              % number of transmitting antennas
%x是1*(Lp+8)的复向量，其中向量元素表示表示经过OFDM解调后的输出
%h是信道特性估计,由信道模块给出,为1*(Lp+8)的复向量，其中向量元素表示OFDM接收端处理的结果
%y为s通过均衡器后的输出信号,是1x(Lp+8) 的复向量，其中向量元素表示经过均衡处理的结果
 %L=length(h)/4;
%q=0;
% for u=1:U%得到对应信道的上的信道函数
%    for s=1:S
%        q=q+abs(h((u-1)*2+s,:)).^2;
%    end
%end
 
%q=abs(h(4,:)).^2;
%均衡权值向量为q,由信道模块给出的h来计算
%q=1./h;
%均衡输出为r
y=x./h;
