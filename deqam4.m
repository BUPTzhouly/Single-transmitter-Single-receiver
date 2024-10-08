%该程序用来完成对输入信号进行4QAM解调，属于《链路级仿真软件设计》程序二的4QAM解调模块
%作者：赵亚利  
%编程日期：2005－3－7


function [llr]=deqam4(x,h,SNR_db)
%[llr]=deqam4(x)
%x是1*(Lp+8)的复向量，其中向量元素表示输入到解调器中的信号
%h为信道估值向量，是1*(Lp+8)的复向量，其中向量元素表示对应采样时刻的信道状态。
%SNR_db为信道噪比。
%llr为x通过解调器后的输出信号,是1*2(Lp+8) 的实向量，表示经过相关解调后的数据“似然比“信息。

%SNR_db为信道信噪比，信号功率为1，噪声功率为n_power
SNR_linr=10^(SNR_db/10);
%噪声方差
n_power=1/SNR_linr;

%len_input为输入信号的长度
len_input=length(x);



%计算低位的LLR信息
temp_1=real(x)>0;
temp_2=2*temp_1-1;
s0_L=(temp_2+j)*sqrt(2)/2;
s1_L=(temp_2-j)*sqrt(2)/2;

%计算信号点与标准信号星座点之间的距离
         d0_L=abs(x-s0_L);
         d1_L=abs(x-s1_L);
%计算似然比信息llr
llr_L=1.*(d1_L.^2-d0_L.^2)./(n_power.*(1./((abs(h)).^2)));

%计算高位的LLR信息
temp_1=imag(x)>0;
temp_2=2*temp_1-1;
s0_H=(1+j*temp_2)*sqrt(2)/2;
s1_H=(-1+j*temp_2)*sqrt(2)/2;


%%计算信号点与标准信号星座点之间的距离
       d0_H=abs(x-s0_H);
       d1_H=abs(x-s1_H);
%计算似然比信息llr
llr_H=1.*(d1_H.^2-d0_H.^2)./(n_power.*(1./((abs(h)).^2)));

%将高低位排序，输出解调结果的似然比序列 
llr_output=[];
for I=1:len_input
       llr_output=[llr_output,llr_L(I),llr_H(I)];
end
llr=llr_output;      
    