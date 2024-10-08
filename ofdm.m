%该程序用来完成对输入信号进行OFDM调制
%作者：赵亚利  
%编程日期：2005－3－7
%修改日期：2006－4－17

function [y]=ofdm(x)
%[y]=ofdm(x)
%x是1x(Lp+8)的复向量，其中向量元素为4QAM调制符号，表示经过4QAM调制的结果数据。
%y是1x(67.5*10e(-6)/Ts) 的复向量，其中向量元素是OFDM发送端处理的结果，表示发送端的基带信号在时域上，采样速率为Ts的采样结果，信号持续时间为 （3 symbols）。

%参数初始化设置
n=nargin;
%输入的频域信号中包含的OFDM符号数目
if n~=1
  error('input argument error');
end

%输入一帧中含有的OFDM符号个数
s=7;
%输入的信号长度为
L=length(x);
%子载波数
subcarrier_num=1024;
%使用的子载波数
used_subcarrier_num=150;

%循环前缀的长度
%前六个符号的循环前缀是cp1=72
cp1=72;
%第七个符号的循环前缀是cp2=80
cp2=80;
%每一个OFDM符号中间应该补‘0’个数zeros_pad
zeros_pad=subcarrier_num-(used_subcarrier_num+1);

%每一个OFDM符号占用的子载波
subcarrier_used_left=floor(used_subcarrier_num/2);
subcarrier_used_right=ceil(used_subcarrier_num/2);

%对输入信号进行分割，分割为s个符号，再对每个符号进行FFT运算，实现OFDM解调,并保证能量不变
time_domain_x_link=[];
for I=1:s
    %对输入进行分割,分割成为7段，每段长度为data_num  
        x_temp_data=x(used_subcarrier_num*(I-1)+1:(used_subcarrier_num*I));
    %对每个分割的部分进行补零操作，使其长为sub_carriers
        x_temp_pad=[0,x_temp_data(1:subcarrier_used_left),zeros(1,zeros_pad),x_temp_data((1+subcarrier_used_left):used_subcarrier_num)];
    %对每个符号进行IFFT运算
        time_domain_x_temp=ifft(x_temp_pad)*sqrt(subcarrier_num);
  
    %对每个符号添加循环前缀，前6个符号CP=72,第7个符号CP=80
       if(I~=s)
             time_domain_x_cp_temp=[time_domain_x_temp((subcarrier_num-cp1+1):subcarrier_num),time_domain_x_temp];
       else
             time_domain_x_cp_temp=[time_domain_x_temp((subcarrier_num-cp2+1):subcarrier_num),time_domain_x_temp];
       end
   
    %将符号连接成为串行数据流
       time_domain_x_link=[time_domain_x_link,time_domain_x_cp_temp];
end
%将信息流输出
y=time_domain_x_link;
