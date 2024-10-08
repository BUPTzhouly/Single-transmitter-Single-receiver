%该程序用来完成对输入信号进行OFDM解调，属于《链路级仿真软件设计》程序二的调制解调模块
%作者：赵亚利  
%编程日期：2005－3－7
              
function [y]=deofdm(x)
%[y]=deofdm(x)
%x中的向量元素表示OFDM信号通过信道后的结果。
%h中向量元素表示OFDM接收端处理的结果，h为均衡用的信道信息

%参数初始化

%子载波数
subcarrier_num=1024;
%使用的子载波数
used_subcarrier_num=150;

%循环前缀的长度
%前六个符号的循环前缀是cp1=73
cp1=72;
%第七个符号的循环前缀是cp2=74
cp2=80;

%输入的数据长度
input_s_len=length(x);

%每个OFDM符号的采样点数目
symbol_sample1=subcarrier_num+cp1;
%输入的符号数目
symbol_num=floor(input_s_len/symbol_sample1);
%未使用的子载波数目
delete_s=subcarrier_num-(used_subcarrier_num+1);

%每一个OFDM符号占用的子载波
subcarrier_used_left=floor(used_subcarrier_num/2);
subcarrier_used_right=ceil(used_subcarrier_num/2);
y=[];
%对输入信号进行分割，分割为symbol_num个符号，再对每个符号分别进行IFFT运算，实现OFDM调制,并保证能量不变
for I=1:symbol_num
    %对每个符号去循环前缀
        if(I~=symbol_num) 
              x_temp=x(((I-1)*symbol_sample1+cp1+1):(I*symbol_sample1));
        else
              x_temp=x(((I-1)*symbol_sample1+cp2+1):(I*symbol_sample1+cp2-cp1));
        end
    %对每个符号进行FFT运算
        fre_domain_x_temp=fft(x_temp)/sqrt(subcarrier_num);
    %去除调制时添加的信息点
        fre_domain_x_del=[fre_domain_x_temp(2:subcarrier_used_left+1),fre_domain_x_temp((subcarrier_num-subcarrier_used_right+1):subcarrier_num)];
   %得到OFDM输出的结果
       y=[y,fre_domain_x_del];
end


