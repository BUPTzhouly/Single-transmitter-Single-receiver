% 该程序用来模拟信号通过AWGN信道，属于《链路级仿真软件设计》程序一的信道模块
function [y]=myawgn(x,SNR_in_dB)
%y= myawgn(x,SNR_in_dB)
% x是长度为1×Lp的已调输入信号；SNR_in_dB为输入信噪比，单位为dB。
%SNR_in_dB=S/N,Eb＝1为输入信号平均符号能量，N=N0/2= σ^2为噪声能量，N0为白噪声的单边功率谱密度 
%SNR=S/N=log2(M)*Eb/N0                             
% y为x通过AWGN信道后的输出信号
%作者：娄文科 日期：2005－2－22
if nargin~=2 
    error('input arguments are not matched ')
end
%计算输入信号序列长度                         
N=length(x);
%利用系统给定的信噪比计算出高斯白噪声的能量
SNR_linr=exp(SNR_in_dB*log(10)/10);
noise_power=1/SNR_linr;
%产生长度为N的复高斯白噪声向量
noise_vector= (sqrt(noise_power/2))*(randn(1,N)+j*randn(1,N));
%对信号叠加噪声
y=x+noise_vector;