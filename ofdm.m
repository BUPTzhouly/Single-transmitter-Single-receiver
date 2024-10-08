%�ó���������ɶ������źŽ���OFDM����
%���ߣ�������  
%������ڣ�2005��3��7
%�޸����ڣ�2006��4��17

function [y]=ofdm(x)
%[y]=ofdm(x)
%x��1x(Lp+8)�ĸ���������������Ԫ��Ϊ4QAM���Ʒ��ţ���ʾ����4QAM���ƵĽ�����ݡ�
%y��1x(67.5*10e(-6)/Ts) �ĸ���������������Ԫ����OFDM���Ͷ˴���Ľ������ʾ���Ͷ˵Ļ����ź���ʱ���ϣ���������ΪTs�Ĳ���������źų���ʱ��Ϊ ��3 symbols����

%������ʼ������
n=nargin;
%�����Ƶ���ź��а�����OFDM������Ŀ
if n~=1
  error('input argument error');
end

%����һ֡�к��е�OFDM���Ÿ���
s=7;
%������źų���Ϊ
L=length(x);
%���ز���
subcarrier_num=1024;
%ʹ�õ����ز���
used_subcarrier_num=150;

%ѭ��ǰ׺�ĳ���
%ǰ�������ŵ�ѭ��ǰ׺��cp1=72
cp1=72;
%���߸����ŵ�ѭ��ǰ׺��cp2=80
cp2=80;
%ÿһ��OFDM�����м�Ӧ�ò���0������zeros_pad
zeros_pad=subcarrier_num-(used_subcarrier_num+1);

%ÿһ��OFDM����ռ�õ����ز�
subcarrier_used_left=floor(used_subcarrier_num/2);
subcarrier_used_right=ceil(used_subcarrier_num/2);

%�������źŽ��зָ�ָ�Ϊs�����ţ��ٶ�ÿ�����Ž���FFT���㣬ʵ��OFDM���,����֤��������
time_domain_x_link=[];
for I=1:s
    %��������зָ�,�ָ��Ϊ7�Σ�ÿ�γ���Ϊdata_num  
        x_temp_data=x(used_subcarrier_num*(I-1)+1:(used_subcarrier_num*I));
    %��ÿ���ָ�Ĳ��ֽ��в��������ʹ�䳤Ϊsub_carriers
        x_temp_pad=[0,x_temp_data(1:subcarrier_used_left),zeros(1,zeros_pad),x_temp_data((1+subcarrier_used_left):used_subcarrier_num)];
    %��ÿ�����Ž���IFFT����
        time_domain_x_temp=ifft(x_temp_pad)*sqrt(subcarrier_num);
  
    %��ÿ���������ѭ��ǰ׺��ǰ6������CP=72,��7������CP=80
       if(I~=s)
             time_domain_x_cp_temp=[time_domain_x_temp((subcarrier_num-cp1+1):subcarrier_num),time_domain_x_temp];
       else
             time_domain_x_cp_temp=[time_domain_x_temp((subcarrier_num-cp2+1):subcarrier_num),time_domain_x_temp];
       end
   
    %���������ӳ�Ϊ����������
       time_domain_x_link=[time_domain_x_link,time_domain_x_cp_temp];
end
%����Ϣ�����
y=time_domain_x_link;
