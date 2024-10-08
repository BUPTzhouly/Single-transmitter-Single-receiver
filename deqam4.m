%�ó���������ɶ������źŽ���4QAM��������ڡ���·�����������ơ��������4QAM���ģ��
%���ߣ�������  
%������ڣ�2005��3��7


function [llr]=deqam4(x,h,SNR_db)
%[llr]=deqam4(x)
%x��1*(Lp+8)�ĸ���������������Ԫ�ر�ʾ���뵽������е��ź�
%hΪ�ŵ���ֵ��������1*(Lp+8)�ĸ���������������Ԫ�ر�ʾ��Ӧ����ʱ�̵��ŵ�״̬��
%SNR_dbΪ�ŵ���ȡ�
%llrΪxͨ��������������ź�,��1*2(Lp+8) ��ʵ��������ʾ������ؽ��������ݡ���Ȼ�ȡ���Ϣ��

%SNR_dbΪ�ŵ�����ȣ��źŹ���Ϊ1����������Ϊn_power
SNR_linr=10^(SNR_db/10);
%��������
n_power=1/SNR_linr;

%len_inputΪ�����źŵĳ���
len_input=length(x);



%�����λ��LLR��Ϣ
temp_1=real(x)>0;
temp_2=2*temp_1-1;
s0_L=(temp_2+j)*sqrt(2)/2;
s1_L=(temp_2-j)*sqrt(2)/2;

%�����źŵ����׼�ź�������֮��ľ���
         d0_L=abs(x-s0_L);
         d1_L=abs(x-s1_L);
%������Ȼ����Ϣllr
llr_L=1.*(d1_L.^2-d0_L.^2)./(n_power.*(1./((abs(h)).^2)));

%�����λ��LLR��Ϣ
temp_1=imag(x)>0;
temp_2=2*temp_1-1;
s0_H=(1+j*temp_2)*sqrt(2)/2;
s1_H=(-1+j*temp_2)*sqrt(2)/2;


%%�����źŵ����׼�ź�������֮��ľ���
       d0_H=abs(x-s0_H);
       d1_H=abs(x-s1_H);
%������Ȼ����Ϣllr
llr_H=1.*(d1_H.^2-d0_H.^2)./(n_power.*(1./((abs(h)).^2)));

%���ߵ�λ�����������������Ȼ������ 
llr_output=[];
for I=1:len_input
       llr_output=[llr_output,llr_L(I),llr_H(I)];
end
llr=llr_output;      
    