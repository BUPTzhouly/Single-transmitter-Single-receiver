%�ó���������ɶ������źŽ���OFDM��������ڡ���·�����������ơ�������ĵ��ƽ��ģ��
%���ߣ�������  
%������ڣ�2005��3��7
              
function [y]=deofdm(x)
%[y]=deofdm(x)
%x�е�����Ԫ�ر�ʾOFDM�ź�ͨ���ŵ���Ľ����
%h������Ԫ�ر�ʾOFDM���ն˴���Ľ����hΪ�����õ��ŵ���Ϣ

%������ʼ��

%���ز���
subcarrier_num=1024;
%ʹ�õ����ز���
used_subcarrier_num=150;

%ѭ��ǰ׺�ĳ���
%ǰ�������ŵ�ѭ��ǰ׺��cp1=73
cp1=72;
%���߸����ŵ�ѭ��ǰ׺��cp2=74
cp2=80;

%��������ݳ���
input_s_len=length(x);

%ÿ��OFDM���ŵĲ�������Ŀ
symbol_sample1=subcarrier_num+cp1;
%����ķ�����Ŀ
symbol_num=floor(input_s_len/symbol_sample1);
%δʹ�õ����ز���Ŀ
delete_s=subcarrier_num-(used_subcarrier_num+1);

%ÿһ��OFDM����ռ�õ����ز�
subcarrier_used_left=floor(used_subcarrier_num/2);
subcarrier_used_right=ceil(used_subcarrier_num/2);
y=[];
%�������źŽ��зָ�ָ�Ϊsymbol_num�����ţ��ٶ�ÿ�����ŷֱ����IFFT���㣬ʵ��OFDM����,����֤��������
for I=1:symbol_num
    %��ÿ������ȥѭ��ǰ׺
        if(I~=symbol_num) 
              x_temp=x(((I-1)*symbol_sample1+cp1+1):(I*symbol_sample1));
        else
              x_temp=x(((I-1)*symbol_sample1+cp2+1):(I*symbol_sample1+cp2-cp1));
        end
    %��ÿ�����Ž���FFT����
        fre_domain_x_temp=fft(x_temp)/sqrt(subcarrier_num);
    %ȥ������ʱ��ӵ���Ϣ��
        fre_domain_x_del=[fre_domain_x_temp(2:subcarrier_used_left+1),fre_domain_x_temp((subcarrier_num-subcarrier_used_right+1):subcarrier_num)];
   %�õ�OFDM����Ľ��
       y=[y,fre_domain_x_del];
end


