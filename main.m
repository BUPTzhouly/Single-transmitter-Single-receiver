 %������ʵ�� ��������С����̿���ģ�顱�Ĺ��ܣ��μ�����·������������_v2.1����
%���ߣ��ֻ�
%�޸ļ�¼��2005��2��25
%��չ������е����б���
clear;
%�����趨:
%SCM����
U=1;% ��������
S=1;%��������
mul_path=6;

%�������е�����ȵ���
M=6;
%���Ƹ�������ȵ�ѡ���������а���
N=1;

%����������
SNs=[ 0  2  4  6  8  10];
T=7680;
%��·�����ݰ��Ĵ�С
Lp=696;%(bits)

%�������еĽ���������
%�������е���������ݰ�����
Ps=[5000,5000,5000,5000,5000,5000,5000,5000,5000]*2;
%�������е���������
Pes=[100,100,100,100,100,100,100,100,100];
%��������
Tb=0.5*10^(-3)/Lp;%(s)
%�������
Ts=1*10^(-3)/15/1024;%(S)
%����ʱ��
T_symbo=0.5*10^(-3)/7;
%%%%%%%%%%%%%%%���ý���洢����:
%������ȵ�������ݰ���Ŀ
Ps_run_s=zeros(N,M);
%������ȵ���ܱ�����
Bits_run_s=zeros(N,M);
%������ȵ���ܴ������ݰ���Ŀ
Ps_err_s=zeros(N,M);
%������ȵ���ܴ��������
Bits_err_s=zeros(N,M);
%������ȵ�������
P_err_rate_s=zeros(N,M);
%������ȵ��������
Bit_err_rate_s=zeros(N,M);

%������֯���
interleave_table=interleav_matrix(3*(Lp+4));

%���ò��������������
rand('state',sum(100*clock));
randn('state',sum(7*clock));

%ʹ��ѭ����������M������ȵ�
for I=1:M
    I
    tic
	%��ʼ��:
	%���͵����ݰ�����
 	Ps_Tx=0;
	%��������ݰ�����
	Ps_Rx_err=0;
	%���͵��ܱ�����
	Bits_Tx=0;
	%����ı�����Ŀ
	Bits_Rx_err=0;
	%�÷����������
	SN=SNs(I);
	%�÷����������е����ݰ���Ŀ
	P_max=Ps(I);
	%�÷����������еĴ������ݰ���Ŀ
	Pe_max=Pes(I);
	%ʹ��ѭ�������������P_max�����ݰ�
	for II=1:P_max
        II
		%����1����ΪLp�����ݰ�:
		P_Tx=(rand(1,Lp)>0.5);
        
        %�ӿ�1��turbo���룩       
	    sub_interleaver=sub_interleave(Lp);
        turbo_out=turbo_encode_1_3(P_Tx,sub_interleaver);
        
        %�ӿ�1.1����֯������
        %interleav_out = interleaving(turbo_out ,interleave_table);
		
		%�ӿ�2��4QAM ���ƣ�
		%Qam_out=qam4(interleav_out);
        Qam_out=qam4(turbo_out);
                   
        %�ӿ�2.1��OFDMģ�飩��
        ofdm_out=ofdm(Qam_out);
   
        %��ʼ�������ŵ����Ƶ�ʱ��弤���У���ÿ��ofdm���ŵ��м�λ�ó���
        %һ�����ݰ������߸����ţ�ǰ�������ŵ�ѭ��ǰ׺��73�����߸����ŵ�ѭ��ǰ׺��74
        train_ofdm=[zeros(1,72),zeros(1,512),1,zeros(1,511),zeros(1,72),zeros(1,512),1,zeros(1,511),zeros(1,72),zeros(1,512),1,zeros(1,511),...
        zeros(1,72),zeros(1,512),1,zeros(1,511),zeros(1,72),zeros(1,512),1,zeros(1,511),zeros(1,72),zeros(1,512),1,zeros(1,511),zeros(1,80),zeros(1,512),1,zeros(1,511),];
       
        %�жϵ�ǰ���ݰ��Ƿ��ǵ�һ�����ݰ�
        if II==1
           pre_train_fere=zeros(U*S,10);
           pre_fere=zeros(U*S,10);
        end     
		 %�ӿ�3���ŵ�ģ�飩��
        [h,delays]=wdy_scm_para(mul_path,T,U,S);%mul_pathΪ�Ӿ�����TΪ֡��
       
        
        %���ŵ������룺pre_train_fereǰ��ĸ���,train_ofdm��temp_train_ofdm�ֱ�Ϊ�����źţ�mul_pathΪ�Ӿ�����NΪ��������X����������
        [channel_out_train,train_fere]=scm1_channel(train_ofdm,pre_train_fere,h,delays,mul_path,U,S);
        
        %���ŵ������룺pre_fereǰ��ĸ���,ofdm_out,temp_ofdm_out�ֱ�Ϊ�����źţ�mul_pathΪ�Ӿ�����NΪ��������X����������
        [scm_out,signal_fere]=scm1_channel(ofdm_out,pre_fere,h,delays,mul_path,U,S);
        
        %�����signal_fere����һ����ĸ���
        pre_train_fere=train_fere;
        pre_fere=signal_fere;
        
        %��ʱ����ŵ�����ֵ�任��Ƶ��
       for i=1:U*S %�����߸������У�ǰ����Ϊ��Ƶ�źţ�����������������ʱ����ͳ�ƣ�ֱ����ȡ��Ҫ�������ź�
		   pre_temp1(i,:)=channel_out_train(i,1:1096);%��һ�����Ŷ�Ӧ��Ƶ���ŵ�����ֵ74:1097
		   pre_temp2(i,:)=channel_out_train(i,1097:2192);%�ڶ������Ŷ�Ӧ��Ƶ���ŵ�����ֵ1171:2194
		   pre_temp3(i,:)=channel_out_train(i,2193:3288);%���������Ŷ�Ӧ��Ƶ���ŵ�����ֵ2268:3291
           pre_temp4(i,:)=channel_out_train(i,3289:4384);%���ĸ����Ŷ�Ӧ��Ƶ���ŵ�����ֵ3365:4388
           pre_temp5(i,:)=channel_out_train(i,4385:5480);%��������Ŷ�Ӧ��Ƶ���ŵ�����ֵ4462:5485
           pre_temp6(i,:)=channel_out_train(i,5481:6576);%���������Ŷ�Ӧ��Ƶ���ŵ�����ֵ5559:6582
           pre_temp7(i,:)=channel_out_train(i,6577:7680);%���߸����Ŷ�Ӧ��Ƶ���ŵ�����ֵ6657:7680
       end
       %Ϊ�õ�׼ȷ���ŵ��弤��Ӧ���������ŵ�����������������,
       for p=1:U*S
           train_out=[pre_temp1(p,1:72),pre_temp1(p,584:1096),pre_temp1(p,73:583),...
                      pre_temp2(p,1:72),pre_temp2(p,584:1096),pre_temp2(p,73:583),...
                      pre_temp3(p,1:72),pre_temp3(p,584:1096),pre_temp3(p,73:583),...
                      pre_temp4(p,1:72),pre_temp4(p,584:1096),pre_temp4(p,73:583),...
                      pre_temp5(p,1:72),pre_temp5(p,584:1096),pre_temp5(p,73:583),...
                      pre_temp6(p,1:72),pre_temp6(p,584:1096),pre_temp6(p,73:583),...
                      pre_temp7(p,1:80),pre_temp7(p,592:1104),pre_temp7(p,81:591)];
           channel_F_domain(p,:)=deofdm(train_out)*32;
       end
       %���߽����ź�(������˹����)��
       awgn_out=myawgn(scm_out,SN);
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
       re_siga=deofdm(awgn_out);%��OFDM�ź�
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
       %���⣬������ֿ��п��ޣ�
       equal_out=equal(re_siga,channel_F_domain);
       
       %�ӿ�4.2�����ģ�飩��
       deqam_out=deqam4(equal_out,channel_F_domain,SN);
       
       %�ӿ�5���⽻֯ģ�飩��
       %deinter_out= de_interleaving( deqam_out ,interleave_table);
        
       %�ӿ�5.1(turbo������ģ��)��
       %P_Rx= turbo_decode(deinter_out,sub_interleaver);
       P_Rx= turbo_decode(deqam_out,sub_interleaver);
       
       %�ӿ�6�������Ǳ�:
        % sigcon(P_Rx);
        
        
       %ͳ�ƴ���ı��������������λ�á����ݰ����������λ��
	   Bit_compare_temp=(P_Rx~=P_Tx);
		%�����
		Bits_Rx_err=Bits_Rx_err+sum(Bit_compare_temp);
        
		%���
		if (sum(Bit_compare_temp))~=0
			Ps_Rx_err=Ps_Rx_err+1;
        end
		%���������е����ݰ���Ŀ�ͱ�����Ŀ
        I
        Ps_Tx=Ps_Tx+1;
		Bits_Tx=Bits_Tx+Lp;
		%���������Ƿ�ﵽ��������
		if Ps_Rx_err>=Pe_max
		    break;
        end
	end
	
	%��¼������ȵ�ķ�����
	%�����ݰ���Ŀ
	Ps_run_s(I)=Ps_Tx;
	%�ܱ�����
	Bits_run_s(I)=Bits_Tx;
	%�������ݰ���Ŀ
	Ps_err_s(I)=Ps_Rx_err;
	%���������Ŀ
	Bits_err_s(I)=Bits_Rx_err;
	%�����
	P_err_rate_s(I)=Ps_Rx_err/Ps_Tx;
	%������
	Bit_err_rate_s(I)=Bits_Rx_err/Bits_Tx;
    save('S1_R1.mat','Bit_err_rate_s','P_err_rate_s','SNs');
end