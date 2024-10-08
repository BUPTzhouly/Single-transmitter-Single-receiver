 %本程序实现 程序二”中“流程控制模块”的功能（参见《链路级仿真软件设计_v2.1》）
%作者：林辉
%修改记录：2005年2月25
%清空工作区中的所有变量
clear;
%参数设定:
%SCM参数
U=1;% 发射天线
S=1;%接收天线
mul_path=6;

%仿真运行的信噪比点数
M=6;
%控制各个信噪比点选择的最大运行包数
N=1;

%各点的信噪比
SNs=[ 0  2  4  6  8  10];
T=7680;
%链路级数据包的大小
Lp=696;%(bits)

%仿真运行的结束条件：
%各点运行的最大总数据包数量
Ps=[5000,5000,5000,5000,5000,5000,5000,5000,5000]*2;
%各点运行的最大误包数
Pes=[100,100,100,100,100,100,100,100,100];
%比特速率
Tb=0.5*10^(-3)/Lp;%(s)
%采样间隔
Ts=1*10^(-3)/15/1024;%(S)
%符号时间
T_symbo=0.5*10^(-3)/7;
%%%%%%%%%%%%%%%设置结果存储参数:
%各信噪比点的总数据包数目
Ps_run_s=zeros(N,M);
%各信噪比点的总比特数
Bits_run_s=zeros(N,M);
%各信噪比点的总错误数据包数目
Ps_err_s=zeros(N,M);
%各信噪比点的总错误比特数
Bits_err_s=zeros(N,M);
%各信噪比点的误包率
P_err_rate_s=zeros(N,M);
%各信噪比点的误码率
Bit_err_rate_s=zeros(N,M);

%产生交织码表
interleave_table=interleav_matrix(3*(Lp+4));

%设置产生随机数的种子
rand('state',sum(100*clock));
randn('state',sum(7*clock));

%使用循环，以运行M个信噪比点
for I=1:M
    I
    tic
	%初始化:
	%发送的数据包数量
 	Ps_Tx=0;
	%错误的数据包数量
	Ps_Rx_err=0;
	%发送的总比特数
	Bits_Tx=0;
	%错误的比特数目
	Bits_Rx_err=0;
	%该仿真点的信噪比
	SN=SNs(I);
	%该仿真点最多运行的数据包数目
	P_max=Ps(I);
	%该仿真点最多运行的错误数据包数目
	Pe_max=Pes(I);
	%使用循环，以运行最多P_max个数据包
	for II=1:P_max
        II
		%产生1个长为Lp的数据包:
		P_Tx=(rand(1,Lp)>0.5);
        
        %接口1（turbo编码）       
	    sub_interleaver=sub_interleave(Lp);
        turbo_out=turbo_encode_1_3(P_Tx,sub_interleaver);
        
        %接口1.1（交织器）：
        %interleav_out = interleaving(turbo_out ,interleave_table);
		
		%接口2（4QAM 调制）
		%Qam_out=qam4(interleav_out);
        Qam_out=qam4(turbo_out);
                   
        %接口2.1（OFDM模块）：
        ofdm_out=ofdm(Qam_out);
   
        %初始化用于信道估计的时域冲激序列，在每个ofdm符号的中间位置抽样
        %一个数据包里有七个符号，前六个符号的循环前缀是73，第七个符号的循环前缀是74
        train_ofdm=[zeros(1,72),zeros(1,512),1,zeros(1,511),zeros(1,72),zeros(1,512),1,zeros(1,511),zeros(1,72),zeros(1,512),1,zeros(1,511),...
        zeros(1,72),zeros(1,512),1,zeros(1,511),zeros(1,72),zeros(1,512),1,zeros(1,511),zeros(1,72),zeros(1,512),1,zeros(1,511),zeros(1,80),zeros(1,512),1,zeros(1,511),];
       
        %判断当前数据包是否是第一个数据包
        if II==1
           pre_train_fere=zeros(U*S,10);
           pre_fere=zeros(U*S,10);
        end     
		 %接口3（信道模块）：
        [h,delays]=wdy_scm_para(mul_path,T,U,S);%mul_path为子径数，T为帧长
       
        
        %过信道，输入：pre_train_fere前块的干扰,train_ofdm、temp_train_ofdm分别为输入信号，mul_path为子径数，N为发射天线X接收天线数
        [channel_out_train,train_fere]=scm1_channel(train_ofdm,pre_train_fere,h,delays,mul_path,U,S);
        
        %过信道，输入：pre_fere前块的干扰,ofdm_out,temp_ofdm_out分别为输入信号，mul_path为子径数，N为发射天线X接收天线数
        [scm_out,signal_fere]=scm1_channel(ofdm_out,pre_fere,h,delays,mul_path,U,S);
        
        %输出：signal_fere对下一个块的干扰
        pre_train_fere=train_fere;
        pre_fere=signal_fere;
        
        %将时域的信道估计值变换到频域
       for i=1:U*S %由于七个符号中，前两个为导频信号，所以在这里我们暂时不作统计，直接提取所要的有用信号
		   pre_temp1(i,:)=channel_out_train(i,1:1096);%第一个符号对应的频域信道函数值74:1097
		   pre_temp2(i,:)=channel_out_train(i,1097:2192);%第二个符号对应的频域信道函数值1171:2194
		   pre_temp3(i,:)=channel_out_train(i,2193:3288);%第三个符号对应的频域信道函数值2268:3291
           pre_temp4(i,:)=channel_out_train(i,3289:4384);%第四个符号对应的频域信道函数值3365:4388
           pre_temp5(i,:)=channel_out_train(i,4385:5480);%第五个符号对应的频域信道函数值4462:5485
           pre_temp6(i,:)=channel_out_train(i,5481:6576);%第六个符号对应的频域信道函数值5559:6582
           pre_temp7(i,:)=channel_out_train(i,6577:7680);%第七个符号对应的频域信道函数值6657:7680
       end
       %为得到准确的信道冲激响应函数而对信道函数进行重新排序,
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
       %天线接收信号(经过高斯白噪)：
       awgn_out=myawgn(scm_out,SN);
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
       re_siga=deofdm(awgn_out);%解OFDM信号
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
       %均衡，这个部分可有可无：
       equal_out=equal(re_siga,channel_F_domain);
       
       %接口4.2（解调模块）：
       deqam_out=deqam4(equal_out,channel_F_domain,SN);
       
       %接口5（解交织模块）：
       %deinter_out= de_interleaving( deqam_out ,interleave_table);
        
       %接口5.1(turbo码译码模块)：
       %P_Rx= turbo_decode(deinter_out,sub_interleaver);
       P_Rx= turbo_decode(deqam_out,sub_interleaver);
       
       %接口6（虚拟仪表）:
        % sigcon(P_Rx);
        
        
       %统计错误的比特数、错误比特位置、数据包数、错误包位置
	   Bit_compare_temp=(P_Rx~=P_Tx);
		%误比特
		Bits_Rx_err=Bits_Rx_err+sum(Bit_compare_temp);
        
		%误包
		if (sum(Bit_compare_temp))~=0
			Ps_Rx_err=Ps_Rx_err+1;
        end
		%更新已运行的数据包数目和比特数目
        I
        Ps_Tx=Ps_Tx+1;
		Bits_Tx=Bits_Tx+Lp;
		%检查误包率是否达到结束条件
		if Ps_Rx_err>=Pe_max
		    break;
        end
	end
	
	%记录该信噪比点的仿真结果
	%总数据包数目
	Ps_run_s(I)=Ps_Tx;
	%总比特数
	Bits_run_s(I)=Bits_Tx;
	%错误数据包数目
	Ps_err_s(I)=Ps_Rx_err;
	%错误比特数目
	Bits_err_s(I)=Bits_Rx_err;
	%误包率
	P_err_rate_s(I)=Ps_Rx_err/Ps_Tx;
	%误码率
	Bit_err_rate_s(I)=Bits_Rx_err/Bits_Tx;
    save('S1_R1.mat','Bit_err_rate_s','P_err_rate_s','SNs');
end