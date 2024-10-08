function[final_sig,signal_fere]=scm_channel(signal,pre_interfere,H,delays,mul_path,U,S)
%mul_path=6;   
%S    =2;               % number of receiving antennas
%U    = 2;              % number of transmitting antennas
%Ts=1*10^(-3)/15/1024;
%Delay=10^(-9)*[0,310,710,1090,1730,2510];
%relative_power=[0,-1,-9,-10,-15,-20];
%Am=sqrt(10.^(0.1*relative_power));
%Am=Am./sqrt(sum(Am.^2));
%m_max=floor(max(Delay)/Ts);
N=length(signal);

%H=H./abs(H);
pre_seq_max=size(pre_interfere,2);
now_seq=floor(delays);
now_seq_max=max(now_seq);
inte_signal=zeros(U*S,N+now_seq_max);
%h=zeros(U*S,mul_path,N+now_seq_max);
for u=1:U%得到对应信道的上的信道函数
    for s=1:S
        for b=1:mul_path
            delay_add=delays(b)+1;
            %h((u-1)*2+s,b,:)=H(u,s,b,:);%*Am(i);
            temp_h(1,:)=H(u,s,b,:);%h((u-1)*2+s,b,:);
            inte_signal((u-1)*2+s,delay_add:delay_add+N-1)=inte_signal((u-1)*2+s,delay_add:delay_add+N-1)+signal.*temp_h(1,delay_add:delay_add+N-1);
        end
    end
end
    %h(1,b,:)=H(1,1,b,:);%*Am(i);
    %h(2,b,:)=H(1,2,b,:);%*Am(i);
    %h(3,b,:)=H(2,1,b,:);%*Am(i);
    %h(4,b,:)=H(2,2,b,:);%*Am(i);

%for n=1:mul_path%多径叠加
%    delay_add=delays(n)+1;
%    for m=1:U*S
%        aa(m,:)=h(m,n,:);
%    end
%    inte_signal(1,delay_add:delay_add+N-1)=inte_signal(1,delay_add:delay_add+N-1)+signal.*aa(1,delay_add:delay_add+N-1);
    %inte_signal(2,delay_add:delay_add+N-1)=inte_signal(2,delay_add:delay_add+N-1)+signal.*aa(2,delay_add:delay_add+N-1);
    %inte_signal(3,delay_add:delay_add+N-1)=inte_signal(3,delay_add:delay_add+N-1)+temp_signal.*aa(3,delay_add:delay_add+N-1);
    %inte_signal(4,delay_add:delay_add+N-1)=inte_signal(4,delay_add:delay_add+N-1)+temp_signal.*aa(4,delay_add:delay_add+N-1);
    %end
for j=1:U*S%上个块对本块的干扰
    inte_signal(j,1:pre_seq_max)=inte_signal(j,1:pre_seq_max)+pre_interfere(j,:);
end
%输出信号
final_sig=inte_signal(:,1:N);
signal_fere=inte_signal(:,N+1:N+now_seq_max);