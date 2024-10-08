%----------------------------------------------------------------------------------------------------
%perform the turbo coder/linhui 2001.2,7
%----------------------------------------------------------------------------------------------------
function post_code=turbo_encode_1_3(before_code,interleave_order)
%----------------------------------------------------------------------------------------------------
code_len=length(before_code);
%----------------------------------------------------------------------------------------------------
%get the length of the data input to the coder
%----------------------------------------------------------------------------------------------------
coder1=zeros(3,1);
Z1=zeros(1,code_len);
for code_i=1:code_len
   Z1(code_i)=xor(before_code(code_i),xor(coder1(1),coder1(2)));
   temp=xor(before_code(code_i),xor(coder1(3),coder1(2)));
   coder1(3)=coder1(2);
   coder1(2)=coder1(1);
   coder1(1)=temp;
end
%----------------------------------------------------------------------------------------------------
%get the second branch "Z1" of the first coder
%-----------------------------------------------------------------------------------------------------
s_order=interleave_order;
%----------------------------------------------------------------------------------------------------
%call the function "interleave" to get the permuted order
%----------------------------------------------------------------------------------------------------
coder2=zeros(3,1);
Z2=zeros(1,code_len);
for code_i=1:code_len
   Z2(code_i)=xor(before_code(1+s_order(code_i)),xor(coder2(1),coder2(2)));
   temp=xor(before_code(1+s_order(code_i)),xor(coder2(3),coder2(2)));
	coder2(3)=coder2(2);
   coder2(2)=coder2(1);
   coder2(1)=temp;
end
%----------------------------------------------------------------------------------------------------
%get the second branch "Z2" of the second coder
%----------------------------------------------------------------------------------------------------
tail1=zeros(3,1);
tail2=zeros(3,1);
tail3=zeros(3,1);
tail4=zeros(3,1);
for tail_i=1:3
   tail1(tail_i)=xor(coder1(3),coder1(2));
   tail2(tail_i)=xor(coder1(3),coder1(1));
   coder1(3)=coder1(2);
   coder1(2)=coder1(1);
   coder1(1)=0;
   tail3(tail_i)=xor(coder2(3),coder2(2));
   tail4(tail_i)=xor(coder2(3),coder2(1));
   coder2(3)=coder2(2);
   coder2(2)=coder2(1);
   coder2(1)=0;
end
%----------------------------------------------------------------------------------------------------
%get 4x3 tail bits of the trelis terminal
%----------------------------------------------------------------------------------------------------
post_len=code_len*3+12;
post_code=zeros(1,post_len);
for post_i=1:code_len
   post_code(post_i*3-2)=before_code(post_i);
   post_code(post_i*3-1)=Z1(post_i);
   post_code(post_i*3)=Z2(post_i);
end
for post_ii=1:3
   post_code(code_len*3+post_ii*2-1)=tail1(post_ii);
   post_code(code_len*3+post_ii*2)=tail2(post_ii);
   post_code(code_len*3+post_ii*2+5)=tail3(post_ii);
   post_code(code_len*3+post_ii*2+6)=tail4(post_ii);
end
%-----------------------------------------------------------------------------------------------------
%put branches and tails together to get the encode result
%-----------------------------------------------------------------------------------------------------


