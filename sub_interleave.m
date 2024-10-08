%interleave for tubor encode and decode /lin hui 2001.2.6 
%-----------------------------------------------------------------------------------------------------
function post_interleave=sub_interleave(interleave_length)
%-----------------------------------------------------------------------------------------------------
%interleave_length is the length of data before tubor code
%-----------------------------------------------------------------------------------------------------
before_interleave=zeros(interleave_length,1);
for before_i=1:interleave_length
   before_interleave(before_i)=before_i-1;
end
%-----------------------------------------------------------------------------------------------------
%make a matrix as the order before interleave 
%-----------------------------------------------------------------------------------------------------
%post_inter is a same size matrix as before_inter such as "nx1". it is the order after interleave 
%-----------------------------------------------------------------------------------------------------
%table
%-----------------------------------------------------------------------------------------------------
p_table=[7;11;13;17;19;23;29;31;37;41;43;...
      47;53;59;61;67;71;73;79;83;89;97;...
      101;103;107;109;113;127;131;137;139;149;151;...
      157;163;167;173;179;181;191;193;197;199;211;...
      223;227;229;233;239;241;251;257];
v_table=[3;2;2;3;2;5;2;3;2;6;3;...
      5;2;2;2;2;7;5;3;2;3;5;...
      2;5;2;6;3;3;2;3;2;2;6;...
      5;2;5;2;2;2;19;5;2;3;2;...
   	3;2;6;3;7;7;6;3];
pat1=[19;9;14;4;0;2;5;7;12;18;10;8;13;17;3;1;16;6;15;11];
pat2=[19;9;14;4;0;2;5;7;12;18;16;13;17;15;3;1;6;11;8;10];
pat3=[9;8;7;6;5;4;3;2;1;0];
pat4=[4;3;2;1;0];
prime_num=[7;11;13;17;19;23;29;31;37;41;43;47;53;59;61;67;71;...
      	73;79;83;89;97;101;103;107;109;113;127;131;137;139;149;...
      	151;157;163;167;173;179;181;191;193;197;199;211;223;227;...
   		229;233;239;241;251;257];

%-----------------------------------------------------------------------------------------------------
%some constant to be used
%-----------------------------------------------------------------------------------------------------
K=size(before_interleave,1);
%-----------------------------------------------------------------------------------------------------
%return the rows of the input matrix
%-----------------------------------------------------------------------------------------------------
if K>=40&K<=159 
   R=5;
elseif (K>=160&K<=200)|(K>=481&K<=530)
   R=10;
else
   R=20;
end
%-----------------------------------------------------------------------------------------------------
%get the rows of the interleave matrix;
%-----------------------------------------------------------------------------------------------------
if K>=481&K<=530
   p=53;
   v=2;
   C=p;
else
   for p_i=1:52
      if p_table(p_i)+1-K/R>=0
         break;
      end      
   end
   p=p_table(p_i);
   v=v_table(p_i);
   if p-K/R>=0 
      if p-1-K/R>=0
         C=p-1;
      else
         C=p;
      end
   else   
      C=p+1;
   end
end
%-----------------------------------------------------------------------------------------------------
%get the column of the interleave matrix and 'p'\'v';
%-----------------------------------------------------------------------------------------------------
interleave_matrix_step1=-ones(R,C);
for R_i=1:R
   for C_i=1:C
      K_i=(R_i-1)*C+C_i;
      if K_i<=K
         interleave_matrix_step1(R_i,C_i)=before_interleave(K_i);
      end
   end
end
%-----------------------------------------------------------------------------------------------------
%put source into interleave matrix
%-----------------------------------------------------------------------------------------------------
s=ones(p-1,1);
for s_i=2:(p-1)
   s(s_i)=mod(v*s(s_i-1),p);
end;
%-----------------------------------------------------------------------------------------------------
%get the intra-row permutation base sequence : 's'
%-----------------------------------------------------------------------------------------------------
q=ones(R,1);
q_i=2;
prime_i=1;
while q_i<R+1
   if mod(p-1,prime_num(prime_i))~=0&prime_num(prime_i)~=p-1
      q(q_i)=prime_num(prime_i);
      q_i=q_i+1;
   end
   prime_i=prime_i+1;
end      
%-----------------------------------------------------------------------------------------------------
%get the base sequent 'q'
%-----------------------------------------------------------------------------------------------------
if K>=40&K<=159
   T=pat4;
elseif K<=200
   T=pat3;
elseif K<=480 
   T=pat1;
elseif K<=530
   T=pat3;
elseif K<=2280
   T=pat1;
elseif K<=2480
   T=pat2;
elseif K<=3160
   T=pat1;
elseif K<=3210
   T=pat2;
else
   T=pat1;
end
r=zeros(R,1);
for r_i=1:R
   r(T(r_i)+1)=q(r_i);
end
%-----------------------------------------------------------------------------------------------------
%permute 'q' to get 'r' 
%-----------------------------------------------------------------------------------------------------
interleave_matrix_step2=zeros(R,C);
if C==p 
   for U_j=1:R
      for U_i=1:p-1
         interleave_matrix_step2(U_j,U_i)=interleave_matrix_step1(U_j,1+s(1+mod((U_i-1)*r(U_j),p-1)));
      end
      interleave_matrix_step2(U_j,p)=interleave_matrix_step1(U_j,1);
   end 
elseif C==p+1
   for U_j=1:R
      if U_j==R&K==C*R
         interleave_matrix_step2(U_j,1)=interleave_matrix_step1(U_j,p+1);
         interleave_matrix_step2(U_j,p+1)=interleave_matrix_step1(U_j,1+s(1));
		else
         interleave_matrix_step2(U_j,1)=interleave_matrix_step1(U_j,1+s(1));
         interleave_matrix_step2(U_j,p+1)=interleave_matrix_step1(U_j,p+1);
      end   
      for U_i=2:p-1
         interleave_matrix_step2(U_j,U_i)=interleave_matrix_step1(U_j,1+s(1+mod((U_i-1)*r(U_j),p-1)));
      end
      interleave_matrix_step2(U_j,p)=interleave_matrix_step1(U_j,1);
   end 
else
   for U_j=1:R
   	for U_i=1:p-1
      	   interleave_matrix_step2(U_j,U_i)=interleave_matrix_step1(U_j,1+s(1+mod((U_i-1)*r(U_j),p-1))-1);
   	end     
   end
end
%-----------------------------------------------------------------------------------------------------
%perform the intra row permutation
%------------------------------------------------------------------------------------------------------
interleave_matrix_step3=zeros(R,C);
for T_i=1:R
   interleave_matrix_step3(T_i,:)=interleave_matrix_step2(T(T_i)+1,:);
end 
%------------------------------------------------------------------------------------------------------
%perform the inter row permutation
%------------------------------------------------------------------------------------------------------
post_ii=1;
for post_i=1:R*C
   if interleave_matrix_step3(post_i)~=-1
      post_interleave(post_ii)=interleave_matrix_step3(post_i);
      post_ii=post_ii+1;
   end
end
post_interleave=post_interleave';
%-----------------------------------------------------------------------------------------------------
%delete unneeded unit and get the result 
%-----------------------------------------------------------------------------------------------------
   
   




       
      
      
      

      
 
