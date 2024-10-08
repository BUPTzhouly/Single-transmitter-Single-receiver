function[H,delay_out]=wdy_scm_para(N,T,U,S)
%这个生成scm信道快衰的函数，根据输入的路径数和抽样点数来生成信道函数、延迟点数
%用户可以通过修改N来调整生成的路径数，修改T来调整信道函数的抽样点数（长度）
%N    = 6;             % number of paths,这里作为输入变量，来控制生成的路径数
%T    = 7680;          % number of time samples，这里作为变量，来控制生成的抽样点数
%S    =2;               % number of receiving antennas
%U    = 2;              % number of transmitting antennas
M    = 20;             % number of subpaths
speed_of_light=2.99792458e8;
CenterFrequency=2e9;          %载波频率
ThetaBs=360*(rand-0.5);       %U(-180,180) degrees, U denotes uniform pdf
ThetaMs=360*(rand-0.5);       %U(-180,180) degrees,即MS的角度
MsDirection=360*(rand-0.5);   %U(-180,180) degrees with respect to broadside
MsVelocity=70/3.6;
LM = 1:M;
MsElementPosition=[0,0.5];    %天线的参考距离
BsElementPosition=[0,0.5];    %天线的参考距离
wavelength=speed_of_light/CenterFrequency;
k_CONST = 2*pi/wavelength;
bulkpar=struct(  'Scenario','urban_macro',...
                 'BsUrbanMacroAS','eight',...            % choices: 'eight' and 'fifteen'. 
                 'DelaySamplingInterval',1*10^(-3)/15/1024/16); %  *5e-9     
Scenario=bulkpar.Scenario;        
switch lower(Scenario)
    case {'suburban_macro','urban_macro'}
         basepar=macro(N,M,ThetaBs,ThetaMs,bulkpar);  %求得功率分配、延迟时间、AOD、AOA及子径的大小
     case {'urban_micro'}        
        basepar=micro(N,M,ThetaBs,ThetaMs,bulkpar);   %求得功率分配、延迟时间、AOD、AOA及子径的大小
end     
delays=basepar.delays;
path_powers_all = basepar.path_powers;
%path_powers_all=ones(1,N);
delta_t=bulkpar.DelaySamplingInterval*16;
%Delay=10^(-9)*[0,310,710,1090,1730,2510];
delay_out=floor(delays/delta_t); %输出六条径延迟的点数（一条径时为0）
max_delay=max(delay_out);        %最大延迟点数
%max_delay=max(delay_out);
H = zeros(U,S,N,T+max_delay);    %快衰信道函数
%t = repmat(delta_t,1,T+max_delay).*repmat([0:T+max_delay-1],1); %抽样时间（点数）
t = repmat(delta_t,1,T+max_delay).*([0:T+max_delay-1]); 
for u = 1:U % cycles (MS) antennas
    du = MsElementPosition(u)*wavelength;
    for s = 1:S % cycles Tx (BS) atennas
        ds = BsElementPosition(s)*wavelength;
        for n = 1:N % cycles paths
            LM_index = 0; % 
            temp                 = zeros(M,T+max_delay);   % oversized, just to keep it always the same size
             for m=1:M % cycles subpaths
                 LM_index = LM_index+1;
                 temp(m,:) = exp(j * (k_CONST * ds * sin((basepar.aods(n,LM(LM_index)))*pi/180) +(basepar.subpath_phases(n,LM(LM_index))*pi/180)+k_CONST * du * sin((basepar.aoas(n,LM(LM_index)))*pi/180)...
                 )) *exp(j * k_CONST * MsVelocity * cos((basepar.aoas(n,LM(LM_index)) - MsDirection)*pi/180) * (t));
             end % subpaths
               H(u,s,n,:) =  sqrt(path_powers_all(n) / M) * sum(temp,1);
         end % paths 
    end % Tx antennas 
end % Rx antennas 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%生成macro环境下的功率分配、延迟时间、AOD、AOA及子径的大小
function basepar=macro(N,M,ThetaBs,ThetaMs,bulkpar)
Scenario              = bulkpar.Scenario;
DelaySamplingInterval = bulkpar.DelaySamplingInterval;  %对延迟时间进行量化的单位时间
Bsq=[0 0 0; 0 0 0; 0 0 0.7071];
% pre-computed value: C=sqrtm(A-B)
C = [0.8997 0.1926 -0.3917; 0.1926 0.8997 -0.3917; -0.3917 -0.3917  0.4395];
switch lower(Scenario)% general environment parameters for suburban macro [1, Table 5.1]
    case {'suburban_macro'}
        r_as      = 1.2;        
        r_ds      = 1.4;        
        sigma_rnd = 3;    % per-path shadowing std in dB, needed in step 5
        mu_as      =  0.69 ;
        epsilon_as =  0.13 ;
        mu_ds      = -6.80 ;
        epsilon_ds =  0.288;
        sigma_sf_ave   =  8    ; % in dB
         % generate alphas, betas and gammas for all sites       
        abc = C*randn(3,1) + Bsq*randn(3,1); 
        %step3:determine as ds sf
        sigma_ds = 10.^(epsilon_ds*abc(1) + mu_ds);       
        sigma_as = 10.^(epsilon_as*abc(2) + mu_as);
        sigma_sf = 10.^(0.1*sigma_sf_ave*abc(3));
    case {'urban_macro'}
        r_as      = 1.3;        
        r_ds      = 1.7;        
        sigma_rnd = 3;    % per-path shadowing std in dB, needed in step 5
        if strcmpi(bulkpar.BsUrbanMacroAS,'fifteen')% general environment parameters for urban macro [1, Table 5.1]
            mu_as      =  1.18 ;
            epsilon_as =  0.210;
        else     % Note: 8 degree angle spread is set automatically if no match to 'fifteen'
            mu_as      =  0.810;
            epsilon_as =  0.34 ;
        end    
            mu_ds      = -6.18 ;         
            sigma_sf_ave   =  8    ; % in dB            
            epsilon_ds =  0.18 ;            
        % generate alphas, betas and gammas for all sites
        abc = C*randn(3,1) + Bsq*randn(3,1); 
        %step3:determine as ds sf
        sigma_ds = 10.^(epsilon_ds*abc(1) + mu_ds);       
        sigma_as = 10.^(epsilon_as*abc(2) + mu_as);
        sigma_sf = 10.^(0.1*sigma_sf_ave*abc(3));
end
% step 4: generate delays in a (NumLinks x N) matrix
taus         = sort(-r_ds*sigma_ds*log(rand(1,N)),2);
taus_sorted  = taus - taus(1);       % normalize min. delay to zero
taus_rounded=DelaySamplingInterval*floor( taus_sorted/DelaySamplingInterval + 0.5);
% step 5: determine random average powers in a (NumLinks x N) matrix
ksi    = randn(1,N)*sigma_rnd;           % per-path shadowing
Pprime = exp((1-r_ds)/r_ds*taus_sorted./sigma_ds ).*10.^(-ksi/10);
P      = Pprime./sum(Pprime);     % power normalization
   % step 6: determine AoDs
sigma_aod     = r_as*sigma_as;                              % AoD angle spreads for all users
deltas        = abs(randn(1,N).*repmat(sigma_aod,1,N));
deltas_sorted = sort(deltas,2);
delta_aod     = sign(rand(1,N)-0.5).*deltas_sorted;  % a (NumLinks x N) matrix of path AoDs
delta_aod     = repmat(delta_aod,M,1);                    % a (M x (NumLinks*N)) matrix
     % The phases are computed in step 13
aod_2deg     = [0.0894 0.2826 0.4984 0.7431 1.0257 1.3594 1.7688 2.2961 3.0389 4.3101];     % [1, Table 5.2]
delta_nm_aod = [aod_2deg; -aod_2deg];
delta_nm_aod = repmat(delta_nm_aod(:),1,N);  % a (M x (NumLinks*N)) matrix
    % step 9: determine the AoAs
sigma_aoa = 104.12*(1-exp(-0.2175*abs(10*log10(P))));
delta_aoa = randn(1,N).*sigma_aoa;     % a (NumLinks x N) matrix of path AoAs
delta_aoa = repmat(delta_aoa,M,1);          % a (M x (NumLinks*N)) matrix
 % step 10: determine the offset AoAs at the MS
aoa_35deg    = [1.5679 4.9447 8.7224 13.0045 17.9492 23.7899 30.9538 40.1824 53.1816 75.4274];      % [1, Table 5.2]
delta_nm_aoa = [aoa_35deg; -aoa_35deg];
delta_nm_aoa = repmat(delta_nm_aoa(:),1,N); % a (M x N*NumLinks) matrix
 % step 11: pair AoA subpaths randomly with AoD subpaths (within a path)
[dummy h]           = sort(rand(M,N),1);       % create N*NumLinks random permutations of integers [1:M]
inds                = h+repmat([1:M:M*N],M,1)-1;
delta_nm_aoa_paired = delta_nm_aoa(inds);    % random permutation of columns, a (M x N*NumLinks) matrix
% Values of theta_nm_aoa and theta_nm_aod may be outside (-180,180).
theta_nm_aoa=prin_value(ThetaMs + delta_aoa+delta_nm_aoa_paired)';
theta_nm_aod=prin_value(ThetaBs+delta_aod+delta_nm_aod)'; 
phi= 360*rand(N,M);        % random phases for all users
basepar=struct( 'delays',taus_rounded,...
                'path_powers',P,...             % before: 'subpath_powers',Psub,...
                'aods',theta_nm_aod,...         % in degrees
                'aoas',theta_nm_aoa,...         % in degrees
                'subpath_phases',phi);           % in degrees
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%生成micro环境下的功率分配、延迟时间、AOD、AOA及子径的大小
function basepar=micro(N,M,ThetaBs,ThetaMs,bulkpar)
DelaySamplingInterval   =bulkpar.DelaySamplingInterval;
% general environment parameters for urban micro [1, Table 5.1]
max_ds      =1.2e-6;  % maximum excess delay in seconds
max_aod     =40;     % maximum AoD angle in degrees
sigma_rnd   =3;    % per-path shadowing std in dB, needed in step 6
taus=sort(max_ds*rand(1,N),2);  % sorted path delays for all users
taus=taus-taus(1);              % normalize min. delay to zero
taus_rounded=DelaySamplingInterval*floor( taus/DelaySamplingInterval + 0.5);
% step 6: determine random average powers in a (NumLinks x N) matrix
z=randn(1,N)*sigma_rnd;       % per-path shadowing
Pprime=10.^( -(taus/1e-6 + z/10) );
P=Pprime./sum(Pprime);    % power normalization
   % step 7: determine AoDs
delta_aod=2*max_aod*(rand(1,N)-0.5);
delta_aod=repmat(delta_aod,M,1);      % a (M x (NumLinks*N)) matrix
% step 9: determine  offset AoDs at the BS
aod_5deg=[0.2236 0.7064 1.2461 1.8578 2.5642 3.3986 4.4220 5.7403 7.5974 10.7753]; % [1, Table 5.2]
delta_nm_aod = [aod_5deg; -aod_5deg];
delta_nm_aod=repmat(delta_nm_aod(:),1,N);  % a (NumLinks x N) matrix
% step 10: determine the AoAs
sigma_aoa = 104.12*(1-exp(-0.265*abs(10*log10(P))));
delta_aoa = randn(1,N).*sigma_aoa;     % a (NumLinks x N) matrix of path AoAs
delta_aoa = repmat(delta_aoa,M,1);          % a (M x (NumLinks*N)) matrix
% step 11: determine the offset AoAs at the MS
aoa_35deg       =[1.5679 4.9447 8.7224 13.0045 17.9492 23.7899 30.9538 40.1824 53.1816 75.4274];      % [1, Table 5.2]
delta_nm_aoa    = [aoa_35deg; -aoa_35deg];
delta_nm_aoa    =repmat(delta_nm_aoa(:),1,N); % a (M x N*NumLinks) matrix
    % step 12: pair AoA subpaths randomly with AoD subpaths (within a path)
[dummy h]           = sort(rand(M,N),1);       % create N*NumLinks random permutations of integers [1:M]
inds                =h+repmat([1:M:M*N],M,1)-1;
delta_nm_aoa_paired =delta_nm_aoa(inds);    % random permutation of columns, a (M x N*NumLinks) matrix
% Values of theta_nm_aoa and theta_nm_aod may be outside (-180,180).
theta_nm_aoa=prin_value(ThetaMs+delta_aoa+delta_nm_aoa_paired)';          
theta_nm_aod=prin_value(ThetaBs+delta_aod+delta_nm_aod)'; 
phi= 360*rand(N,M);        % random phases for all users
basepar=struct( 'delays',taus_rounded,...
                'path_powers',P,...             % before: 'subpath_powers',Psub,...
                'aods',theta_nm_aod,...
                'aoas',theta_nm_aoa,...
                'subpath_phases',phi); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A function that maps inputs from (-inf,inf) to (-180,180)
function y=prin_value(x)
y=mod(x,360);
y=y-360*floor(y/180);