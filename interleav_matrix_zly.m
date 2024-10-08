function interleaver = interleav_matrix_zly(N)
    % 生成Turbo交织码表
    % 输入N是Turbo码内部交织器的位数

    % 检查N是否为正整数
    if ~isscalar(N) || N < 1 || floor(N) ~= N
        error('N必须是一个正整数。');
    end

    % 确定参数
    H_prime = 30;
    R = [0, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1];
    Q_prime = floor(N / H_prime);
    Q = Q_prime + R(mod(N, H_prime) + 1);
    
    % 初始化交织器索引
    interleaver = zeros(1, N);
    
    % 生成交织码表
    for i = 0:N-1
        f1 = 1; % 这些值(f1, f2)根据标准进行选择
        f2 = 1; % 通常是固定的或者基于某些条件的选择值
        % 根据3GPP TS 36.212计算交织索引
        % 注意: 实际应用中f1和f2的值应根据具体标准确定
        interleaver(i + 1) = mod(f1 * i + f2 * i^2, N) + 1;
    end
    
    % 返回交织码表
end