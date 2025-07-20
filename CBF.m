%% CBF空間スペクトル重ね合わせ実験
% サイドローブを推定値と判定してしまうケースを可視化
% Issue #42対応

rng('default')
clear; clc; close all;

theta_0 = deg2rad(0);        % 真のDOA [rad]
SNR_dB = -5;                  % SNR [dB]
N = 4;                       % スナップショット数
K = 8;                       % アンテナ素子数
ant2ant = 1/2;               % アンテナ間隔 (λ/2)
num_trials = 100;            % 試行回数

% アンテナアレイ設定
z = (0:K-1).';               % 素子位置インデックス
d = z*ant2ant;               % 素子位置 [λ単位]
a = @(theta) exp(-1i*2*pi*d*sin(theta));  % ステアリングベクトル
a0 = a(theta_0);             % 目標方向のステアリングベクトル

% 電力設定
Pn = 0;                      % ノイズ電力 [dB]
Ps = SNR_dB + Pn;            % 信号電力 [dB]

% 角度走査設定
scan_theta = linspace(-pi/2, pi/2, 2^9);  % 角度走査範囲 [rad]
scan_a = a(scan_theta);      % 走査用ステアリングベクトル


% 理論的なアレーファクタ（真のDOA基準）
AF = zeros(size(scan_theta));
for i = 1:length(scan_theta)
    theta_i = scan_theta(i);
    u = 2*pi*ant2ant * (sin(theta_i) - sin(theta_0));
    if abs(u) < 1e-10
        AF(i) = K;  % sinc関数のu=0での値
    else
        AF(i) = abs(sin(K * u / 2) / sin(u / 2));
    end
end
AF_dB = 20 * log10(AF / max(AF));

P_CBF_all = zeros(num_trials, length(scan_theta));

for trial = 1:num_trials
    % 受信信号生成
    s = wgn(1, N, Ps, 'complex');      % 送信信号
    V = wgn(K, N, Pn, 'complex');      % ノイズ
    X = a0 * s + V;                    % 受信信号
    R_hat = (X * X') / N;              % サンプル相関行列
    
    % CBF空間スペクトル計算
    tmp = R_hat * scan_a;
    P_CBF = abs(sum(conj(scan_a) .* tmp, 1)).';
    
    % 正規化（dB表示）
    P_CBF_dB = 10 * log10(P_CBF / max(P_CBF));
    P_CBF_all(trial, :) = P_CBF_dB;
end

%% プロット
figure('Position', [100, 100, 800, 600]);

hold on;

yyaxis left;
% 100試行の空間スペクトルを薄い青で重ね合わせ
for trial = 1:num_trials
    plot(rad2deg(scan_theta), P_CBF_all(trial, :),'-','Color', [0.7, 0.9, 1.0], 'LineWidth', 0.5);
end
hold off;

% グラフ設定
ylabel('CBF空間スペクトル [dB]');
ylim([-20, 0]);

% 右軸（アレーファクタ）
yyaxis right;

% アレーファクタを赤い太線で重ね合わせ
plot(rad2deg(scan_theta), AF_dB, 'r-', 'LineWidth', 3, 'DisplayName', 'Array Factor');
hold off;

% 右軸の設定
ylabel('アレーファクタ [dB]');
ylim([-80, 0]);

% 共通設定
xlabel('角度 [degree]');
title(['CBF空間スペクトル重ね合わせ (N=' num2str(N) ', K=' num2str(K) ', SNR=' num2str(SNR_dB) 'dB, ' num2str(num_trials) '試行)']);
grid on;
xlim([-90, 90]);

% 真のDOAにマーカー
xline(rad2deg(theta_0), 'k--', 'LineWidth', 2, 'DisplayName', ['True DOA: ' num2str(rad2deg(theta_0)) '°']);

% 凡例の設定
%legend('CBF Spectrum (100 trials)', 'Array Factor', ['True DOA: ' num2str(rad2deg(theta_0)) '°'], 'Location', 'northeast');

%% 結果の保存
%条件設定
conditions = struct();
conditions.theta_0_deg = rad2deg(theta_0);     % 真のDOA [degree]
conditions.SNR_dB = SNR_dB;                    % SNR [dB]
conditions.N = N;                              % スナップショット数
conditions.K = K;                              % アンテナ素子数
conditions.ant2ant = ant2ant;                  % アンテナ間隔 [λ/2]
conditions.num_trials = num_trials;            % 試行回数
conditions.Pn = Pn;                            % ノイズ電力 [dB]
conditions.Ps = Ps;                            % 信号電力 [dB]
conditions.scan_theta_deg = rad2deg(scan_theta); % 角度走査範囲 [degree]
conditions.timestamp = datetime('now');

% 結果設定
results = struct();
results.P_CBF_all = P_CBF_all;                 % CBF空間スペクトル（全試行）
results.AF_dB = AF_dB;                         % アレーファクタ [dB]
results.scan_theta_deg = rad2deg(scan_theta);  % 角度走査範囲 [degree]

% 保存実行
try
    saved_files = save_simulation(...
        'conditions', conditions, ...
        'results', results, ...
        'experiment_name', 'spatial_spectrum_overlay', ...
        'filename_tag', 'cbf_sidelobe_analysis', ...
        'figure_visible', true, ...
        'auto_commit', true, ...
        'auto_push', true);

    fprintf('実験結果を保存しました:\n');
    for i = 1:length(saved_files)
        fprintf('  %s\n', saved_files{i});
    end
catch ME
    fprintf('保存エラー: %s\n', ME.message);
    % 手動保存
    save('spatial_spectrum_overlay_results.mat', 'conditions', 'results');
end

fprintf('実験完了: CBF空間スペクトル重ね合わせ\n');
fprintf('パラメータ: N=%d, K=%d, SNR=%ddB, 試行回数=%d\n', N, K, SNR_dB, num_trials);
