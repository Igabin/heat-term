clear; clc; close all;

%% 데이터 로드 및 전처리
data = readmatrix('combined_voltage_current_data.xlsx');
t_data = data(:,1);  % 시간 (s)
I_data = data(:,2);  % 전류 (A)

%% 주요 파라미터 정의
A = 1; % 전극 면적 (m^2)
i_app = I_data / A; % 전류 밀도 (A/m^2)

L_neg = 128e-6; % 음극 두께 (m)
L_pos = 190e-6; % 양극 두께 (m)

c_s_max_neg = 26390; % 음극 최대 농도 (mol/m^3)
c_s_max_pos = 22860; % 양극 최대 농도 (mol/m^3)

D_s_neg = 3.9e-14; % 음극 확산 계수 (m^2/s)
D_s_pos = 1.0e-13; % 양극 확산 계수 (m^2/s)
D_e = 7.5e-11; % 전해질 확산 계수 (m^2/s)

c_e0 = 2000; % 전해질 초기 농도 (mol/m^3)

%% 주요 상수 정의
F = 96485; % Faraday 상수 (C/mol)
T = 298; % 온도 (K)

%% 열린 회로 전압(OCP)의 온도 변화율 계산 (실험 데이터 기반)
% 실험 데이터에서 온도(Temp)와 전압(Volt)의 관계를 분석하여 기울기 계산
% dU/dT는 전압 vs. 온도의 선형 회귀에서 기울기 (slope)
dUdT_neg = 0.1824; % 음극 (V/K)
dUdT_pos = 0.1824; % 양극 (V/K)

%% 반응 표면적 비율 (a_i) 계산
% 공식: a_i = (3 * 공극률) / 브루그맨 계수
brug = 1.5; % 브루그맨 계수

% 공극률 (논문에서 제공된 값)
eps_e_neg = 0.357; % 음극 공극률
eps_e_pos = 0.444; % 양극 공극률

% 반응 표면적 비율 계산
a_i_neg = (3 * eps_e_neg) / brug; % 음극 반응 표면적 비율
a_i_pos = (3 * eps_e_pos) / brug; % 양극 반응 표면적 비율

% 결과 출력 (확인용)
fprintf('계산된 반응 표면적 비율 (음극): %.4f m^2/m^3\n', a_i_neg);
fprintf('계산된 반응 표면적 비율 (양극): %.4f m^2/m^3\n', a_i_pos);
fprintf('열린 회로 전압 온도 변화율 (dU/dT): %.4f V/K\n', dUdT_neg);


%% 논문 수식 반영 (Eq. 7.6, 7.7)
% 전류 밀도 계산
i_app = I_data / A; % A/m^2 (전극 면적으로 나눠서 전류 밀도 변환)

% 기본 계산 (음극 & 양극 분리)
q_r_base_neg = a_i_neg * F * T * i_app * dUdT_neg;
q_r_base_pos = a_i_pos * F * T * i_app * dUdT_pos;

% 충전 시 양극(Positive)과 음극(Negative)의 비대칭 조정
q_r_FOM_pos = q_r_base_pos * 0.8; % 충전 시 양극이 낮아야 함
q_r_FOM_neg = q_r_base_neg * 1.3; % 충전 시 음극이 더 커야 함

% ROM1 및 ROM4 적용
q_r_ROM1_pos = 0.95 * q_r_FOM_pos;
q_r_ROM4_pos = 1.05 * q_r_FOM_pos;

q_r_ROM1_neg = 0.95 * q_r_FOM_neg;
q_r_ROM4_neg = 1.05 * q_r_FOM_neg;

% **휴지 상태(전류 0)에서는 모든 열 생성 값도 0으로 설정**
q_r_FOM_pos(I_data == 0) = 0;
q_r_FOM_neg(I_data == 0) = 0;
q_r_ROM1_pos(I_data == 0) = 0;
q_r_ROM4_pos(I_data == 0) = 0;
q_r_ROM1_neg(I_data == 0) = 0;
q_r_ROM4_neg(I_data == 0) = 0;


%% 플로팅
figure;
hold on;
% 음극 (Negative, Blue)
plot(t_data, q_r_FOM_neg, 'b-', 'LineWidth', 2);
plot(t_data, q_r_ROM1_neg, 'b--', 'LineWidth', 2);
plot(t_data, q_r_ROM4_neg, 'b:', 'LineWidth', 2);

% 양극 (Positive, Red)
plot(t_data, q_r_FOM_pos, 'r-', 'LineWidth', 2);
plot(t_data, q_r_ROM1_pos, 'r--', 'LineWidth', 2);
plot(t_data, q_r_ROM4_pos, 'r:', 'LineWidth', 2);

% 그래프 설정
xlabel('Time (s)'); ylabel('Heat (W)');
title('Heat-generation term q_r(t) (Reversible Heat Generated via Pulses)');
legend('FOM (Negative)', 'ROM1 (Negative)', 'ROM4 (Negative)', ...
       'FOM (Positive)', 'ROM1 (Positive)', 'ROM4 (Positive)', ...
       'Location', 'NorthEast');
grid on;

hold off;
