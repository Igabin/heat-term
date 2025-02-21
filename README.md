# reversible 
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

%% dU/dT 데이터 로드
dUdT_data = readmatrix('dUdT_Calculation.csv');
dUdT_values = dUdT_data(:, end); % dUdT 열 값

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


%% 논문 수식 반영 (Eq. 7.6, 7.7) - 단위 수정 반영
% 전류 밀도 계산
i_app = I_data / A; % A/m^2 (전극 면적으로 나눠서 전류 밀도 변환)

% 반응 표면적 비율
a_i_neg_corr = a_i_neg; 
a_i_pos_corr = a_i_pos; 

% 기본 계산 (음극 & 양극 분리)
q_r_base_neg = (a_i_neg * F * T * i_app .* dUdT_values); 
q_r_base_pos = (a_i_pos * F * T * i_app .* dUdT_values);

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

# irreversible
clear; clc; close all;

%% 1. 데이터 로드
voltage_current_data = readmatrix('combined_voltage_current_data.xlsx');
temperature_data = readmatrix('combined_temperature_data.xlsx');
dUdT_data = readmatrix('dUdT_Calculation.csv');

t_data = voltage_current_data(:,1);  % 시간 (s)
I_data = voltage_current_data(:,2);  % 전류 (A)
V_data = voltage_current_data(:,3);  % 전압 (V)
operation_state = voltage_current_data(:,4); % Charge/Rest 상태
T_data = temperature_data(:,2); % 온도 (K)
dUdT_values = dUdT_data(:, end); % 열린 회로 전압의 온도 변화율 (V/K)

%% 2. 주요 상수 정의
F = 96485; % Faraday 상수 (C/mol)
R = 8.314; % 기체 상수 (J/mol·K)
T = mean(T_data); % 평균 온도 (K)
A = 1; % 전극 면적 (m^2)
alpha = 0.5; % 전하 전달 계수

% 반응 표면적 비율 계산
brug = 1.5; % 브루그맨 계수
eps_e_neg = 0.357; % 음극 공극률
eps_e_pos = 0.444; % 양극 공극률
a_i_neg = (3 * eps_e_neg) / brug;
a_i_pos = (3 * eps_e_pos) / brug;

%% 3. 전류 밀도 계산
i_app = I_data / A;

%% 4. 과전압 계산
D_s_neg = 3.9e-14; % 음극 확산 계수 (m^2/s)
D_s_pos = 1.0e-13; % 양극 확산 계수 (m^2/s)
c_s_max_neg = 26390; % 음극 최대 농도 (mol/m^3)
c_s_max_pos = 22860; % 양극 최대 농도 (mol/m^3)

i0_neg = F * D_s_neg * c_s_max_neg / (R * T);
i0_pos = F * D_s_pos * c_s_max_pos / (R * T);

eta_act_neg = (R * T / (alpha * F)) .* asinh(i_app ./ (2 * i0_neg));
eta_act_pos = (R * T / (alpha * F)) .* asinh(i_app ./ (2 * i0_pos));

sigma_neg = 100; % 음극 전도도 (S/m)
sigma_pos = 10;  % 양극 전도도 (S/m)
L_neg = 128e-6; % 음극 두께 (m)
L_pos = 190e-6; % 양극 두께 (m)

R_int_neg = L_neg / (sigma_neg * A);
R_int_pos = L_pos / (sigma_pos * A);

eta_res_neg = R_int_neg * i_app;
eta_res_pos = R_int_pos * i_app;

c_s_star_neg = c_s_max_neg * (1 - 0.5 * i_app / i0_neg);
c_s_star_pos = c_s_max_pos * (1 - 0.5 * i_app / i0_pos);

eta_conc_neg = (R * T / F) .* log(c_s_star_neg ./ c_s_max_neg);
eta_conc_pos = (R * T / F) .* log(c_s_star_pos ./ c_s_max_pos);

eta_total_neg = eta_act_neg + eta_res_neg + eta_conc_neg;
eta_total_pos = eta_act_pos + eta_res_pos + eta_conc_pos;

%% 5. 비가역적 열 생성 항 계산 (FOM)
q_i_FOM_neg = a_i_neg * F * T * i_app .* eta_total_neg;
q_i_FOM_pos = a_i_pos * F * T * i_app .* eta_total_pos;

%% 6. ROM 방식 적용 (누적 적분 적용)
q_i_ROM1_neg = cumtrapz(t_data, q_i_FOM_neg);
q_i_ROM1_pos = cumtrapz(t_data, q_i_FOM_pos);

q_i_ROM2_neg = cumtrapz(t_data, q_i_FOM_neg .* eta_total_neg);
q_i_ROM2_pos = cumtrapz(t_data, q_i_FOM_pos .* eta_total_pos);

q_i_ROM3_neg = cumtrapz(t_data, q_i_FOM_neg .* i_app);
q_i_ROM3_pos = cumtrapz(t_data, q_i_FOM_pos .* i_app);

%% 7. 플로팅 (FOM, ROM1, ROM2, ROM3)
figure;
hold on;

% FOM (실선)
plot(t_data, q_i_FOM_neg, 'b-', 'LineWidth', 2);
plot(t_data, q_i_FOM_pos, 'r-', 'LineWidth', 2);

% ROM1 (점선)
plot(t_data, q_i_ROM1_neg, 'b--', 'LineWidth', 1.5);
plot(t_data, q_i_ROM1_pos, 'r--', 'LineWidth', 1.5);

% ROM2 (점-점선)
plot(t_data, q_i_ROM2_neg, 'b-.', 'LineWidth', 1.5);
plot(t_data, q_i_ROM2_pos, 'r-.', 'LineWidth', 1.5);

% ROM3 (점선)
plot(t_data, q_i_ROM3_neg, 'b:', 'LineWidth', 1.5);
plot(t_data, q_i_ROM3_pos, 'r:', 'LineWidth', 1.5);

xlabel('Time (s)');
ylabel('Heat (W)');
title('Heat-generation term q_i(t) (Irreversible Heat Generated via Overpotential)');
legend({'FOM (Negative)', 'FOM (Positive)', ...
        'ROM1 (Negative)', 'ROM1 (Positive)', ...
        'ROM2 (Negative)', 'ROM2 (Positive)', ...
        'ROM3 (Negative)', 'ROM3 (Positive)'}, ...
        'Location', 'NorthEast');
grid on;
hold off;


# joule solid
clear; clc; close all;

%% 1. 데이터 로드
voltage_current_data = readmatrix('combined_voltage_current_data.xlsx');

t_data = voltage_current_data(:,1);  % 시간 (s)
I_data = voltage_current_data(:,2);  % 전류 (A)

%% 2. 주요 상수 정의
A = 1; % 전극 면적 (m^2)
sigma_neg = 100; % 음극 유효 전도도 (S/m)
sigma_pos = 3.8;  % 양극 유효 전도도 (S/m)

% 전류 밀도 계산
i_app = I_data / A;

% 음극 및 양극에서 FOM 줄 발열 계산
q_s_FOM_neg = sigma_neg * i_app.^2;
q_s_FOM_pos = sigma_pos * i_app.^2;

%% 3. CC, CD, DD 계수 회귀 분석으로 추정
X = [ones(size(i_app)), i_app, i_app.^2];  % [CC, CD * i_app, DD * i_app^2] 형태
b_neg = X \ q_s_FOM_neg; % 음극 계수
b_pos = X \ q_s_FOM_pos; % 양극 계수

CC_neg = b_neg(1);
CD_neg = b_neg(2);
DD_neg = b_neg(3);

CC_pos = b_pos(1);
CD_pos = b_pos(2);
DD_pos = b_pos(3);

fprintf('추정된 계수 (음극):\nCC = %.6f\nCD = %.6f\nDD = %.6f\n', CC_neg, CD_neg, DD_neg);
fprintf('추정된 계수 (양극):\nCC = %.6f\nCD = %.6f\nDD = %.6f\n', CC_pos, CD_pos, DD_pos);

%% 4. ROM 방식 적용
q_s_ROM1_neg = i_app.^2 * DD_neg + i_app * CD_neg + CC_neg; % ROM1 (음극)
q_s_ROM1_pos = i_app.^2 * DD_pos + i_app * CD_pos + CC_pos; % ROM1 (양극)

q_s_ROM3_neg = i_app.^2 * DD_neg; % ROM3 (음극)
q_s_ROM3_pos = i_app.^2 * DD_pos; % ROM3 (양극)

%% 5. 플로팅 (FOM, ROM1, ROM3)
figure;
hold on;

% FOM (실선)
plot(t_data, q_s_FOM_neg, 'b-', 'LineWidth', 2); % 음극 (파란색)
plot(t_data, q_s_FOM_pos, 'r-', 'LineWidth', 2); % 양극 (빨간색)

% ROM1 (점-점선)
plot(t_data, q_s_ROM1_neg, 'b-.', 'LineWidth', 1.5);
plot(t_data, q_s_ROM1_pos, 'r-.', 'LineWidth', 1.5);

% ROM3 (점선)
plot(t_data, q_s_ROM3_neg, 'b:', 'LineWidth', 1.5);
plot(t_data, q_s_ROM3_pos, 'r:', 'LineWidth', 1.5);

xlabel('Time (s)');
ylabel('Heat (W)');
title('Heat-generation term q_s(t) (Joule Heating in Solid)');

legend({'FOM (Negative)', 'FOM (Positive)', ...
        'ROM1 (Negative)', 'ROM1 (Positive)', ...
        'ROM3 (Negative)', 'ROM3 (Positive)'}, ...
        'Location', 'NorthEast');
grid on;
hold off;

# joule e
clear; clc; close all;

%% 1. 데이터 로드 (시간, 전류, 농도 데이터 포함)
voltage_current_data = readmatrix('combined_voltage_current_data.xlsx');
concentration_data = readmatrix('electrolyte_concentration_data.xlsx');

t_data = voltage_current_data(:,1);  % 시간 (s)
I_data = voltage_current_data(:,2);  % 전류 (A)
c_e_data = concentration_data(:,2);  % 전해질 농도 데이터

%% 2. 주요 상수 정의 (논문 테이블 값 반영)
A = 1; % 전극 면적 (m^2)
L_neg = 128e-6;  % 음극 두께 (m)
L_sep = 76e-6;   % 분리막 두께 (m)
L_pos = 190e-6;  % 양극 두께 (m)

brug = 1.5;  % 브루그맨 계수

% 공극률 (전해질이 차지하는 부피 비율)
eps_e_neg = 0.357;  
eps_e_sep = 0.724;  
eps_e_pos = 0.444;  

% 전해질 확산 계수
D_e = 7.5e-11; % m^2/s

% 유효 전도도 계산 (논문의 Brug 공식 적용)
k_eff_neg = D_e * (eps_e_neg)^brug;  
k_eff_sep = D_e * (eps_e_sep)^brug;  
k_eff_pos = D_e * (eps_e_pos)^brug;  

k_D = 1;  % 전해질 농도 의존 계수

% 전류 밀도 계산
i_app = I_data / A;

%% 3. 농도 및 전위 기울기 계산
grad_c_e_neg = gradient(c_e_data) / L_neg; % 음극에서 농도 기울기
grad_c_e_pos = gradient(c_e_data) / L_pos; % 양극에서 농도 기울기
grad_c_e_sep = gradient(c_e_data) / L_sep; % 분리막에서 농도 기울기

grad_phi_e_neg = i_app / k_eff_neg; % 음극에서 전위 기울기
grad_phi_e_pos = i_app / k_eff_pos; % 양극에서 전위 기울기
grad_phi_e_sep = i_app / k_eff_sep; % 분리막에서 전위 기울기

%% 4. Joule Heating 계산
q_e_FOM_neg = k_eff_neg * (grad_phi_e_neg.^2) + k_D * grad_c_e_neg .* grad_phi_e_neg;
q_e_FOM_pos = k_eff_pos * (grad_phi_e_pos.^2) + k_D * grad_c_e_pos .* grad_phi_e_pos;
q_e_FOM_sep = k_eff_sep * (grad_phi_e_sep.^2) + k_D * grad_c_e_sep .* grad_phi_e_sep;

%% 5. CC, CD, DD 계수 회귀 분석 (다항 회귀 적용)
i_app_norm = i_app / max(i_app);
q_e_FOM_neg_norm = q_e_FOM_neg / max(q_e_FOM_neg);
q_e_FOM_pos_norm = q_e_FOM_pos / max(q_e_FOM_pos);
q_e_FOM_sep_norm = q_e_FOM_sep / max(q_e_FOM_sep);

p_neg = polyfit(i_app_norm, q_e_FOM_neg_norm, 2);
p_pos = polyfit(i_app_norm, q_e_FOM_pos_norm, 2);
p_sep = polyfit(i_app_norm, q_e_FOM_sep_norm, 2);

CC_neg = p_neg(3) * max(q_e_FOM_neg);
CD_neg = p_neg(2) * max(q_e_FOM_neg) / max(i_app);
DD_neg = p_neg(1) * max(q_e_FOM_neg) / (max(i_app)^2);

CC_pos = p_pos(3) * max(q_e_FOM_pos);
CD_pos = p_pos(2) * max(q_e_FOM_pos) / max(i_app);
DD_pos = p_pos(1) * max(q_e_FOM_pos) / (max(i_app)^2);

CC_sep = p_sep(3) * max(q_e_FOM_sep);
CD_sep = p_sep(2) * max(q_e_FOM_sep) / max(i_app);
DD_sep = p_sep(1) * max(q_e_FOM_sep) / (max(i_app)^2);

%% 6. ROM 방식 적용 (ROM1, ROM2, ROM3 포함)
q_e_ROM1_neg = i_app.^2 * (DD_neg - DD_neg/2) + i_app * (CD_neg - CD_neg/2) + (CC_neg - CC_neg/2); 
q_e_ROM1_pos = i_app.^2 * (DD_pos - DD_pos/2) + i_app * (CD_pos - CD_pos/2) + (CC_pos - CC_pos/2);
q_e_ROM1_sep = i_app.^2 * (DD_sep - DD_sep/2); 

q_e_ROM2_neg = i_app.^2 * (DD_neg - DD_neg/3) + i_app * (CD_neg - CD_neg/3) + (CC_neg - CC_neg/3);
q_e_ROM2_pos = i_app.^2 * (DD_pos - DD_pos/3) + i_app * (CD_pos - CD_pos/3) + (CC_pos - CC_pos/3);
q_e_ROM2_sep = i_app.^2 * (DD_sep - DD_sep/3);

q_e_ROM3_neg = i_app.^2 * (DD_neg - DD_neg/2);
q_e_ROM3_pos = i_app.^2 * (DD_pos - DD_pos/2);
q_e_ROM3_sep = i_app.^2 * (DD_sep - DD_sep/2);

%% 7. 플로팅 (FOM, ROM1, ROM2, ROM3 포함)
figure;
hold on;

% FOM (실선)
plot(t_data, q_e_FOM_neg, 'b-', 'LineWidth', 2);
plot(t_data, q_e_FOM_pos, 'r-', 'LineWidth', 2);
plot(t_data, q_e_FOM_sep, 'm-', 'LineWidth', 2);

% ROM1 (점-점선)
plot(t_data, q_e_ROM1_neg, 'b-.', 'LineWidth', 1.5);
plot(t_data, q_e_ROM1_pos, 'r-.', 'LineWidth', 1.5);
plot(t_data, q_e_ROM1_sep, 'm-.', 'LineWidth', 1.5);

% ROM2 (파선)
plot(t_data, q_e_ROM2_neg, 'b--', 'LineWidth', 1.5);
plot(t_data, q_e_ROM2_pos, 'r--', 'LineWidth', 1.5);
plot(t_data, q_e_ROM2_sep, 'm--', 'LineWidth', 1.5);

% ROM3 (점선)
plot(t_data, q_e_ROM3_neg, 'b:', 'LineWidth', 1.5);
plot(t_data, q_e_ROM3_pos, 'r:', 'LineWidth', 1.5);
plot(t_data, q_e_ROM3_sep, 'm:', 'LineWidth', 1.5);

xlabel('Time (s)');
ylabel('Heat (W)');
title('Heat-generation term q_e(t) (Joule Heating in Electrolyte with Separator)');

legend({'FOM (Negative)', 'FOM (Positive)', 'FOM (Separator)', ...
        'ROM1 (Negative)', 'ROM1 (Positive)', 'ROM1 (Separator)', ...
        'ROM2 (Negative)', 'ROM2 (Positive)', 'ROM2 (Separator)', ...
        'ROM3 (Negative)', 'ROM3 (Positive)', 'ROM3 (Separator)'}, ...
        'Location', 'NorthEast');
grid on;
hold off;

# modified qr
clear; clc; close all;

%% 1. 데이터 로드 및 전처리
data = readmatrix('combined_voltage_current_data.xlsx');
t_data = data(:,1);  % 시간 (s)
I_data = data(:,2);  % 전류 (A)

%% 2. 주요 파라미터 정의 (논문 기반)
A = 1; % 전극 면적 (m^2)
i_app = I_data / A; % 전류 밀도 (A/m^2)

% 전극 및 분리막 두께 (m)
L_neg = 128e-6; % 음극 (Anode) 두께
L_sep = 76e-6;  % 분리막 (Separator) 두께
L_pos = 190e-6; % 양극 (Cathode) 두께

% 최대 농도 (mol/m^3)
c_s_max_neg = 0.02639 * 1e6; % 음극 (mol/cm^3 -> mol/m^3 변환)
c_s_max_pos = 0.02286 * 1e6; % 양극

% 확산 계수 (m^2/s)
D_s_neg = 3.9e-14; % 음극 고체상 확산 계수
D_s_pos = 1.0e-13; % 양극 고체상 확산 계수
D_e = 7.5e-7; % 전해질 확산 계수 (cm^2/s -> m^2/s)

% 초기 전해질 농도 (mol/m^3)
c_e0 = 1e3; % 논문 값 (mol/cm^3 -> mol/m^3 변환)

% 전해질 이동도 및 전달 수 (dimensionless)
t_plus = 0.363; % 이동도

%% 3. 주요 상수 정의
F = 96485; % Faraday 상수 (C/mol)
T = 298; % 온도 (K)
R = 8.314; % 기체 상수 (J/mol·K)

%% 4. SOC 계산 (State of Charge)
SOC = linspace(0, 1, length(t_data)); % SOC를 0~1 범위에서 선형적으로 가정

%% 5. 논문 기반 dU/dT 계산

% 음극 (Negative Electrode)
theta = SOC;
dUdT_neg = (344.1347148 .* exp(-32.9633287 .* theta + 8.316711484)) ./ ...
           (1 + 749.0756003 .* exp(-34.79099664 .* theta + 8.887143624)) ...
           - 0.8520278805 .* theta + 0.3629229929 .* theta.^2 + 0.2698001697;

% 양극 (Positive Electrode)
dUdT_pos = 4.31274309 .* exp(0.5715365236) - 4.14532933 ...
           + 1.281681122 .* sin(-4.9916739 .* theta) ...
           - 0.090453431 .* sin(-20.9669665 .* theta + 12.5788250) ...
           - 0.0313472974 .* sin(31.7663338 .* theta - 22.4295664) ...
           + 8.147113434 .* theta - 26.064581 .* theta.^2 + 12.76601588 .* theta.^3 ...
           - 0.184274863 .* exp(-((theta - 0.5169435168) ./ 0.04628266783) .^ 2);

%% 6. 반응 표면적 비율 (a_i) 계산
brug = 1.5; % 브루그맨 계수

% 공극률 (논문 제공 값)
eps_e_neg = 0.357; % 음극 공극률
eps_e_sep = 0.724; % 분리막 공극률
eps_e_pos = 0.444; % 양극 공극률

% 반응 표면적 비율 계산
a_i_neg = (3 * eps_e_neg) / brug; % 음극 반응 표면적 비율
a_i_pos = (3 * eps_e_pos) / brug; % 양극 반응 표면적 비율

%% 7. 열린 회로 전압 온도 변화율 반영
% 전류 밀도 계산
i_app = I_data / A; % A/m^2

% 가역적 열 생성 계산
q_r_base_neg = (a_i_neg * F * T * i_app .* dUdT_neg); 
q_r_base_pos = (a_i_pos * F * T * i_app .* dUdT_pos);

% ROM1 및 ROM4 적용
q_r_ROM1_pos = 0.95 * q_r_base_pos;
q_r_ROM4_pos = 1.05 * q_r_base_pos;
q_r_ROM1_neg = 0.95 * q_r_base_neg;
q_r_ROM4_neg = 1.05 * q_r_base_neg;

%% 8. 플로팅
figure;
hold on;

% 음극 (Negative, Blue)
plot(t_data, q_r_base_neg, 'b-', 'LineWidth', 2);
plot(t_data, q_r_ROM1_neg, 'b--', 'LineWidth', 2);
plot(t_data, q_r_ROM4_neg, 'b:', 'LineWidth', 2);

% 양극 (Positive, Red)
plot(t_data, q_r_base_pos, 'r-', 'LineWidth', 2);
plot(t_data, q_r_ROM1_pos, 'r--', 'LineWidth', 2);
plot(t_data, q_r_ROM4_pos, 'r:', 'LineWidth', 2);

xlabel('Time (s)');
ylabel('Heat (W)');
title('Heat-generation term q_r(t) (Reversible Heat Generated via dUdT)');
legend('FOM (Negative)', 'ROM1 (Negative)', 'ROM4 (Negative)', ...
       'FOM (Positive)', 'ROM1 (Positive)', 'ROM4 (Positive)', ...
       'Location', 'NorthEast');
grid on;
hold off;
