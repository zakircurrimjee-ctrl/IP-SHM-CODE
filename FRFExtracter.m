%% ============================================================
%  PLOT ALL ACCELEROMETERS FOR ONE RUN IN dB
%% ============================================================

clear; clc; close all;

filename = "/Users/zakir/Desktop/Thesis/Expiremental data/Section1_FRF_only/Section1_FRF_combined_long.csv";
T = readtable(filename);

T.run = string(T.run);
T.sensor = string(T.sensor);

selectedRun = "EIBS2M2S1";

sensors = unique(T.sensor);

figure;
hold on;

for i = 1:length(sensors)

    sensorName = sensors(i);

    idx = T.run == selectedRun & T.sensor == sensorName;
    FRF = T(idx, :);

    if isempty(FRF)
        continue;
    end

    FRF = sortrows(FRF, "x_value");

    f = FRF.x_value;
    H = FRF.real + 1i * FRF.imaginary;

    mag_dB = 20 * log10(abs(H));

    plot(f, mag_dB, "LineWidth", 1.3, "DisplayName", sensorName);

end

hold off;
grid on;

xlabel("Frequency [Hz]");
ylabel("FRF Magnitude [dB]");
title("FRF Magnitude in dB for All Sensors - " + selectedRun);
legend("Location", "best");

xlim([0 400]);