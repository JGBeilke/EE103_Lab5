close all
clear all
clc

m_1 = repmat(0.9,1,1250);
m_2 = repmat(0.99,1,1250);
C   = repmat(133e-9,1,1250);
R   = repmat(10000,1,1250);

f = 0:0.16:199.84;

w = f*(2 * pi);

H_1 = arrayfun(@transfer, m_1, R, C, w);
H_2 = arrayfun(@transfer, m_2, R, C, w);


figure('name', 'J.G.B. Lab 5 Question 2');
subplot(211);
plot(f, abs(H_1));
xlabel('freq (Hz)');
ylabel('|H(f)|');
title('m = 0.9')

subplot(212);
plot(f, abs(H_2));
xlabel('freq (Hz)');
ylabel('|H(f)|');
title('m = 0.99')

load ecg_signal.mat;

x = ecg;
t = linspace(0, 2.5, length(x));
dt = 0.002;
X = fft(x, length(x))*dt;
W = linspace(-pi, pi, length(x));
w_2 = W/dt;
f_2 = w_2/2/pi;
H_3 = arrayfun(@transfer, m_1, R, C, w_2);

figure('name', 'J.G.B. Lab 5 Question 2 Part B')
subplot(411)
plot(t, ecg);
xlabel('t');
ylabel('x(t)');
title('x(t) vs t')


subplot(412)
plot(f_2, abs(fftshift(X)));
xlabel('f');
ylabel('|X(\omega)|');
title('X(\omega) vs f')

Z = H_3.*fftshift(X);
subplot(413);
plot(f_2, abs(Z), 'b');
xlabel('freq (Hz)');
ylabel('|Z(\omega)|');
title('Z(\omega) vs f');

z = ifft(ifftshift(Z))/dt;
subplot(414);
plot(t, z);
xlabel('t (s)');
ylabel('z(t)');
title('z(t) vs t');

function y = transfer(m, R, C, w)

  y = ((1+m)*((2*j*w*R*C)^2+1))/((2*j*w*R*C)^2+4*(1-m)*j*w*R*C+1);

end
