close all
clear all
clc

dt = 0.1;
dw = 0.1;
w  = -31.4:dw:31.4;
t  = linspace(-100, 100, length(w));

G = arrayfun(@Q1_func, w);

g = invFourier(G, w, t);

figure('Name', 'J.G.B. Lab 5 Question 1');
subplot(311);
plot(w, G);
xlabel(' \omega ');
ylabel('G(\omega)');
title('G(\omega) vs \omega');

subplot(312);
plot(t, real(g));
xlabel('t');
ylabel('Real Parts of g(t)');
title('Real Parts of g(t) vs t');

subplot(313);
plot(t, imag(g));
xlabel('t');
ylabel('Imaginary Parts of g(t)')
title('Imaginary Parts of g(t) vs t')

function y = Q1_func(w)

  w = abs(w);
  if (w >= 5) & (w <= 10)

    y = 2;

  else

    y = 0;

  end
end

function f = invFourier(F, w, t)

  for counter = 1:length(w)

    f(counter) = trapz(t, F.*exp(j*w(counter)*t));

  end
end
