clear all;
close all;
clc;

f = 0:10:500e3;
d = 15;
L = 0.6e-6;
C = 0.05e-9;
k = 0.2;
delta_t = 0.1e-6;
tmax = 10e-3;
t = 0:delta_t:tmax;

gammaf = 1j .* 2 .* pi .* f .* sqrt(L .* C .* (1 + ((1-1j) .* k) ./ (L .* sqrt(2 .* pi .* f))));
Hf = exp(-d .* gammaf);
Hfdb = 20 .* log(Hf);

figure,
plot(f, Hfdb);
ylabel("|H(f)|db");
xlabel("f(Hz)");
grid();

argH = 0;
Htmp = Hf .* exp(1j * argH);
Htmp(1) = 0;
Hf = [Htmp  conj(fliplr(Htmp(2:end)))];
h_t = ifft(Hf);

figure,
plot(t(1:550), h_t(1:550), 'r');
ylabel("h(t)");
xlabel("t(s)");
grid();


fm = 500;
Bg = 400;
f = -5e6:100:5e6;

mi_t = zeros(32, length(t));

for i = 1:32
    mi_t(i,:) = generisiSignal(t,fm);
end

c_t = zeros(32, length(t));
x_f = zeros(1, length(f));

for k = 0:31
    fmc = (fm + Bg) .* k;
    c_t(k+1, :) = mi_t(k+1, :) .* cos(2 .* pi .* fmc .* t);
    signal = ussbFilter(c_t(k+1, :), f, fm, fmc);
    x_f = x_f + signal;
end

figure,
plot(f, abs(x_f), 'r');
ylabel("|X(F)|");
xlabel("f(Hz)");
xlim([0 50e3]);
grid();

%argX = 0;
%Xtmp = x_f .* exp(1j * argX);
%x_f = [Xtmp 0 conj(fliplr(Xtmp(2:end)))];
x_t = ifft(x_f);

figure,
plot(t, x_t);
ylabel("x(t)");
xlabel("t(s)");
grid();


fcdsb = 200e3;
y_t = x_t .* cos(2 .* pi .* fcdsb .* t);
y_f = fft(y_t);

figure,
plot(t, y_t);
ylabel("y(t)");
xlabel("t(s)");
grid();

figure,
plot(f, abs(y_f), 'r');
ylabel("|Y(f)|");
xlabel("f(Hz)");
xlim([0 500e3]);
grid();

z_f = y_f .* Hf;

figure,
plot(f, abs(z_f), 'r');
ylabel("|Z(f)|");
xlabel("f(Hz)");
grid();
xlim([0 500e3]);

%{
argZ = 0;
Ztmp = z_f .* exp(1j * argZ);
z_f = [Ztmp 0 conj(fliplr(Ztmp(2:end)))];
z_t = ifft(z_f);
%}

z_t = conv(y_t, h_t);

center = length(z_t)/2;

figure,
plot(t, z_t(center:end));
ylabel("z(t)");
xlabel("t(s)");
grid();

w_t = z_t(center:end) .* cos(2 .* pi .* fcdsb .* t);
w_f = fft(w_t);

figure,
plot(t, w_t);
ylabel("w(t)");
xlabel("t(s)");
grid();

figure,
plot(f, abs(w_f), 'r');
ylabel("|W(f)|");
xlabel("|f(Hz)|");
xlim([0 500e3]);
grid();

r_t = zeros(32, length(t));
r_f = zeros(32, length(t));
fs = 2*10e6;

for k = 0:31
    fmc = (fm + Bg) .* k;
    r_t(k+1, :) = w_t .* cos(2 .* pi .* fmc .* t);
    [b,a] = butter(5, fm/(fs/2));
    r_t(k+1, :) = filter(b, a, r_t(k+1, :));
end

figure,
hold on
plot(t, mi_t(5, :));
plot(t, r_t(5, :));
%xlim([0.004 0.006]);


