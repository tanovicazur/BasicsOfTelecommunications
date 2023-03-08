clear all;
close all;
clc;

fm = 4000;
f = 0:10:fm;
fi = pi/2;
f0 = 100000;
omega0 = 2 * pi * f0;

x1f = f .* 2/fm .* (stepfun(f, 0) - stepfun(f, fm/2)) + (f-fm) .* -2/fm .* (stepfun(f, fm/2) - stepfun(f, fm));
x2f = -2/fm .* (f - fm/2) .* (stepfun(f, 0.1) - stepfun(f, fm/2)) + 2/fm .* (f - fm/2) .* (stepfun(f, fm/2) - stepfun(f, fm)); 

f = -50*f0:10:50*f0;
t = linspace(0, 0.1, 1000000);


extension = zeros(1,((length(f)-1)/2) - length(x1f));
x1f = [extension x1f];
x2f = [extension x2f];

argX = 0;
xtmp1 = x1f .* exp(1j * argX);
xtmp2 = x2f .* exp(1j * argX);
x1f = [xtmp1 0 conj(fliplr(xtmp1(2:end)))];
x2f = [xtmp2 0 conj(fliplr(xtmp2(2:end)))];
x1t = ifft(x1f);
x2t = ifft(x2f);

intFun1 = @(f) abs((f .* 2/fm .* (stepfun(f, 0) - stepfun(f, fm/2)) + (f-fm) .* -2/fm .* (stepfun(f, fm/2) - stepfun(f, fm)))).^2;
intFun2 = @(f) abs(-2/fm .* (f - fm/2) .* (stepfun(f, 0.1) - stepfun(f, fm/2)) + 2/fm .* (f - fm/2) .* (stepfun(f, fm/2) - stepfun(f, fm))).^2;

E1 = integral(intFun1, 0, fm); % E1 = E2 = 1333.33
E2 = integral(intFun2, 0, fm);

center = round(length(f)/2);

figure(1)
subplot(2, 1, 1);
plot(f(center-600:center+600), x1f(center-600:center+600));
title("Signali X1 i X2 u frekventnom domenu");
xlabel("f[Hz]");
ylabel("|X1(f)|");
grid on;

subplot(2, 1, 2);
plot(f(center-600:center+600), x2f(center-600:center+600));
xlabel("f[Hz]");
ylabel("|X2(f)|");
grid on;

figure(2)
subplot(2, 1, 1);
plot(t(1:40000), x1t(1:40000));
title("Signali X1 i X2 u vremenskom domenu");
xlabel("t[s]");
ylabel("X1(t)");
grid on;

subplot(2, 1, 2);
plot(t(1:40000), x2t(1:40000), 'r');
xlabel("t[s]");
ylabel("X2(t)");
grid on;

un1t = sin(omega0.*t);
un12t = sin(omega0.*t - fi);


figure(3);
subplot(2, 1, 1);
plot(t(1:40000), un1t(1:40000));
xlabel("t[s]");
ylabel("un1(t)");
grid on;

subplot(2, 1, 2);
plot(t(1:40000), un12t(1:40000));
xlabel("t[s]");
ylabel("un12(t)");
grid on;

C1 = x1t .* un1t;
C1f = fft(C1);

%Tacke C1 i C2 u vremenskom domenu
%===================================
figure(4);
subplot(2, 1, 1);
plot(t(1:40000), C1(1:40000));
title("Signali C1 i C2 u vremenskom domenu");
xlabel("t[s]");
ylabel("C1(t)");
grid on;

C2 = x2t .* un12t;
C2f = fft(C2);

subplot(2, 1, 2);
plot(t(1:40000), C2(1:40000), 'r');
xlabel("t[s]");
ylabel("C2(t)");
grid on;
%====================================
figure(5);
subplot(2, 1, 1);
plot(f(center+9000:center+11000), abs(C1f(center+9000:center+11000)));
title("Signali C1 i C2 u frekventnom domenu");
xlabel("f[Hz]");
ylabel("C1(f)");
grid on;

subplot(2, 1, 2);
plot(f(center+9000:center+11000), abs(C2f(center+9000:center+11000)));
xlabel("f[Hz]");
ylabel("C2(f)");
grid on;

D = C1 + C2;
Df = fft(D);

%H(jomega) = 20log(Ay/Ax)
%Ay = Ae ; Ax = Ad
%Ae/Ad = 10^-5 => Ae = 10^-0.5 * Ad
%S obzirom da je fazna karakteristika 0 => Ae = Ef

Df = [0 Df];
Ef = 10.^-0.5 * (stepfun(f, f0-fm) - stepfun(f, f0 + fm)).* Df;

Etmp = Ef .* exp(1j * argX);
Ef = [conj(fliplr(Etmp(round(length(Ef)/2):end))) Etmp(round(length(Ef)/2):end-1)];
E = ifft(Ef);

%Tacke D i E u vremenskom domenu
%====================================
figure(6);
subplot(2, 1, 1);
plot(f(center+9000:center+11000), abs(Df(center+9000:center+11000)));
title("Signali D i E u frekventnom domenu");
xlabel("f[Hz]");
ylabel("D(f)");
grid on;

subplot(2, 1, 2);
plot(f(center+9000:center+11000), abs(Ef(center+9000:center+11000)));
xlabel("f[Hz]");
ylabel("E(f)");
grid on;

figure(7);
subplot(2, 1, 1);
plot(t(1:40000), D(1:40000));
title("Signali D i E u vremenskom domenu");
xlabel("t[s]");
ylabel("D(t)");
grid on;

subplot(2, 1, 2);
plot(t(1:40000), E(1:40000), 'r');
xlabel("t[s]");
ylabel("E(t)");
grid on;
%====================================
un1t = [un1t 0];
un12t = [un12t 0];

F1t = E .* un1t;
F2t = E .* un12t;

figure(8);
subplot(2, 1, 1);
plot(t(1:40000), F1t(1:40000));
title("Signali F1 i F2 u vremenskom domenu");
xlabel("t[s]");
ylabel("F1(t)");
grid on;

subplot(2, 1, 2);
plot(t(1:40000), F2t(1:40000), 'r');
xlabel("t[s]");
ylabel("F2(t)");
grid on;

F1f = fft(F1t);
F2f = fft(F2t);

figure(9);
subplot(2, 1, 1);
plot(f(center-600:center+600), abs(F1f(center-600:center+600)));
title("Signali F1 i F2 u frekventnom domenu");
xlabel("f[Hz]");
ylabel("F1(f)");
grid on;

subplot(2, 1, 2);
plot(f(center-600:center+600), abs(F2f(center-600:center+600)));
xlabel("f[Hz]");
ylabel("F2(f)");
grid on;

B1f =  F1f .* (stepfun(f, -fm) - stepfun(f, fm));
B2f = F2f .* (stepfun(f, -fm) - stepfun(f, fm));

figure(10);
subplot(2, 1, 1);
plot(f(center-600:center+600), abs(B1f(center-600:center+600)));
title("Signali B1 i B2 u frekventnom domenu");
xlabel("f[Hz]");
ylabel("B1(f)");
grid on;

subplot(2, 1, 2);
plot(f(center-600:center+600), abs(B2f(center-600:center+600)));
xlabel("f[Hz]");
ylabel("B2(f)");
grid on;

B1tmp = B1f .* exp(1j * argX);
B1f = [conj(fliplr(B1tmp(round(length(B1f)/2):end))) B1tmp(round(length(B1f)/2):end-1)];
B1 = ifft(B1f);

B2tmp = B2f .* exp(1j * argX);
B2f = [conj(fliplr(B1tmp(round(length(B2f)/2):end))) B2tmp(round(length(B2f)/2):end-1)];
B2 = ifft(B2f);

figure(11);
subplot(2, 1, 1);
plot(t(1:40000), B1(1:40000));
title("Signali B1 i B2 u vremenskom domenu");
xlabel("t[s]");
ylabel("B1(t)");
grid on;

subplot(2, 1, 2);
plot(t(1:40000), B2(1:40000), 'r');
xlabel("t[s]");
ylabel("B2(t)");
grid on;




