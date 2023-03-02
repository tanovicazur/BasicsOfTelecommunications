function x = generisiSignal(t,fm)
    Nsf = 0;
while(Nsf < 1)
    Nsf = round(rand*32);
end
    x = zeros(1,size(t,2));
for k = 1:Nsf
    znak = 1;
    if(rand < 0.5);
        znak = -1;
    end
    r_fm = round(rand*fm);
    r_theta = rand*pi;
    r_A = rand;
    x = x + znak*r_A*cos(2*pi*r_fm*t + r_theta);
end