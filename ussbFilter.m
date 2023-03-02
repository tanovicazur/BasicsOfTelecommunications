function signal = ussbFilter(x,f, fm, fmc)
     %{ 
    step_f = stepfun(f, -fmc - fm) - stepfun(f, -fmc) + stepfun(f, fmc) - stepfun(f, fmc + fm);
    step_t = ifft(step_f);
    signal = conv(x, step_t, "same");
     %}
    
    step = stepfun(f, -fmc - fm) - stepfun(f, -fmc) + stepfun(f, fmc) - stepfun(f, fmc + fm);
    temp = fftshift(fft(x));
    signal = temp .* step;
    
    %{
    argSignal = 0;
    signalTmp = signalf .* exp(1j * argSignal);
    signalf = [signalTmp 0 conj(fliplr(signalTmp(2:end)))];
    signal = ifft(signalf);
    %}
end