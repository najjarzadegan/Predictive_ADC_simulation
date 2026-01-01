function noise = generate_noise(num_levels, power, size)

% if power = 0  --> distribution = uniform with num_level levels
% if power > 0  --> distribution = normal with sigma = sqrt(power)

% num_levels: number of noise levels
% signal_power: power of signal
% filename: name of the output text file

x = floor(num_levels/2);
% Generate random uniform integer noise

if(power == 0)
    noise = randi([x-num_levels+1,x], [size,1]);
    p = 10*log10(rms(noise)^2)
    filename = sprintf('noise_%d_%s_%.2fdB_m%d_%d.txt',size,'uniform',p,-min(noise),max(noise))

else
    noise = round(10+sqrt(power)*randn(size,1))-10;
    p = 10*log10(rms(noise)^2)
    filename = sprintf('noise_%d_%s_%.2fdB_m%d_%d.txt',size,'normal',p,-min(noise),max(noise))

end

% Save noise to text file

dlmwrite(filename, noise, 'delimiter', '\n');

hist(noise)

end
