syms fc t slope dt
deltaT = t-dt;
Tx = expand(cos( 2 * pi * (fc * t + slope * t * t/2)));    
Rx = expand(cos( 2 * pi * (fc * deltaT + slope * deltaT * deltaT/2)));

mix = expand(simplify(Tx*Rx));


syms range_resolution speed_of_light max_range
% Bandwidth
Bsweep = speed_of_light/(2 * range_resolution);

% Chirp Time
chirp_factor = 5.5;
Tchirp = chirp_factor * 2 * max_range/speed_of_light;

% Slope
slope = Bsweep/Tchirp;


n = 10000;
t = linspace(0,10,n);
ks = 1:1:10000;
f = 1/10;
x = zeros(1,n);
Y = zeros(1,n);
B = 2;
for i=1:n
    x(i) = B*(-1/pi*sum(sin(2*pi*ks*f*t(i))./ks)+0.5);
    Y(i) = B*(1/(2*pi*pi*f)*sum(cos(2*pi*ks*f*t(i))./(ks.*ks))+0.5*t(i));
end

figure()
hold on;
plot(t,x)
plot(t,Y)