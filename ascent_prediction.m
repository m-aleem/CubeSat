% Modified from: http://zmatt.net/weather-balloon-physics (Accessed 2014)

% balloon parameters
r0=0.54; % uninflated radius (m)
initial=0.95; % initial radius (m)
burst=3.5; % burst radius (m)
mb=0.8; % balloon mass (kg)
mp=1.52; % payload mass (kg)
M=2; % molecular mass of gas

% range of altitudes to sweep (m)
h=[0:100:50000];
[rho,a,T,p]=stdatmo(h,0,'SI',true);
rubberrho=1100; % density of rubber (kgm-3)
d0=mb/(4*pi*r0ˆ2*rubberrho); % uninflated thickness
maxstretch=(burst/r0)*1.08; % theoretical max stretch

% (1.08 fudge factor over burst point, need to compare to experiment)
Jm=2*maxstretchˆ2+maxstretchˆ-4-3; % Gent parameter
m=mb+mp;

% determine n corresponding to initial radius
n=moles(p(1),T(1),initial,r0,d0,Jm);
r2=radius_mooneyrivlin(n,p,T,r0,d0,initial);
i2=find(r2>burst,1)-1;
h2=h(i2)
p2=gasp(n,r2,T)-p;
l2=lift(n,M,r2,rho,m);
v2=terminalvelocity(n,M,r2,rho,m);
function n = moles(p,T,r,r0,d0,Jm)

% Determine number of moles of gas in balloon for given radius r
R = 8.3144621; % gas constant
pin = p+gent(r,r0,d0,Jm);
V = 4/3*pi*r.ˆ3;
n = pin.*V./(R*T);

function r = radius_mooneyrivlin(n,p,T,r0,d0,initial)
% Determine balloon radius using Mooney-Rivlin model
mu = 300000; % Pa
R = 8.3144621; % gas constant
p0 = 2*mu*d0/r0;
k = 0.1; % set k=0 for neo-Hookean
X = 3*n*
R/(4*pi);
for i=1:length(p)
  % Upon rearrangement of the equations we get a 7th order polynomial. Since
  % we are stepping from a known initial value, an approach like gradient
  % descent would be more sensible. But who can resist just calling 'roots'?
  poly = [ (p0*k/r0) p(i) (p0*r0) 0 (-X*T(i)) 0 (-p0*k*r0ˆ5) 0 ...(-p0*r0ˆ7) ];rtemp = roots(poly);
  valid = (imag(rtemp) == 0) & (rtemp > 0);
  rtemp = rtemp(valid);
  [w,iw] = min(abs(rtemp-initial));
  r(i) = rtemp(iw);
  initial = rtemp(iw);
end
