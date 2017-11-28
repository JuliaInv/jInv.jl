function velToConductMod(v,mid,a,b)
d = (b-a)./2.0;
dinv = 10;
tt = dinv.*(mid - v);
t = (d.*(tanh.(tt)+1) + a);
dt = -(dinv*d)*(sech.(tt)).^2;
dt = (2.0-v./mid).*dt + (-1./mid).*t;
t = t.*(2.0-v./mid);
return vec(t),spdiagm(vec(dt))
end
