function [z, dz, dzz] = disc(t, opts)

arguments
    t
    opts.radius = 1
    opts.center = [0 0]
end

r = opts.radius;
x0 = opts.center(1);
y0 = opts.center(2);

x   =  r*cos(t) + x0;
y   =  r*sin(t) + y0;
dx  = -r*sin(t);
dy  =  r*cos(t);
dxx = -r*cos(t);
dyy = -r*sin(t);

z   = [ x(:).'   ; y(:).'   ];
dz  = [ dx(:).'  ; dy(:).'  ];
dzz = [ dxx(:).' ; dyy(:).' ];

end
