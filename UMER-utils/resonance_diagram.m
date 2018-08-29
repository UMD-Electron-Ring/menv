% resonance_diagram.m
%
% Kiersten Ruisard
% 10/18/16
%
% plot resonance diagram up to specified order
% mx + ny = k
% x = (k-ny)/m
% x = 1 where y = (k-m)/n
%
function lineobj = resonance_diagram(order)

offset = 0;
linestyle = '-';
color = 'k';
linewidth = .5;
y_limits = ylim; x_limits = xlim;
xoff = ceil(x_limits(1));
yoff = ceil(y_limits(1));

linex=[]; liney=[];
for m = 0:order
    %for n = -order:order
    for n = -(order-m):(order-m)
        for k = -(order+2):(order+2)
            x1 = 0; y1 = k/n; x2 = 1; y2 = (k-m)/n;
            if ((0<=x1 && x1<=1) && (0<=x2 && x2<=1))
                % lines of x projection shorter than 1
                if not(0<=y2 && y2<=1)
                    if y2 < 0 
                        y2 = 0; x2 = k/m;
                    elseif y2 > 1
                        y2 = 1; x2 = (k-n)/m;
                    end
                end
                if not(0<=y1 && y1<=1)
                    if y1 < 0 
                        y1 = 0; x1 = k/m;
                    elseif y1 > 1
                        y1 = 1; x1 = (k-n)/m;
                    end
                end
                linex = [linex,[x1;x2]+xoff]; liney = [liney,[y1;y2]+yoff];
            end
        end
    end
end

% -- get horizontal integer line (nux = k)
x1 = 1; y1 = 0; x2 = 1; y2 = 1;
linex = [linex,[x1;x2]+xoff]; liney = [liney,[y1;y2]+yoff];


npermutations = floor(x_limits(2))-ceil(x_limits(1));
for n=2:npermutations
   linex = [linex,linex+1,linex,linex+1];
   liney = [liney,liney,liney+1,liney+1];
end

lineobj = line(linex,liney,'Color',color,'LineStyle',linestyle,'LineWidth',linewidth);
ylim(y_limits); xlim(x_limits);
end