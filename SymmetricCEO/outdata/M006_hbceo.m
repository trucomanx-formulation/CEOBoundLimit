U=[1.296620e-02	2.462477e-02	4.054135e-02	6.075012e-02	8.515390e-02	1.135561e-01	1.456844e-01	1.812086e-01	2.197551e-01	2.609185e-01	3.042710e-01	3.493702e-01	3.957657e-01	4.430050e-01	4.906381e-01	5.382220e-01	5.853235e-01	6.315238e-01	6.764199e-01	7.196280e-01	7.607849e-01	7.995506e-01	8.356095e-01	8.686720e-01	8.984759e-01	9.247871e-01	9.474013e-01	9.661438e-01	9.808714e-01	9.914715e-01	9.978638e-01	1.000000e+00	];
V=[3.017758e-03	6.588194e-03	1.221930e-02	2.033003e-02	3.129106e-02	4.541684e-02	6.295969e-02	8.410544e-02	1.089703e-01	1.375991e-01	1.699638e-01	2.059637e-01	2.454256e-01	2.881050e-01	3.336878e-01	3.817927e-01	4.319735e-01	4.837233e-01	5.364778e-01	5.896209e-01	6.424897e-01	6.943817e-01	7.445630e-01	7.922786e-01	8.367643e-01	8.772618e-01	9.130371e-01	9.434015e-01	9.677367e-01	9.855216e-01	9.963592e-01	1.000000e+00	];
X=[5.000000e-02	6.451613e-02	7.903226e-02	9.354839e-02	1.080645e-01	1.225806e-01	1.370968e-01	1.516129e-01	1.661290e-01	1.806452e-01	1.951613e-01	2.096774e-01	2.241935e-01	2.387097e-01	2.532258e-01	2.677419e-01	2.822581e-01	2.967742e-01	3.112903e-01	3.258064e-01	3.403226e-01	3.548387e-01	3.693548e-01	3.838710e-01	3.983871e-01	4.129032e-01	4.274194e-01	4.419355e-01	4.564516e-01	4.709677e-01	4.854839e-01	5.000000e-01	];
h=figure(1);
hp=semilogy(X,U,'o-',X,V,'>-');
set(hp,'markersize',14);
grid on;
hx=xlabel('rho');
hy=ylabel('H, M=6');
xlim([min(X),max(X)]);
hl = legend('hb(BER)','H(U0|Omega)');
legend (hl, 'location', 'northoutside', 'orientation', 'horizontal');
ha = gca(); 
FONTSIZE=16;
set(hx,'fontsize',FONTSIZE);
set(hy,'fontsize',FONTSIZE);
set(ha,'fontsize',FONTSIZE);
set(hl,'fontsize',FONTSIZE);
refresh();
saveas (h,'outdata/M006_hbceo.eps');


