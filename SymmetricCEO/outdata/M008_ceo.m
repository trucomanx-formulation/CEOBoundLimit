U=[1.935781e-04	5.174394e-04	1.123153e-03	2.124334e-03	3.643211e-03	5.806686e-03	8.742779e-03	1.257744e-02	1.743172e-02	2.341930e-02	3.064430e-02	3.919945e-02	4.916457e-02	6.060528e-02	7.357205e-02	8.809943e-02	1.042056e-01	1.218922e-01	1.411442e-01	1.619301e-01	1.842023e-01	2.078979e-01	2.329388e-01	2.592334e-01	2.866768e-01	3.151525e-01	3.445334e-01	3.746830e-01	4.054573e-01	4.367056e-01	4.682727e-01	5.000000e-01	];
V=[1.850660e-04	5.352222e-04	1.118778e-03	2.117724e-03	3.816396e-03	5.849958e-03	8.497643e-03	1.198726e-02	1.673203e-02	2.350457e-02	3.160689e-02	4.111789e-02	5.139014e-02	6.187311e-02	7.238796e-02	8.893521e-02	1.052632e-01	1.352708e-01	1.428173e-01	1.663418e-01	1.751625e-01	1.876833e-01	2.403756e-01	2.567703e-01	2.865137e-01	2.990654e-01	3.331165e-01	3.781388e-01	3.855422e-01	4.238411e-01	4.688645e-01	5.044335e-01	];
X=[5.000000e-02	6.451613e-02	7.903226e-02	9.354839e-02	1.080645e-01	1.225806e-01	1.370968e-01	1.516129e-01	1.661290e-01	1.806452e-01	1.951613e-01	2.096774e-01	2.241935e-01	2.387097e-01	2.532258e-01	2.677419e-01	2.822581e-01	2.967742e-01	3.112903e-01	3.258064e-01	3.403226e-01	3.548387e-01	3.693548e-01	3.838710e-01	3.983871e-01	4.129032e-01	4.274194e-01	4.419355e-01	4.564516e-01	4.709677e-01	4.854839e-01	5.000000e-01	];
h=figure(1);
hp=semilogy(X,U,'o-',X,V,'>-');
set(hp,'markersize',14);
grid on;
hx=xlabel('rho');
hy=ylabel('BER, M=8');
xlim([min(X),max(X)]);
hl = legend('Theoretical','Experimental');
legend (hl, 'location', 'northoutside', 'orientation', 'horizontal');
ha = gca(); 
FONTSIZE=16;
set(hx,'fontsize',FONTSIZE);
set(hy,'fontsize',FONTSIZE);
set(ha,'fontsize',FONTSIZE);
set(hl,'fontsize',FONTSIZE);
refresh();
saveas (h,'outdata/M008_ceo.eps');

