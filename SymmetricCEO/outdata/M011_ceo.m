U=[5.801346e-06	2.508624e-05	7.936574e-05	2.042030e-04	4.535558e-04	9.023896e-04	1.648135e-03	2.810878e-03	4.532265e-03	6.973191e-03	1.031036e-02	1.473193e-02	2.043235e-02	2.760668e-02	3.644462e-02	4.712435e-02	5.980654e-02	7.462877e-02	9.170016e-02	1.110971e-01	1.328595e-01	1.569877e-01	1.834413e-01	2.121379e-01	2.429538e-01	2.757251e-01	3.102503e-01	3.462938e-01	3.835901e-01	4.218486e-01	4.607596e-01	5.000000e-01	];
V=[5.732737e-05	1.867886e-04	4.787740e-04	1.048368e-03	1.972295e-03	3.281105e-03	5.338004e-03	7.999625e-03	1.261053e-02	2.021558e-02	2.356406e-02	3.244819e-02	3.990336e-02	5.606044e-02	6.392808e-02	7.509533e-02	9.430835e-02	1.087742e-01	1.363515e-01	1.568628e-01	1.847042e-01	2.048820e-01	2.215491e-01	2.700422e-01	2.932417e-01	3.141104e-01	3.411059e-01	3.459460e-01	4.207067e-01	4.331641e-01	4.625113e-01	5.094528e-01	];
X=[5.000000e-02	6.451613e-02	7.903226e-02	9.354839e-02	1.080645e-01	1.225806e-01	1.370968e-01	1.516129e-01	1.661290e-01	1.806452e-01	1.951613e-01	2.096774e-01	2.241935e-01	2.387097e-01	2.532258e-01	2.677419e-01	2.822581e-01	2.967742e-01	3.112903e-01	3.258064e-01	3.403226e-01	3.548387e-01	3.693548e-01	3.838710e-01	3.983871e-01	4.129032e-01	4.274194e-01	4.419355e-01	4.564516e-01	4.709677e-01	4.854839e-01	5.000000e-01	];
h=figure(1);
hp=semilogy(X,U,'o-',X,V,'>-');
set(hp,'markersize',14);
grid on;
hx=xlabel('rho');
hy=ylabel('BER, M=11');
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
saveas (h,'outdata/M011_ceo.eps');

