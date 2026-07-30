clc; clear;
prob.n=5; prob.q=2; prob.lb=zeros(5,1); prob.ub=ones(5,1);
prob.f=@(x)test_problems('ZDT1',x);
b.N=10; b.max_fevals=500; b.Delta0=.3; b.Dmin=1e-7; b.Dmax=1;
b.eta1=.1; b.eta2=.75; b.gamma=2; b.mu_max=5; b.c_mu=1;
b.p=.9; b.mu_mode='dynamic'; b.ref_pt=[1.1 1.1];
cfg={}; labels={};
o=b; o.mu_mode='dynamic'; cfg{end+1}=o; labels{end+1}='mu dynamic';
o=b; o.mu_mode='fixed'; o.mu_fixed=1; cfg{end+1}=o; labels{end+1}='mu fixed 1';
o=b; o.mu_mode='fixed'; o.mu_fixed=1e-8; cfg{end+1}=o; labels{end+1}='mu fixed 1e-8';
for p=[.3 .6 .9], o=b; o.p=p; cfg{end+1}=o; labels{end+1}=sprintf('p %.1f',p); end
for c=[.01 .1 1 10], o=b; o.c_mu=c; cfg{end+1}=o; labels{end+1}=sprintf('c %.2g',c); end
for i=1:numel(cfg)
    val=zeros(20,3);
    for r=1:20
        rng(r*1000); res=pMODFTR_STCH(prob,cfg{i}); m=evaluate_front(res.F,'ZDT1');
        val(r,:)=[hv2d(res.F,b.ref_pt),m.igd,m.spacing_cv];
    end
    fprintf('%-14s HV %.8f IGD %.8f CV %.8f\n',labels{i},mean(val,1));
end
