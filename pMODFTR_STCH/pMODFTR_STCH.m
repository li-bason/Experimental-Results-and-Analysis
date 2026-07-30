function results = pMODFTR_STCH(prob, opts)
% p-MODFTR with STCH models and MOEA/D-style decomposition.
% The acceptance and radius rules follow Algorithm 1 in 多目标.txt.
opts = fill_opts(opts, prob);
n = prob.n; q = prob.q; N = opts.N;
Lambda = make_weights(N,q);

cap = opts.max_fevals + N*n + 10;
aX = zeros(cap,n); aF = zeros(cap,q); nfev = 0;
for i = 1:N                         % plain random initialization (no LHS)
    x = prob.lb(:)+rand(n,1).*(prob.ub(:)-prob.lb(:));
    f = prob.f(x); nfev=nfev+1; aX(nfev,:)=x'; aF(nfev,:)=f(:)';
end
z = min(aF(1:nfev,:),[],1)';
sX=aX(1:nfev,:); sF=aF(1:nfev,:); nsol=nfev;
X=zeros(N,n); FC=zeros(N,q); Delta=opts.Delta0*ones(N,1);
for i=1:N
    [X(i,:),FC(i,:)] = archive_best(aX(1:nfev,:),aF(1:nfev,:),Lambda(i,:)',z);
end

rec_gap=max(1,round(opts.max_fevals/200)); last_rec=nfev;
hv_h=hv2d(FC,opts.ref_pt); fev_h=nfev;
while nfev < opts.max_fevals
    z=min(aF(1:nfev,:),[],1)';
    for ii=randperm(N)
        if nfev>=opts.max_fevals, break; end
        x=X(ii,:)'; f0=FC(ii,:)'; D=Delta(ii); lam=Lambda(ii,:)';
        mu=get_mu(opts,D,q);

        % A model is controllably accurate with conditional probability p.
        % Both events use the shared archive, so p is not confounded with
        % a different function-evaluation budget.
        accurate = rand <= opts.p;
        % Equal-cost geometry set for every p configuration. The random
        % event controls model quality, never the evaluation budget.
        xgeom=zeros(n,n); fgeom=zeros(n,q); hgeom=zeros(n,1); ngeom=0;
        for j=1:n
            if nfev>=opts.max_fevals, break; end
            xp=x; roomp=prob.ub(j)-x(j); roomm=x(j)-prob.lb(j);
            if roomp>=roomm, h=min(D,roomp); else, h=-min(D,roomm); end
            if abs(h)<=eps, continue; end
            xp(j)=xp(j)+h; fp=prob.f(xp); nfev=nfev+1;
            aX(nfev,:)=xp'; aF(nfev,:)=fp(:)';
            ngeom=ngeom+1; xgeom(ngeom,:)=xp'; fgeom(ngeom,:)=fp(:)';
            hgeom(ngeom)=h;
        end
        Xa=aX(1:nfev,:); Fa=aF(1:nfev,:);
        dist=sqrt(sum((Xa-x').^2,2)); ids=find(dist<=2*D+1e-14);
        if numel(ids)>opts.model_max_points
            [~,ord]=sort(dist(ids)); ids=ids(ord(1:opts.model_max_points));
        end
        gm=zeros(numel(ids),1);
        for k=1:numel(ids), gm(k)=stch_val(Fa(ids(k),:)',lam,z,mu); end
        if accurate
            g0=stch_val(f0,lam,z,mu); grd=zeros(n,1);
            for k=1:ngeom
                j=find(abs(xgeom(k,:)'-x)>1e-14,1);
                gp=stch_val(fgeom(k,:)',lam,z,mu);
                grd(j)=(gp-g0)/hgeom(k);
            end
        else
            % Outside the controllable-accuracy event no descent direction
            % is certified. Treat it as an unsuccessful TR iteration.
            grd=zeros(n,1); g0=stch_val(f0,lam,z,mu);
        end
        ng=norm(grd);
        if ng<=1e-12
            if accurate, Delta(ii)=max(D/opts.gamma,opts.Dmin); end
            continue
        end
        s=-D*grd/ng; xn=min(max(x+s,prob.lb(:)),prob.ub(:)); s=xn-x;
        pred=g0-(g0+grd'*s);
        if norm(s)<=1e-14 || pred<=1e-14
            Delta(ii)=max(D/opts.gamma,opts.Dmin); continue
        end
        if nfev>=opts.max_fevals, break; end
        fn=prob.f(xn); nfev=nfev+1; aX(nfev,:)=xn'; aF(nfev,:)=fn(:)';
        actual=tch_val(f0,lam,z)-tch_val(fn,lam,z);
        rho=actual/pred;

        if rho>=opts.eta1
            X(ii,:)=xn'; FC(ii,:)=fn(:)';
            nsol=nsol+1; sX(nsol,:)=xn'; sF(nsol,:)=fn(:)';
            if ng < opts.eta2*D
                Delta(ii)=max(D/opts.gamma,opts.Dmin);
            else
                Delta(ii)=min(opts.gamma*D,opts.Dmax);
            end
        else
            Delta(ii)=max(D/opts.gamma,opts.Dmin);
        end

        % Decomposition collaboration: every evaluated point is available
        % to every subproblem; retain the best Phi value for each weight.
        z=min(aF(1:nfev,:),[],1)';
        for jj=1:N
            [xb,fb]=archive_best(sX(1:nsol,:),sF(1:nsol,:),Lambda(jj,:)',z);
            if tch_val(fb',Lambda(jj,:)',z)<tch_val(FC(jj,:)',Lambda(jj,:)',z)
                X(jj,:)=xb; FC(jj,:)=fb;
            end
        end
        if nfev-last_rec>=rec_gap
            hv_h(end+1)=hv2d(FC,opts.ref_pt); fev_h(end+1)=nfev; last_rec=nfev;
        end
    end
end

% Report the decomposition representatives, not arbitrary initial points.
z=min(aF(1:nfev,:),[],1)';
for i=1:N
    [X(i,:),FC(i,:)]=archive_best(sX(1:nsol,:),sF(1:nsol,:),Lambda(i,:)',z);
end
[FC,ia]=unique(FC,'rows','stable'); X=X(ia,:); nd=nd_idx(FC);
hv_h(end+1)=hv2d(FC,opts.ref_pt); fev_h(end+1)=nfev;
results.hv_hist=hv_h; results.fev_hist=fev_h;
results.F=FC(nd,:); results.X=X(nd,:); results.F_all=FC; results.fevals=nfev;
results.archive_F=aF(1:nfev,:); results.archive_X=aX(1:nfev,:);
results.solution_archive_F=sF(1:nsol,:);
end

function v=tch_val(f,lambda,z)
    v=max(lambda(:).*abs(f(:)-z(:)));
end

function [x,f]=archive_best(X,F,lambda,z)
    v=max(abs(F-z').*lambda',[],2); [~,i]=min(v); x=X(i,:); f=F(i,:);
end
