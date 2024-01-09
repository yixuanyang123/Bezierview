   clf;

% -------- functionality
   shw=1;
   do_prt = 0; % matlab plot:  only makes sense when shw>0
   bvout=1;
   bbnet=1;

% -------- sizes
   cfs = 6;  % coefficients per patch (total degree quadratic)
   idx1 =1:cfs;  % orientation
   idx0 =[3 2 1 5 4 6];
   pats = 6; % patches per triangle
   n = 32;   % evaluation density
   dim = 3;  % dimension of range (3-space)
   tet=0;    % use tetrahedra data

% --------- evaluated BB basis functions
   u0 = linspace(0,1,n);
   [u, v] = meshgrid(u0);
   bbb{1} = (1-u-v).^2;
   bbb{2} = 2*(1-u-v).*u;
   bbb{3} = u.*u;
   bbb{4} = 2*(1-u-v).*v;
   bbb{5} = 2*v.*u;
   bbb{6} = v.*v;
   mask = ones(n,n);  % suppress half of the 4-sided patch
   bidx = find(flipud(triu(ones(n,n),1)));
   mask(bidx) = nan;

% (3 - 3*c0)/(c0 + 1),   3/(2*(c0 + 1))
%  yields for c0 = -1/2     9, 3
% 9 
% 3  1
% 3  1  1
   nn = 3:4;   % allowed valences
   c0 = cos(2*pi./nn);
% rational weights
      wni = 3*(1 - c0)./(c0 + 1);
      wti = 3./(2*(c0 + 1));

% -------- GEOMETRY+CONNECTIVITY:  double simplex
   % vertices 
   V = [...
      cos(2*pi*[0 1 2]/3) 0 0;...
      sin(2*pi*[0 1 2]/3) 0 0;...
      0      0     0      -1 1];
   val = [ 4 4 4 3 3];

   nbr = [ 1 4 3 5 4 1;...
	   2 2 2 2 3 3;...
	   4 3 5 1 1 5]';

   if tet==1,
      % -------- tet
      V =3*[-1  1  1 -1;...
	    -1  1 -1  1;...
	    -1 -1  1  1];
      val = [3 3 3 3];
      nbr = [...
	  2     3     4;...
	   1     4     3;...
	   4     1     2;...
	   3     2     1];
   end;

   nbl = nbr(:,[3 2 1]);
   [dim,vts] = size(V);
   [fcs,vfc] = size(nbr);

%----draw funnel
   if shw>0,
      clr = ['y','c','r'];
      fidx = [1,2,5];
      for jj=1:3,
	 ii = fidx(jj);
      end;
      view(-V(:,vts))
      axis equal
   end;

% --- complete Euclidean part of a single quadratic

   if bvout==1, fp = fopen('mixW.bv','w'); end;

for orient=1:2,
   if orient==1, 
      nbs = nbr; 
      idx = idx1;
      if bvout==1, fprintf(fp,"group %d odd \n",orient); end;
   else
      nbs = nbl;
      idx = idx0;
      if bvout==1, fprintf(fp,"group %d even \n",orient); end;
   end;
   for ff=1:fcs, 
      fc{ff} = V(:,nbs(ff,:));
      ctr(:,ff) = fc{ff}*ones(vfc,1)/vfc;
      for kk=1:vfc,
	 dual{ff}(:,kk) = (3*ctr(:,ff)+fc{ff}(:,kk))/4; % 
      end;
   end;
   % --- project vtx neighbors MISSING  (not needed for specific geometry)
   for ii=1:vts,
      % project duals   -- currently left out since valence 3 or 4
   end;
   for ff=1:fcs, 

      for kk=1:vfc, % k = index inside face if rem(pat,2)==1,
	 pat = 2*kk+orient;
	 km = kk-1; if km < 1, km= kk+vfc-1; end;
	 kp = kk+1; if kp > vfc, kp= kk-vfc+1; end;
	 top = nbs(ff,kk); % top global vtx index
	 nxt = nbs(ff,kp); % bottom left vtx index
	 prv = nbs(ff,km); % bottom right vtx index
	 %  kk
	 %  kp  km
	 % nxt  prv
	 ntp =0;
	 for tt=1:fcs, % neighbor triangle nxt
	    for ii =1:vfc,
		ip = ii+1; if ip>vfc, ip=1; end;
		if (nbs(tt,ii) == nxt) && nbs(tt,ip) == top,
		    ntp = tt;  break; end;
	     end;
	     if ntp ~= 0, break; end;
	 end;
	 ntm =0;
	 for jj=1:fcs, % neighbor triangle prv-top edge
	    for nn =1:vfc,
		np = nn+1; if np>vfc, np=1; end;
		if (nbs(jj,nn) == top) && (nbs(jj,np) == prv),
		    ntm = jj; break; end;
	     end;
	     if ntm ~= 0, break; end;
	 end;
	 % v10  v00 ff
	 % v11  v01
	 vf0 = dual{ff}(:,kk);  
	 vfp = dual{ff}(:,kp); 
	 vfm = dual{ff}(:,km); 
	 v0p = dual{ntp}(:,ip);  
	 v1p = dual{ntp}(:,ii); 
	 v0m = dual{ntm}(:,nn);
	 if (pat==1 && ff==1),
	    plot3(vf0(1),vf0(2),vf0(3),'ro'); hold on;
	    plot3(v0p(1),v0p(2),v0p(3),'g*'); hold on;
	    plot3(vfp(1),vfp(2),vfp(3),'b*'); hold on;
	    plot3(v1p(1),v1p(2),v1p(3),'k+'); hold on;
	 end;

	 wti_top = wti(val(top)-2);
	 wti_bot = wti(val(nxt)-2);
	 wni_top = wni(val(top)-2);
         %w_top = 2*(1+c0(val(top)-2))/3;
         %w_bot = 2*(1+c0(val(nxt)-2))/3;
         %w_prv = 2*(1+c0(val(prv)-2))/3;
         w_top = 1;
         w_bot = 1;
         w_prv = 1;
	 w1 = sqrt(wti_top*wti_bot);
	 w5 = w_top;
	 w4 = w_top*wti_top;
	 wgt = [ w1,  w5,  (w_top+w_bot+w_prv)/3, ...
                 w4,  w5, ...
                 w_top*wni_top ];

	 % --- assemble by averaging in coeff_i weight_i
	 qE(:,5) = wgt(5)*vf0;
	 qE(:,2) = wgt(2)*(vf0+vfp)/2;  
	 midopp  = wgt(2)*(v0p+v1p)/2;
	 qE(:,4) = wgt(4)*(vf0+v0p)/2;
	 qE(:,1) = wgt(1)*(vf0+vfp+v0p+v1p)/4;
	 qE(:,3) = (w_top*vf0+w_bot*vfp+w_prv*vfm)/3; 

	 % top is average:  v0p-o + v0m-o = (vf0-o)*2c0 
	 %  v0p+v0m-vf0(2c0) = 2o-2c0o = 2(1-c0)o
	 cc = c0(val(top)-2);
	 tv =  (v0p+v0m-2*cc*vf0)/(2*(1-cc));  % top vertex  Euclidian
	 qE(:,6) =  wgt(6)*tv; 
	 if (pat==1 && ff==1 &&kk==1),
	    plot3(v0m(1),v0m(2),v0m(3),'co'); hold on;
	    plot3(tv(1),tv(2),tv(3),'r+'); hold on;
	 end;

	 % --- assemble
	 bbase = [qE; wgt]';
	 if shw==4, if ff==1, show(bbase,dim,cfs,bbb,mask,'b'); end; end;

	 % -----export --------------------
	 bez{ff}{pat} = bbase(idx,:);
	 if bvout==1, 
	    fprintf(fp,"%d\n%d\n",11,2);
            for cf=1:cfs, % BB-coeffs of patch
               fprintf(fp,"%f %f %f %f\n", bez{ff}{pat}(cf,:)); end;
	 end;
	 if shw==3, show(bez{ff}{pat},dim,cfs,bbb,mask,'r'); end;
	 if bbnet==1,
	    of = 0.01;
	    for ii=1:cfs,
	       ww = bbase(ii,4);
	       xx = bbase(ii,1)/ww;
	       yy = bbase(ii,2)/ww;
	       zz = bbase(ii,3)/ww; 
	       text(xx+of,yy+of,zz+of, num2str(ww)); hold on;
	    end;
	    ids = [1 2 4; 2 3 5; 4 5 6];
	    for ii=1:3,% three subtriangles of bb-net of one quadratic
	       for jj=1:3,% each corner of a subtriangle
		  jp = jj+1; if jp==4; jp=1; end;
		  ll = [ids(ii,jj), ids(ii,jp)];
		  if rem(ff,2)==1,
		  p = plot3(bbase(ll,1)./bbase(ll,4), bbase(ll,2)./bbase(ll,4), ...
			bbase(ll,3)./bbase(ll,4),'k'); hold on;
		  else
		  p = plot3(bbase(ll,1)./bbase(ll,4), bbase(ll,2)./bbase(ll,4), ...
			bbase(ll,3)./bbase(ll,4),'r'); hold on;
		  end
		  set(p,'Linewidth',3');
	       end;
	    end;
	 end;
      end;
   end;
end;
if bvout==1, fclose(fp); end;

if do_prt==1,
   print('tst','-dpng'); 
end




