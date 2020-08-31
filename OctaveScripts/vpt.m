% External iputs:
n0=2; %toroidal number
m0=3; %poloidal number
ns=3000 % number of flux surfaces
s=[0:1:ns]/ns; %normalized toroidal flux
iota=1./(1+2*s.^2); %iota
% End external iputs

nsbar=ns  %number of points (flux surfaces) for rescaled label
iota_res=n0/m0; %resonant iota value
ns_vpt=8000;
nphi_vpt=300;

nphi_vpt = make_even(nphi_vpt); %nphi_vpt must be even

s_tor=[0:1:ns_vpt]/ns_vpt;
iota_vpt = interp1(s,iota,s_tor,'linear','extrap');
phi_vpt=[0:1:nphi_vpt]/nphi_vpt*2*pi;
cosphi=cos(phi_vpt);
s_res=interp1(iota_vpt-iota_res,s_tor,0); %resonant radius

del_s=0.01 %approximate separatrix half-width
% Igichine model: del_s = 4 sqrt(s_res) ka / e,   ka=0.06
ka=0.06
del_s = Igichine_model_del_s(s_res, ka)

iota_pr=interp1(s_tor,gradient(iota_vpt)./gradient(s_tor),s_res) %derivative of iota at resonant radius
s_pol=cumtrapz(s_tor,iota_vpt); %poloidal flux normalized to toroidal flux at the edge
s_hel=n0*s_tor-m0*s_pol; %unperturbed helical flux

ampl=m0*iota_pr*del_s^2/4; %amplitude of perturbed potential
pow=(s_tor/s_res).^(m0/2);
wid=del_s*m0/2;
ex=exp((pow-1)/wid);
ex0=exp(-1/wid);
pertshel=ampl*(ex-ex0+pow)./(ex.*pow-ex0+1); %amplitude times shape function of perturbed potetial

%pertshel(s>=s_res)=ampl*(s_res./s(s>=s_res)).^(m0/2);
%pertshel(s<s_res)=ampl*(s(s<s_res)./s_res).^(m0/2);

s_hel_pert=s_hel'*ones(size(cosphi))+pertshel'*cosphi; %perturbed helical flux (2D)

s_hel_xpoint=s_hel_pert(:,nphi_vpt/2+1); %perturbed helical flux at phi=pi (section crossing X-point)

s_hel_sep=interp1(s_tor,s_hel_xpoint,s_res); %helical flux at the separatrix

s_in=s_tor(s_tor<s_res);                     %values of array of unperturbed toroidal flux for the inner region
s_pol_in=s_pol(s_tor<s_res);                 %values of array of unperturbed poloidal flux for the inner region
s_hel_xpoint_in=s_hel_xpoint(s_tor<s_res);   %values of array of perturbed flux at phi=pi for the inner region

s_out=s_tor(s_tor>=s_res);                   %values of array of unperturbed toroidal flux for the outer region
s_pol_out=s_pol(s_tor>=s_res);               %values of array of unperturbed poloidal flux for the outer region
s_hel_xpoint_out=s_hel_xpoint(s_tor>=s_res); %values of array of perturbed flux at phi=pi for the outner region

%Initialize the arrays of unperturbed toroidal and poloidal fluxes for inner and outer regions:
s_of_shel_in=zeros(size(s_in,2),size(phi_vpt,2));        % s(s_helical,phi) inner region
s_of_shel_out=zeros(size(s_out,2),size(phi_vpt,2));      % s(s_helical,phi) outer region
s_pol_of_shel_in=zeros(size(s_in,2),size(phi_vpt,2));    % s_pol(s_helical,phi) inner region
s_pol_of_shel_out=zeros(size(s_out,2),size(phi_vpt,2));  % s_pol(s_helical,phi) outner region

% Invert the depedence s_helical(s,phi). Here we obtain s=s(s_helical(s_x,pi),phi) for the equidistant
% grid of s_x where s_helical(s_x,pi/2) is a helical flux at the cross-section phi=pi which contains X-point
% We obtain also such a dependence for the unperturbed the poloidal flux, s_pol=s_pol(s_helical(s_x,pi),phi)
for i=1:1:size(phi_vpt,2)
  s_hel_in=s_hel_pert(s_tor<s_res,i);
  s_of_shel_in(:,i)=interp1(s_hel_in,s_in,s_hel_xpoint_in);
  s_pol_of_shel_in(:,i)=interp1(s_hel_in,s_pol_in,s_hel_xpoint_in);
  s_hel_out=s_hel_pert(s_tor>=s_res,i);
  s_of_shel_out(:,i)=interp1(s_hel_out,s_out,s_hel_xpoint_out,'linear','extrap');
  s_pol_of_shel_out(:,i)=interp1(s_hel_out,s_pol_out,s_hel_xpoint_out,'linear','extrap');
end

% Averages over the helical angle of unperturbed toroidal (s) and unperturbed poloidal (s_pol) fluxes as functions of helical flux:
avs_of_shel_in=sum(s_of_shel_in(:,2:end)')/(size(phi_vpt,2)-1);
avs_of_shel_out=sum(s_of_shel_out(:,2:end)')/(size(phi_vpt,2)-1);
avs_pol_of_shel_in=sum(s_pol_of_shel_in(:,2:end)')/(size(phi_vpt,2)-1);
avs_pol_of_shel_out=sum(s_pol_of_shel_out(:,2:end)')/(size(phi_vpt,2)-1);

% Derivative of toroidal flux (as function of helical flux and helical angle) over the average toroidal flux (which is a function of helical flux only):
ds_dsbar_in=zeros(size(s_in,2),size(phi_vpt,2));
ds_dsbar_out=zeros(size(s_out,2),size(phi_vpt,2));
for i=1:1:size(phi_vpt,2)
  ds_dsbar_in(:,i)=gradient(s_of_shel_in(:,i))./gradient(avs_of_shel_in');
  ds_dsbar_out(:,i)=gradient(s_of_shel_out(:,i))./gradient(avs_of_shel_out');
end

% Computation of:
% 1) perturbed helical Boozer angle phi_B as function of helical angle phi on perturbed flux surfaces (result - phib_in and phib_out),
% 2) of helical angle phi as function of perturbed helical Boozer angle phi_B (inverse to the above dependence, result - phi_of_phib_in and phi_of_phib_out)
% 3) old radial Boozer coordinate s as function of perturbed flux surface label sx and perturbed Boozer angle:
phib_in=zeros(size(s_in,2),size(phi_vpt,2));
phi_of_phib_in=zeros(size(s_in,2),size(phi_vpt,2));
s_of_phib_in=zeros(size(s_in,2),size(phi_vpt,2));
for i=1:1:size(ds_dsbar_in,1)
  phib_in(i,:)=cumtrapz(phi_vpt,ds_dsbar_in(i,:));
  phi_of_phib_in(i,:)=interp1(phib_in(i,:),phi_vpt,phi_vpt);
  s_of_phib_in(i,:)=interp1(phi_vpt,s_of_shel_in(i,:),phi_of_phib_in(i,:));
end
phib_out=zeros(size(s_out,2),size(phi_vpt,2));
phi_of_phib_out=zeros(size(s_out,2),size(phi_vpt,2));
s_of_phib_out=zeros(size(s_out,2),size(phi_vpt,2));
for i=1:1:size(ds_dsbar_out,1)
  phib_out(i,:)=cumtrapz(phi_vpt,ds_dsbar_out(i,:));
  phi_of_phib_out(i,:)=interp1(phib_out(i,:),phi_vpt,phi_vpt);
  s_of_phib_out(i,:)=interp1(phi_vpt,s_of_shel_out(i,:),phi_of_phib_out(i,:));
end

% compute iota_b (perturbed iota):
iotab_in=gradient(avs_pol_of_shel_in)./gradient(avs_of_shel_in);
iotab_out=gradient(avs_pol_of_shel_out)./gradient(avs_of_shel_out);

% Plot old and new iota:
plot(s_tor,iota_vpt,'b',s,s*0+iota_res,':',avs_of_shel_in,iotab_in,'r',avs_of_shel_out,iotab_out,'r'),xlim([0 1])
legend('\iota(s)','\iota_{res}','\iota_b(s_b)')

% Redefine perturbed normalized toroidal flux - close the gap at island location, rescale to the range [0:1]:
sbarrange_embedded_in=avs_of_shel_in(end)-avs_of_shel_in(1) %range of average toroidal flux for inner region
sbarrange_embedded_out=avs_of_shel_out(end)-avs_of_shel_out(1) %range of average toroidal flux for outer region
sbarrange=sbarrange_embedded_in+sbarrange_embedded_out %total range (will be rescaled to 1)

sbar_resc=[0:1:nsbar-1]/(nsbar-1); %rescaled label (changes between 0 and 1)

hsbar=sbarrange/(nsbar-1)                                 % step size over non-rescaled label without the gap
nsbar_in=fix(sbarrange_embedded_in*(nsbar-1)/sbarrange)+1 %number of points in the inner region
nsbar_out=nsbar-nsbar_in                                  %number of points in the outer region
sbar_equi=([1:1:nsbar]-1)*hsbar;
avs_equi_in=sbar_equi(1:nsbar_in);                                  %grid on original averged s in the inner region
avs_equi_out=sbar_equi(nsbar_in+1:end)+avs_of_shel_out(1)-avs_of_shel_in(end); %grid on original averged s in the outer region

% iota_b for the rescaled label:
iota_resc=zeros(size(sbar_resc));
iota_resc(1:nsbar_in)=interp1(avs_of_shel_in,iotab_in,avs_equi_in);
iota_resc(nsbar_in+1:end)=interp1(avs_of_shel_out,iotab_out,avs_equi_out);
figure(2)
plot(s_tor,iota_vpt,sbar_resc,iota_resc)
legend('\iota(s)','\iota_b(s_{rescaled})')

% Unperturbed s and phi as functions of rescaled (perturbed) label and perturbed Boozer phi:
s_resc=zeros(size(sbar_resc,2),size(phi_vpt,2));
phi_resc=zeros(size(sbar_resc,2),size(phi_vpt,2));
for i=1:1:size(phi_vpt,2)
  s_resc(1:nsbar_in,i)=interp1(avs_of_shel_in,s_of_phib_in(:,i),avs_equi_in);
  s_resc(nsbar_in+1:end,i)=interp1(avs_of_shel_out,s_of_phib_out(:,i),avs_equi_out);
  phi_resc(1:nsbar_in,i)=interp1(avs_of_shel_in,phi_of_phib_in(:,i),avs_equi_in);
  phi_resc(nsbar_in+1:end,i)=interp1(avs_of_shel_out,phi_of_phib_out(:,i),avs_equi_out);
end


% plot perturbed flux surfaces

figure(3)
%plot(phi_of_phib_in(end-3,:),s_of_phib_in(end-3,:),':')
%hold on
%for i=100:100:size(s_in,2)
%plot(phi_of_phib_in(i,:),s_of_phib_in(i,:),':')
%end
%plot(phi_of_phib_out(3,:),s_of_phib_out(3,:),':')
%for i=100:100:size(s_out,2)
%plot(phi_of_phib_out(i,:),s_of_phib_out(i,:),':')
%end

plot(phi_resc(1,:),s_resc(1,:))
hold on
nplot=300;
dnplot=nsbar/nplot;
for i=1:dnplot:nsbar
  plot(phi_resc(i,:),s_resc(i,:))
end
hold off

% \brief Make input(array) even.
%
% At the moment this is done by 'rounding' towards zero (fix).
function output = make_even(input)
  output = 2*fix(input/2);
end

% Igichine model: del_s = 4 sqrt(s_res) ka / e
% ka=0.06
function del_s = Igichine_model_del_s(s_res, ka)
  del_s = 4*sqrt(s_res)*ka/exp(1);
end
