clear; close all

data_type='AMSR2';

if data_type=='SSMI2' % SSMI NT2 25 km
   year_start=1992;
   year_stop=2007;
   T=40;
   max_gap = 100; 
elseif data_type=='AMSRH' % NT2 analysis 12.5 km
   year_start=2002;
   year_stop=2011;
   T=80;
   max_gap = 100; 
elseif data_type=='AMSEB' % ASI analysis 6.25 km
   grid_res=6.25;
   year_start=2002;
   year_stop=2011;
   T=160;
   max_gap = 100; 
elseif data_type=='AMSR2' % ASI analysis 3.125 km
   year_start=2012;
   year_stop=2018;
   T=320; % smoothing parameters  
   max_gap = 100; 
elseif data_type=='5km5d' % Mike's model output
   year_start=2009;
   year_stop=2009;
   T=160; % smoothing parameters
   max_gap = 100; 
else 
   error('Unknown data type');
end


SDtime=datenum(year_start,1,1):1:datenum(year_stop,12,31);
%nsect=NaN*ones(length(SDtime),1); int_dd=nsect; int_ddsm=nsect;

% smoothing parameters
Wn=[1/T];   % lowpass
[b,a]=butter(2,Wn,'low');

SLAT=70; SLON=0; HEMI='s';


for i=1:length(SDtime)
    xsm=[];ysm=[];
    AE(i).SDtime=SDtime(i);
    AE(i).nsect=[]; AE(i).x=[]; AE(i).y=[];
    AE(i).lon=[]; AE(i).lat=[];
    AE(i).xsm=[]; AE(i).ysm=[];
    %[IE,S]=find_main_ice_edge(SDtime1);
    [IE]=find_main_ice_edge(SDtime(i),data_type,max_gap);
%    if(~isempty(IE.nsect)); nsect(i)=IE.nsect; end

    if(~isempty(IE.x));
        xsm=filtfilt(b,a,IE.x); % smooth the iceedge
        ysm=filtfilt(b,a,IE.y);

        dx   = IE.x(2:end)-IE.x(1:(end-1)); dxsm=xsm(2:end)-xsm(1:(end-1));
        dy   = IE.y(2:end)-IE.y(1:(end-1)); dysm=ysm(2:end)-ysm(1:(end-1));
        dd   = sqrt(dx.^2+dy.^2);
        ddsm = sqrt(dxsm.^2+dysm.^2);
        int_dd(i)=sum(dd); int_ddsm(i)=sum(ddsm);
        AE(i).nsect=IE.nsect;
        AE(i).x=IE.x; AE(i).y=IE.y;
        AE(i).xsm=xsm; AE(i).ysm=ysm;
        %AE(i).lon=IE.lon; AE(i).lat=IE.lat;
        AE(i).rawlen=int_dd(i);
        AE(i).filtlen=int_ddsm(i);

    end
end

META.SLAT=SLAT; META.SLON=SLON; META.HEMI=HEMI; META.filtlen=T; META.data_type=data_type;
META.year_start=year_start;META.year_stop=year_stop;

data_dir=['/net/esrdata1/springer/data/IceEdge/data_',META.data_type,'/'];
data_file = fullfile(data_dir,[META.data_type,'_iceedge','_f',num2str(META.filtlen),'.mat'])
save(data_file,'AE','META');

make_plots=0;
plot_iceedge_finder(META.data_type,META.filtlen,make_plots);

