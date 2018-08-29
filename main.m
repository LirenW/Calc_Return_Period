%%
% 
%  This Program is for calculating the return period at TWO specific
%    return level. Using PWMs method. Here we take 30-year return 
%    period as an example.
%  
%  Created by Liren
%  May 2018
% 
clc;
clear;

%% SECTION 1: Setting module
% Set for the input information
flag      = true
fiFutureA = './';  % Two input field
fiFutureB = './';
fiRef     = './';
outFile   = './example.nc';
varName   = 'var';

lon     = ncread(fiRef, 'lon');
lat     = ncread(fiRef, 'lat');

FutureA   = ncread(fiFutureA, varName);
FutureB   = ncread(fiFutureB, varName);
FutureRef = ncread(fiRef    , varName);
nlon      = length(lon);
nlat      = length(lat);

%% Calculation Module
% reyears = [10 20 30 50 100];
reyears        = [30];
nreyrs         = length(reyears);

[alpha, zeta]  = PWMs(FutureB, flag);
for iyr = 1: nreyrs;
    X          = calcReturnPeriod(FutureRef, reyears(iyr));
    TMP        = -((X - zeta) ./ alpha);
    F          = exp(-exp(TMP));

    RP_A(iyr,: ,: ) = 1 ./ (1 - F);
end

[alpha, zeta]  = PWMs(FutureA, flag);
for iyr = 1: nreyrs;
    X         = calcReturnPeriod(FutureRef, reyears(iyr));
    TMP       = -((X - zeta) ./ alpha);
    F         = exp(-exp(TMP));

    RP_B(iyr,: ,: ) = 1 ./ (1 - F);
end

ind       = find(isnan(RP_B) == 1);
RP_B(ind) = 1.e20;
ind       = find(isnan(RP_A) == 1);
RP_A(ind) = 1.e20;

%% Output Module
ncid   = netcdf.create(outFile, 'NETCDF4');
% Global Attributions setting
varidGlo = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncid , varidGlo , 'creation_date' , datestr(now));
if flag == true
    netcdf.putAtt(ncid , varidGlo , 'estimate_method' , 'PWMs');
else
    netcdf.putAtt(ncid , varidGlo , 'estimate_method' , 'Regular moments estimate method');
end

dimidx = netcdf.defDim(ncid, 'lat'   , nlat  );
dimidy = netcdf.defDim(ncid, 'lon'   , nlon  );
dimidn = netcdf.defDim(ncid, 'return', nreyrs);

varid1 = netcdf.defVar(ncid, 'return' ,'double' ,dimidn);
varid2 = netcdf.defVar(ncid, 'lat'    ,'double' ,dimidx);
varid3 = netcdf.defVar(ncid, 'lon'    ,'double' ,dimidy);

netcdf.putAtt(ncid, varid1, 'units', 'year'         );
netcdf.putAtt(ncid, varid2, 'units', 'degrees_north');
netcdf.putAtt(ncid, varid3, 'units', 'degrees_east' );

varidA = netcdf.defVar(ncid, 'returnA', 'double', [dimidn dimidy dimidx]);
varidB = netcdf.defVar(ncid, 'returnB', 'double', [dimidn dimidy dimidx]);
netcdf.defVarFill(ncid , varidA , false , 1.e20); % Must under netCDF4
netcdf.defVarFill(ncid , varidB , false , 1.e20); % Must under netCDF4

netcdf.endDef(ncid);

netcdf.putVar(ncid, varid1, reyears);
netcdf.putVar(ncid, varid2, lat    );
netcdf.putVar(ncid, varid3, lon    );
netcdf.putVar(ncid, varidA, RP_B   );
netcdf.putVar(ncid, varidB, RP_A   );

netcdf.close(ncid);
