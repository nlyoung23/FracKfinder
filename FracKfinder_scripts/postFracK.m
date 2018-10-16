clear all
close all

% =======================================================================
% = FracKFinder: A 3-D K Tensor-Determining Program for HydroGeoSphere ==
% =============== By Nathan Young and Jacqueline Reber ==================
% ========================Postprocessor Notes============================
% 
% This script reads the HSPLOT-extracted model output from the preprocessor
% script, and produces a 6-component hydraulic conductivity (K) tensor for 
% both the matrix and the fracture domain. Also provided are single values
% of K for the matrix and fracture domains. 
% The K tensors are computed using an adaptation of the methods 
% presented in Wang et al, (2005). We regress a vector of the mean 
% groundwater velocity for a given model run against the magnitude of 
% the gradient rotation for that model run. Symmetry of the off diagonals in
% the resulting 9-component K tensor is assumed, but in practice the
% magnitudes of the off-diagonal components are averaged in order to impose
% symmetry (Zhang et al., 2006). This script outputs a vector with the
% principle components of hydraulic conductivity (Kx, Ky, and Kz), as well
% as a second vector with the values of the off-diagonal components (Kxy,
% Kxz, Kyz). 

%%______user input variables______

%%Enter location of preprocessor folder
cd 'E:\Big_HGS_runs\K_tensor_in_prog925'

ne=0.298; %effective porosity 

% ===================================================
% ====NO USER INPUT REQUIRED AFTER THIS LINE=========
% ===================================================
% =================THE CODE==========================
% ===================================================

tic
PMfilename = dir('*o.pm.dat'); %list all the vtk files in your current directory.
frfilename = dir('*o.frac.dat'); %list all the vtk files in your current directory.
grads1 = csvread('headrotationmagnitudes.csv'); %gradient magnitudes

%Loading XYZ runs

%%create output filename list file
for a = 1:26
        if a<10
        textFileName = [ '0' num2str(a) 'o.frac' '.dat'];
	  	fid = fopen(textFileName, 'rt');
		data = textscan(fid, '%s', 'delimiter', '\n');
        else
        textFileName2 = [ num2str(a) 'o.frac' '.dat'];
		fid2 = fopen(textFileName2, 'rt');
		data = textscan(fid2, '%s', 'delimiter', '\n');
		end 
 
 %%find headers
 idx1 = find(strcmp(data{1}, '# x linear velocity (cell centred)'), 1, 'first');
 idx2 = find(strcmp(data{1}, '# y linear velocity (cell centred)'), 1, 'first');
 idx3 = find(strcmp(data{1}, '# z linear velocity (cell centred)'), 1, 'first');
 idx4 = find(strcmp(data{1}, '# br'), 1, 'first');
 %%create read startpoints
 idx1s = idx1+1;
 idx2s = idx2+1;
 idx3s = idx3+1;
 %%create read endpoints
 idx1e = idx2-2;
 idx2e = idx3-2;
 idx3e = idx4-2;
 
 Afrx{a} = mean(dlmread (frfilename(a).name,'',[idx1s 0 idx1e 0])); %matrix a containing values from vtk file_in_loadpath
 Afry{a} = mean(dlmread (frfilename(a).name,'',[idx2s 0 idx2e 0]));
 Afrz{a} = mean(dlmread (frfilename(a).name,'',[idx3s 0 idx3e 0]));
 
 fclose(fid);
end


for aa = 1:26
        if aa<10
        textFileName = [ '0' num2str(aa) 'o.pm' '.dat'];
	  	fid = fopen(textFileName, 'rt');
		data = textscan(fid, '%s', 'delimiter', '\n');
        else
        textFileName2 = [ num2str(aa) 'o.pm' '.dat'];
		fid = fopen(textFileName2, 'rt');
		data = textscan(fid, '%s', 'delimiter', '\n');
		end 
 
 %%find headers
 idx4 = find(strcmp(data{1}, '# x linear velocity (cell centred)'), 1, 'first');
 idx5 = find(strcmp(data{1}, '# y linear velocity (cell centred)'), 1, 'first');
 idx6 = find(strcmp(data{1}, '# z linear velocity (cell centred)'), 1, 'first');
 idx7 = find(strcmp(data{1}, '# br'), 1, 'first');
 %%create read startpoints
 idx4s = idx4+1;
 idx5s = idx5+1;
 idx6s = idx6+1;
 %%create read endpoints
 idx4e = idx5-2;
 idx5e = idx6-2;
 idx6e = idx7-2;
 
Apmx{aa} = mean(dlmread (PMfilename(aa).name,'',[idx4s 0 idx4e 0])); %matrix a containing values from vtk file_in_loadpath
Apmy{aa} = mean(dlmread (PMfilename(aa).name,'',[idx5s 0 idx5e 0]));
Apmz{aa} = mean(dlmread (PMfilename(aa).name,'',[idx6s 0 idx6e 0]));
 
fclose(fid);
end


for b = 27:52
        if b<10
        textFileName = [ '0' num2str(b) 'o.frac' '.dat'];
	  	fid = fopen(textFileName, 'rt');
		data = textscan(fid, '%s', 'delimiter', '\n');
        else
        textFileName2 = [ num2str(b) 'o.frac' '.dat'];
		fid2 = fopen(textFileName2, 'rt');
		data = textscan(fid2, '%s', 'delimiter', '\n');
		end 
 
 %%find headers
 idx1 = find(strcmp(data{1}, '# x linear velocity (cell centred)'), 1, 'first');
 idx2 = find(strcmp(data{1}, '# y linear velocity (cell centred)'), 1, 'first');
 idx3 = find(strcmp(data{1}, '# z linear velocity (cell centred)'), 1, 'first');
 idx4 = find(strcmp(data{1}, '# br'), 1, 'first');
 %%create read startpoints
 idx1s = idx1+1;
 idx2s = idx2+1;
 idx3s = idx3+1;
 %%create read endpoints
 idx1e = idx2-2;
 idx2e = idx3-2;
 idx3e = idx4-2;
 
 Bfrx{b} = mean(dlmread (frfilename(b).name,'',[idx1s 0 idx1e 0])); %matrix a containing values from vtk file_in_loadpath
 Bfry{b} = mean(dlmread (frfilename(b).name,'',[idx2s 0 idx2e 0]));
 Bfrz{b} = mean(dlmread (frfilename(b).name,'',[idx3s 0 idx3e 0]));
 
 fclose(fid);
end

for bb = 27:52
        if bb<10
        textFileName = [ '0' num2str(bb) 'o.pm' '.dat'];
	  	fid = fopen(textFileName, 'rt');
		data = textscan(fid, '%s', 'delimiter', '\n');
        else
        textFileName2 = [ num2str(bb) 'o.pm' '.dat'];
		fid2 = fopen(textFileName2, 'rt');
		data = textscan(fid2, '%s', 'delimiter', '\n');
		end 
 
%  find headers
 idx4 = find(strcmp(data{1}, '# x linear velocity (cell centred)'), 1, 'first');
 idx5 = find(strcmp(data{1}, '# y linear velocity (cell centred)'), 1, 'first');
 idx6 = find(strcmp(data{1}, '# z linear velocity (cell centred)'), 1, 'first');
 idx7 = find(strcmp(data{1}, '# br'), 1, 'first');
%  create read startpoints
 idx4s = idx4+1;
 idx5s = idx5+1;
 idx6s = idx6+1;
%  create read endpoints
 idx4e = idx5-2;
 idx5e = idx6-2;
 idx6e = idx7-2;
 
Bpmx{bb} = mean(dlmread (PMfilename(bb).name,'',[idx4s 0 idx4e 0])); %matrix a containing values from vtk file_in_loadpath
Bpmy{bb} = mean(dlmread (PMfilename(bb).name,'',[idx5s 0 idx5e 0]));
Bpmz{bb} = mean(dlmread (PMfilename(bb).name,'',[idx6s 0 idx6e 0]));
 
 fclose(fid);
end

for c = 53:78
        if c<10
        textFileName = [ '0' num2str(c) 'o.frac' '.dat'];
	  	fid = fopen(textFileName, 'rt');
		data = textscan(fid, '%s', 'delimiter', '\n');
        else
        textFileName2 = [ num2str(c) 'o.frac' '.dat'];
		fid2 = fopen(textFileName2, 'rt');
		data = textscan(fid2, '%s', 'delimiter', '\n');
		end 
 
%  find headers
 idx1 = find(strcmp(data{1}, '# x linear velocity (cell centred)'), 1, 'first');
 idx2 = find(strcmp(data{1}, '# y linear velocity (cell centred)'), 1, 'first');
 idx3 = find(strcmp(data{1}, '# z linear velocity (cell centred)'), 1, 'first');
 idx4 = find(strcmp(data{1}, '# br'), 1, 'first');
%  create read startpoints
 idx1s = idx1+1;
 idx2s = idx2+1;
 idx3s = idx3+1;
%  create read endpoints
 idx1e = idx2-2;
 idx2e = idx3-2;
 idx3e = idx4-2;
 
 Cfrx{c} = mean(dlmread (frfilename(c).name,'',[idx1s 0 idx1e 0])); %matrix a containing values from vtk file_in_loadpath
 Cfry{c} = mean(dlmread (frfilename(c).name,'',[idx2s 0 idx2e 0]));
 Cfrz{c} = mean(dlmread (frfilename(c).name,'',[idx3s 0 idx3e 0]));
 
 fclose(fid);
end

for cc = 53:78
        if cc<10
        textFileName = [ '0' num2str(cc) 'o.pm' '.dat'];
	  	fid = fopen(textFileName, 'rt');
		data = textscan(fid, '%s', 'delimiter', '\n');
        else
        textFileName2 = [ num2str(cc) 'o.pm' '.dat'];
		fid2 = fopen(textFileName2, 'rt');
		data = textscan(fid2, '%s', 'delimiter', '\n');
		end 
 
%  find headers
 idx4 = find(strcmp(data{1}, '# x linear velocity (cell centred)'), 1, 'first');
 idx5 = find(strcmp(data{1}, '# y linear velocity (cell centred)'), 1, 'first');
 idx6 = find(strcmp(data{1}, '# z linear velocity (cell centred)'), 1, 'first');
 idx7 = find(strcmp(data{1}, '# br'), 1, 'first');
%  create read startpoints
 idx4s = idx4+1;
 idx5s = idx5+1;
 idx6s = idx6+1;
%  create read endpoints
 idx4e = idx5-2;
 idx5e = idx6-2;
 idx6e = idx7-2;
 
Cpmx{cc} = mean(dlmread (PMfilename(cc).name,'',[idx4s 0 idx4e 0])); %matrix a containing values from vtk file_in_loadpath
Cpmy{cc} = mean(dlmread (PMfilename(cc).name,'',[idx5s 0 idx5e 0]));
Cpmz{cc} = mean(dlmread (PMfilename(cc).name,'',[idx6s 0 idx6e 0]));
 
 fclose(fid);
end

% Editing and reshaping output

% extract porous media velocities in the x direction
vxi= cell2mat(Apmx); 
vxii= cell2mat(Bpmx);
vxiii = cell2mat(Cpmx);

% extract porous media velocities in the y direction
vyi= cell2mat(Apmy); 
vyii= cell2mat(Bpmy);
vyiii = cell2mat(Cpmy);

%extract porous media velocities in the z direction
vzi= cell2mat(Apmz); 
vzii= cell2mat(Bpmz);
vziii = cell2mat(Cpmz);


%extract fracture velocities in the x direction
fx1= cell2mat(Afrx); 
fx2= cell2mat(Bfrx);
fx3 = cell2mat(Cfrx);

%extract fracture velocities in the y direction
fy1= cell2mat(Afry); 
fy2= cell2mat(Bfry);
fy3 = cell2mat(Cfry);

%extract fracture velocities in the z direction
fz1= cell2mat(Afrz); 
fz2= cell2mat(Bfrz);
fz3 = cell2mat(Cfrz);

vx = horzcat(vxi,vxii,vxiii); %combine x value porous media output into a single vector
vy = horzcat(vyi,vyii,vyiii); %as above for y
vz = horzcat(vzi,vzii,vziii); %as above for z

vfx = horzcat(fx1,fx2,fx3); %combine x value fracture output into a single vector
vfy = horzcat(fy1,fy2,fy3); %as above for y
vfz = horzcat(fz1,fz2,fz3); %as above for z

v = vertcat(vx,vy,vz); %vertically stack above vectors to show average x y and z velocity
vf = vertcat(vfx,vfy,vfz);

%_______K tensor computation____________

v1=v.*ne; %Multiply matrix groundwater velocities by porosity to get K values. 
v2 = v1.'; %Transpose matrix K values for OLS calculation.
vf2=vf.'; %Transpose fracture K values for OLS calculation. Porosity of 
          %fractures in HGS is 1. Fracture aperture is used to convert 
          %from fracture velocity to fracture K from the equations
          %presented in the HydroGeoSphere Manual (Aquanty 2015).
grads11 = grads1.'; %transpose gradients for OLS calculation.
%_____________OLS computation_______________


[b,se_b,mse] = lscov(grads11,v2); 
pm_K_tensor = b;
frac_K_tensor = lscov(grads11,vf2);
single_pm_k=norm(pm_K_tensor)
single_frac_K=norm(frac_K_tensor);

kxy = (pm_K_tensor(2,1)+pm_K_tensor(1,2))/2;
kyz = (pm_K_tensor(3,1)+pm_K_tensor(1,3))/2;
kxz = (pm_K_tensor(3,2)+pm_K_tensor(2,3))/2;

kx = pm_K_tensor(1,1);
ky = pm_K_tensor(2,2);
kz = pm_K_tensor(3,3);


matrix_component_values = [kx,ky,kz]'
matrix_off_diagonals = [kxy,kyz,kxz]'


pm_K_tensor(1,2) = kxy;
pm_K_tensor(2,1) = kxy;
pm_K_tensor(1,3) = kyz;
pm_K_tensor(3,1) = kyz;
pm_K_tensor(3,2) = kxz;
pm_K_tensor(2,3) = kxz;

kfxy = (frac_K_tensor(2,1)+frac_K_tensor(1,2))/2;
kfyz = (frac_K_tensor(3,1)+frac_K_tensor(1,3))/2;
kfxz = (frac_K_tensor(3,2)+frac_K_tensor(2,3))/2;

kfx = frac_K_tensor(1,1);
kfy = frac_K_tensor(2,2);
kfz = frac_K_tensor(3,3);

frac_K_tensor(1,2) = kfxy;
frac_K_tensor(2,1) = kfxy;
frac_K_tensor(1,3) = kfyz;
frac_K_tensor(3,1) = kfyz;
frac_K_tensor(3,2) = kfxz;
frac_K_tensor(2,3) = kfxz;

fracture_component_values = [kfx,kfy,kfz]'
fracture_diagonals = [kfxy,kfyz,kfxz]'

toc
%%%write an output file that collates tensors and single k values