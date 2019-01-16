%% part 1
%% read phase shifted holograms
I1 = double(imread('holo1.bmp'));
I2 = double(imread('holo2.bmp'));
I3 = double(imread('holo3.bmp'));
I4 = double(imread('holo4.bmp'));
I5 = double(imread('holo5.bmp'));

pha = atan2(2.*(I2-I4), (2.*I3-I1-I5));
amp = sqrt((2.*(I2-I4)).^2+(2.*I3-I1-I5).^2);

%% display obtained phase and amplitude
% f1 = figure('Name','Input phase');
% imagesc(pha)
% saveas(gcf, 'input_pha.png')
% f2 = figure('Name','Input amplitude');
% imagesc(amp)
% saveas(gcf, 'input_amp.png')

%% find plane
%definition of parameters for propagation and further stpes
inp_field = amp.*exp(1i*pha); % input field
lamb =  0.6328; % used wavelength
pix_size = 3.45/10;
refr_ind = 1;

%% iterations for manual analysis
% %definition of parameters
% no_it = 10; %number of iterations
% step = 5; %step size [um]
% minimal = 95; %starting value for scanning [um]
% f3 = figure('Name','Manual analysis of planes');
% 
% for i=1:no_it
%     sol = propagation(inp_field, minimal+(i*step),lamb,refr_ind, pix_size); %propagate to selected dostance
%     %display consecutive planes in single window (used for ten planes only for clarity)
%     subplot(2,5,i)
%     imagesc(abs(sol))
%     title(strcat(num2str(minimal+(i*step)), ' um'))
% end
% saveas(gcf, 'compare.png')


%% iterations for standard deviation and search for AOI
%definition of parameters
no_it = 20; %number of iterations
step = 5; %step size [um]
minimal = 70; %starting value for scanning [um]
stdAoi = zeros(1,no_it); %vector of standard veriations
aoiX = 270:600;
aoiY = 360:690;

% %find AOI
% f4 = figure('Name', 'AOI');
% sol = propagation(inp_field, 0,lamb,refr_ind, pix_size); %propagate to selected dostance 
% aoi = sol(aoiX, aoiY);
% imagesc(abs(aoi))
% saveas(gcf, 'aoi.png')

for i=1:no_it
    sol = propagation(inp_field, minimal+(i*step),lamb,refr_ind, pix_size); %propagate to selected dostance 
    %calculate standard deviation of single lens
    aoi = sol(aoiX, aoiY);
    stdAoi(i) = std2(abs(aoi));
end

% %display plot of std dev
% f5 = figure('Name', 'Plot of std dev');
% plot(minimal+step:step:(minimal+no_it*step),stdAoi)
% xlabel('dist[um]')
% ylabel('std. dev. [a.u.]')
% title('Standard deviation for different planes')
% saveas(gcf, 'std_dev_plot.png')

%% display in-focus image
% f6 = figure('Name', 'In-focus field (amplitude)');
focus_dist = minimal + step*find(stdAoi == min(stdAoi));
sol = propagation(amp.*exp(1i*pha), focus_dist, 0.6328,1, 3.45/10);
% imagesc(abs(sol))
% saveas(gcf, 'output_amp.png')
% 
% imagesc(angle(sol))
% saveas(gcf, 'output_pha.png')


%% removal of the background
%figure; imagesc(pha); title('phase');
maxDistance = 0.01; %maximal deviation from plane
[X,Y] = meshgrid(1:size(pha,1), 1:size(pha,2));
X = X.';
Y = Y.';
pha = angle(sol);
plane = fit([X(:), Y(:)],pha(:) ,'poly11');
meanPhase = mean(mean(plane(X,Y)));
pha = pha - plane(X, Y);
%figure; imagesc(pha); title('phase without angled background');

%phase unwrapping
dx = 1;
lambda = 1;
[unph, dIn] = TIE_unwrap2DCdct2(pha, dx,lambda);
figure; imagesc(unph); title('Unwrapped phase');
unph = unph - min(unph);
unph = uint8(unph);
imwrite(unph, 'unph_focused.bmp', 'bmp')

%% part 2
% calls c++ code for circle detection
system('circ_detection\tomo_template.exe');

%% part 3
%% variables for table
Name = string.empty();%vector with name of the parameter
Unit = string.empty();%vector with unit of the parameter
Value = double.empty();%vector with the parameter value
StdDev = double.empty();%vector with standar deviation of the parameter

%% draw obtained circles
figure; imagesc(unph);
load('save.txt') % file containing detected circles
ds = 3.45/10; %pixel size in reality
viscircles(save(:,1:2), save(:,3)); title('Detected circles');

%% binarization
% BW = imbinarize(unph);
% figure; imagesc(BW); title('Binarization');
% figure; imshowpair(unph,BW,'montage');

%% segmentation
% L = bwlabel(BW);
% figure; imagesc(L);title('Segmentation');
 
%% find and remove circles that are on the edges
edgePercent = 15;
%find elements that are in the vicinity of edges
removeX = find(save(:,1) <= size(pha,2)*edgePercent/100 |  save(:,1) >= size(pha,2)*(100-edgePercent)/100);
removeY = find(save(:,2) <= size(pha,1)*edgePercent/100 |  save(:,2) >= size(pha,1)*(100-edgePercent)/100);
%remove selected circles
save(removeX, :) = [];
save(removeY, :) = [];
figure; imagesc(unph);
viscircles(save(:,1:2), save(:,3)); title('After removing edge elements');

%% calculate diameter and its std dev
r_mat = save(:,3)*ds; % temporary variable storing all radius'
r = mean(r_mat);
r_std = std(r_mat);
d = 2*r;
d_std = 2* r_std;
Name = [Name 'diameter'];
Unit = [Unit 'um'];
Value = [Value d];
StdDev = [StdDev d_std];

%% calculate h and its std dev (TEA)
n0 = 1; % refractive index of environment
ns = 1.45701; % refractive index of material
h_mat = zeros(size(save(:,3))); % vector storing final values of height
tmp_h_part = lamb/(2*pi*(ns-n0)); % repeating constant part of calculation
% calculations done for each radius
for i = 1:size(save(:,3),1)
    tmp_r = save(i,3);
    h_cent = double(unph(save(i,2),save(i,1)))*tmp_h_part;
    h_edge1 = double(unph(save(i,2)-tmp_r,save(i,1)))*tmp_h_part;
    h_edge2 = double(unph(save(i,2)+tmp_r,save(i,1)))*tmp_h_part;
    h_edge3 = double(unph(save(i,2),save(i,1)+tmp_r))*tmp_h_part;
    h_edge4 = double(unph(save(i,2),save(i,1)-tmp_r))*tmp_h_part;
    h_edge = mean([h_edge1, h_edge2, h_edge3, h_edge4]);
    h_mat(i) = h_cent; - h_edge;
end
h = mean(h_mat);
h_std = std(h_mat);
Name = [Name 'lens sag'];
Unit = [Unit 'um'];
Value = [Value h];
StdDev = [StdDev h_std];

%% calculate ROC and its std dev
ROC_mat = h_mat./2 + r_mat.*r_mat./(2*h_mat);
ROC = mean(ROC_mat);
ROC_std = std(ROC_mat);
Name = [Name 'ROC'];
Unit = [Unit 'um'];
Value = [Value ROC];
StdDev = [StdDev ROC_std];

%% calculate focal length and its std dev
f_mat = ROC_mat/(ns-1);
f = mean(f_mat);
f_std = std(f_mat);
Name = [Name 'focal length'];
Unit = [Unit 'um'];
Value = [Value f];
StdDev = [StdDev f_std];

%% calculate NA and its std dev
NA_mat = h_mat./r_mat;
NA = mean(NA_mat);
NA_std = std(NA_mat);
Name = [Name 'NA'];
Unit = [Unit 'unitless'];
Value = [Value NA];
StdDev = [StdDev NA_std];

%% calculate spherical error and its std dev
% NA_mat = h_mat./r_mat;
% NA = mean(NA_mat);
% NA_std = std(NA_mat);

%% create table with output data
outputParams = table(Name', Unit', Value', StdDev')
outputParams.Properties.VariableNames = {'Name' 'Unit' 'Value' 'StdDev'}