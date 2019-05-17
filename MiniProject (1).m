clear all;
clc;

%DirImages = '\\fs2\18234602\MATLAB\MiniProject\30Kph';
%DirImages = '\\fs2\18234602\MATLAB\MiniProject\50Kph';
%DirImages = '\\fs2\18234602\MATLAB\MiniProject\80Kph';
DirImages = '\\fs2\18234602\MATLAB\MiniProject\20Kph';
%DirImages = '\\fs2\18234602\MATLAB\MiniProject\stress';
%DirImages = '\\fs2\18234602\MATLAB\MiniProject\100Kph';
%DirImages = '\\fs2\18234602\MATLAB\MiniProject\GoldStandard2';


% just run the code by giving the accuarcy function various directories in
% the development dataset. To run the stress dataset accuracy just call the
% function below. accuarcyStress()
accuarcy(DirImages,'*.jpg',20)
accuarcyStress()


% calculates the accuracy of the classification of the test folders, The
% directory of the images is given as a parameter along with the labels in
% the directory each prediction is compared to the coressponding label and
% the equivalent accuracy of that folder is returned.
function acc = accuarcy(dir,extension,label)
 file = GetFileFromDirectory(dir,extension);
 tempArray = zeros(1,length(file));
for i = 1:length(file)
    % Get file name and file path
    fileName = file(i).name;
    filePath = fullfile(dir, fileName);
    image = imread(filePath);
    cropped = ExtractSign(image);
    prediction = similarity(cropped);
    
    if prediction == label
         tempArray(i) = 1;
    else
        tempArray(i) = 0;
    end
acc = (sum(tempArray)/length(tempArray))*100;
end
end

% This function uses the histogram of oriented gradients to convert the
% unkown image to a feature vector representation, This vector is then
% compared to the features of each of the images in the gold standard
% datset by using the Euclidian distance metric. The image with the lowest
% distance is then classified as the output.
function ans = similarity(image)
GoldStandardImages = '\\fs2\18234602\MATLAB\MiniProject\GoldStandard';
file = GetFileFromDirectory(GoldStandardImages,'*.jpg');
reference = HOG(image);
Array = zeros(1,length(file));
x = {};
x{1} = [0,0.0132649352752722,0.322612585935337,0.211607450608346,0.612226920875000,0.0633013464361293,0.210024859656076,0.608068751628225,0,0,0.0352604997051313,0.496366109665864,0.334248538834376,0.173714260457325,0.102898634069911,0.444770996453633,0.571849313245941,0,0,0.239452657020994,0.530365977259560,0.271015696432369,0.194115003123940,0.123726364051834,0.517182443066565,0.421572372979435,0,0,0.220108484100393,0.0934504648286556,0.407273145141005,0.722778619955329,0.0890963358709809,0.176682555445222,0.299130526119750,0,0,0.0984202185035984,0.201222042529487,0.673609982182292,0.380425760435909,0.146967445064113,0.292218766223402,0.305259787891604,0,0,0.245735887336606,0.365344568452201,0.436771637818230,0.295076437799749,0.144856577184978,0.270356365575536,0.541531878005713,0,0,0.505150729342836,0.125307576277766,0.242827921054228,0.651609016655234,0.101315253761357,0.364468686391628,0.0335069160047839,0,0,0.481980079334707,0.297436471626755,0.392355518307668,0.221588016544370,0.0919124735869692,0.588713563720344,0.175031294605209,0,0,0.331284341572268,0.490891217703830,0.214193896345819,0.116052485928374,0.0581721796785302,0.592507384757162,0.371596962239038,0];
x{2} = [0,0.0214888865015625,0.327189647364463,0.203904268387166,0.613046997612779,0.0547509755857816,0.178112441958145,0.609109444009954,0,0,0.0685540597017751,0.511875467583866,0.293702280915149,0.186640179720033,0.109548919385609,0.395732574225778,0.594377782953173,0,0,0.259115774225366,0.525578253697792,0.233584707625602,0.198955308511149,0.121185389971296,0.509304848463727,0.429754929325028,0,0,0.285509115871384,0.176584155451684,0.226848965312167,0.726693084183203,0.146085953000382,0.172963781982451,0.357320894920174,0,0,0.241958273585207,0.359762014672484,0.381077776006337,0.337013964400631,0.291250280351433,0.353452055009581,0.221027574101947,0,0,0.375328367065462,0.294969379504912,0.318942397280141,0.371689928537498,0.274897883264464,0.290751883896775,0.303291145792960,0,0,0.589265563344597,0.175375162357956,0.0667491167831127,0.615543441335458,0.179312595789679,0.330301741841141,0.0304722759184263,0,0,0.598205776349224,0.366929556178447,0.117134101695348,0.156166875754405,0.258438502829900,0.528929079699469,0.0611636979876581,0,0,0.444415528101808,0.440370818653580,0.148936989166881,0.223204255434406,0.226938250307594,0.519812263030372,0.243460889566548,0];
x{3} = [0,0.0267136042252374,0.310694995434817,0.125306112961451,0.740916426894493,0.0513638546068019,0.103033764498136,0.507712397021321,0,0,0.0753304687947961,0.560645023462291,0.158367177581075,0.313672502887396,0.0746352836754753,0.356350176029029,0.562650472261738,0,0,0.184805736004884,0.631190058305390,0.104527513946766,0.0971033034210594,0.0714640824611557,0.547680094320548,0.356368614343334,0,0,0.272061137368978,0.108215194783124,0.0975436315459885,0.840802218865059,0.135015122341374,0.123105584483288,0.247450293826557,0,0,0.284822311142510,0.351893781234145,0.176683986919886,0.471196790868575,0.230895835270461,0.365303499763615,0.202172718767287,0,0,0.272554862478667,0.423304092404950,0.174047274562868,0.252170255699107,0.194306622327607,0.460699625453000,0.255442523387194,0,0,0.583466273906949,0.131106048036589,0.0464141952495906,0.645294629960077,0.166217015262897,0.351983805744437,0.0286874795241638,0,0,0.643665776203313,0.332152252504785,0.0731681435707297,0.131024306469534,0.227452656007440,0.534047683559770,0.0629014484145491,0,0,0.412903480909392,0.475688706913965,0.138648075555780,0.199021557227216,0.195668428154240,0.568781263412198,0.185635751642401,0];
x{4} = [0,0.0307802576128280,0.326204724535627,0.164697475481563,0.653457536300507,0.0489050206962803,0.128389909355892,0.576373420517018,0,0,0.0801059032146968,0.550826732546074,0.184172174935608,0.161069844616763,0.0734854519375262,0.362212195366567,0.585840606873895,0,0,0.179338309736127,0.565532650284533,0.112255153747367,0.0833571894264070,0.0733795864549287,0.544560258424840,0.402871203330310,0,0,0.250302123514192,0.0850089726377821,0.155598913163374,0.746421539369893,0.107318334835508,0.117419108178983,0.358860393245529,0,0,0.240878363737532,0.210774031156247,0.213260319076135,0.293655322269182,0.171420084957369,0.227732357026231,0.192239060072635,0,0,0.251622476556370,0.248099366761062,0.153102318136176,0.168561688850088,0.144109547943989,0.267059476809337,0.239111154396481,0,0,0.551543375578086,0.112626421409700,0.0613263561953992,0.645679608460098,0.134053497016255,0.322003649112805,0.0498755637080813,0,0,0.575762921325632,0.285646701454495,0.0798053277699548,0.138869252181898,0.138062798084978,0.513497878366070,0.0935427066014250,0,0,0.398215851068013,0.441424424015784,0.0788376218238366,0.0948947911067221,0.0996864974669069,0.495055780502107,0.176134778697409,0];
x{5} = [0,0.0370368236705146,0.350403978084771,0.101477582034587,0.674240979379953,0.0693979833381442,0.161565821543376,0.547276484457606,0,0,0.0885366456176351,0.483631015015018,0.110963595928570,0.257047434156162,0.0856386089643338,0.340434340381124,0.616098253557075,0,0,0.107158349685964,0.487320535293412,0.151853491155268,0.252489889752135,0.0614729282835442,0.542046527499783,0.417465785164464,0,0,0.311189590050118,0.0957640184356492,0.0714474826509781,0.808744362060331,0.0981247513115214,0.126429714562314,0.260818272762198,0,0,0.252780771084058,0.0984878955610213,0.102573929544526,0.473083240295084,0.138198830395342,0.129405009359785,0.105207556624937,0,0,0.170068460614886,0.150890220340726,0.141033595095778,0.406877537641205,0.129448143603869,0.146752895032281,0.128141714526605,0,0,0.576269437030294,0.100399287759448,0.0183721219325931,0.692267102443776,0.0932058366181916,0.297132318120605,0.0175003709852942,0,0,0.665752340090600,0.251118812274126,0.0296976841759001,0.313516050671929,0.104020459274848,0.417716030458119,0.0479821786328498,0,0,0.473044520351064,0.494925912827014,0.0431717057745832,0.247476109216594,0.156884965061332,0.482297395812816,0.0732395881154570,0];

for i = 1:length(file)
    fileName = file(i).name;
    filePath = fullfile(GoldStandardImages, fileName);
    tempimage = imread(filePath);
    Gold = ExtractSign(tempimage);
    %GoldHog = HOG(Gold);
    GoldHog = x{i};
    Array(i) = sqrt(sum((GoldHog(:)-reference(:)).^ 2));
end
[value4, idx] = min(Array);
if idx == 1
    %ans = 100
    ans = 20;
elseif idx == 2
    %ans = 20
    ans = 30;
elseif idx == 3
    %ans = 30;
    ans = 50;
elseif idx == 4
    %ans = 50
    ans = 80;
elseif idx == 5
    %ans = 80
    ans = 100;
end
end

% calculates the accuracy of the classification on detected images in the
% stress datset. Using again the HOG feature represtation and the Euclidian
% distance metric to find the similarity.
function ACC = accuarcyStress()
DirImages = '\\fs2\18234602\MATLAB\MiniProject\stress';
file = GetFileFromDirectory(DirImages,'*.TIF');
predictionArray = zeros(1,length(file));
labelsOFImage = [40 40 30 30 20 20 20 40 40 40 40 40 40 40 40 20];
for i = 1:length(file)
    % Get file name and file path
    fileName = file(i).name;
    filePath = fullfile(DirImages, fileName);
    [X,cmap] = imread(filePath);
    cropped = redcars(X);
    prediction = similarity2(cropped);

    if prediction == labelsOFImage(i)
        predictionArray(i) = 1;
    end
end
ACC = (sum(predictionArray)/length(predictionArray))*100
end

% This function also uses the histogram of oriented gradients to convert the
% unkown image to a feature vector representation, This vector is then
% compared to the features of each of the images in the gold standard
% datset by using the Euclidian distance metric. The image with the lowest
% distance is then classified as the output. However it uses a smaller gold
% standard of only the three images found in the stress dataset 20,30,40 to
% classify the sign which of course increased the accuracy of the model.
function ans = similarity2(image)
GoldStandardImages = '\\fs2\18234602\MATLAB\MiniProject\GoldStandard2';
file = GetFileFromDirectory(GoldStandardImages,'*.jpg');
reference = HOG(image);
x = {};
x{1} = [0,0.0132649352752722,0.322612585935337,0.211607450608346,0.612226920875000,0.0633013464361293,0.210024859656076,0.608068751628225,0,0,0.0352604997051313,0.496366109665864,0.334248538834376,0.173714260457325,0.102898634069911,0.444770996453633,0.571849313245941,0,0,0.239452657020994,0.530365977259560,0.271015696432369,0.194115003123940,0.123726364051834,0.517182443066565,0.421572372979435,0,0,0.220108484100393,0.0934504648286556,0.407273145141005,0.722778619955329,0.0890963358709809,0.176682555445222,0.299130526119750,0,0,0.0984202185035984,0.201222042529487,0.673609982182292,0.380425760435909,0.146967445064113,0.292218766223402,0.305259787891604,0,0,0.245735887336606,0.365344568452201,0.436771637818230,0.295076437799749,0.144856577184978,0.270356365575536,0.541531878005713,0,0,0.505150729342836,0.125307576277766,0.242827921054228,0.651609016655234,0.101315253761357,0.364468686391628,0.0335069160047839,0,0,0.481980079334707,0.297436471626755,0.392355518307668,0.221588016544370,0.0919124735869692,0.588713563720344,0.175031294605209,0,0,0.331284341572268,0.490891217703830,0.214193896345819,0.116052485928374,0.0581721796785302,0.592507384757162,0.371596962239038,0];
x{2} = [0,0.0214888865015625,0.327189647364463,0.203904268387166,0.613046997612779,0.0547509755857816,0.178112441958145,0.609109444009954,0,0,0.0685540597017751,0.511875467583866,0.293702280915149,0.186640179720033,0.109548919385609,0.395732574225778,0.594377782953173,0,0,0.259115774225366,0.525578253697792,0.233584707625602,0.198955308511149,0.121185389971296,0.509304848463727,0.429754929325028,0,0,0.285509115871384,0.176584155451684,0.226848965312167,0.726693084183203,0.146085953000382,0.172963781982451,0.357320894920174,0,0,0.241958273585207,0.359762014672484,0.381077776006337,0.337013964400631,0.291250280351433,0.353452055009581,0.221027574101947,0,0,0.375328367065462,0.294969379504912,0.318942397280141,0.371689928537498,0.274897883264464,0.290751883896775,0.303291145792960,0,0,0.589265563344597,0.175375162357956,0.0667491167831127,0.615543441335458,0.179312595789679,0.330301741841141,0.0304722759184263,0,0,0.598205776349224,0.366929556178447,0.117134101695348,0.156166875754405,0.258438502829900,0.528929079699469,0.0611636979876581,0,0,0.444415528101808,0.440370818653580,0.148936989166881,0.223204255434406,0.226938250307594,0.519812263030372,0.243460889566548,0];
x{3} = [0,0,0.375878029964870,0.254229995022539,0.490275691258525,0.115559088646609,0.163425230419509,0.716466349608973,0,0,0.0459240634962492,0.665700942715707,0.252582349229370,0.113656258512438,0.206658285733121,0.324732167178394,0.574050793703115,0,0,0.359751478148277,0.624393923908575,0.196228078989969,0.0925028035420110,0.196228078989969,0.462514017710055,0.425160837811600,0,0,0.279711192606115,0.359610329210110,0.355996063316874,0.539415493815165,0.406852643790713,0.251727230447077,0.381424353553794,0,0,0.137015613871562,0.565160573232613,0.388210905969427,0.226064229293045,0.433882777259948,0.468275903535594,0.228359356452604,0,0,0.470592147695102,0.395151192334913,0.323532101540383,0.187176880579696,0.323532101540383,0.395151192334913,0.470592147695102,0,0,0.659896529007694,0.276513961803144,0.0977624487418805,0.483899433155502,0.293287346225642,0.397488820092019,0,0,0,0.478930572646912,0.451540074181584,0.136837306470546,0.112885018545396,0.228062177450910,0.693436542493146,0.0456124354901821,0,0,0.418389570756643,0.477905242758059,0.160919065675632,0.113786962561443,0.160919065675632,0.637206990344079,0.354021944486390,0];
Array = zeros(1,length(file));
for i = 1:length(file)
    fileName = file(i).name;
    filePath = fullfile(GoldStandardImages, fileName);
    tempimage = imread(filePath);
    %Gold = ExtractSign(tempimage);
    %GoldHog = HOG(Gold)
    GoldHog = x{i};
    Array(i) = sqrt(sum((GoldHog(:)-reference(:)).^ 2)); 
end

[value4, idx] = min(Array);

if idx == 1
   ans = 20;
elseif idx == 2
   ans = 40;
elseif idx == 3
    ans = 30;
end

end

% This function extracts a semi circle from the sign containing the
% first Digit at it's centre.
function croppedImageBW = ExtractSign(image)
    % enhance image contrast with histogram equalisation
    image = histeq(image);
    % resize all the images to a standardized size
    image = imresize(image, [100 100]);
    % convert to greyscale
    greyimage = rgb2gray(image);
    % if the average image pixel or brightness is below a 
    % certain threhold than adjust the image to increase brightness.
    if mean(greyimage(:)) < 120
     greyimage = imadjust(greyimage);
    end
    
    % If the motion blur of the image is above a certain threhold than use
    % a wiener filter to reduce the blur.
    if blurMetric(greyimage) > 0.55
        greyimage = wiener2(greyimage,[3 3]);
    end
    
    % find the centre and radius of every circle in the image.
    [centersBright, radiiBright] = imfindcircles(greyimage,[15 45],'ObjectPolarity','bright','Sensitivity',0.99);
    
    % get the image dimensions.
    [height,width,numberOfColorChannels] = size(image);
    % find the image centre.
    p2 = [height/2,width/2];
    X = zeros(1,length(centersBright));
    % extract the circle who is closet to the centre of the image as this
    % will refer to the signs boundary.
    for i = 1:length(centersBright)
        p1 = centersBright(i,:);
        X(i) =  norm(p1 - p2);
    end
    
    [value, idx] = min(X);
        
    %Circurly crop the image.
    [xx,yy] = ndgrid((1:height)-centersBright(idx,1),(1:width)-centersBright(idx,2));
    mask = uint8((xx.^2 + yy.^2)< radiiBright(1)^2);
    croppedImage = uint8(zeros(size(greyimage)));
    croppedImage = greyimage.*mask;
    croppedImage = imsharpen(croppedImage);

    % convert the cropped image to binary.
    croppedImageBW = 1-im2bw(croppedImage);
    % remove any imperfections in the image.
    croppedImageBW = bwareaopen(croppedImageBW, 10);
    
    % crop this new image using the geometry of the circle found to extract
    % the semi circle of the image.
    croppedImageBW = imcrop(croppedImageBW,[centersBright(idx,1)-radiiBright(1) centersBright(idx,2)-radiiBright(1) radiiBright(1) 2*radiiBright(1)]);
    % resize the image to a standard for comparison
    croppedImageBW = imresize(croppedImageBW, [60 30]);
end

% gets the files from a specific directory. 
function file_ = GetFileFromDirectory(Path,extension)
if ~isdir(Path)
    disp('no such directory exists');
    return;
end
filePath = fullfile(Path, extension);
file_ = dir(filePath);
end

% this function detects signs in the stress dataset by converting the image
% to HSV format an thresholding the image to extract the red features.
% These 
function cropped_image = redcars(image)
    %image1 = decorrstretch(image);
    image1 = histeq(image);
    imageHSV = rgb2hsv(image1);
    % Define thresholds for 'H' channel
    hMin = 0.9;
    hMax = 0.15;
    % Define thresholds for 'S' channel
    sMin = 0.4;
    sMax = 1.0;
    % Threshold image
    red = (imageHSV(:, :, 1) >= hMin | imageHSV(:, :, 1) <= hMax) & (imageHSV(:, :, 2) >= sMin & imageHSV(:, :, 2) <= sMax);
	% group all red object within 5 pixels of each other as one object
	% delete all objects smaller than 35 pixels in area
    red = bwareaopen(red, 40);
    % Remove connected components touching the image border
    red = imclearborder(red);
    % Construct a disk-shaped structuring element
    se = strel('disk', 10);
    % Dilation followed by erosion (worked better than 'imclose')
    red = imdilate(red, se);
    red = imerode(red, se);
    
    % try and catch used when the program fails to find any region of
    % interest. 
    try
        % find all the circles in the threholded image.
    [centers, radii] = imfindcircles(red,[15 45],'ObjectPolarity','bright','Sensitivity',0.8);
    
    % filter the cirlces found in the image by extracting the region which
    % contains the circle with the largest radius.
    if length(centers) > 0
        %figure();
        %imshow(image);
        [max_value, max_idx] = max(radii);
        %viscircles(centers(max_idx,:),radii(max_idx));
        new_radius = radii(1)+25;
        cropped_image = imcrop(image, [centers(max_idx,1)-new_radius centers(max_idx,2)-new_radius 2*new_radius 2*new_radius]);
    else
        % if the function fails to find any circles move to extracting
        % sqaures. from the image.
        red = imfill(red,'holes');
        stats = regionprops('table',red,'Centroid','MajorAxisLength','MinorAxisLength','BoundingBox','Area');
        % to classify a sqaure get the ratio of the minor and major axises
        % of a blob the ratio of the area and predicted area using the
        % minor and major axis. If both these values are below a certain
        % threshold then the red blob is a sqaure.
        square = {};
        for i = 1:length(stats.Centroid)
            temp = stats.MajorAxisLength(i)/stats.MinorAxisLength(i);
            tempArea = stats.MajorAxisLength(i)*stats.MinorAxisLength(i);
            areaRatio = tempArea/stats.Area(i);
            if temp <= 1.15 && temp >= 0.8 && areaRatio <= 1.4 && areaRatio >=0.7
                square = [square, i];
           end
        end
         % as there may be multiple squares in the image, excratct the
         % image closest to the centre of the image.
        [height,width,numberOfColorChannels] = size(red);
        p2 = [width/2,height/2];
        X = zeros(1,length(square));
        for i = 1:length(square)
            p1 = stats.Centroid(square{i},:);
            X(i) =  norm(p1 - p2);
        end
        
        [min_V, min_id] = min(X);
        theSquare = stats.BoundingBox(square{min_id},:);
        cropped_image =  imcrop(image,theSquare);      
    end
    catch
        % if no regions of interest are found just return the uncropped
        % image.
        warning('No Circles Found');
        [height,width,numberOfColorChannels] = size(red);
        p2 = [width/2,height/2];
        X = zeros(1,length(stats.Centroid));
        for i = 1:length(stats.Centroid)
            p1 = stats.Centroid(i,:);
            X(i) =  norm(p1 - p2);
        end
        
        [min_V, min_id] = min(X);
        theSquare = stats.BoundingBox(min_id,:);
        cropped_image =  imcrop(image,theSquare);
    end
end


function blur = blurMetric(original)
% Written by DO Quoc Bao, PhD Student in L2TI Laboratory, Paris 13 University, France
I = double(original);
[y x] = size(I);
Hv = [1 1 1 1 1 1 1 1 1]/9;
Hh = Hv';
B_Ver = imfilter(I,Hv);%blur the input image in vertical direction
B_Hor = imfilter(I,Hh);%blur the input image in horizontal direction
D_F_Ver = abs(I(:,1:x-1) - I(:,2:x));%variation of the input image (vertical direction)
D_F_Hor = abs(I(1:y-1,:) - I(2:y,:));%variation of the input image (horizontal direction)
D_B_Ver = abs(B_Ver(:,1:x-1)-B_Ver(:,2:x));%variation of the blured image (vertical direction)
D_B_Hor = abs(B_Hor(1:y-1,:)-B_Hor(2:y,:));%variation of the blured image (horizontal direction)
T_Ver = D_F_Ver - D_B_Ver;%difference between two vertical variations of 2 image (input and blured)
T_Hor = D_F_Hor - D_B_Hor;%difference between two horizontal variations of 2 image (input and blured)
V_Ver = max(0,T_Ver);
V_Hor = max(0,T_Hor);
S_D_Ver = sum(sum(D_F_Ver(2:y-1,2:x-1)));
S_D_Hor = sum(sum(D_F_Hor(2:y-1,2:x-1)));
S_V_Ver = sum(sum(V_Ver(2:y-1,2:x-1)));
S_V_Hor = sum(sum(V_Hor(2:y-1,2:x-1)));
blur_F_Ver = (S_D_Ver-S_V_Ver)/S_D_Ver;
blur_F_Hor = (S_D_Hor-S_V_Hor)/S_D_Hor;
blur = max(blur_F_Ver,blur_F_Hor);
end


function H=HOG(Im)
% created by Oswaldo Ludwig.
%Image descriptor based on Histogram of Orientated Gradients for gray-level images. This code 
%was developed for the work: O. Ludwig, D. Delgado, V. Goncalves, and U. Nunes, 'Trainable 
%Classifier-Fusion Schemes: An Application To Pedestrian Detection,' In: 12th International IEEE 
%Conference On Intelligent Transportation Systems, 2009, St. Louis, 2009. V. 1. P. 432-437. In 
%case of publication with this code, please cite the paper above.
nwin_x=3;%set here the number of HOG windows per bound box
nwin_y=3;
B=9;%set here the number of histogram bins
[L,C]=size(Im); % L num of lines ; C num of columns
H=zeros(nwin_x*nwin_y*B,1); % column vector with zeros
m=sqrt(L/2);
if C==1 % if num of columns==1
    Im=im_recover(Im,m,2*m);%verify the size of image, e.g. 25x50
    L=2*m;
    C=m;
end
Im=double(Im);
step_x=floor(C/(nwin_x+1));
step_y=floor(L/(nwin_y+1));
cont=0;
hx = [-1,0,1];
hy = -hx';
grad_xr = imfilter(double(Im),hx);
grad_yu = imfilter(double(Im),hy);
angles=atan2(grad_yu,grad_xr);
magnit=((grad_yu.^2)+(grad_xr.^2)).^.5;
for n=0:nwin_y-1
    for m=0:nwin_x-1
        cont=cont+1;
        angles2=angles(n*step_y+1:(n+2)*step_y,m*step_x+1:(m+2)*step_x); 
        magnit2=magnit(n*step_y+1:(n+2)*step_y,m*step_x+1:(m+2)*step_x);
        v_angles=angles2(:);    
        v_magnit=magnit2(:);
        K=max(size(v_angles));
        %assembling the histogram with 9 bins (range of 20 degrees per bin)
        bin=0;
        H2=zeros(B,1);
        for ang_lim=-pi+2*pi/B:2*pi/B:pi
            bin=bin+1;
            for k=1:K
                if v_angles(k)<ang_lim
                    v_angles(k)=100;
                    H2(bin)=H2(bin)+v_magnit(k);
                end
            end
        end
                
        H2=H2/(norm(H2)+0.01);        
        H((cont-1)*B+1:cont*B,1)=H2;
    end
end

end
