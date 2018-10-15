# DEBLUR
%Image processing assginment 2
%GLADIS ZACHARIA
%17307R004
%Codes except the inbuilt functions for FFT , IFFT, SSIM has been written by myself.
function varargout = DEBLUR(varargin)
% DEBLUR MATLAB code for DEBLUR.fig
%      DEBLUR, by itself, creates a new DEBLUR or raises the existing
%      singleton*.
%
%      H = DEBLUR returns the handle to a new DEBLUR or the handle to
%      the existing singleton*.
%
%      DEBLUR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DEBLUR.M with the given input arguments.
%
%      DEBLUR('Property','Value',...) creates a new DEBLUR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DEBLUR_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DEBLUR_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DEBLUR

% Last Modified by GUIDE v2.5 14-Oct-2018 19:56:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DEBLUR_OpeningFcn, ...
                   'gui_OutputFcn',  @DEBLUR_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before DEBLUR is made visible.
function DEBLUR_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DEBLUR (see VARARGIN)

% Choose default command line output for DEBLUR
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DEBLUR wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DEBLUR_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in loadimage.
function loadimage_Callback(hObject, eventdata, handles)
% hObject    handle to loadimage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global initial_image  Finalimage K path
[file , path]=uigetfile('*.jpg','All files');%opening the directory
initial_image=imread([path,file]);%reading the image
imag=initial_image;
%imag = imread('Gladis .jpg');
ker=imread('kernel1.png');% reading the kernel
imag=rgb2hsv(imag);%converting image to hsv
img=imag(:,:,3);%taking the v value of the loaded image
ker=rgb2hsv(ker);%converting kernel to hsv
ker=ker(:,:,3);%taking the v value of the kernel
[m,n]=size(img);% finding size of image
[r,c]=size(ker);%finding size of kernel
ker=padarray(ker,[(m-r)/2,(n-c)/2],0,'both');%padding kernel to the size of image
I=fft2((img));%FFT of image
K=fft2(ker);%FFT of kernel
J=I.*K; % blurring the image
D=ifft2(J);% inverse FFT to obtain blurred image
D=ifftshift(D);
%%now we have to normalize the output on an intensity scale of [0,1]
maxm=max(D(:));
minm=min(D(:));
[m,n]=size(J);
for i= 1:m
    for j=1:n
         D(i,j)= (D(i,j))/(maxm-minm);
    end
end
imag(:,:,3)=D;% assigning the obtained value to v value of the loaded image
Finalimage=hsv2rgb(imag);% converting back to rgb
axes(handles.axes1);
imshow(Finalimage);% displaying the blurred image

% --- Executes on button press in fullinverse.
function fullinverse_Callback(~, eventdata, handles)
% hObject    handle to fullinverse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Finalimage K
Blurred=Finalimage;
Blurred=rgb2hsv(Blurred);%converting to hsv
Blur=Blurred(:,:,3);% obtaining the v value of the blurred image
B=fft2(Blur);% FFT of the blurred image
S=B./K; % inverse filter
L=ifft2(S);% taking inverse FFT
L=ifftshift(L);
%now we have to normalize the output on an intensity scale of [0,1]
maxm=max(L(:));
minm=min(L(:));
[m,n]=size(K);
MSE=0;
for i= 1:m
    for j=1:n
         L(i,j)= (L(i,j))/(maxm-minm);
         MSE = MSE+((Blur(i,j)-L(i,j)).^2)/(m*n);
    end
end
Blurred(:,:,3)=L;% assigning the obtained value to v value of the blurred image
PSNR = 10*log10(256*256/MSE);
ssimval = ssim(L,Blur);
disp(ssimval);
disp(MSE);
disp(PSNR);
Final=hsv2rgb(Blurred);% converting the deblurred image to rgb
axes(handles.axes2);
imshow(Final);% displaying the deblurred image


% --- Executes on slider movement.
function truncated_Callback(hObject, eventdata, handles)
% hObject    handle to truncated (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global Finalimage K
Blurred=Finalimage;
n=get(hObject,'value')%Order of butterworth filter

Blurred=rgb2hsv(Blurred);%converting to hsv
Blur=Blurred(:,:,3);% obtaining the v value of the blurred image
B=fft2(Blur);% FFT of the blurred image
% creating a butterworth lowpass filter
[M N] = size(K);
filter=zeros(M,N);
d0=10;				%Cutoff Frequency		
	for i=1:M
    	for j=1:N
        	dist=(i-M/2)^2 + (j-N/2)^2;
        	filter(i,j)= ( 1 + (dist/d0)^(2*n))^(-1);
        end
    end
z=K.*filter;%truncating the kernel
S=B./z;%truncated inverse filter
L=ifft2(S);% IFFT
L=ifftshift(L);
%now we have to normalize the output on an intensity scale of [0,1]
maxm=max(L(:));
minm=min(L(:));
[m,n]=size(K);
MSE=0;
for i= 1:m
    for j=1:n
         L(i,j)= (L(i,j))/(maxm-minm);
         MSE = MSE+((Blur(i,j)-abs(L(i,j))).^2)/(m*n);
    end
end
Blurred(:,:,3)=abs(L);% assigning the obtained value to v value of the blurred image
PSNR = 10*log10(256*256/MSE);
ssimval = ssim(abs(L),Blur);
disp(ssimval);
disp(MSE);
disp(PSNR);
Final=hsv2rgb(Blurred);%converting the deblurred image to rgb
axes(handles.axes2)
imshow(Final);% displaying the deblurred image



% --- Executes during object creation, after setting all properties.
function truncated_CreateFcn(hObject, eventdata, handles)
% hObject    handle to truncated (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function weiner_Callback(hObject, eventdata, handles)
% hObject    handle to weiner (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global Finalimage K
a=get(hObject,'value');
a=a*50;
disp(a);
Blurred=Finalimage;
Blurred=rgb2hsv(Blurred);%converting blurred image to hsv
Blur=Blurred(:,:,3);%taking v value of blurred image
B=fft2(Blur);% FFT of blurred image
%implementing weiner filter
[m,n]=size(K);
for i= 1:m
    for j=1:n
     S(i,j)=((conj(K(i,j))/((abs(K(i,j)).^2)+a))).*B(i,j);
    end
end
L=ifft2(S); % IFFT
L=ifftshift(L);
%now we have to normalize the output on an intensity scale of [0,1]
maxm=max(L(:));
minm=min(L(:));
MSE=0;
for i= 1:m
    for j=1:n
         L(i,j)= (L(i,j))/(maxm-minm);
         MSE = MSE+(((Blur(i,j)-L(i,j)).^2)/(m*n));
    end
end
Blurred(:,:,3)=L;% assigning the obtained value to v value of the blurred image
PSNR = 10*log10(256*256/MSE);
ssimval = ssim(L,Blur);
disp(ssimval);
disp(MSE);
disp(PSNR);
Final=hsv2rgb(Blurred);% converting deblurred image to rgb
axes(handles.axes2);
imshow(Final);%displaying the deblurred image

% --- Executes during object creation, after setting all properties.
function weiner_CreateFcn(hObject, eventdata, handles)
% hObject    handle to weiner (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function constrained_Callback(hObject, eventdata, handles)
% hObject    handle to constrained (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global Finalimage K
a=get(hObject,'value');
a=50*a;
disp(a);
Blurred=Finalimage;
Blurred=rgb2hsv(Blurred);%converting blurred image to hsv
Blur=Blurred(:,:,3);%taking v value of blurred image
B=fft2(Blur);%FFT 
[m,n]=size(K);
p=[0 -1 0 , -1 4 -1 , 0 -1 0];% p(x,y)
% padding p to size of image
p=padarray(p,[(m)/2,(n-16)/2],0,'pre');
p=padarray(p,[(m-2)/2,(n-2)/2],0,'post');
P=fft2(p);% FFT 
%constrained least square filter
[m,n]=size(K);
for i= 1:m
    for j=1:n
        S(i,j)=((conj(K(i,j))/((abs(K(i,j))^2)+(a*(abs(P(i,j))^2))))).*B(i,j);
    end
end
L=ifft2(S);%IFFT
L=ifftshift(L);
%now we have to normalize the output on an intensity scale of [0,1]
maxm=max(L(:));
minm=min(L(:));
[m,n]=size(K);
MSE=0;
for i= 1:m
    for j=1:n
         L(i,j)= (L(i,j))/(maxm-minm);
         MSE = MSE+((Blur(i,j)-L(i,j)).^2)/(m*n);
    end
end
Blurred(:,:,3)=L;% assigning the obtained value to v value of the blurred image
PSNR = 10*log10(256*256/MSE);
ssimval = ssim(L,Blur);
disp(ssimval);
disp(MSE);
disp(PSNR);
Final=hsv2rgb(Blurred);%converting deblurred image to rgb
axes(handles.axes2);
imshow(Final);%diplaying deblurred image



% --- Executes during object creation, after setting all properties.
function constrained_CreateFcn(hObject, eventdata, handles)
% hObject    handle to constrained (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
