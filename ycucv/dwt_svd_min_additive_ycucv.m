

clear all;
N = 64;
key = 100;
k1 = .16;
band = 'hl';


%loading the host image
rgb_image=imread('lena.jpg');

%colorTransform = makecform('srgb2lab');
%lab_image = applycform(rgb_image, colorTransform);

% RGB to Lab conversion

R = rgb_image(:,:,1);
G = rgb_image(:,:,2);
B = rgb_image(:,:,3);

[Y_img, Cu_img, Cv_img] = rgb2ycucv(R,G,B);

% L_img = lab_image(:,:,1);
% A_img = lab_image(:,:,2);
% B_img = lab_image(:,:,3);

[a,b]=size(Cv_img);
I=imresize(Cv_img,[1024,1024]);
figure(1),imshow(rgb_image);

%loading watermark
o=imread('isi.jpg');
o=rgb2gray(o);
o=imresize(o,[64,64]);
level=graythresh(o);
o=im2bw(o,level);
figure(2),imshow(o);


%perform arnold trnsform key times
in = o; 
out = zeros(N);

for s=1:key
    
    for y=0:N-1
        for x=0:N-1
            p = [1 1; 1 2] * [ x; y ];
            out( mod(p(2),N)+1 , mod(p(1),N)+1 ) = in(y+1, x+1);
        end
    end
    
    in = out;
    
end
    
w = out;


% w = o;

%perform 4 level dwt on I 
[ll1,lh1,hl1,hh1] = dwt2(I,'db1','d');
[ll2,lh2,hl2,hh2] = dwt2(ll1,'db1','d');
[ll3,lh3,hl3,hh3] = dwt2(ll2,'db1','d');
[ll4,lh4,hl4,hh4] = dwt2(ll3,'db1','d');

%image adaptive pixelwise masking parameters
l = 3;

% frequency sensetivity

    if (band == 'hh')
        F_band = 1.414;
    else
        F_band = 1;
    end
    
    switch l
        case ( 0 )
            F_level = 1;
        case ( 1 )
            F_level = 0.32;
        case ( 2 )
            F_level = 0.16;
        case ( 3 )
            F_level = 0.10;
    end
    
    F = F_band * F_level;
    
    
% luminance matrix

for i = 0:63
    for j = 0:63
        L(i+1,j+1) = (1/256) * ll4( 1 + floor( i / 2^(3-l) ), 1 + floor( j / 2^(3-l) ) );

        if ( L(i+1,j+1) < 0.5 )
            L_mod(i+1,j+1) = 1 - L(i+1,j+1);
        else
            L_mod(i+1,j+1) = L(i+1,j+1);
        end
  
        Lum(i+1,j+1) = 1 + L_mod(i+1,j+1);
    end
end



% contrast masking

C = hl4;        % last row and colom is added as hl4

for i = 1:63
    for j = 1:63
        
            for k = 0:(3-l)
                for x = 0:1
                    for y = 0:1
                       % for m = l:3
                        
                        p = hl4( y+ i/2^k ,x+ j/2^k );
                        q = lh4( y+ i/2^k ,x+ j/2^k );
                        r = hh4( y+ i/2^k ,x+ j/2^k );
                        
                        c_edge = (1/16)^k * [ p^2 + q^2 + r^2 ];
                        %end
                    end
                end
            end
            
            
                % generally var( ll4( x+i/2^(3-l) , y+i/2^(3-l) ))
                
                     p = ll4( i/2^(3-l) , j/2^(3-l) );       % x=0,y=0;
                     q = ll4( i/2^(3-l) , 1+j/2^(3-l) );     % x=0,y=1;
                     r = ll4( 1+i/2^(3-l) , j/2^(3-l) );     % x=1,y=0;
                     s = ll4( 1+i/2^(3-l) , 1+j/2^(3-l) );   % x=1,y=1;
                
                c_texture = var([p q r s]);
        
        C(i,j) = ( c_edge * c_texture )^0.2 ;
    end
end

% considering frequency ,luminance and contrast together

for i = 1:64
    for j = 1:64
        
        Q(i,j) = F * Lum(i,j) * C(i,j);
        
        %if( Q(i,j) <= 0)
        %    Q(i,j) = 0.0001;   % nullify the coefficients which are negative in value
        %end
    end
end



% apply svd to hl4
[u_s ,s_s ,v_s] = svd(hl4);

%apply svd to watermark

[ u_w, s_w, v_w] = svd(w);

%calculating watermarking strength

jnd = 0.5 * Q ;


alpha = abs(min(min( jnd /det(u_s * v_s') ))/s_w(1,1));

%watermark embedding
s = zeros(64,64);
for i = 1:64
    
    s(i,i) = s_s(i,i) + alpha * s_w(i,i) ;
    
end

%[uu,ss,vv] = svd(s);

new_hl4 = u_s * s * v_s' ;

ll3 = idwt2(ll4,lh4,new_hl4,hh4,'db1','d');
ll2 = idwt2(ll3,lh3,hl3,hh3,'db1','d');
ll1 = idwt2(ll2,lh2,hl2,hh2,'db1','d');
I   = idwt2(ll1,lh1,hl1,hh1,'db1','d');

%embedded watermark
Cv_img =imresize(I,[a,b]);

% [ R, G, B] = Lab2rgb(L_img, A_img, B_img);
% y = cat(3,R,G,B);
  
ycucv_image = cat(3,Y_img,Cu_img,Cv_img);
%figure(10), imshow(ycocg_image);
y = ycucv2rgb(ycucv_image);
figure(3), imshow(y/max(max(max(y))));

% PSNR measurment

[row column k] = size(rgb_image);

mse_R_image = (double(rgb_image(:,:,1)) - double(y(:,:,1))).^2;
mse_G_image = (double(rgb_image(:,:,2)) - double(y(:,:,2))).^2;
mse_B_image = (double(rgb_image(:,:,3)) - double(y(:,:,3))).^2;

mse_R = (sum(sum(mse_R_image)))/(row*column);
mse_G = (sum(sum(mse_G_image)))/(row*column);
mse_B = (sum(sum(mse_B_image)))/(row*column);

mse = (mse_R + mse_G + mse_B)/3 ;

PSNR = 10*log10( 255^2/mse )

% Storing the image to lossy file formats.jpeg
q=input('Quality Factor q = ');
imwrite(y,'svd_min_aditive.jpg','jpg','quality',q);
y=imread('svd_min_aditive.jpg');
figure(4),imshow((y));                  
title ('compressed image');



% watermark extraction

[p,q]=size(y);

% colorTransform = makecform('srgb2lab');
% lab_image = applycform(y, colorTransform);

R = y(:,:,1);
G = y(:,:,2);
B = y(:,:,3);

[Y_img, Cu_img, Cv_img] = rgb2ycucv( R, G, B);

% L_img = lab_image(:,:,1);
% A_img = lab_image(:,:,2);
% B_img = lab_image(:,:,3);

I=imresize(Cv_img,[1024,1024]);

[ll1,lh1,hl1,hh1]=dwt2(I,'db1','d');
[ll2,lh2,hl2,hh2]=dwt2(ll1,'db1','d');
[ll3,lh3,hl3,hh3]=dwt2(ll2,'db1','d');
[ll4,lh4,hl4,hh4]=dwt2(ll3,'db1','d');

[u1,s1,v1] = svd(hl4);

s_new = zeros(64,64);
for i = 1:64
        
        s_new(i,i) = abs( s1(i,i) - s_s(i,i))/alpha;
   
end

extracted_watermark = u_w * s_new * v_w' ;
w1 = extracted_watermark;


%perform inverse arnold tx key times
in = w1; 
out = zeros(N);

for k=1:key
    
    for y=0:N-1
        for x=0:N-1
            p = [2 -1 ; -1 1] * [ x; y ];
            out( mod(p(2),N)+1 , mod(p(1),N)+1 ) = in(y+1, x+1);
        end
    end
    
    in = out;
    
end

w1 = out;


%gray to binary convertion
level = graythresh(w1);
w1 = im2bw(w1,level);

figure(5),imshow(w1);
title('extracted image');
        
%Normalised Coefficient (NC) Measurement 
nc_num=0; nc_den=0;a=0; b=0;
for  i=1:64
    for j=1:64
        nc_num=nc_num+(w1(i,j)*o(i,j));
        a=a+(o(i,j))^2;
        b=b+(w1(i,j))^2;
           
    end
end
nc_den=sqrt(nc_den+a*b); 
 
nc=(nc_num/nc_den)

