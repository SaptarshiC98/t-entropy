row=256;
col=256;
img=rand(row,col);
img=round(img);
imshow(img);
img = imread('lena.png');
imshow(img)
size(img)
im=img(:,:,1);
[X,Y]=Enc(BW);
Z=Dec(X,Y);
BW=im2bw(img,0.5);
tic
[X,Y]=Enc1(im);
toc
tic
Z=Dec1(X,Y);
toc
n=100
m=100
Y=zeros(m,n);
for i=1:m
    for j=1:n
        r1=randsample(256,1)-1;
        Y(i,j)=120;
    end
end
Y=uint8(Y);

X1=img(:,:,1);
X2=img(:,:,2);
X3=img(:,:,3);
[Y(:,:,1),Z(:,:,1)]=Enc1(X1);
[Y(:,:,2),Z(:,:,2)]=Enc1(X2);
[Y(:,:,3),Z(:,:,3)]=Enc1(X3);
Y=uint8(Y);
Z=uint8(Z);
rec=zeros(r,c,3);
rec(:,:,1)=Dec1(Y(:,:,1),Z(:,:,1));
rec(:,:,2)=Dec1(Y(:,:,2),Z(:,:,2));
rec(:,:,3)=Dec1(Y(:,:,3),Z(:,:,3));
rec=uint8(rec);
[r,c]=size(X1);
Y=zeros(r,c,3);
Z=zeros(r,c,3);





Dom::GaloisField(q)


m = 4;
x = [3 2 9];
y = gf(x,m)
y=gf(1:256,8);
z=int(y(1))
Z=int16(im);
y=gf(im,8);

%%%%%%%%%%%%%________DEFINE MULTIPLICATION_________%%%%%%%%%%

M=zeros(256,256);
y=gf(0:255,8);
dummy=gf(M,8);
for i= 1:256
    for j=i:256
        dummy(i,j)=y(i).*y(j);
        disp([i,j])
    end
end

for i= 1:256
    for j=i:256
    M(i,j)=gf2dec(dummy(i,j));
    disp([i,j])
    end
end   

for i=2:256
    for j=1:(i-1)
        M(i,j)=M(j,i);
        disp([i,j])
    end
end

%%%%%%%%%%%%%%%%%%%__________DEFINE ADDITION_________%%%%%%%%%%%%%%%

A=zeros(256,256);
y=gf(0:255,8);
dummy=gf(A,8);
for i= 1:256
    for j=i:256
        dummy(i,j)=y(i)+y(j);
        disp([i,j])
    end
end

for i= 1:256
    for j=i:256
    A(i,j)=gf2dec(dummy(i,j));
    disp([i,j])
    end
end   

for i=2:256
    for j=1:(i-1)
        A(i,j)=A(j,i);
        disp([i,j])
    end
end

%%%%%%%%%%%%%%%%%%_______DIVISION TABLE______%%%%%%%%%%%%%%

D=zeros(256,255);
y=gf(0:255,8);
dummy=gf(D,8);
for i= 1:256
    for j=1:255
        dummy(i,j)=y(i)/y(1+j);
        disp([i,j])
    end
end

for i= 1:256
    for j=i:256
    A(i,j)=gf2dec(dummy(i,j));
    disp([i,j])
    end
end   

for i=2:256
    for j=1:(i-1)
        A(i,j)=A(j,i);
        disp([i,j])
    end
end

for i= 1:256
    for j=1:255
    D(i,j)=gf2dec(dummy(i,j));
    disp([i,j])
    end
end   
%%%%%%%%%%__________SUBSTRACTION TABLE_______%%%%%%%%%%%%%%
S=zeros(256,256);
y=gf(0:255,8);
dummy=gf(S,8);
for i= 1:256
    for j=1:256
        dummy(i,j)=y(i)-y(j);
        disp([i,j])
    end
end
tic
for i= 1:256
    for j=1:256
    S(i,j)=gf2dec(dummy(i,j));
    disp([i,j])
    end
end   
toc

im=[0,50,255
    30,255,0
    150,6,230]


im=uint8(im);
imshow(im)



A=load('Add.csv');
D=load('Div.csv');
M=load('Mult.csv');
S=load('Sub.csv');
B=Shamir_n_k_grey(6,4,im,A,M);
Y=B(:,:,1);
Y=uint8(Y);
imshow(Y);
R=Combine_S_N_k(B,2,[1,2],A,M,D,S);
imshow(R)

B=Shamir_n_k_color(4,2,img,A,M);

im=100;
Y=B(:,:,:,1);
Y=uint8(Y);
imshow(Y);
R=Combine_S_N_k(B,2,[1,2],A,M,D,S);

R=Combine_S_n_k_color(B,2,[1,2],A,M,D,S);
imshow(R);

A=load('Add.csv');
D=load('Div.csv');
M=load('Mult.csv');
S=load('Sub.csv');
img=rgb2gray(img);
img = imread('peppers.jpg');
img = imread('kodim23.png');
img=imread('dollar.tif');
imshow(img);
tic
B=Shamir_n_k_color(5,4,img,A,M);
toc
imshow(uint8(B(:,:,:,5)));
tic
R=Combine_S_n_k_color(B,2,[4,5],A,M,D,S);
toc
imshow(R);
subplot(3,3,1)       % add first plot in 2 x 1 grid
imshow(uint8(B(:,:,:,1)));
title('Share 1')
subplot(3,3,2)       % add first plot in 2 x 1 grid
imshow(uint8(B(:,:,:,2)));
title('Share 2')
subplot(3,3,3)       % add first plot in 2 x 1 grid
imshow(uint8(B(:,:,:,3)));
title('Share 3')
subplot(3,3,4)       % add first plot in 2 x 1 grid
imshow(uint8(B(:,:,:,4)));
title('Share 4')
subplot(3,3,5)       % add first plot in 2 x 1 grid
imshow(uint8(B(:,:,:,5)));
title('Share 5')
subplot(3,3,6)       % add first plot in 2 x 1 grid
imshow(R);
title('Reconstruction from shares 1,2,3')
subplot(3,3,7)       % add first plot in 2 x 1 grid
imshow(R);
title('Reconstruction from shares 3,4,5')
subplot(3,3,8)       % add first plot in 2 x 1 grid
imshow(R);
title('Reconstruction from shares 4,1,5')
subplot(3,3,9)       % add first plot in 2 x 1 grid
imshow(R);
title('Reconstruction from shares 2,3,4')


tic
B=Shamir_n_k_grey(6,4,img,A,M);
toc
imshow(uint8(B(:,:,6)));
tic
R=Combine_S_N_k(B,3,[1,3,5],A,M,D,S);
toc
imshow(R);
