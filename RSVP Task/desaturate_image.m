
function image = desaturate_image(img)

mask=img;
mask=rgb2hsv(mask);
[h,s,v]=imsplit(mask);
s=s*0.4; %scale saturation
v=v*0.8; %scale brightness
h=h*0.3; %rotate hue
mask=cat(3,h,s,v);
mask=hsv2rgb(mask);
mask=mask * 255;
image = mask;

end 