
for i=1:1:10
str = int2str(i);
s=strcat(str,'.tif')
a=10;
fileID = fopen('airfoil.txt','a+');
fprintf(fileID,'%6s %12s\r\n','x','exp(x)');
fprintf(fileID,'%s ',s);
fprintf(fileID,'%d\r\n',a);
fclose(fileID);
end