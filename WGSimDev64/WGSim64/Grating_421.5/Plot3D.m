InputFile="WGArgsSpectrum.m";
OutputFile1="TestWL4.59e-07RayPoint.3D"
OutputFile2="TestWL4.59e-07RayPoints.3D"


run(InputFile);

figure(10);
clf;
hold on
x=[-BoundaryWidth, BoundaryWidth, BoundaryWidth,-BoundaryWidth,-BoundaryWidth];
y=[BoundaryHeight,BoundaryHeight,-BoundaryHeight,-BoundaryHeight,BoundaryHeight];
z=[0,0,0,0,0];
plot3(x,y,z,'b');
x=[-InputRadius, InputRadius, InputRadius,-InputRadius,-InputRadius];
y=[InputRadius,InputRadius,-InputRadius,-InputRadius,InputRadius];
x=x+InputGratingOffsetX;
plot3(x,y,z,'b');

fid=fopen(OutputFile1,"r","l");
x=fread(fid,Inf,"single",8);
fseek(fid,4);
y=fread(fid,Inf,"single",8);
fseek(fid,8);
z=fread(fid,Inf,"single",8);
fclose(fid);

fid=fopen(OutputFile2,"r","l");
n=fread(fid,Inf,"int32");
fclose(fid);

i = 1;
o=1;
l=length(n);
while i < l
  t=abs(n(i))+o-1;
  a=x(o:t);
  b=y(o:t);
  c=z(o:t);
  plot3(a,b,c,"r");
  o=t+1;
  i++;
endwhile
