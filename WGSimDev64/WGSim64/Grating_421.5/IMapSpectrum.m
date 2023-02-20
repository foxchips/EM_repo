close all;
clear all;
PhotonPlotUpperLimit=50;
InputFile1="WGArgsSpectrumBase.m";
InputFile2="WGArgsSpectrum.m";
run(InputFile1);
run(InputFile2);

a = InputAngleLimitX(2); % half angle
b = InputAngleLimitY(2);
SolidAngle = 4*asin( sind(a) * sind(b)); % solid angle
DetectorRes = FOVXSamples * FOVYSamples; % detector resolution
if DetectorCircular ==1
  DetectorSize = pi*DetectorXYZW'(4)^2;
else
  DetectorSize = (2*DetectorXYZW'(4))^2;
endif


if(isunix())
  graphics_toolkit fltk;
endif

if 1
  for WL=1:Wavelengths
    InputFile=strcat("WGArgsSpectrumWork",num2str(WL),".m");
    if(isunix())
    system(["cp " InputFile1 " " InputFile]);
    else
    system(["copy " InputFile1 " " InputFile]);
    endif
    fid=fopen(InputFile,"a");
    fdisp (fid, "SpectrumWavelengths=1;");
    fdisp (fid, "SpectrumWavelengthAndReleativeIntensity=[");
    fdisp (fid, strcat(num2str(WavelengthAndReleativeIntensity(WL,1)), ",",num2str(WavelengthAndReleativeIntensity(WL,2)),";"));
    fdisp (fid,"];");
    fdisp (fid, strcat("WavelengthFolder=",num2str(WavelengthFolder(WL)),";"));
    fdisp (fid, strcat("WavelengthRefractiveIndex=",num2str(WavelengthRefractiveIndex(WL)),";"));
    fclose(fid);
    OutputFile=strcat("TestWL",num2str(WavelengthAndReleativeIntensity(WL,1)));
    if(isunix())
      [st, outp]=system(["../WGSimDev64 " InputFile " " OutputFile]);
      #[st, outp]=system(["wine ../WGSim64.exe " InputFile " " OutputFile]);
    else
      [st, outp]=system(["..\\WGSim64.exe " InputFile " " OutputFile]);
    endif
    Status(WL)=st;
    printf("Wavelength Index:%d Status:%d\n",WL,st);
    fflush(stdout);
    Output{WL}=outp;
 endfor
endif

if (AllRayDetectors == 1)
   for d=0:(AllRayDetectors-1) 
   for WL=1:Wavelengths
     OutputFile=strcat("Results/TestWL",num2str(WavelengthAndReleativeIntensity(WL,1)),"AllRay");
     fid=fopen(strcat(OutputFile,".",num2str(d)),"r","l");
     WM=fread(fid,[2*AllRayDetectorWHPixels(1),2*AllRayDetectorWHPixels(2)],"float");
     WM=WM';
     fclose(fid);
     %{
     figure('name',strcat("All Ray Detector:",num2str(d)));
     clf;
     tx = linspace (-90,90,2*AllRayDetectorWHPixels(1));
     ty = linspace (-90,90,2*AllRayDetectorWHPixels(2));
     mesh(tx,ty,WM);
     shading("faceted");
     %}
     fname=strcat("Results/AllRayDetectorMatrix",num2str(WavelengthAndReleativeIntensity(WL,1)),"-",num2str(d),".csv");
     csvwrite (fname, WM);
    if(WL==1)
      IntensityMap=WM .* WavelengthPhotonToLm(WL);
    else
      IntensityMap+=WM .* WavelengthPhotonToLm(WL);
    endif
   endfor
     figure('name',strcat("All Ray Angular Detector, All Wavelength"));
     clf;
     tx = linspace (-90,90,2*AllRayDetectorWHPixels(1));
     ty = linspace (-90,90,2*AllRayDetectorWHPixels(2));
     mesh(tx,ty,IntensityMap);
     shading("faceted");
     fname=strcat("Results/AllRayAngularDetectorMatrixAllWavelength.csv");
     csvwrite (fname, IntensityMap);
   endfor

   for WL=1:Wavelengths
     OutputFile=strcat("Results/TestWL",num2str(WavelengthAndReleativeIntensity(WL,1)),"AllRaySpatial");
     fid=fopen(strcat(OutputFile,".All"),"r","l");
     WM=fread(fid,[2*AllRayDetectorSpatialWHPixels(1),2*AllRayDetectorSpatialWHPixels(2)],"float");
     WM=WM';
     fclose(fid);
     %{
     figure('name',strcat("All Ray Spatial Detector:",num2str(WavelengthAndReleativeIntensity(WL,1))));
     clf;
     tx = linspace (-AllRayDetectorSpatialWH(1),AllRayDetectorSpatialWH(1),2*AllRayDetectorSpatialWHPixels(1));
     ty = linspace (-AllRayDetectorSpatialWH(2),AllRayDetectorSpatialWH(2),2*AllRayDetectorSpatialWHPixels(2));
     mesh(tx,ty,WM);
     shading("faceted");
     %}
     fname=strcat("Results/AllRaySpatialDetectorMatrix",num2str(WavelengthAndReleativeIntensity(WL,1)),"-All",".csv");
     csvwrite (fname, WM);
     if(WL==1)
      IntensityMap=WM .* WavelengthPhotonToLm(WL);
    else
      IntensityMap+=WM .* WavelengthPhotonToLm(WL);
    endif    
   endfor
     figure('name',strcat("All Ray Spatial Detector, All Wavelength"));
     clf;
     tx = linspace (-AllRayDetectorSpatialWH(1),AllRayDetectorSpatialWH(1),2*AllRayDetectorSpatialWHPixels(1));
     ty = linspace (-AllRayDetectorSpatialWH(2),AllRayDetectorSpatialWH(2),2*AllRayDetectorSpatialWHPixels(2));
     mesh(tx,ty,IntensityMap);
     shading("faceted");
     fname=strcat("Results/AllRaySpatialDetectorMatrixAllWavelength.csv");
     csvwrite (fname, IntensityMap);
   
endif


   for sp=0:1
   if (sp==0) sps="S"; else sps="P";endif
   for d=0:(Detectors-1)
       C=rand(FOVYSamples,FOVXSamples,3);
       PhotonsIn=0;
       PhotonsOut=0;
   for WL=1:Wavelengths
     OutputFile=strcat("Results/TestWL",num2str(WavelengthAndReleativeIntensity(WL,1)),"Counts");
     fid=fopen(strcat(OutputFile,".",num2str(d)),"r","l");
     photons=fread(fid,1,"uint64");
     PhotonsIn+=photons/2*WavelengthPhotonToLm(WL);
     fclose(fid);
     OutputFile=strcat("Results/TestWL",num2str(WavelengthAndReleativeIntensity(WL,1)),sps,".",num2str(d));
     fid=fopen(OutputFile,"r","l");
     WM=fread(fid,[FOVXSamples,FOVYSamples],"float");
     WM=WM';
    if(WL==1)
      IntensityMap=WM .* WavelengthPhotonToLm(WL);
      C(:,:,1)=IntensityMap.*(WavelengthRGB(WL,1)/255.0);
      C(:,:,2)=IntensityMap.*(WavelengthRGB(WL,2)/255.0);
      C(:,:,3)=IntensityMap.*(WavelengthRGB(WL,3)/255.0);
    else
      IntensityMap+=WM .* WavelengthPhotonToLm(WL);
      C(:,:,1)+=IntensityMap.*(WavelengthRGB(WL,1)/255.0);
      C(:,:,2)+=IntensityMap.*(WavelengthRGB(WL,2)/255.0);
      C(:,:,3)+=IntensityMap.*(WavelengthRGB(WL,3)/255.0);
     endif
     fclose(fid);
    endfor
    PhotonsOut = sum(sum(IntensityMap));
    %figure('name',strcat(sps, " All Spectrum Detector:",num2str(d+1)));
    %clf;
    tx = linspace (FOVX(1),FOVX(2), FOVXSamples);
    ty = linspace (FOVY(1),FOVY(2), FOVYSamples);
    IntensityMap = IntensityMap ./ PhotonsIn ./ DetectorSize ./ SolidAngle .* DetectorRes;
    %{
    C/=PhotonPlotUpperLimit; 
    mesh(tx,ty,IntensityMap,C);
    shading("faceted");
    pbaspect([1,2,1]);
    caxis([0,PhotonPlotUpperLimit]);
    zlim([0,PhotonPlotUpperLimit]);
    title (strcat("Detector:",num2str(d+1)," Lm In:",num2str(PhotonsIn)," Lm Out:",num2str(PhotonsOut)));
    %}
    figure('name',strcat(sps, " All Spectrum Image Detector:",num2str(d+1)));
    clf;
    IntMapRot = rot90 (IntensityMap);
    imagesc(tx,ty,IntMapRot);
    pbaspect([InputAngleLimitY(2) InputAngleLimitX(2) 1])
    colormap('jet');
    title (strcat("Detector:",num2str(d+1)," Lm In:",num2str(PhotonsIn)," Lm Out:",num2str(PhotonsOut), " Ratio: ", num2str(PhotonsOut/PhotonsIn)));
    fname=strcat("Results/DetectorMatrix",num2str(d+1),sps,".csv");
    csvwrite (fname, IntensityMap);
   endfor
   endfor
   
if 1   
  for d=0:(Detectors-1)

  C=rand(FOVYSamples,FOVXSamples,3);
  PhotonsIn=0;
  PhotonsOut=0;
  for WL=1:Wavelengths
     OutputFile=strcat("Results/TestWL",num2str(WavelengthAndReleativeIntensity(WL,1)),"Counts");
     fid=fopen(strcat(OutputFile,".",num2str(d)),"r","l");
     photons=fread(fid,1,"uint64");
     printf("Detector %d:  Wavelength Index:%d PhotonsIn:%d",d,WL,photons);
     PhotonsIn+=photons * WavelengthPhotonToLm(WL);
     photons=fread(fid,1,"uint64");
     printf(" PhotonsOut:%d\n",photons);
     PhotonsOut+=photons * WavelengthPhotonToLm(WL);
     fclose(fid);
     
     OutputFile=strcat("Results/TestWL",num2str(WavelengthAndReleativeIntensity(WL,1)));
     fid=fopen(strcat(OutputFile,".",num2str(d)),"r","l");
     WM=fread(fid,[FOVXSamples,FOVYSamples],"uint64");
     WM=WM';
     if(WL==1)
      IntensityMap=WM .* WavelengthPhotonToLm(WL);
      C(:,:,1)=IntensityMap.*(WavelengthRGB(WL,1)/255.0);
      C(:,:,2)=IntensityMap.*(WavelengthRGB(WL,2)/255.0);
      C(:,:,3)=IntensityMap.*(WavelengthRGB(WL,3)/255.0);
    else
      IntensityMap+=WM .* WavelengthPhotonToLm(WL);
      C(:,:,1)+=IntensityMap.*(WavelengthRGB(WL,1)/255.0);
      C(:,:,2)+=IntensityMap.*(WavelengthRGB(WL,2)/255.0);
      C(:,:,3)+=IntensityMap.*(WavelengthRGB(WL,3)/255.0);
     endif
     fclose(fid);
    endfor

    %figure('name',strcat("All Spectrum Detector:",num2str(d+1)));
    %clf;
    tx = linspace (FOVX(1),FOVX(2), FOVXSamples);
    ty = linspace (FOVY(1),FOVY(2), FOVYSamples);  
    IntensityMap = IntensityMap ./ PhotonsIn ./ DetectorSize ./ SolidAngle .* DetectorRes;
    %{
    C/=PhotonPlotUpperLimit;
    mesh(tx,ty,IntensityMap,C);
    shading("faceted");
    pbaspect([1,2,1]);
    caxis([0,PhotonPlotUpperLimit]);
    zlim([0,PhotonPlotUpperLimit]);
    title (strcat("Detector:",num2str(d+1)," Lm In:",num2str(PhotonsIn)," Lm Out:",num2str(PhotonsOut)));
    %}
    figure('name',strcat(" All Spectrum Detector:",num2str(d+1)));
    clf;
    IntMapRot = rot90 (IntensityMap);
    imagesc(tx,ty,IntMapRot);
    pbaspect([InputAngleLimitY(2) InputAngleLimitX(2) 1])
    colormap('jet');
    title (strcat("Detector:",num2str(d+1)," Lm In:",num2str(PhotonsIn)," Lm Out:",num2str(PhotonsOut), " Ratio: ", num2str(PhotonsOut/PhotonsIn)));
    fname=strcat("Results/DetectorMatrix",num2str(d+1),".csv");
    csvwrite (fname, IntensityMap);
    endfor
   
   if  (SpatialDetector == 1)
     
   SpatialPhotonPlotUpperLimit=2500;
   C=rand(2*SpatialBoundaryWidthPixels,2*SpatialBoundaryHeightPixels,3);
   for WL=1:Wavelengths
     OutputFile=strcat("Results/TestWL",num2str(WavelengthAndReleativeIntensity(WL,1)),"SpatialReflective.All");
     fid=fopen(OutputFile,"r","l");
     WM=fread(fid,[2*SpatialBoundaryWidthPixels,2*SpatialBoundaryHeightPixels],"uint64");
     WM=WM';
     if(WL==1)
      IntensityMap=WM .* WavelengthPhotonToLm(WL);
      C(:,:,1)=IntensityMap.*(WavelengthRGB(WL,1)/255.0);
      C(:,:,2)=IntensityMap.*(WavelengthRGB(WL,2)/255.0);
      C(:,:,3)=IntensityMap.*(WavelengthRGB(WL,3)/255.0);
    else
      IntensityMap+=WM .* WavelengthPhotonToLm(WL);
      C(:,:,1)+=IntensityMap.*(WavelengthRGB(WL,1)/255.0);
      C(:,:,2)+=IntensityMap.*(WavelengthRGB(WL,2)/255.0);
      C(:,:,3)+=IntensityMap.*(WavelengthRGB(WL,3)/255.0);
     endif
     fclose(fid);
     fname=strcat("Results/SpatialReflectiveMatrix",num2str(WavelengthAndReleativeIntensity(WL,1)),".csv");
     csvwrite (fname, WM);
    endfor
     figure('name',strcat("Spatial Reflective All Wavelength"));
     clf;
     tx = linspace (-BoundaryWidth,BoundaryWidth,2*SpatialBoundaryWidthPixels);
     ty = linspace (-BoundaryHeight,BoundaryHeight,2*SpatialBoundaryHeightPixels);
     C/=SpatialPhotonPlotUpperLimit;
     mesh(tx,ty,IntensityMap,C);
     shading("faceted");
     pbaspect([1,2,1]);
     caxis([0,SpatialPhotonPlotUpperLimit]);
     zlim([0,SpatialPhotonPlotUpperLimit]);
     fname=strcat("Results/SpatialReflectiveMatrixALL.csv");
     csvwrite (fname, IntensityMap);

   for WL=1:Wavelengths
     OutputFile=strcat("Results/TestWL",num2str(WavelengthAndReleativeIntensity(WL,1)),"SpatialTransmissive.All");
     fid=fopen(OutputFile,"r","l");
     WM=fread(fid,[2*SpatialBoundaryWidthPixels,2*SpatialBoundaryHeightPixels],"uint64");
     WM=WM';
     if(WL==1)
      IntensityMap=WM .* WavelengthPhotonToLm(WL);
      C(:,:,1)=IntensityMap.*(WavelengthRGB(WL,1)/255.0);
      C(:,:,2)=IntensityMap.*(WavelengthRGB(WL,2)/255.0);
      C(:,:,3)=IntensityMap.*(WavelengthRGB(WL,3)/255.0);
    else
      IntensityMap+=WM .* WavelengthPhotonToLm(WL);
      C(:,:,1)+=IntensityMap.*(WavelengthRGB(WL,1)/255.0);
      C(:,:,2)+=IntensityMap.*(WavelengthRGB(WL,2)/255.0);
      C(:,:,3)+=IntensityMap.*(WavelengthRGB(WL,3)/255.0);
     endif
     fclose(fid);
     fname=strcat("Results/SpatialTransmissiveMatrix",num2str(WavelengthAndReleativeIntensity(WL,1)),".csv");
     csvwrite (fname, WM);
    endfor
     figure('name',strcat("Spatial Transmissive All Wavelength"));
     clf;
     tx = linspace (-BoundaryWidth,BoundaryWidth,2*SpatialBoundaryWidthPixels);
     ty = linspace (-BoundaryHeight,BoundaryHeight,2*SpatialBoundaryHeightPixels);
     C/=SpatialPhotonPlotUpperLimit;
     mesh(tx,ty,IntensityMap,C);
     shading("faceted");
     pbaspect([1,2,1]);
     caxis([0,SpatialPhotonPlotUpperLimit]);
     zlim([0,SpatialPhotonPlotUpperLimit]);
     fname=strcat("Results/SpatialTransmissiveMatrixAll.csv");
     csvwrite (fname, IntensityMap);

   for WL=1:Wavelengths
     OutputFile=strcat("Results/TestWL",num2str(WavelengthAndReleativeIntensity(WL,1)),"SpatialSTransmissive.All");
     fid=fopen(OutputFile,"r","l");
     WM=fread(fid,[2*SpatialBoundaryWidthPixels,2*SpatialBoundaryHeightPixels],"float");
     WM=WM';
     if(WL==1)
      IntensityMap=WM .* WavelengthPhotonToLm(WL);
      C(:,:,1)=IntensityMap.*(WavelengthRGB(WL,1)/255.0);
      C(:,:,2)=IntensityMap.*(WavelengthRGB(WL,2)/255.0);
      C(:,:,3)=IntensityMap.*(WavelengthRGB(WL,3)/255.0);
    else
      IntensityMap+=WM .* WavelengthPhotonToLm(WL);
      C(:,:,1)+=IntensityMap.*(WavelengthRGB(WL,1)/255.0);
      C(:,:,2)+=IntensityMap.*(WavelengthRGB(WL,2)/255.0);
      C(:,:,3)+=IntensityMap.*(WavelengthRGB(WL,3)/255.0);
     endif
     fclose(fid);
     fname=strcat("Results/SpatialSTransmissiveMatrix",num2str(WavelengthAndReleativeIntensity(WL,1)),".csv");
     csvwrite (fname, WM);
    endfor
     figure('name',strcat("Spatial S Transmissive All Wavelength"));
     clf;
     tx = linspace (-BoundaryWidth,BoundaryWidth,2*SpatialBoundaryWidthPixels);
     ty = linspace (-BoundaryHeight,BoundaryHeight,2*SpatialBoundaryHeightPixels);
     mesh(tx,ty,IntensityMap);
     shading("faceted");
     fname=strcat("Results/SpatialSTransmissiveMatrixAll.csv");
     csvwrite (fname, IntensityMap);

   for WL=1:Wavelengths
     OutputFile=strcat("Results/TestWL",num2str(WavelengthAndReleativeIntensity(WL,1)),"SpatialPTransmissive.All");
     fid=fopen(OutputFile,"r","l");
     WM=fread(fid,[2*SpatialBoundaryWidthPixels,2*SpatialBoundaryHeightPixels],"float");
     WM=WM';
     if(WL==1)
      IntensityMap=WM .* WavelengthPhotonToLm(WL);
      C(:,:,1)=IntensityMap.*(WavelengthRGB(WL,1)/255.0);
      C(:,:,2)=IntensityMap.*(WavelengthRGB(WL,2)/255.0);
      C(:,:,3)=IntensityMap.*(WavelengthRGB(WL,3)/255.0);
    else
      IntensityMap+=WM .* WavelengthPhotonToLm(WL);
      C(:,:,1)+=IntensityMap.*(WavelengthRGB(WL,1)/255.0);
      C(:,:,2)+=IntensityMap.*(WavelengthRGB(WL,2)/255.0);
      C(:,:,3)+=IntensityMap.*(WavelengthRGB(WL,3)/255.0);
     endif
     fclose(fid);
     fname=strcat("Results/SpatialPTransmissiveMatrix",num2str(WavelengthAndReleativeIntensity(WL,1)),".csv");
     csvwrite (fname, WM);
    endfor
     figure('name',strcat("Spatial P Transmissive All Wavelength"));
     clf;
     tx = linspace (-BoundaryWidth,BoundaryWidth,2*SpatialBoundaryWidthPixels);
     ty = linspace (-BoundaryHeight,BoundaryHeight,2*SpatialBoundaryHeightPixels);
     mesh(tx,ty,IntensityMap);
     shading("faceted");
     fname=strcat("Results/SpatialPTransmissiveMatrixAll.csv");
     csvwrite (fname, IntensityMap);

   for WL=1:Wavelengths
     OutputFile=strcat("Results/TestWL",num2str(WavelengthAndReleativeIntensity(WL,1)),"SpatialSReflective.All");
     fid=fopen(OutputFile,"r","l");
     WM=fread(fid,[2*SpatialBoundaryWidthPixels,2*SpatialBoundaryHeightPixels],"float");
     WM=WM';
     if(WL==1)
      IntensityMap=WM .* WavelengthPhotonToLm(WL);
      C(:,:,1)=IntensityMap.*(WavelengthRGB(WL,1)/255.0);
      C(:,:,2)=IntensityMap.*(WavelengthRGB(WL,2)/255.0);
      C(:,:,3)=IntensityMap.*(WavelengthRGB(WL,3)/255.0);
    else
      IntensityMap+=WM .* WavelengthPhotonToLm(WL);
      C(:,:,1)+=IntensityMap.*(WavelengthRGB(WL,1)/255.0);
      C(:,:,2)+=IntensityMap.*(WavelengthRGB(WL,2)/255.0);
      C(:,:,3)+=IntensityMap.*(WavelengthRGB(WL,3)/255.0);
     endif
     fclose(fid);
     fname=strcat("Results/SpatialSReflectiveMatrix",num2str(WavelengthAndReleativeIntensity(WL,1)),".csv");
     csvwrite (fname, WM);
    endfor
     figure('name',strcat("Spatial S Reflective All Wavelength"));
     clf;
     tx = linspace (-BoundaryWidth,BoundaryWidth,2*SpatialBoundaryWidthPixels);
     ty = linspace (-BoundaryHeight,BoundaryHeight,2*SpatialBoundaryHeightPixels);
     mesh(tx,ty,IntensityMap);
     shading("faceted");
     fname=strcat("Results/SpatialSReflectiveMatrixAll.csv");
     csvwrite (fname, IntensityMap);

   for WL=1:Wavelengths
     OutputFile=strcat("Results/TestWL",num2str(WavelengthAndReleativeIntensity(WL,1)),"SpatialPReflective.All");
     fid=fopen(OutputFile,"r","l");
     WM=fread(fid,[2*SpatialBoundaryWidthPixels,2*SpatialBoundaryHeightPixels],"float");
     WM=WM';
     if(WL==1)
      IntensityMap=WM .* WavelengthPhotonToLm(WL);
      C(:,:,1)=IntensityMap.*(WavelengthRGB(WL,1)/255.0);
      C(:,:,2)=IntensityMap.*(WavelengthRGB(WL,2)/255.0);
      C(:,:,3)=IntensityMap.*(WavelengthRGB(WL,3)/255.0);
    else
      IntensityMap+=WM .* WavelengthPhotonToLm(WL);
      C(:,:,1)+=IntensityMap.*(WavelengthRGB(WL,1)/255.0);
      C(:,:,2)+=IntensityMap.*(WavelengthRGB(WL,2)/255.0);
      C(:,:,3)+=IntensityMap.*(WavelengthRGB(WL,3)/255.0);
     endif
     fclose(fid);
     fname=strcat("Results/SpatialPReflectiveMatrix",num2str(WavelengthAndReleativeIntensity(WL,1)),".csv");
     csvwrite (fname, WM);
    endfor
     figure('name',strcat("Spatial P Reflective All Wavelength"));
     clf;
     tx = linspace (-BoundaryWidth,BoundaryWidth,2*SpatialBoundaryWidthPixels);
     ty = linspace (-BoundaryHeight,BoundaryHeight,2*SpatialBoundaryHeightPixels);
     mesh(tx,ty,IntensityMap);
     shading("faceted");
     fname=strcat("Results/SpatialPReflectiveMatrixAll.csv");
     csvwrite (fname, IntensityMap);

    endif
   endif
