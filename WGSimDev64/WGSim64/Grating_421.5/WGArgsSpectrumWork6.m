    OutputGratingPeriod=375.3e-9;
    WaveGuideDepth=0.5e-3;
    InputGratingHalfWidth=3e-3;
    InputGratingHalfHeight=3.5e-3;
    InputStandOff=1e-3;
    InputRadius=2e-3;
    FOVX=[-13.5,13.5];
    FOVXSamples=54;
    FOVY=[-24,24];
    FOVYSamples=96;
    SamplesMC=512;
    BoundaryWidth=20e-3;
    BoundaryHeight=20e-3;
    BoundaryRoofEnd=-10e-3;    
    BoundaryRoofIndent=12e-3;    
    InputGratingAngle=-30;
    InputGratingOffsetX=-25.0e-3;
    InputGratingRadiusSamples=16;
    InputAngleLimitX=[-13.5;13.5];
    InputAngleLimitY=[-24,24];

    Detectors=30; % count from the top left row by row
    DetectorCircular=0;
    DetectorReflectedZIsFromSurface=1;
    DetectorXYZW=[
    5e-3, 0,     23e-3,1.5e-3;
    10e-3,10e-3, 23e-3,1.5e-3;
    10e-3,5e-3,  23e-3,1.5e-3;
    10e-3,0,     23e-3,1.5e-3;
    10e-3,-5e-3, 23e-3,1.5e-3;
    10e-3,-10e-3,23e-3,1.5e-3;
    5e-3, 10e-3, 23e-3,1.5e-3;
    5e-3, 5e-3,  23e-3,1.5e-3;
    5e-3, 0,     23e-3,1.5e-3;
    5e-3, -5e-3, 23e-3,1.5e-3;
    5e-3, -10e-3,23e-3,1.5e-3;
    0,    10e-3, 23e-3,1.5e-3;
    0,    5e-3,  23e-3,1.5e-3;
    0,    0,     23e-3,1.5e-3;
    0,    -5e-3, 23e-3,1.5e-3;
    0,    -10e-3,23e-3,1.5e-3;
    10e-3,10e-3, -23e-3,1.5e-3;
    10e-3,5e-3,  -23e-3,1.5e-3;
    10e-3,0,     -23e-3,1.5e-3;
    10e-3,-5e-3, -23e-3,1.5e-3;
    10e-3,-10e-3,-23e-3,1.5e-3;
    5e-3, 10e-3, -23e-3,1.5e-3;
    5e-3, 5e-3,  -23e-3,1.5e-3;
    5e-3, 0,     -23e-3,1.5e-3;
    5e-3, -5e-3, -23e-3,1.5e-3;
    5e-3, -10e-3,-23e-3,1.5e-3;
    0,    10e-3, -23e-3,1.5e-3;
    0,    5e-3,  -23e-3,1.5e-3;
    0,    0,     -23e-3,1.5e-3;
    0,    -5e-3, -23e-3,1.5e-3;
    0,    -10e-3,-23e-3,1.5e-3;
    ];
    
    AllRayDetectors=0;
    AllRayDetectorCircular=1;
    AllRayDetectorXYW=[
    0,0,17.5e-3;
    ];
    AllRayDetectorAbsorbing=[
    0
    ];
    AllRayDetectorWHPixels=[
    128,128
    ];
    AllRayDetectorSpatialWH=[
    20e-3,20e-3
    ];
    AllRayDetectorSpatialWHPixels=[
    128,128
    ];

    
    OutputSPAngles=[34,88;0,180];    
    OutputGratingAngle=30;
    OutputGratingAngle2=-30;   
    InputGratingEfficiencyWavelength=[
    430e-9,440e-9,450e-9,460e-9,470e-9,480e-9,490e-9,500e-9,510e-9,520e-9,530e-9,540e-9,550e-9,560e-9,570e-9,580e-9,590e-9,600e-9,610e-9,620e-9,630e-9,640e-9
    ];
    InputGratingEfficiencyWavelengths=22;
    InputGratingEfficiencyAngles=[177,-88,1];
    RegionsInX=1;
    RegionsInY=1;
    RegionEfficiency=[
    100;
    ];
    OrderLimit=1;
    
    SpatialDetector = 0; 
    SpatialBoundaryWidthPixels=32;
    SpatialBoundaryHeightPixels=32;
        
    Ray3DOn=0;
    Ray3DX=0;
    Ray3DY=0;

SpectrumWavelengths=1;
SpectrumWavelengthAndReleativeIntensity=[
4.55e-007,0.89142;
];
WavelengthFolder=459;
WavelengthRefractiveIndex=1.74;
