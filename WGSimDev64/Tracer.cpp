#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <chrono>
#ifdef __linux__
#include <unistd.h>
#else
#ifndef WG_SINGLETHREAD
#undef _GLIBCXX_HAS_GTHREADS
#include "mingw.thread.h"
#include <mutex>
#include "mingw.mutex.h"
#include "mingw.condition_variable.h"
#include <atomic>
#endif
#include <assert.h>
#endif
#include <vector>
#include <stdio.h>
#include <string>
#include <sstream>
#include <iomanip>

#include "../glm/glm/glm.hpp"

std::uint64_t TotalSamples;

typedef float wgFloat;
typedef glm::vec2 wgVec2;
typedef glm::vec3 wgVec3;

#define Pi 3.14159265358979323846264338327950288419716939937510f
#define HalfPi (Pi*0.5f)
#define Pi2  (2 * Pi)
#define AngleTolerance (Pi2 / (360 * 32))


int RegionsInX,RegionsInY,AllRayDetectorCircular,DetectorCircular,DetectorReflectedZIsFromSurface,AllRayDetectorWHPixels[2],AllRayDetectorSpatialWHPixels[2];
wgFloat AllRayDetectorSpatialWH[2];

#define MaxOrdersEfficiencies 64
struct OrderEfficiencyS
{
    int Orders;
    int Order[MaxOrdersEfficiencies][2];
    wgFloat e[MaxOrdersEfficiencies][8];
};

struct OrderEfficiencyST
{
    OrderEfficiencyS *oe;
    OrderEfficiencyST()
    {
        oe = new OrderEfficiencyS[RegionsInX*RegionsInY*360*360];
    }
    OrderEfficiencyS& Get(int o,int p,int t)
    {
        return oe[o*360*360+p*360+t];
    }
    void Save(const char *EfficiencyCache)
    {
        std::string outf(EfficiencyCache);
        outf="Results/"+outf+".cch";
        std::ofstream outfile(outf,std::ofstream::binary);
        outfile.write((char *)&oe[0],RegionsInX*RegionsInY*360*360*sizeof(OrderEfficiencyS));
        outfile.close();
    }
    bool Load(const char *EfficiencyCache)
    {
        std::string inf(EfficiencyCache);
        inf="Results/"+inf+".cch";
        std::ifstream infile(inf,std::ofstream::binary);
        if(!infile.is_open()) return false;
        infile.read((char *)&oe[0],RegionsInX*RegionsInY*360*360*sizeof(OrderEfficiencyS));
        infile.close();
        return true;
    }
};

OrderEfficiencyST *OrderEfficiency=0;
OrderEfficiencyST *OrderTEfficiency=0;


//1r first interaction subsequent are 0r TIR1
wgFloat InputEfficiency[2][64*180];
wgFloat InputEfficiencyS[2][64*180];
wgFloat InputEfficiencyP[2][64*180];


struct Angle
{
    wgFloat theta, phi;
};

struct DoubleDiffractionS
{
    int Index;
    wgFloat RelIntensity;
};

const int MAXREGIONS=8*8;

struct AngleInfoS
{
    Angle Incident;
    bool InputGratingAngle;
    unsigned char Exits;
    wgVec2 Offset;
    wgFloat InputEfficiency[2];
    wgFloat InputEfficiencyS[2];
    wgFloat InputEfficiencyP[2];
    int Angles1;
    Angle DiffractedAngle1[32];
    wgFloat DiffractedAngleEfficiency1[MAXREGIONS][32];
    wgFloat DiffractedAngleSPEfficiency1[MAXREGIONS][32][4];
    wgFloat DiffractedAngleTEfficiency1[MAXREGIONS][32];
    wgFloat DiffractedAngleTSPEfficiency1[MAXREGIONS][32][4];
    int DiffractedAngleInfoIndex1[32];
    int DiffractedAngleOrder1[32][2];
};

enum { MaxAngleInfos = 32, MaxDetectors=64, MaxWavelengths=64  };

struct RayTracerInfoS
{
    int AngleInfos;
    AngleInfoS AngleInfo[MaxAngleInfos];
    wgFloat TotalInputIntensity[MaxDetectors],TotalS[MaxDetectors],TotalP[MaxDetectors], WaveLength, TotalRays, TotalExitRays;
    wgFloat S,P;
    int MaxDepth, ExitRayIndex;
    Angle InputGratingAngle, ExitAngle, AirExitAngle;
    wgVec2 AirExitAngleOffset,ExitAngleOffset;
    std::uint64_t DetectorOn;
    wgVec2 DetectorProjected[MaxDetectors][4];
    bool DetectorProjectedTransmissive[MaxDetectors];
    int ThreadID;
    unsigned int RandomIndex[4];
    std::uint64_t SamplesMC;
    wgFloat InputEfficiency,TIREfficiency;
    volatile int Progress;
    volatile bool Finished;
};

wgFloat WavelengthRefractiveIndex,CriticalAngle;
wgFloat MinDistance, MaxDistance;
wgFloat InputStandOff, InputRadius, OutputGratingPeriod, InputGratingOffsetX, InputGratingAngle, OutputGratingAngle, OutputGratingAngle2;
wgFloat WaveGuideDepth;
wgFloat InputAngleLimitX[2], InputAngleLimitY[2];
wgFloat FOVX[2];
wgFloat FOVY[2];
wgFloat InputGratingHalfWidth, InputGratingHalfHeight;
bool Halt;
int SpatialBoundaryWidthPixels,SpatialBoundaryHeightPixels;
int FOVXSamples;
int FOVYSamples;
int InputGratingRadiusSamples;
int SceneWidth, SceneHeight;
int MaxUseableThreads;
int SamplesMC;
int Detectors;
wgFloat DetectorXYZW[MaxDetectors*4];
int AllRayDetectors;
wgFloat AllRayDetectorXYW[MaxDetectors*3];
int SpectrumWavelengths;
wgFloat SpectrumWavelengthAndReleativeIntensity[MaxWavelengths*2];
int Ray3DOn;
wgFloat Ray3DX,Ray3DY;
int OrderLimit;
int InputGratingEfficiencyWavelengths;
wgFloat InputGratingEfficiencyWavelength[128];
int InputGratingEfficiencyAngles[3];
int OutputSPAngles[4];
int WavelengthFolder,RegionEfficiency[64*64];

int AllRayDetectorAbsorbing[MaxDetectors]= {0};
wgFloat *AllRayDetectorPixel[MaxDetectors]= {0};
wgFloat *AllRayDetectorSpatialPixel=0;
std::uint64_t *Pixel[MaxDetectors]= {0};
std::uint64_t *SpatialPixel[2]= {0,0};
wgFloat *SpatialSP[2]= {0,0};
wgFloat *PixelS[MaxDetectors]= {0};
wgFloat *PixelP[MaxDetectors]= {0};

wgFloat BoundaryWidth, BoundaryHeight;

wgFloat BoundaryRoofEnd,BoundaryRoofIndent;


uint32_t xorshift128plus(uint64_t *s)
{
    uint64_t x = s[0];
    uint64_t const y = s[1];
    s[0] = y;
    x ^= x << 23;
    s[1] = x ^ y ^ (x >> 17) ^ (y >> 26);
    return (s[1]+y);
}

unsigned int Random32(unsigned int *rndx)
{
    rndx[0]=rndx[0]*1664525+1013904223;
    return rndx[0];
}

wgFloat MCRandom(unsigned int *rndx)
{
    const wgFloat f=1.0f/wgFloat(0xffffffffu);
    //wgFloat t = Random32(rndx)*f;
    wgFloat t = xorshift128plus((uint64_t *)rndx)*f;
    return t;
}

bool AllRayIntensityCheck(wgFloat Intensity,wgFloat Theta, wgFloat Phi, const wgVec2 *Distance, RayTracerInfoS *p)
{
    if(!AllRayDetectors) return false;
    int Absorbed=0;
    Theta/=Pi*0.5f;
    wgFloat sx=Distance->x/AllRayDetectorSpatialWH[0],sy=Distance->y/AllRayDetectorSpatialWH[1];
    if((fabs(sx)<=1)&&(fabs(sy)<=1))
    {
        int x=int(sx*(AllRayDetectorSpatialWHPixels[0]-1)+AllRayDetectorSpatialWHPixels[0]);
        int y=int(sy*(AllRayDetectorSpatialWHPixels[1]-1)+AllRayDetectorSpatialWHPixels[1]);
        AllRayDetectorSpatialPixel[x+y*AllRayDetectorSpatialWHPixels[0]*2] += Intensity;
    }
    for(int i=0; i<AllRayDetectors; i++)
    {
        if(DetectorCircular)
        {
            wgFloat dx=(Distance->x-AllRayDetectorXYW[i*3+0]);
            wgFloat dy=(Distance->y-AllRayDetectorXYW[i*3+1]);
            wgFloat dr=AllRayDetectorXYW[i*3+2];
            if((dx*dx+dy*dy)>(dr*dr)) continue;
        }
        else
        {
            if (Distance->y < (AllRayDetectorXYW[i*3+1]-AllRayDetectorXYW[i*3+2])) continue;
            if (Distance->y > (AllRayDetectorXYW[i*3+1]+AllRayDetectorXYW[i*3+2])) continue;
            if (Distance->x < (AllRayDetectorXYW[i*3+0]-AllRayDetectorXYW[i*3+2])) continue;
            if (Distance->x > (AllRayDetectorXYW[i*3+0]+AllRayDetectorXYW[i*3+2])) continue;
        }
        int x=int(Theta*cos(Phi)*(AllRayDetectorWHPixels[0]-1)+AllRayDetectorWHPixels[0]);
        int y=int(Theta*sin(Phi)*(AllRayDetectorWHPixels[1]-1)+AllRayDetectorWHPixels[1]);
        /*
        if((x<0)||(x>(AllRayDetectorWHPixels[0]*2))||(y<0)||(y>(AllRayDetectorWHPixels[0]*2)))
        {
            std::cout<<"theta:"<<Theta<<" Phi:"<<Phi<<" x:"<<x<<" y:"<<y<<"\n";
        }
        else
        */
        AllRayDetectorPixel[i][x+y*AllRayDetectorWHPixels[0]*2] += Intensity;


        if(AllRayDetectorAbsorbing[i]) Absorbed++;
    }
    bool r=(Absorbed==AllRayDetectors);
    //static int t=0;
    //if(r) std::cout<<"Absorbed:"<<t++<<"\n";
    return r;
}

void IntensitySPCheck(wgFloat Intensity,wgFloat S, wgFloat P, const wgVec2 *Distance, RayTracerInfoS *p,bool Transmissive=false)
{
    for(int i=0; i<Detectors; i++)
        if(p->DetectorOn&(1<<i))
        {
            if(DetectorCircular)
            {
                if(Distance->x> BoundaryWidth) continue;
                if(Distance->x<-BoundaryWidth) continue;
                if(Distance->y> BoundaryHeight) continue;
                if(Distance->y<-BoundaryHeight) continue;
                wgFloat dx=(Distance->x-p->DetectorProjected[i][2].x);
                wgFloat dy=(Distance->y-p->DetectorProjected[i][2].y);
                wgFloat dr=p->DetectorProjected[i][3].x;
                if((dx*dx+dy*dy)>(dr*dr)) continue;
            }
            else
            {
                if (Distance->y < p->DetectorProjected[i][0].y) continue;
                if (Distance->y > p->DetectorProjected[i][1].y) continue;
                if (Distance->x < p->DetectorProjected[i][0].x) continue;
                if (Distance->x > p->DetectorProjected[i][1].x) continue;
            }
            if(Transmissive==p->DetectorProjectedTransmissive[i]) p->TotalInputIntensity[i] += Intensity;
            if(Transmissive==p->DetectorProjectedTransmissive[i]) p->TotalS[i] += S;
            if(Transmissive==p->DetectorProjectedTransmissive[i]) p->TotalP[i] += P;
        }
}

inline wgFloat ToDeg(wgFloat a)
{
    return a * 90 / HalfPi;
}
inline wgFloat ToRad(wgFloat a)
{
    return a*HalfPi / 90;
}

bool AngleLess(const Angle *a, const Angle *b)
{
    wgFloat t[2] = { a->theta, b->theta };
    wgFloat p[2] = { a->phi, b->phi };
    if (t[0] < t[1]) return true;
    if (t[0] == t[1])
    {
        if (p[0] < p[1]) return true;
    }
    return false;

}

bool AngleEqual(const Angle *a, const Angle *b)
{
    if (fabs(a->theta - b->theta)>AngleTolerance) return false;
    if (fabs(a->theta) <= AngleTolerance) return true;
    if (fabs(a->phi - b->phi)>AngleTolerance) return false;
    return true;
}

bool AngleNotEqual(const Angle *a, const Angle *b)
{
    return (!AngleEqual(a, b));
}
#define AngleEvanescent 1e20f

Angle DiffractedAngleOfOrder2(int OrderM, int OrderN,Angle IncidentAngle, wgFloat GratingPeriod, wgFloat GratingSkewAngle, wgFloat WaveLength, wgFloat RefractiveIndex,bool Cross=true)
{
    Angle r;
    wgFloat lambda=WaveLength;
    wgFloat nmaterial=RefractiveIndex;

    wgFloat k=2*Pi*nmaterial/lambda;

    wgFloat l1=GratingPeriod;
    wgFloat l2;
    if(Cross) l2=GratingPeriod;
    else l2=std::numeric_limits<float>::infinity();

    wgFloat K1=2*Pi/l1;
    wgFloat K2=2*Pi/l2;

    wgFloat skew=GratingSkewAngle;

    wgFloat phiin=IncidentAngle.phi;
    wgFloat thetain=IncidentAngle.theta;

    wgFloat alpha0=k*sin(thetain)*cos(phiin);
    wgFloat beta0=k*sin(thetain)*sin(phiin+skew);
    //wgFloat gamma00=k*cos(thetain);

    wgFloat alpham=alpha0+OrderM*K1;
    wgFloat betan=beta0+OrderN*K2;

    wgFloat gammamn2=pow(k,2)-pow((1/cos(skew)),2) * (pow(alpham,2) + pow(betan,2) - 2*alpham*betan*sin(skew));
    if(gammamn2<0)
    {
        r.theta = AngleEvanescent;
        return r;
    }
    wgFloat gammamn=sqrt(gammamn2);

    r.theta=acos(gammamn/k);

    r.phi=atan2(betan-alpham*sin(skew),alpham*cos(skew));

    return r;
}



Angle Refract(const Angle *Incident, wgFloat FromRefractiveIndex, wgFloat ToRefractiveIndex)
{
    Angle r = *Incident;
    r.theta = asinf(FromRefractiveIndex*sinf(Incident->theta) / ToRefractiveIndex);
    return r;
}

void ProjectDetector(const Angle *InputGIncident, RayTracerInfoS *p)
{
    p->DetectorOn=0;
    for(int i=0; i<Detectors; i++)
    {
        wgFloat dx=DetectorXYZW[i*4+0],dy=DetectorXYZW[i*4+1],dz=DetectorXYZW[i*4+2],dw=DetectorXYZW[i*4+3];
        p->DetectorProjectedTransmissive[i]=(dz<0);
        wgFloat eyel = tanf(InputGIncident->theta)*std::abs(dz);
        p->DetectorProjected[i][2] = wgVec2(dx - eyel*cos(InputGIncident->phi), dy - eyel*sin(InputGIncident->phi));
        p->DetectorProjected[i][3] = wgVec2(dw, dw);
        p->DetectorProjected[i][0] = p->DetectorProjected[i][2] - p->DetectorProjected[i][3];
        p->DetectorProjected[i][1] = p->DetectorProjected[i][0] + 2.0f*p->DetectorProjected[i][3];
        if (p->DetectorProjected[i][0].x > BoundaryWidth) continue;
        if (p->DetectorProjected[i][1].x < -BoundaryWidth) continue;
        if (p->DetectorProjected[i][0].y > BoundaryHeight) continue;
        if (p->DetectorProjected[i][1].y < -BoundaryHeight) continue;
        if (p->DetectorProjected[i][0].x < -BoundaryWidth) p->DetectorProjected[i][0].x = -BoundaryWidth;
        if (p->DetectorProjected[i][1].x > BoundaryWidth) p->DetectorProjected[i][1].x = BoundaryWidth;
        if (p->DetectorProjected[i][0].y < -BoundaryHeight) p->DetectorProjected[i][0].y = -BoundaryHeight;
        if (p->DetectorProjected[i][1].y > BoundaryHeight) p->DetectorProjected[i][1].y = BoundaryHeight;
        p->DetectorOn|=(1<<i);
    }
}

void CalcInputEfficiency(AngleInfoS &ai, RayTracerInfoS *p)
{
    Angle a=ai.Incident;
    int d=int(ToDeg(atan2(ai.Offset.x,2*WaveGuideDepth)));
    int d2=int(ToDeg(a.theta*cos(a.phi)));
    const int DisplayDBG=0;//ai.Exits;
    if(DisplayDBG) std::cout<<"Angle In:("<<ToDeg(a.theta)<<","<<ToDeg(a.phi)<<") Angle Out:"<<d<<"("<<d2<<")";
    d-=InputGratingEfficiencyAngles[1];
    d/=InputGratingEfficiencyAngles[2];
    if(DisplayDBG) std::cout<<" Index:"<<d;
    if(d<0) d=0;
    if(d>=InputGratingEfficiencyAngles[0]) d=InputGratingEfficiencyAngles[0]-1;
    int fw=0;
    wgFloat w=p->WaveLength;
    for(int i=1; i<InputGratingEfficiencyWavelengths; i++)
        if(w<=InputGratingEfficiencyWavelength[i])
        {
            if((w-InputGratingEfficiencyWavelength[i-1])<(InputGratingEfficiencyWavelength[i]-w))
                fw=i-1;
            else
                fw=i;
            break;
        }
    if(DisplayDBG) std::cout<<" Wavelength:("<<w<<") Index:"<<fw;
    ai.InputEfficiency[0]=InputEfficiency[0][fw*InputGratingEfficiencyAngles[0]+d];
    ai.InputEfficiency[1]=InputEfficiency[1][fw*InputGratingEfficiencyAngles[0]+d];
    ai.InputEfficiencyS[0]=InputEfficiencyS[0][fw*InputGratingEfficiencyAngles[0]+d];
    ai.InputEfficiencyS[1]=InputEfficiencyS[1][fw*InputGratingEfficiencyAngles[0]+d];
    ai.InputEfficiencyP[0]=InputEfficiencyS[0][fw*InputGratingEfficiencyAngles[0]+d];
    ai.InputEfficiencyP[1]=InputEfficiencyS[1][fw*InputGratingEfficiencyAngles[0]+d];
    if(DisplayDBG) std::cout<<" Efficiency:("<<ai.InputEfficiency[0]<<","<<ai.InputEfficiency[1]<<")\n";
}

int PreRayTraceAdd1(const Angle *AngleIncident, RayTracerInfoS *p)
{
    assert(ToDeg(AngleIncident->theta>=0));
    assert(ToDeg(AngleIncident->theta<90));
    for (int i = 0; i < p->AngleInfos; i++) if (AngleEqual(&p->AngleInfo[i].Incident, AngleIncident)) return i;
    assert(p->AngleInfos < MaxAngleInfos);
    int f = p->AngleInfos;
    p->AngleInfos++;
    p->AngleInfo[f].Angles1 = 0;
    p->AngleInfo[f].Incident = *AngleIncident;
    p->AngleInfo[f].Exits = (fabs(AngleIncident->theta) < CriticalAngle);

    wgVec2 DistanceTravelled;
    DistanceTravelled.x = (WaveGuideDepth)*tan(AngleIncident->theta)*cos(AngleIncident->phi);
    DistanceTravelled.y = (WaveGuideDepth)*tan(AngleIncident->theta)*sin(AngleIncident->phi);

    p->AngleInfo[f].Offset = 2.0f*DistanceTravelled;
    CalcInputEfficiency(p->AngleInfo[f],p);
    if (p->AngleInfo[f].Exits)
    {
        p->InputEfficiency=p->AngleInfo[f].InputEfficiency[1];
        p->S=p->AngleInfo[f].InputEfficiencyS[1];
        p->P=p->AngleInfo[f].InputEfficiencyP[1];
        if (fabs(AngleIncident->theta - p->ExitAngle.theta) > AngleTolerance)
            p->AngleInfo[f].Exits |= 2;
        else
            p->ExitRayIndex |= (1 << f);
    }
    return f;
}

inline int TOrderIndex(int x,int y)
{
    int OrderBound=OrderLimit;
    int n=((x+OrderBound)*((2*OrderBound)+1))+(y+OrderBound);
    return n;
}

const int SkewOffset=30;
void PreRayTraceMC1(Angle AngleIncident, RayTracerInfoS *p)
{
    const int aoffset=SkewOffset;
    int SIndex = 0;
    PreRayTraceAdd1(&AngleIncident, p);
    while (SIndex < p->AngleInfos)
    {
        int f = SIndex;
        AngleIncident = p->AngleInfo[f].Incident;
        p->AngleInfo[f].Angles1=0;

        int phi=ToDeg(AngleIncident.phi),theta=ToDeg(AngleIncident.theta);
        phi=(phi+aoffset)%360;
        if(phi<0) phi+=360;
        theta=(theta%360);
        if(theta<0) theta+=360;
        assert(phi>=0);
        assert(phi<360);
        assert(theta>=0);
        assert(theta<360);

        Angle ra=AngleIncident;
        ra.phi=ra.phi+ToRad(aoffset);

        OrderEfficiencyS &OE=OrderEfficiency->Get(0,phi,theta);
        for(int i=0; i<OE.Orders; i++)
        {
            //if(e<=1e-12) continue;
            //if(std::abs(OE.Order[i][0])>1) continue;
            //if(std::abs(OE.Order[i][1])>1) continue;
            Angle a = DiffractedAngleOfOrder2(OE.Order[i][0], OE.Order[i][1],ra, OutputGratingPeriod, OutputGratingAngle, p->WaveLength, WavelengthRefractiveIndex);
            if (a.theta < AngleEvanescent)
            {
                a.phi=a.phi-ToRad(aoffset);
                for(int j=0; j<(RegionsInX*RegionsInY); j++)
                {
                    p->AngleInfo[f].DiffractedAngleEfficiency1[j][p->AngleInfo[f].Angles1]=OrderEfficiency->Get(j,phi,theta).e[i][4];
                    for(int k=0; k<4; k++) p->AngleInfo[f].DiffractedAngleSPEfficiency1[j][p->AngleInfo[f].Angles1][k]=OrderEfficiency->Get(j,phi,theta).e[i][k];
                    OrderEfficiencyS &oe=OrderEfficiency->Get(j,phi,theta);
                    int oi=TOrderIndex(oe.Order[i][0],oe.Order[i][1]);
                    p->AngleInfo[f].DiffractedAngleTEfficiency1[j][p->AngleInfo[f].Angles1]=OrderTEfficiency->Get(j,phi,theta).e[oi][4];
                    for(int k=0; k<4; k++) p->AngleInfo[f].DiffractedAngleTSPEfficiency1[j][p->AngleInfo[f].Angles1][k]=OrderTEfficiency->Get(j,phi,theta).e[i][k];
                }
                p->AngleInfo[f].DiffractedAngle1[p->AngleInfo[f].Angles1] = a;
                p->AngleInfo[f].DiffractedAngleInfoIndex1[p->AngleInfo[f].Angles1] = PreRayTraceAdd1(&a, p);
                p->AngleInfo[f].DiffractedAngleOrder1[p->AngleInfo[f].Angles1][0] = OE.Order[i][0];
                p->AngleInfo[f].DiffractedAngleOrder1[p->AngleInfo[f].Angles1][1] = OE.Order[i][1];
                p->AngleInfo[f].Angles1++;
            }
        }

        SIndex++;
    }

    if(Ray3DOn)
        for(int f=0; f<p->AngleInfos; f++)
        {
            const int aoffset3=aoffset;
            Angle ra=p->AngleInfo[f].Incident;
            ra.phi=ra.phi+ToRad(aoffset3);
            std::cout<<"\n"<<f<<":Input Angle "<<ToDeg(p->AngleInfo[f].Incident.theta)<<","<<ToDeg(p->AngleInfo[f].Incident.phi)<<"\n";
            std::cout<<"    Exit Angle "<<ToDeg(p->ExitAngle.theta)<<","<<ToDeg(p->ExitAngle.phi)<<"\n";
            std::cout<<"    Exit Flags:"<<int(p->AngleInfo[f].Exits)<<"\n";

            int phi=ToDeg(p->AngleInfo[f].Incident.phi),theta=ToDeg(p->AngleInfo[f].Incident.theta);
            phi=(phi+aoffset)%360;
            if(phi<0) phi+=360;
            theta=(theta%360);
            if(theta<0) theta+=360;
            assert(phi>=0);
            assert(phi<360);
            assert(theta>=0);
            assert(theta<360);

            OrderEfficiencyS &OE=OrderEfficiency->Get(0,phi,theta);
            std::cout<<theta<<","<<phi<<" Efficiency Orders:"<<OE.Orders<<"\n";
            for(int i=0; i<OE.Orders; i++)
                std::cout<<"Order:"<<OE.Order[i][0]<<","<<OE.Order[i][1]<<" e:"<<OE.e[i]<<"\n";
            for(int i=0; i<p->AngleInfo[f].Angles1; i++)
            {
                std::cout<<i<<": "<<ToDeg(p->AngleInfo[f].DiffractedAngle1[i].theta)<<","<<ToDeg(p->AngleInfo[f].DiffractedAngle1[i].phi)<<" Index:"<<p->AngleInfo[f].DiffractedAngleInfoIndex1[i];
                std::cout<<" Order:"<<p->AngleInfo[f].DiffractedAngleOrder1[i][0]<<","<<p->AngleInfo[f].DiffractedAngleOrder1[i][1]<<"\n";
            }
        }


}

void LineNull(wgVec3 a,int LineType)
{

}

std::vector<int> Points3D;
std::vector<wgVec3> Point3D;

void Line3D(wgVec3 a,int LineType)
{
    static int n=0;
    n++;
    switch(LineType)
    {
    case -1:
        Points3D.push_back(-n);
        n=0;
        break;
    case  1:
        Points3D.push_back(n);
        n=0;
        break;
    default:
        break;
    }
    Point3D.push_back(a);
}

wgFloat NormaliserSP(wgFloat S,wgFloat P,wgFloat Scale=1)
{
    //return 1;
    return Scale/(S+P);
}

template <typename T>
bool RayTraceTIRMC1(int AngleIndex, wgFloat Intensity, int &Depth, wgVec2 &Distance, RayTracerInfoS *p,T Line)
{
    wgVec2 c(InputGratingOffsetX,0);
    wgFloat e=p->InputEfficiency,PolS=p->S,PolP=p->P;
    wgVec2 DistanceIncidentTravelled = p->AngleInfo[AngleIndex].Offset;
    if (DistanceIncidentTravelled.x <= 0) return false;
    while(Distance.x<-BoundaryWidth)
    {
        if(AllRayIntensityCheck(e*SamplesMC,p->AngleInfo[AngleIndex].Incident.theta,p->AngleInfo[AngleIndex].Incident.phi,&Distance,p)) return false;
        Line(wgVec3(Distance,0),0);
        Line(wgVec3(Distance + p->AngleInfo[AngleIndex].Offset, 2.0f*WaveGuideDepth),0);
        Distance+=2.0f*DistanceIncidentTravelled;
        wgVec2 l=(Distance-c);
        if(std::abs(l.x)<=InputGratingHalfWidth)
            if(std::abs(l.y)<=InputGratingHalfHeight)
            {
                e*=p->AngleInfo[AngleIndex].InputEfficiency[0];
                wgFloat ts=(PolS+PolP)*p->AngleInfo[AngleIndex].InputEfficiencyS[0];
                wgFloat tp=(PolS+PolP)*p->AngleInfo[AngleIndex].InputEfficiencyP[0];
                PolS=ts;
                PolP=tp;
            }
    }
    if(AllRayIntensityCheck(e*SamplesMC,p->AngleInfo[AngleIndex].Incident.theta,p->AngleInfo[AngleIndex].Incident.phi,&Distance,p)) return false;
    p->TIREfficiency=e;
    wgFloat spl=NormaliserSP(PolS,PolP);
    p->S=PolS*spl;
    p->P=PolP*spl;
    //std::cout<<"Input Efficiency:"<<p->InputEfficiency<<","<<p->TIREfficiency<<"\n";
    Line(wgVec3(Distance,0), 1);
    if (fabs(Distance.y)>BoundaryHeight)
    {
        return false;
    }

    return true;
}

template <typename T>
void RayTraceMC1(int AngleIndex, wgFloat Intensity, int Depth, wgVec2 Distance, RayTracerInfoS *p,T Line)
{

    Line(wgVec3(Distance,0), 0);
    while (true)
    {
        if ((Distance.x > (BoundaryWidth))||(Distance.x < (-BoundaryWidth))||(Distance.y < -BoundaryHeight)||(Distance.y > BoundaryHeight))
        {
            Line(wgVec3(Distance + p->AngleInfo[AngleIndex].Offset, 2.0f*WaveGuideDepth),1);
            return;
        }
        if(Distance.x<BoundaryRoofEnd)
        {
            if(std::abs(Distance.y)>((BoundaryHeight-BoundaryRoofIndent)+((BoundaryRoofIndent/(BoundaryRoofEnd+BoundaryWidth))*(Distance.x+BoundaryWidth))))
            {
                //std::cout<<Distance.x<<","<<Distance.y<<":"<<int((BoundaryWidth-Distance.x)*1e3)<<","<<int((BoundaryHeight-Distance.y)*1e3)<<"\n";
                //assert((std::abs(Distance.y)>(BoundaryHeight-BoundaryRoofIndent)));
                Line(wgVec3(Distance + p->AngleInfo[AngleIndex].Offset, 2.0f*WaveGuideDepth),1);
                return;
            }
        }

        wgFloat r=MCRandom(&p->RandomIndex[0])/**p->AngleInfo[AngleIndex].TotalEfficiency1*/,tr=0;
        int f=0,Region=0;

        int rx=int((Distance.x-(-BoundaryWidth))/(2*BoundaryWidth/RegionsInX));
        int ry=int((Distance.y-(-BoundaryHeight))/(2*BoundaryHeight/RegionsInY));
        Region=rx*RegionsInY+ry;

        for(int i=0; i<p->AngleInfo[AngleIndex].Angles1; i++)
        {
            tr+=p->AngleInfo[AngleIndex].DiffractedAngleEfficiency1[Region][i];
            if(r<=tr)
            {
                wgFloat S=p->S*p->AngleInfo[AngleIndex].DiffractedAngleSPEfficiency1[Region][i][0]+p->P*p->AngleInfo[AngleIndex].DiffractedAngleSPEfficiency1[Region][i][2];
                wgFloat P=p->S*p->AngleInfo[AngleIndex].DiffractedAngleSPEfficiency1[Region][i][1]+p->P*p->AngleInfo[AngleIndex].DiffractedAngleSPEfficiency1[Region][i][3];
                wgFloat spl=NormaliserSP(S,P);
                p->S=S*spl;
                p->P=P*spl;
                AngleIndex=p->AngleInfo[AngleIndex].DiffractedAngleInfoIndex1[i];
                f=1;
                break;
            }

            int na=p->AngleInfo[AngleIndex].DiffractedAngleInfoIndex1[i];
            if (p->AngleInfo[na].Exits)
                if (!(p->AngleInfo[na].Exits&2))
                {
                    tr+=p->AngleInfo[AngleIndex].DiffractedAngleTEfficiency1[Region][i];
                    if(r<=tr)
                    {
                        wgFloat S=p->S*p->AngleInfo[AngleIndex].DiffractedAngleTSPEfficiency1[Region][i][0]+p->P*p->AngleInfo[AngleIndex].DiffractedAngleTSPEfficiency1[Region][i][2];
                        wgFloat P=p->S*p->AngleInfo[AngleIndex].DiffractedAngleTSPEfficiency1[Region][i][1]+p->P*p->AngleInfo[AngleIndex].DiffractedAngleTSPEfficiency1[Region][i][3];
                        wgFloat spl=NormaliserSP(S,P);
                        p->S=S*spl;
                        p->P=P*spl;
                        AngleIndex=na;
                        f=-1;
                        break;
                    }
                }



        }
        if(AllRayIntensityCheck(Intensity,p->AngleInfo[AngleIndex].Incident.theta,p->AngleInfo[AngleIndex].Incident.phi,&Distance,p)) return;;
        if(!f)
        {
            Line(wgVec3(Distance + p->AngleInfo[AngleIndex].Offset, 2.0f*WaveGuideDepth),1);
            return;
        }

        if (p->AngleInfo[AngleIndex].Exits)
        {
            if (p->AngleInfo[AngleIndex].Exits&2) return;
            if(DetectorReflectedZIsFromSurface)
                if(f==1) Distance += p->AngleInfo[AngleIndex].Offset;
            int x=(((Distance.x+BoundaryWidth)*(2*SpatialBoundaryWidthPixels))/(2*BoundaryWidth));
            int y=(((Distance.y+BoundaryHeight)*(2*SpatialBoundaryHeightPixels))/(2*BoundaryHeight));
            SpatialPixel[(f==-1)?1:0][x+(y*2*SpatialBoundaryWidthPixels)]++;
            wgFloat spl=NormaliserSP(p->S,p->P);
            SpatialSP[(f==-1)?1:0][(x+(y*2*SpatialBoundaryWidthPixels))]+=p->S*spl;
            SpatialSP[(f==-1)?1:0][(SpatialBoundaryWidthPixels*2*SpatialBoundaryHeightPixels*2)+(x+(y*2*SpatialBoundaryWidthPixels))]+=p->P*spl;
            wgFloat l=-16*WaveGuideDepth;
            wgVec2 dt = float(l)*p->AirExitAngleOffset;
            Line(wgVec3(Distance + dt, l),-1);
            IntensitySPCheck(Intensity,p->S,p->P,&Distance,p,f==-1);
            return;
        }

        Line(wgVec3(Distance + p->AngleInfo[AngleIndex].Offset, 2.0f*WaveGuideDepth),0);
        Line(wgVec3(Distance+2.0f*p->AngleInfo[AngleIndex].Offset,0), 0);
        Distance += 2.0f*p->AngleInfo[AngleIndex].Offset;

    }
}

int OneTraceSamples;
void InitialiseRayTrace(void)
{
    wgFloat rs;
    OneTraceSamples=0;
    if (InputGratingRadiusSamples) rs = InputRadius / InputGratingRadiusSamples;
    else rs = 0;
    for (int igy = -InputGratingRadiusSamples; igy <= InputGratingRadiusSamples; igy++)
        for (int igx = -InputGratingRadiusSamples; igx <= InputGratingRadiusSamples; igx++)
        {
            wgFloat digx = igx;
            wgFloat digy = igy;

            wgFloat dx = digx*rs, dy = digy*rs, rad = sqrtf(dx*dx + dy*dy);
            if (rad <= InputRadius)
            {
                OneTraceSamples++;
            }
        }
}

void FinaliseRayTrace(void)
{
}



template <typename T>
void OneTrace(wgFloat ai, wgFloat aj, RayTracerInfoS *p,T Line)
{

    TotalSamples+=OneTraceSamples*p->SamplesMC;

    for(int i=0; i<Detectors; i++)
    {
        p->TotalInputIntensity[i] = 0;
        p->TotalS[i] = 0;
        p->TotalP[i] = 0;
    };
    wgFloat a;
    wgFloat phi = atan2f(aj, ai);
    if (fabs(ai) >= fabs(aj)) a = ai / cosf(phi);
    else a = aj / sinf(phi);
    Angle InputGIncident = { a, phi };
    ProjectDetector(&InputGIncident, p);
    if (!p->DetectorOn)
    {
        std::cout << "No Projected Detector At:" << ToDeg(ai) << "," << ToDeg(aj) << " For:"<<(p->DetectorOn^((1<<Detectors)-1))<< std::endl;
        return;
    }
    p->AirExitAngle = InputGIncident;
    wgVec2 DistanceTravelled;
    DistanceTravelled.x = tanf(p->AirExitAngle.theta)*cosf(p->AirExitAngle.phi);
    DistanceTravelled.y = tanf(p->AirExitAngle.theta)*sinf(p->AirExitAngle.phi);
    p->AirExitAngleOffset = DistanceTravelled;
    InputGIncident = Refract(&InputGIncident,1, WavelengthRefractiveIndex);
    p->ExitAngle = InputGIncident;
    DistanceTravelled.x = tanf(p->ExitAngle.theta)*cosf(p->ExitAngle.phi);
    DistanceTravelled.y = tanf(p->ExitAngle.theta)*sinf(p->ExitAngle.phi);
    p->ExitAngleOffset = DistanceTravelled;

    p->AngleInfos = 0;
#if 0
    p->InputGratingAngle = DiffractedAngleOfOrder(1, InputGIncident, OutputGratingPeriod, InputGratingAngle, p->WaveLength, WavelengthRefractiveIndex);
#else
    InputGIncident.phi+=ToRad(SkewOffset);
    p->InputGratingAngle = DiffractedAngleOfOrder2(1, 0,InputGIncident, OutputGratingPeriod, InputGratingAngle, p->WaveLength, WavelengthRefractiveIndex,false);
    InputGIncident.phi-=ToRad(SkewOffset);
    p->InputGratingAngle.phi-=ToRad(SkewOffset);
#endif
    if (p->InputGratingAngle.theta >= AngleEvanescent)
    {
        //std::cout << "Evanescent ray at:" << ToDeg(ai) << "," << ToDeg(aj) << "  "<<ToDeg(a)<<","<<ToDeg(phi)<<std::endl;
        return;
    }
    if(fabs(p->InputGratingAngle.theta) < CriticalAngle)
    {
        //std::cout << "Critical ray at:" << ToDeg(ai) << "," << ToDeg(aj) << "  "<<ToDeg(a)<<","<<ToDeg(phi)<<std::endl;
        return;
    }
    p->ExitRayIndex = 0;

    PreRayTraceMC1(p->InputGratingAngle, p);


    if (!p->ExitRayIndex)
    {
        std::cout << "No exit ray at:" << ToDeg(ai) << "," << ToDeg(aj) << "  "<<ToDeg(a)<<","<<ToDeg(phi)<<std::endl;
        return;
    }

    wgVec2 InputStandOffOffset=-(p->ExitAngleOffset*(WaveGuideDepth*2)+p->AirExitAngleOffset*InputStandOff);

    wgFloat rs;
    if (InputGratingRadiusSamples) rs = InputRadius / InputGratingRadiusSamples;
    else rs = 0;
    for (int igy = -InputGratingRadiusSamples; igy <= InputGratingRadiusSamples; igy++)
        for (int igx = -InputGratingRadiusSamples; igx <= InputGratingRadiusSamples; igx++)
        {
            wgFloat digx = igx;
            wgFloat digy = igy;

            wgFloat dx = digx*rs, dy = digy*rs, rad = sqrtf(dx*dx + dy*dy);
            if (rad <= InputRadius)
            {
                wgFloat t = atan2f(dy, dx);
                wgVec2 g(rad*cos(t), rad*sinf(t));
                g+=InputStandOffOffset;
                if(fabs(g.x)>InputGratingHalfWidth) continue;
                if(fabs(g.y)>InputGratingHalfHeight) continue;
                g.x+=InputGratingOffsetX;
                int d=0;
                if(RayTraceTIRMC1(0, 1, d, g, p,Line))
                    for(std::uint64_t i=0; i<(p->SamplesMC*p->TIREfficiency); i++) RayTraceMC1(0, 1, d, g, p,Line);
                //else
                //std::cout << "No Input TIR out ray at:" << ToDeg(ai) << "," << ToDeg(aj) << "  "<<ToDeg(a)<<","<<ToDeg(phi)<<std::endl;


            }
        }
}

wgFloat WavelengthIntegration(wgFloat Wavelength0, wgFloat Intensity0,wgFloat Wavelength1, wgFloat Intensity1)
{
    wgFloat m=(Intensity1-Intensity0)/(Wavelength1-Wavelength0),c=Intensity1-m*Wavelength1;
    return (m*Wavelength1*Wavelength1/2)-(m*Wavelength0*Wavelength0/2)+c*Wavelength1-c*Wavelength0;
}

#ifdef __linux__
void *ExecuteRayTrace(void *pv)
{
    RayTracerInfoS *p=(RayTracerInfoS *)pv;
#else
void ExecuteRayTrace(RayTracerInfoS *p)
{
#endif
    //wgFloat LastTotalInputIntensity[MaxDetectors];
    p->RandomIndex[0] = p->RandomIndex[1] = p->RandomIndex[2] = p->RandomIndex[3] = (p->ThreadID+1)*43+11;
    for (int SampleBlocks = p->ThreadID; SampleBlocks < (SceneWidth*SceneHeight); SampleBlocks += MaxUseableThreads)
    {
        p->Progress++;
        int ts = SampleBlocks;
        int j = (ts / SceneWidth);
        int i = (ts % SceneWidth);
        {
            wgFloat aj = ToRad(FOVY[0]+(((FOVY[1]-FOVY[0])*j) / (FOVYSamples-1)));
            if (aj < ToRad(InputAngleLimitY[0])) continue;
            if (aj > ToRad(InputAngleLimitY[1])) continue;
            wgFloat ai = ToRad(FOVX[0]+(((FOVX[1]-FOVX[0])*i) / (FOVXSamples-1)));
            if (ai < ToRad(InputAngleLimitX[0])) continue;
            if (ai > ToRad(InputAngleLimitX[1])) continue;

            int oy = (FOVYSamples-1 - j);
            int ox = i;
            if (Halt)
            {
                p->Finished = true;
#ifdef __linux
                return 0;
#else
                return;
#endif
            };
            for(int i=0; i<SpectrumWavelengths; i++)
            {
                p->WaveLength = SpectrumWavelengthAndReleativeIntensity[i*2];
                p->SamplesMC=std::uint64_t(SamplesMC*SpectrumWavelengthAndReleativeIntensity[i*2+1]);
                OneTrace(ai, aj, p,LineNull);
                if(SpectrumWavelengths>1)
                    for(int d=0; d<Detectors; d++)
                    {
                        Pixel[d][ox + oy*SceneWidth]+=p->TotalInputIntensity[d];
                        PixelS[d][ox + oy*SceneWidth]+=p->TotalS[d];
                        PixelP[d][ox + oy*SceneWidth]+=p->TotalP[d];
                        //if(i) Pixel[d][ox + oy*SceneWidth]+=WavelengthIntegration(WavelengthAndReleativeIntensity[(i-1)*2],LastTotalInputIntensity[d],WavelengthAndReleativeIntensity[i*2],p->TotalInputIntensity[d]);
                        //LastTotalInputIntensity[d]=p->TotalInputIntensity[d];
                    }
                else
                    for(int d=0; d<Detectors; d++)
                    {
                        Pixel[d][ox + oy*SceneWidth]=p->TotalInputIntensity[d];
                        PixelS[d][ox + oy*SceneWidth]=p->TotalS[d];
                        PixelP[d][ox + oy*SceneWidth]=p->TotalP[d];
                    }

            }
        }
    }
    p->Finished = true;
#ifdef __linux__
    return 0;
#endif
}


std::chrono::time_point<std::chrono::system_clock> T[2];
std::chrono::duration<double> Elapsed(void)
{
    std::chrono::duration<double> Seconds = std::chrono::system_clock::now() - T[0];
    return Seconds;
}

template<typename T>
void RXFind(std::string &s,T &a,std::string as)
{
    std::size_t found = s.find(as);
    if(found==std::string::npos)
    {
        std::cout<<as<<" not found or duplicated in argument file!\n";
        exit(1);
    }
    std::istringstream istr(s.substr(found+as.length()));
    char eq;
    istr>>std::skipws>>eq>>a;
    std::cout<<as<<"="<<a<<";\n";
}

template<typename T>
void RXFindA(std::string &s,T &a,std::string as,int n)
{
    std::size_t found = s.find(as);
    if(found==std::string::npos)
    {
        std::cout<<as<<" not found or duplicated in argument file!\n";
        exit(1);
    }
    std::istringstream istr(s.substr(found+as.length()));
    char w;
    istr>>std::skipws>>w>>w;
    std::cout<<as<<"=[";
    for(int i=0; i<n; i++)
    {
        istr>>std::skipws>>a[i]>>w;
        std::cout<<a[i]<<w;
    }
    std::cout<<";\n";
}

#define RXFind(A) RXFind(str,A,#A)
#define RXFindA(A,B) RXFindA(str,A,#A,B)
void ReOrderWavelengths(void)
{
    int f=1;
    while(f)
    {
        f=0;
        for(int i=1; i<SpectrumWavelengths; i++)
            if(SpectrumWavelengthAndReleativeIntensity[i*2+0]<SpectrumWavelengthAndReleativeIntensity[(i-1)*2+0])
            {
                f=1;
                std::swap(SpectrumWavelengthAndReleativeIntensity[i*2+0],SpectrumWavelengthAndReleativeIntensity[(i-1)*2+0]);
                std::swap(SpectrumWavelengthAndReleativeIntensity[i*2+1],SpectrumWavelengthAndReleativeIntensity[(i-1)*2+1]);
            }
    }
}

#define AssertSimilar(a,b,e) assert(std::abs( ((a)-(b)) < (e) ))


void Phox3TableExtract(bool Transmissive=false)
{
    OrderEfficiencyST *oe;
    if(Transmissive)
    {
        std::cout<<"Extracting Phox3 Table Transmission Efficiencies...\n";
        oe=OrderTEfficiency;
    }
    else
    {
        std::cout<<"Extracting Phox3 Table Efficiencies...\n";
        oe=OrderEfficiency;
    }


    for(int rt=0; rt<(RegionsInX*RegionsInY); rt++)
        for(int i=0; i<360; i++)
            for(int j=0; j<360; j++)
                oe->Get(rt,i,j).Orders=0;

    int MaxOrders=0,OrderBound=OrderLimit;
    for(int rt=0; rt<(RegionsInX*RegionsInY); rt++)
    {
        for(int n=0; n<2; n++)
            for(int i=-OrderBound; i<=OrderBound; i++)
            {

                for(int k=-OrderBound; k<=OrderBound; k++)
                {
                    FILE *fp[5]= {0};
                    for(int fn=0; fn<=4; fn++)
                    {
                        std::stringstream f;
                        f<<"OGWaveEfficiency_"<<WavelengthFolder<<"/";
                        if(fn<4) f<<"SP";
                        f<<RegionEfficiency[rt];
                        if(Transmissive)
                        {
                            if(n) f<<"/Tn180/";
                            else f<<"/T180/";
                            switch(fn)
                            {
                            case 0:
                                f<<"ss/T_table_ss_"<<i;
                                break;
                            case 1:
                                f<<"sp/T_table_sp_"<<i;
                                break;
                            case 2:
                                f<<"ps/T_table_ps_"<<i;
                                break;
                            case 3:
                                f<<"pp/T_table_pp_"<<i;
                                break;
                            case 4:
                                f<<"T_table_"<<i;
                                break;
                            }
                        }
                        else
                        {

                            if(n) f<<"/n180/";
                            else f<<"/180/";
                            switch(fn)
                            {
                            case 0:
                                f<<"ss/table_ss_"<<i;
                                break;
                            case 1:
                                f<<"sp/table_sp_"<<i;
                                break;
                            case 2:
                                f<<"ps/table_ps_"<<i;
                                break;
                            case 3:
                                f<<"pp/table_pp_"<<i;
                                break;
                            case 4:
                                f<<"table_"<<i;
                                break;
                            }
                        }
                        std::stringstream f3;

                        f3<<f.str()<<"_"<<k<<".csv";

                        std::cout<<rt<<":"<<f3.str()<<"\n";
                        fp[fn]=fopen(f3.str().c_str(),"r");
                        if(fn<4)
                        {
                            if(!fp[fn]) printf("Bad file %s\n",f3.str().c_str());
                            fflush(stdout);
                            assert(fp[fn]);
                        }
                    }
                    int maxfn;
                    if(!fp[4]) maxfn=3;
                    else maxfn=4;
                    for(int d=OutputSPAngles[0]; d<=OutputSPAngles[1]; d++)
                    {
                        for(int p=OutputSPAngles[2]; p<=OutputSPAngles[3]; p++)
                        {
                            wgFloat e[5];
                            for(int fn=0; fn<=maxfn; fn++)
                            {
                                char w;
                                int r=fscanf(fp[fn],"%f",&e[fn]);
                                if(p<180) r=fscanf(fp[fn],"%c",&w);
                                else w=' ';
                                //printf("%f%c",e[fn],w);
                                //fflush(stdout);

                            }
                            if(p>=180) continue;
                            if(n) if(!p) continue;
                            int p2=p+0;
                            if(n) p2=-p2;
                            else p2=p2;
                            p2=(p2%360);
                            if(p2<0) p2+=360;
                            assert(p2>=0);
                            assert(p2<=359);
                            OrderEfficiencyS &OE=oe->Get(rt,p2,d);
                            OE.Order[OE.Orders][0]=i;
                            OE.Order[OE.Orders][1]=k;
                            for(int fn=0; fn<=maxfn; fn++) OE.e[OE.Orders][fn]=e[fn];
                            if(maxfn==3)
                            {
                                e[4]=(OE.e[OE.Orders][0]+OE.e[OE.Orders][1]+OE.e[OE.Orders][2]+OE.e[OE.Orders][3])*0.5f;
                                OE.e[OE.Orders][4]=e[4];
                            }
                            AssertSimilar(OE.e[OE.Orders][4],(OE.e[OE.Orders][0]+OE.e[OE.Orders][1]+OE.e[OE.Orders][2]+OE.e[OE.Orders][3])*0.5f,0.01f);
                            if(std::abs(OE.e[OE.Orders][4]-(OE.e[OE.Orders][0]+OE.e[OE.Orders][1]+OE.e[OE.Orders][2]+OE.e[OE.Orders][3])*0.5f)>0.1f)
                                std::cout<<d<<","<<p<<":"<<OE.e[OE.Orders][4]<<"("<<(OE.e[OE.Orders][0]+OE.e[OE.Orders][1]+OE.e[OE.Orders][2]+OE.e[OE.Orders][3])*0.5f<<"),";
                            OE.Orders++;
                            assert(OE.Orders<MaxOrdersEfficiencies);
                            if(OE.Orders>MaxOrders) MaxOrders=OE.Orders;

                        }
                    }
                    for(int fn=0; fn<=4; fn++)
                    {
                        if(fp[fn]) fclose(fp[fn]);
                    }

                }
            }
    }
    std::cout<<"Max Orders:"<<MaxOrders<<"\n";
}


void Phox3TTableExtract(void)
{
    Phox3TableExtract(true);
}


void InputTableExtract(void)
{

    std::cout<<"\n\nExtracting Input Table Efficiencies...\n";

    for(int rt=0; rt<2; rt++)
        for(int a=0; a<180; a++)
            InputEfficiency[rt][a]=0;
    char w;
    std::cout<<"Angles:"<<InputGratingEfficiencyAngles[1]<<" To "<<
             InputGratingEfficiencyAngles[1]+(InputGratingEfficiencyAngles[0]-1)*InputGratingEfficiencyAngles[2]<<" Steps Of "<<InputGratingEfficiencyAngles[2]<<"\n";
    std::cout<<"IGWavesEfficiency/R0.csv\n";
    std::ifstream t("IGWavesEfficiency/R0.csv");
    if(t.good())
    {
        std::string str((std::istreambuf_iterator<char>(t)),std::istreambuf_iterator<char>());
        std::istringstream istr(str);
        std::string ln;
        for(int i=0; i<InputGratingEfficiencyAngles[0]; i++)
        {
            std::getline(istr, ln);
            std::istringstream istr(ln);
            int a;
            for(int j=0; j<InputGratingEfficiencyWavelengths; j++)
            {
                istr>>std::skipws>>InputEfficiency[0][j*InputGratingEfficiencyAngles[0]+i]>>w;
            }

        }
        t.close();
    }
    else
        for(int j=0; j<InputGratingEfficiencyWavelengths; j++)
        {
            for(int i=0; i<InputGratingEfficiencyAngles[0]; i++)
            {
                InputEfficiency[0][j*InputGratingEfficiencyAngles[0]+i]=((InputEfficiencyS[0][j*InputGratingEfficiencyAngles[0]+i]+InputEfficiencyP[0][j*InputGratingEfficiencyAngles[0]+i])/2);
            }
        }


    for(int j=0; j<InputGratingEfficiencyWavelengths; j++)
    {
        std::cout<<InputGratingEfficiencyWavelength[j]<<":\n";
        for(int i=0; i<InputGratingEfficiencyAngles[0]; i++)
        {
            std::cout<<InputEfficiency[0][j*InputGratingEfficiencyAngles[0]+i]<<"("<<(InputEfficiencyS[0][j*InputGratingEfficiencyAngles[0]+i]+InputEfficiencyP[0][j*InputGratingEfficiencyAngles[0]+i])/2<<"),";
            AssertSimilar(InputEfficiency[0][j*InputGratingEfficiencyAngles[0]+i],((InputEfficiencyS[0][j*InputGratingEfficiencyAngles[0]+i]+InputEfficiencyP[0][j*InputGratingEfficiencyAngles[0]+i])/2),0.01f);
            AssertSimilar(InputEfficiencyS[0][j*InputGratingEfficiencyAngles[0]+i]+InputEfficiencyP[0][j*InputGratingEfficiencyAngles[0]+i],2,0.01f);
        }
        std::cout<<"\n";
    }
    std::cout<<";\n";

    std::cout<<"IGWavesEfficiency/R1.csv\n";
    std::ifstream t1("IGWavesEfficiency/R1.csv");
    if(t1.good())
    {
        std::string str1((std::istreambuf_iterator<char>(t1)),std::istreambuf_iterator<char>());
        std::istringstream istr1(str1);
        std::string ln;

        for(int i=0; i<InputGratingEfficiencyAngles[0]; i++)
        {
            int a;
            std::getline(istr1, ln);
            std::istringstream istr(ln);
            for(int j=0; j<InputGratingEfficiencyWavelengths; j++)
            {
                istr>>std::skipws>>InputEfficiency[1][j*InputGratingEfficiencyAngles[0]+i]>>w;
            }

        }
        t1.close();
    }
    else
        for(int j=0; j<InputGratingEfficiencyWavelengths; j++)
        {
            for(int i=0; i<InputGratingEfficiencyAngles[0]; i++)
            {
                InputEfficiency[1][j*InputGratingEfficiencyAngles[0]+i]=((InputEfficiencyS[1][j*InputGratingEfficiencyAngles[0]+i]+InputEfficiencyP[1][j*InputGratingEfficiencyAngles[0]+i])/2);
            }
        }


    for(int j=0; j<InputGratingEfficiencyWavelengths; j++)
    {
        std::cout<<InputGratingEfficiencyWavelength[j]<<":\n";
        for(int i=0; i<InputGratingEfficiencyAngles[0]; i++)
        {
            std::cout<<InputEfficiency[1][j*InputGratingEfficiencyAngles[0]+i]<<"("<<(InputEfficiencyS[1][j*InputGratingEfficiencyAngles[0]+i]+InputEfficiencyP[1][j*InputGratingEfficiencyAngles[0]+i])/2<<"),";
            AssertSimilar(InputEfficiency[1][j*InputGratingEfficiencyAngles[0]+i],((InputEfficiencyS[1][j*InputGratingEfficiencyAngles[0]+i]+InputEfficiencyP[1][j*InputGratingEfficiencyAngles[0]+i])/2),0.01f);
            AssertSimilar(InputEfficiencyS[1][j*InputGratingEfficiencyAngles[0]+i]+InputEfficiencyP[1][j*InputGratingEfficiencyAngles[0]+i],2,0.01f);
        }
        std::cout<<"\n";
    }
    std::cout<<";\n";
}

void InputTableExtractSP(void)
{

    for(int rt=0; rt<2; rt++)
        for(int a=0; a<180; a++)
        {
            InputEfficiencyS[rt][a]=0;
            InputEfficiencyP[rt][a]=0;
        }

    const std::string SP[2]= {"S","P"};
    for(int SorP=0; SorP<2; SorP++)
    {
        wgFloat *p;
        std::cout<<"\n\nExtracting Input Table "<<SP[SorP]<<" Efficiencies...\n";

        char w;
        std::cout<<"Angles:"<<InputGratingEfficiencyAngles[1]<<" To "<<
                 InputGratingEfficiencyAngles[1]+(InputGratingEfficiencyAngles[0]-1)*InputGratingEfficiencyAngles[2]<<" Steps Of "<<InputGratingEfficiencyAngles[2]<<"\n";
        std::cout<<"IGWavesEfficiency/R0"<<SP[SorP]<<".csv\n";
        std::ifstream t("IGWavesEfficiency/R0"+SP[SorP]+".csv");
        assert(t.good());
        std::string str((std::istreambuf_iterator<char>(t)),std::istreambuf_iterator<char>());
        std::istringstream istr(str);
        std::string ln;
        if(SorP) p=InputEfficiencyP[0];
        else p=InputEfficiencyS[0];
        for(int i=0; i<InputGratingEfficiencyAngles[0]; i++)
        {
            std::getline(istr, ln);
            std::istringstream istr(ln);
            int a;
            for(int j=0; j<InputGratingEfficiencyWavelengths; j++)
            {
                istr>>std::skipws>>p[j*InputGratingEfficiencyAngles[0]+i]>>w;
            }

        }
        t.close();
        for(int j=0; j<InputGratingEfficiencyWavelengths; j++)
        {
            std::cout<<InputGratingEfficiencyWavelength[j]<<":\n";
            for(int i=0; i<InputGratingEfficiencyAngles[0]; i++)
            {
                std::cout<<p[j*InputGratingEfficiencyAngles[0]+i]<<",";
            }
            std::cout<<"\n";
        }
        std::cout<<";\n";

        std::cout<<"IGWavesEfficiency/R1"<<SP[SorP]<<".csv\n";
        std::ifstream t1("IGWavesEfficiency/R1"+SP[SorP]+".csv");
        assert(t1.good());
        std::string str1((std::istreambuf_iterator<char>(t1)),std::istreambuf_iterator<char>());
        std::istringstream istr1(str1);

        if(SorP) p=InputEfficiencyP[1];
        else p=InputEfficiencyS[1];
        for(int i=0; i<InputGratingEfficiencyAngles[0]; i++)
        {
            int a;
            std::getline(istr1, ln);
            std::istringstream istr(ln);
            for(int j=0; j<InputGratingEfficiencyWavelengths; j++)
            {
                istr>>std::skipws>>p[j*InputGratingEfficiencyAngles[0]+i]>>w;
            }

        }
        t1.close();
        for(int j=0; j<InputGratingEfficiencyWavelengths; j++)
        {
            std::cout<<InputGratingEfficiencyWavelength[j]<<":\n";
            for(int i=0; i<InputGratingEfficiencyAngles[0]; i++)
            {
                std::cout<<p[j*InputGratingEfficiencyAngles[0]+i]<<",";
            }
            std::cout<<"\n";
        }
        std::cout<<";\n";
    }
}

void Initialise(const char *ArgFile)
{
    std::ifstream t(ArgFile);
    std::string str((std::istreambuf_iterator<char>(t)),std::istreambuf_iterator<char>());
    //std::cout<<str<<"\n";

    RXFind(InputGratingEfficiencyWavelengths);
    RXFindA(InputGratingEfficiencyWavelength,InputGratingEfficiencyWavelengths);
    RXFindA(InputGratingEfficiencyAngles,3);
    RXFindA(OutputSPAngles,4);

    RXFind(RegionsInX);
    RXFind(RegionsInY);
    assert((RegionsInX*RegionsInY)<=MAXREGIONS);
    OrderEfficiency=new OrderEfficiencyST;
    OrderTEfficiency=new OrderEfficiencyST;

    RXFindA(RegionEfficiency,RegionsInX*RegionsInY);

    RXFind(WavelengthFolder);
    RXFind(InputGratingHalfWidth);
    RXFind(InputGratingHalfHeight);
    RXFind(WavelengthRefractiveIndex);
    RXFind(SpectrumWavelengths);
    assert(SpectrumWavelengths<=MaxWavelengths);
    RXFindA(SpectrumWavelengthAndReleativeIntensity,SpectrumWavelengths*2);
    ReOrderWavelengths();
    RXFind(OutputGratingPeriod);
    RXFind(WaveGuideDepth);
    RXFind(InputStandOff);
    RXFind(InputRadius);

    RXFind(Detectors);
    assert(Detectors<=MaxDetectors);
    RXFindA(DetectorXYZW,Detectors*4);
    RXFind(DetectorCircular);
    RXFind(DetectorReflectedZIsFromSurface);

    RXFind(AllRayDetectors);
    assert(AllRayDetectors<=MaxDetectors);
    RXFindA(AllRayDetectorXYW,Detectors*3);
    RXFindA(AllRayDetectorAbsorbing,Detectors);
    RXFind(AllRayDetectorCircular);
    RXFindA(AllRayDetectorWHPixels,2);
    RXFindA(AllRayDetectorSpatialWHPixels,2);
    RXFindA(AllRayDetectorSpatialWH,2);
    for(int i=0; i<AllRayDetectors; i++)
    {
        AllRayDetectorPixel[i]=new wgFloat [AllRayDetectorWHPixels[0]*2*AllRayDetectorWHPixels[1]*2];
        for (int o = 0; o < (AllRayDetectorWHPixels[0]*2*AllRayDetectorWHPixels[1]*2); o++)
        {
            AllRayDetectorPixel[i][o]=0;
        }
    }
    AllRayDetectorSpatialPixel=new wgFloat [AllRayDetectorSpatialWHPixels[0]*2*AllRayDetectorSpatialWHPixels[1]*2];
    for (int o = 0; o < (AllRayDetectorSpatialWHPixels[0]*2*AllRayDetectorSpatialWHPixels[1]*2); o++)
    {
        AllRayDetectorSpatialPixel[o]=0;
    }

    RXFindA(FOVX,2);
    RXFind(FOVXSamples);
    RXFindA(FOVY,2);
    RXFind(FOVYSamples);
    RXFind(SamplesMC);
    RXFind(BoundaryWidth);
    RXFind(BoundaryHeight);

    RXFind(BoundaryRoofEnd);
    RXFind(BoundaryRoofIndent);

    RXFind(SpatialBoundaryWidthPixels);
    RXFind(SpatialBoundaryHeightPixels);
    SpatialPixel[0]=new std::uint64_t [SpatialBoundaryWidthPixels*2*SpatialBoundaryHeightPixels*2];
    SpatialPixel[1]=new std::uint64_t [SpatialBoundaryWidthPixels*2*SpatialBoundaryHeightPixels*2];
    SpatialSP[0]=new wgFloat [SpatialBoundaryWidthPixels*2*SpatialBoundaryHeightPixels*2*2];
    SpatialSP[1]=new wgFloat [SpatialBoundaryWidthPixels*2*SpatialBoundaryHeightPixels*2*2];
    for (int o = 0; o < (SpatialBoundaryWidthPixels*2*SpatialBoundaryHeightPixels*2); o++)
    {
        SpatialPixel[0][o] = SpatialPixel[1][o] = 0;
    }
    for (int o = 0; o < (SpatialBoundaryWidthPixels*2*SpatialBoundaryHeightPixels*2*2); o++)
    {
        SpatialSP[0][o] = SpatialSP[1][o] = 0;
    }

    RXFind(InputGratingAngle);
    RXFind(OutputGratingAngle);
    RXFind(OutputGratingAngle2);
    RXFind(InputGratingOffsetX);
    RXFind(InputGratingRadiusSamples);
    RXFindA(InputAngleLimitX,2);
    RXFindA(InputAngleLimitY,2);
    RXFind(Ray3DOn);
    RXFind(Ray3DX);
    RXFind(Ray3DY);
    RXFind(OrderLimit);

    MaxUseableThreads=1;

#ifndef WG_SINGLETHREAD
#ifdef __linux__
    MaxUseableThreads = 64;
    std::cout<<"Linux:"<<MaxUseableThreads<<" Threads\n";
#else
    MaxUseableThreads = std::thread::hardware_concurrency();
    std::cout<<"Using Hardware Threads\n";

    if (MaxUseableThreads < 1) MaxUseableThreads = 1;
    std::cout << "Hardware Threads:" << MaxUseableThreads <<"("<<std::thread::hardware_concurrency()<<")"<<std::endl;
#endif
#else
    std::cout<<"No Threads\n";
#endif
    //exit(1);

    CriticalAngle = asin(1.0 / WavelengthRefractiveIndex);
    InputGratingAngle=ToRad(InputGratingAngle);
    OutputGratingAngle=ToRad(OutputGratingAngle);
    OutputGratingAngle2=ToRad(OutputGratingAngle2);

    InputTableExtractSP();
    InputTableExtract();
    //exit(1);
    Phox3TTableExtract();
    //exit(1);
    Phox3TableExtract();
    //exit(1);
}

bool Ray3D(const char *DetectorFile)
{
    if(!Ray3DOn) return false;
    InputGratingRadiusSamples=0;
    wgFloat r[2] = { Ray3DX, Ray3DY };

    struct RayTracerInfoS RayTracerInfo;

    InitialiseRayTrace();

    RayTracerInfo.ThreadID = 0;
    RayTracerInfo.Progress = 0;
    RayTracerInfo.Finished = false;

    RayTracerInfo.RandomIndex[0] = RayTracerInfo.RandomIndex[1] = RayTracerInfo.RandomIndex[2] = RayTracerInfo.RandomIndex[3] = (RayTracerInfo.ThreadID+1)*43+11;
    RayTracerInfo.WaveLength = SpectrumWavelengthAndReleativeIntensity[0];
    RayTracerInfo.SamplesMC=SamplesMC;
    OneTrace(ToRad(r[0]), ToRad(r[1]), &RayTracerInfo,Line3D);

    std::string outf(DetectorFile);
    std::ofstream outfile(outf+"RayPoint.3D",std::ofstream::binary);
    outfile.write((char *)&Point3D[0].x,Point3D.size()*3*sizeof(wgFloat));
    outfile.close();
    std::ofstream outfile2(outf+"RayPoints.3D",std::ofstream::binary);
    outfile2.write((char *)&Points3D[0],Points3D.size()*sizeof(int));
    outfile2.close();
    return true;
}

void WaveGuideRayTracer(const char *ArgFile,const char *DetectorFile)
{
    Initialise(ArgFile);

    TotalSamples=0;

    if (Ray3D(DetectorFile)) return;

    Halt = false;

    SceneWidth = FOVXSamples;
    SceneHeight = FOVYSamples;

    T[0] = std::chrono::system_clock::now();

    for(int i=0; i<Detectors; i++)
    {
        Pixel[i] = new std::uint64_t[SceneWidth*SceneHeight];
        PixelS[i] = new wgFloat[SceneWidth*SceneHeight];
        PixelP[i] = new wgFloat[SceneWidth*SceneHeight];
        for (int o = 0; o < (SceneWidth*SceneHeight); o++) Pixel[i][o] = PixelS[i][o] = PixelP[i][o] = 0;
    }

    struct RayTracerInfoS *RayTracerInfo = new RayTracerInfoS[MaxUseableThreads];
    InitialiseRayTrace();

#ifdef WG_SINGLETHREAD
    std::cout<<"No Threads\n";
    ExecuteRayTrace(&RayTracerInfo[0]);
#else
#ifdef __linux__
    pthread_t threads[MaxUseableThreads];

    for (int i = 0; i < MaxUseableThreads; i++)
    {
        RayTracerInfo[i].ThreadID = i;
        RayTracerInfo[i].Progress = 0;
        RayTracerInfo[i].Finished = false;
        int rc = pthread_create(&threads[i], 0, ExecuteRayTrace, &RayTracerInfo[i]);
        assert(!rc);
    }


    int f = 0;
    while (f != MaxUseableThreads)
    {

        std::uint64_t t = 0;
        f = 0;
        for (int i = 0; i < MaxUseableThreads; i++)
        {
            f += RayTracerInfo[i].Finished;
            t += RayTracerInfo[i].Progress;
        }
        float tf = (t * 100.0f) / (SceneWidth*SceneHeight);
        std::cout <<" "<<tf<<"%        \r";
        usleep(100000);
    }

    for (int i = 0; i < MaxUseableThreads; i++)
    {
        void *status;
        int rc = pthread_join(threads[i], &status);
        assert(!rc);
    }
#else

    std::vector<std::thread> Thread;
    for (int i = 0; i < MaxUseableThreads; i++)
    {
        RayTracerInfo[i].ThreadID = i;
        RayTracerInfo[i].Progress = 0;
        RayTracerInfo[i].Finished = false;
        Thread.push_back(std::thread(ExecuteRayTrace, &RayTracerInfo[i]));
    }


    int f = 0;
    while (f != MaxUseableThreads)
    {

        std::uint64_t t = 0;
        f = 0;
        for (int i = 0; i < MaxUseableThreads; i++)
        {
            f += RayTracerInfo[i].Finished;
            t += RayTracerInfo[i].Progress;
        }
        float tf = (t * 100.0f) / (SceneWidth*SceneHeight);
        std::cout <<" "<<tf<<"%        \r";
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }

    for (int i = 0; i < MaxUseableThreads; i++) Thread[i].join();
#endif
#endif

    FinaliseRayTrace();
    delete [] RayTracerInfo;

    std::chrono::duration<double> Secs = Elapsed();
    std::cout << " Took:" << Secs.count() << " Seconds" << std::endl;

    for(int d=0; d<Detectors; d++)
    {
        wgFloat mintp = std::numeric_limits<wgFloat>::max(), maxtp = -std::numeric_limits<wgFloat>::max();
        std::uint64_t tout=0;
        for (int i = 0; i < SceneWidth*SceneHeight; i++)
        {
            //if(!(i%SceneWidth)) std::cout << std::endl;
            wgFloat t = Pixel[d][i];
            tout+=t;
            //std::cout << t<<" ";
            if (t>maxtp) maxtp = t;
            if (t < mintp) mintp = t;
        }
        std::cout << std::endl;

        std::string outf(DetectorFile);
        std::ostringstream stm ;
        stm << d ;
        outf="Results/"+outf+"."+stm.str();
        std::ofstream outfile(outf,std::ofstream::binary);
        outfile.write((char *)Pixel[d],SceneWidth*SceneHeight*sizeof(std::uint64_t));
        outfile.close();

        outf="Results/"+std::string(DetectorFile)+"Counts."+stm.str();
        outfile.open(outf,std::ofstream::binary);
        outfile.write((char *)&TotalSamples,sizeof(std::uint64_t));
        outfile.write((char *)&tout,sizeof(std::uint64_t));
        outfile.write((char *)&mintp,sizeof(wgFloat));
        outfile.write((char *)&maxtp,sizeof(wgFloat));
        outfile.close();

        std::cout <<"Detector:"<<d+1<< " Total Samples In (Photons):" << TotalSamples << " Total Samples Out (Photons):" << tout << " Min:" << mintp << " Max:" << maxtp << std::endl;

    }

    std::string outf(DetectorFile);
    outf="Results/"+outf+"SpatialReflective.All";
    std::ofstream outfile(outf,std::ofstream::binary);
    outfile.write((char *)SpatialPixel[0],2*SpatialBoundaryWidthPixels*2*SpatialBoundaryHeightPixels*sizeof(std::uint64_t));
    outfile.close();
    outf="Results/"+std::string(DetectorFile)+"SpatialTransmissive.All";
    outfile.open(outf,std::ofstream::binary);
    outfile.write((char *)SpatialPixel[1],2*SpatialBoundaryWidthPixels*2*SpatialBoundaryHeightPixels*sizeof(std::uint64_t));
    outfile.close();
    outf="Results/"+std::string(DetectorFile)+"SpatialSReflective.All";
    outfile.open(outf,std::ofstream::binary);
    outfile.write((char *)&SpatialSP[0][0],2*SpatialBoundaryWidthPixels*2*SpatialBoundaryHeightPixels*sizeof(wgFloat));
    outfile.close();
    outf="Results/"+std::string(DetectorFile)+"SpatialPReflective.All";
    outfile.open(outf,std::ofstream::binary);
    outfile.write((char *)&SpatialSP[0][2*SpatialBoundaryWidthPixels*2*SpatialBoundaryHeightPixels],2*SpatialBoundaryWidthPixels*2*SpatialBoundaryHeightPixels*sizeof(wgFloat));
    outfile.close();
    outf="Results/"+std::string(DetectorFile)+"SpatialSTransmissive.All";
    outfile.open(outf,std::ofstream::binary);
    outfile.write((char *)&SpatialSP[1][0],2*SpatialBoundaryWidthPixels*2*SpatialBoundaryHeightPixels*sizeof(wgFloat));
    outfile.close();
    outf="Results/"+std::string(DetectorFile)+"SpatialPTransmissive.All";
    outfile.open(outf,std::ofstream::binary);
    outfile.write((char *)&SpatialSP[1][2*SpatialBoundaryWidthPixels*2*SpatialBoundaryHeightPixels],2*SpatialBoundaryWidthPixels*2*SpatialBoundaryHeightPixels*sizeof(wgFloat));
    outfile.close();

    for(int d=0; d<Detectors; d++)
    {
        std::ostringstream stm ;
        stm << d ;
        outf="Results/"+std::string(DetectorFile)+"S."+stm.str();
        outfile.open(outf,std::ofstream::binary);
        outfile.write((char *)&PixelS[d][0],SceneWidth*SceneHeight*sizeof(wgFloat));
        outfile.close();
        outf="Results/"+std::string(DetectorFile)+"P."+stm.str();
        outfile.open(outf,std::ofstream::binary);
        outfile.write((char *)&PixelP[d][0],SceneWidth*SceneHeight*sizeof(wgFloat));
        outfile.close();
    }

    for(int d=0; d<AllRayDetectors; d++)
    {
        std::ostringstream stm ;
        stm << d ;
        outf="Results/"+std::string(DetectorFile)+"AllRay."+stm.str();
        outfile.open(outf,std::ofstream::binary);
        outfile.write((char *)&AllRayDetectorPixel[d][0],AllRayDetectorWHPixels[0]*AllRayDetectorWHPixels[0]*2*2*sizeof(wgFloat));
        outfile.close();
    }

    std::ostringstream stm ;
    outf="Results/"+std::string(DetectorFile)+"AllRaySpatial.All";
    outfile.open(outf,std::ofstream::binary);
    outfile.write((char *)&AllRayDetectorSpatialPixel[0],AllRayDetectorSpatialWHPixels[0]*AllRayDetectorSpatialWHPixels[0]*2*2*sizeof(wgFloat));
    outfile.close();

}

int main(int argc, char* argv[])
{
    std::cout<<"Wave Guide Simulator Copyright Wave Optics 2016 Coding:G.Symons\n\n";
    if(argc!=3)
    {
        std::cout<<"Usage: <In arguments file name> <Out raw detector(s) intensity map(s)>\n";
        return 0;
    }
    WaveGuideRayTracer(argv[1],argv[2]);
    return 0;
}

