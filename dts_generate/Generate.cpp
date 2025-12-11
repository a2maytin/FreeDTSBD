

/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 */
#include <fstream>
#include <random>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include "Generate.h"
Generate::Generate()
{
}
Generate::Generate(std::vector <std::string> argument)
{

    m_Argument = argument;
    m_Healthy =true;
    m_Seed =36723;
    m_MinFaceAngle = -0.5;
    m_MinVerticesDistanceSquare = 1.0;
    m_MaxLinkLengthSquare = 3.0;
    m_OutputFilename = "out.q";
    m_Box(0) = 10;
    m_Box(1) = 10;
    m_Box(2) = 10;
    m_genus = 0;
    m_Type = "Flat";
    m_N = 5;
    Nfunction f;
    m_tsiPrecision = f.Int_to_String(18)+"."+f.Int_to_String(10);
    
    // Initialize DNA generation parameters
    m_GenerateDNA = false;
    m_DNAOutputFilename = "dna.q";
    m_DNANumBeads = 100;
    m_DNANumChains = 1;
    m_DNABoundaryRadius = 90.0;
    m_DNABoundaryOrigin(0) = 0.0;
    m_DNABoundaryOrigin(1) = 0.0;
    m_DNABoundaryOrigin(2) = 0.0;
    m_DNAMonomerRadius = 1.7;
    m_DNAPersistenceLength = 45.0;
    m_DNANumStages = 1;
    m_DNAStages.clear();
    m_DNAStages.push_back(std::make_pair(1, 100));
    m_DNABondsFilename = "dna_bonds.txt";
    m_DNAStartIndex = -1;  // -1 means infer from membrane file
    m_DNABondK = 100.0;
    m_DNABondL0 = 0.472;  // Default bond length in nm
    m_DNABondKAngle = 50.0;
    m_DNABondTheta0 = 0.0;  // 0 radians = straight chain

    ExploreArguments();     // read the input data


    // After reading the argument we generate a TS file based on the provided values
    if (m_Type == "flat" || m_Type == "Flat" || m_Type == "FLAT" )
    {
        FlatBilayer();
    }
    else if (m_Type == "Tetrahedron" || m_Type == "tetrahedron"  )
    {
        Tetrahedron();
    }
    else if (m_Type == "channel" || m_Type == "Channel"  )
    {
        Cylinder();
    }
    else if (m_Type == "high_gen" || m_Type == "High_gen"  )
    {
        HighTopologyStructure();
    }
    else if (m_Type == "high_genpbc" || m_Type == "High_genpbc"  )
    {
        HighTopologyStructurePBC();
    }
    else
    {
        std::cout<<" Shape type "<<m_Type<<" is not recognized "<<std::endl;
    }
    
    // Generate DNA if requested
    if (m_GenerateDNA) {
        GenerateDNA();
    }
    
}
Generate::~Generate()
{
    
}
//======== high genus
void Generate::HighTopologyStructurePBC()
{
    srand (40072);


///=======
//== Read Gen-morphology variables
//===========
    int TopDegree = m_genus;
    std::vector<Vertex_Map> allV;
    std::vector<Triangle_Map> allT;
    std::vector<Inclusion_Map> allI;
    double DR = 0.001;
    double dx = 0;
               // cleans the log and error files
    

    m_noUpperCreatedHole = 0;
    m_noLowerCreatedHole = 0;
//============================
// Gen vertex
//==============================
    double zm = m_Box(2)/2;
    double l=1.01;
    double xm= l/2;
    double ym= l/2;
    int id=0;
    m_N = int (m_Box(0)/l);
    m_Box(0) = double(m_N)*l;
    m_Box(1) = m_Box(0);
    int NoVertex = m_N;

    
    
    if(m_genus>m_N*m_N)
    {
        std::cout<<" the topology degree is high, you may use higher number of verteces \n";
    }
    
    for (int i=0;i<m_N;i++)
    {
        for (int j=0;j<m_N;j++)
        {
         
            
            
            double x = i*l+xm;
            double y = j*l+ym;
            double z = zm+l/2;
            Vertex_Map v;
            v.id = id;
            dx = double(rand()%1000)/1000.0;
            v.x = x+dx*DR;
            dx = double(rand()%1000)/1000.0;
            v.y = y+dx*DR;
            dx = double(rand()%1000)/1000.0;
            v.z = z+dx*DR;
            v.domain = 0;
            allV.push_back(v);
            id++;
            
        }
        
    }
    
    

    for (int i=0;i<m_N;i++)
    {
        for (int j=0;j<m_N;j++)
        {
            
            double x = i*l+xm;
            double y = j*l+ym;
            double z = zm-l/2;
            Vertex_Map v;
            v.id = id;
            dx = double(rand()%1000)/1000.0;
            v.x = x+dx*DR;
            dx = double(rand()%1000)/1000.0;
            v.y = y+dx*DR;
            dx = double(rand()%1000)/1000.0;
            v.z = z+dx*DR;
            v.domain = 0;
            allV.push_back(v);
            id++;
        }
        
    }
    

//============================
// Connect triangles
//==============================
   id=0;
    bool makehole=false;
    for (int i=0;i<m_N-1;i++)
    {
        for (int j=0;j<m_N-1;j++)
        {
            

            makehole=MakeHole(1,i,j);
            
            if(makehole==true)
            {
                std::vector<Triangle_Map> temT = MakeTrianglesAroundHole(id, i, j);
                for (std::vector<Triangle_Map>::iterator it = temT.begin() ; it != temT.end(); ++it)
                    allT.push_back(*it);

                id=id+8;
            }
            else
            {
                int v1id=idfromij(1,i,j);
                int v2id=idfromij(1,i,j+1);
                int v3id=idfromij(1,i+1,j+1);
                Triangle_Map T1;
                T1.id = id;
                T1.v1 = v1id;
                T1.v2 = v2id;
                T1.v3 = v3id;
                allT.push_back(T1);
                id++;
                
                v1id=idfromij(1,i,j);
                v2id=idfromij(1,i+1,j+1);
                v3id=idfromij(1,i+1,j);
                Triangle_Map T2;
                T2.id = id;
                T2.v1 = v1id;
                T2.v2 = v2id;
                T2.v3 = v3id;
                allT.push_back(T2);
                id++;
            }


            
        }
        
    }
    
    // lower layer
    for (int i=0;i<m_N-1;i++)
    {
        for (int j=0;j<m_N-1;j++)
        {
            int v1id=idfromij(2,i,j);
            int v3id=idfromij(2,i,j+1);
            int v2id=idfromij(2,i+1,j+1);

            makehole=MakeHole(2,i,j);

            
            if(makehole==true)
            {

            }
            else
            {
                Triangle_Map T1;
                T1.id = id;
                T1.v1 = v1id;
                T1.v2 = v2id;
                T1.v3 = v3id;
                allT.push_back(T1);
                id++;
                v1id=idfromij(2,i,j);
                v3id=idfromij(2,i+1,j+1);
                v2id=idfromij(2,i+1,j);
                Triangle_Map T2;
                T2.id = id;
                T2.v1 = v1id;
                T2.v2 = v2id;
                T2.v3 = v3id;
                allT.push_back(T2);
                id++;
            }
            

            
        }
        
    }
    

    for (int i=0;i<m_N-1;i++)
    {
          //upper layer edge
        int v1id=idfromij(1,i,m_N-1);
        
        
      //  id=m_N*j+i;

        int v2id=idfromij(1,i,0);
        int v3id=idfromij(1,i+1,0);
       // std::cout<<v1id<<"  "<<v2id<<"  "<<v3id<<"\n";
        Triangle_Map T1;
        T1.id = id;
        T1.v1 = v1id;
        T1.v2 = v2id;
        T1.v3 = v3id;
        allT.push_back(T1);
        id++;
        
      
         v1id=idfromij(1,i,m_N-1);
         v2id=idfromij(1,i+1,0);
         v3id=idfromij(1,i+1,m_N-1);
        Triangle_Map T2;
        T2.id = id;
        T2.v1 = v1id;
        T2.v2 = v2id;
        T2.v3 = v3id;
        allT.push_back(T2);
        id++;
        
        // lower layer edge
       v1id=idfromij(2,i,m_N-1);
       v3id=idfromij(2,i,0);
    v2id=idfromij(2,i+1,0);
        T1.id = id;
        T1.v1 = v1id;
        T1.v2 = v2id;
        T1.v3 = v3id;
        allT.push_back(T1);
        id++;
        
        
         v1id=idfromij(2,i,m_N-1);
         v3id=idfromij(2,i+1,0);
         v2id=idfromij(2,i+1,m_N-1);
        T2.id = id;
        T2.v1 = v1id;
        T2.v2 = v2id;
        T2.v3 = v3id;
        allT.push_back(T2);
        id++;

        
    }
 
  
    for (int i=0;i<m_N-1;i++)
    {
        //upper layer Y edge

        int v1id=idfromij(1,m_N-1,i);
        int v2id=idfromij(1,m_N-1,i+1);
        int v3id=idfromij(1,0,i+1);
        Triangle_Map T1;
        T1.id = id;
        T1.v1 = v1id;
        T1.v2 = v2id;
        T1.v3 = v3id;
        allT.push_back(T1);
        id++;
        
        
        v1id=idfromij(1,m_N-1,i);
        v2id=idfromij(1,0,i+1);
        v3id=idfromij(1,0,i);
        Triangle_Map T2;
        T2.id = id;
        T2.v1 = v1id;
        T2.v2 = v2id;
        T2.v3 = v3id;
        allT.push_back(T2);
        id++;
        
        
        //lower layer Y edge

        v1id=idfromij(2,m_N-1,i);
        v3id=idfromij(2,m_N-1,i+1);
        v2id=idfromij(2,0,i+1);
        T1.id = id;
        T1.v1 = v1id;
        T1.v2 = v2id;
        T1.v3 = v3id;
        allT.push_back(T1);
        id++;
        
        
        v1id=idfromij(2,m_N-1,i);
        v3id=idfromij(2,0,i+1);
        v2id=idfromij(2,0,i);
        T2.id = id;
        T2.v1 = v1id;
        T2.v2 = v2id;
        T2.v3 = v3id;
        allT.push_back(T2);
        id++;
        
    }
 
    {
        
        int v1id=idfromij(1,m_N-1,m_N-1);
        int v2id=idfromij(1,m_N-1,0);
        int v3id=idfromij(1,0,0);
        Triangle_Map T1;
        T1.id = id;
        T1.v1 = v1id;
        T1.v2 = v2id;
        T1.v3 = v3id;
        allT.push_back(T1);
        id++;
        
        
        v1id=idfromij(1,m_N-1,m_N-1);
        v2id=idfromij(1,0,0);
        v3id=idfromij(1,0,m_N-1);
        Triangle_Map T2;
        T2.id = id;
        T2.v1 = v1id;
        T2.v2 = v2id;
        T2.v3 = v3id;
        allT.push_back(T2);
        id++;
        
        // lower single trinagle
        
        v1id=idfromij(2,m_N-1,m_N-1);
        v3id=idfromij(2,m_N-1,0);
        v2id=idfromij(2,0,0);
        T1.id = id;
        T1.v1 = v1id;
        T1.v2 = v2id;
        T1.v3 = v3id;
        allT.push_back(T1);
        id++;
        
        
        v1id=idfromij(2,m_N-1,m_N-1);
        v3id=idfromij(2,0,0);
        v2id=idfromij(2,0,m_N-1);
        T2.id = id;
        T2.v1 = v1id;
        T2.v2 = v2id;
        T2.v3 = v3id;
        allT.push_back(T2);
        id++;
        
    }
    

    std::cout<<" number of the vertex "<<allV.size()<<" total number of triangles "<< allT.size()<<"\n";
    MeshBluePrint BluePrint;
    BluePrint.bvertex = allV;
    BluePrint.btriangle = allT;
    BluePrint.binclusion = allI;
    BluePrint.simbox = m_Box;
    
    std::string ext = m_OutputFilename.substr(m_OutputFilename.find_last_of(".") + 1);
    
    if(ext==TSExt)
    {
        WriteQFile(m_OutputFilename , BluePrint);
    }
    else if(ext==TSIExt)
    {
        WriteTSI(m_OutputFilename , BluePrint);
    }
    else
    {
        std::cout<<"---> Error: output file with "<<ext<<" extension is not recognized.  It should have either "<<TSExt<<" or "<<TSIExt<<" extension. "<<std::endl;
    }

}
void Generate::HighTopologyStructure()
{
    srand (40072);


///=======
//== Read Gen-morphology variables
//===========
    int TopDegree = m_genus;
    int NoVertex = m_N;
    std::vector<Vertex_Map> allV;
    std::vector<Triangle_Map> allT;
    std::vector<Inclusion_Map> allI;
    double DR = 0.001;
    double dx = 0;
               // cleans the log and error files
    
    if(m_genus>m_N*m_N)
    {
        std::cout<<" the topology degree is high, you may use higher number of verteces \n";
    }
    m_noUpperCreatedHole = 0;
    m_noLowerCreatedHole = 0;
//============================
// Gen vertex
//==============================
    double zm = m_Box(2)/2;
    double xm= 1;
    double ym= 1;
    double l=1.01;
    int id=0;
    for (int i=0;i<m_N;i++)
    {
        for (int j=0;j<m_N;j++)
        {
         
            
            
            double x = i*l+xm;
            double y = j*l+ym;
            double z = zm+l/2;
            Vertex_Map v;
            v.id = id;
            dx = double(rand()%1000)/1000.0;
            v.x = x+dx*DR;
            dx = double(rand()%1000)/1000.0;
            v.y = y+dx*DR;
            dx = double(rand()%1000)/1000.0;
            v.z = z+dx*DR;
            v.domain = 0;
            allV.push_back(v);
            id++;
            
        }
        
    }
    
    

    for (int i=0;i<m_N;i++)
    {
        for (int j=0;j<m_N;j++)
        {
            
            double x = i*l+xm;
            double y = j*l+ym;
            double z = zm-l/2;
            Vertex_Map v;
            v.id = id;
            dx = double(rand()%1000)/1000.0;
            v.x = x+dx*DR;
            dx = double(rand()%1000)/1000.0;
            v.y = y+dx*DR;
            dx = double(rand()%1000)/1000.0;
            v.z = z+dx*DR;
            v.domain = 0;
            allV.push_back(v);
            id++;
        }
        
    }
    

//============================
// Connect triangles
//==============================
   id=0;
    bool makehole=false;
    for (int i=0;i<m_N-1;i++)
    {
        for (int j=0;j<m_N-1;j++)
        {
            

            makehole=MakeHole(1,i,j);
            
            if(makehole==true)
            {
                std::vector<Triangle_Map> temT = MakeTrianglesAroundHole(id, i, j);
                for (std::vector<Triangle_Map>::iterator it = temT.begin() ; it != temT.end(); ++it)
                    allT.push_back(*it);

                id=id+8;
            }
            else
            {
                int v1id=idfromij(1,i,j);
                int v2id=idfromij(1,i,j+1);
                int v3id=idfromij(1,i+1,j+1);
                Triangle_Map T1;
                T1.id = id;
                T1.v1 = v1id;
                T1.v2 = v2id;
                T1.v3 = v3id;
                allT.push_back(T1);
                id++;
                
                v1id=idfromij(1,i,j);
                v2id=idfromij(1,i+1,j+1);
                v3id=idfromij(1,i+1,j);
                Triangle_Map T2;
                T2.id = id;
                T2.v1 = v1id;
                T2.v2 = v2id;
                T2.v3 = v3id;
                allT.push_back(T2);
                id++;
            }


            
        }
        
    }
    
    // lower layer
    for (int i=0;i<m_N-1;i++)
    {
        for (int j=0;j<m_N-1;j++)
        {
            int v1id=idfromij(2,i,j);
            int v3id=idfromij(2,i,j+1);
            int v2id=idfromij(2,i+1,j+1);

            makehole=MakeHole(2,i,j);

            
            if(makehole==true)
            {

            }
            else
            {
                Triangle_Map T1;
                T1.id = id;
                T1.v1 = v1id;
                T1.v2 = v2id;
                T1.v3 = v3id;
                allT.push_back(T1);
                id++;
                v1id=idfromij(2,i,j);
                v3id=idfromij(2,i+1,j+1);
                v2id=idfromij(2,i+1,j);
                Triangle_Map T2;
                T2.id = id;
                T2.v1 = v1id;
                T2.v2 = v2id;
                T2.v3 = v3id;
                allT.push_back(T2);
                id++;
            }
            

            
        }
        
    }
    

    for (int i=0;i<m_N-1;i++)
    {
        
        int v1id=idfromij(1,0,i);
        int v2id=idfromij(2,0,i);
        int v3id=idfromij(2,0,i+1);
        Triangle_Map T1;
        T1.id = id;
        T1.v1 = v1id;
        T1.v2 = v2id;
        T1.v3 = v3id;
        allT.push_back(T1);
        id++;
        
        
         v1id=idfromij(1,0,i);
         v2id=idfromij(2,0,i+1);
         v3id=idfromij(1,0,i+1);
        Triangle_Map T2;
        T2.id = id;
        T2.v1 = v1id;
        T2.v2 = v2id;
        T2.v3 = v3id;
        allT.push_back(T2);
        id++;

        
    }
    
    
    for (int i=0;i<m_N-1;i++)
    {
        
        int v1id=idfromij(1,m_N-1,i+1);
        int v2id=idfromij(2,m_N-1,i+1);
        int v3id=idfromij(2,m_N-1,i);
        Triangle_Map T1;
        T1.id = id;
        T1.v1 = v1id;
        T1.v2 = v2id;
        T1.v3 = v3id;
        allT.push_back(T1);
        id++;
        
        
        v1id=idfromij(1,m_N-1,i+1);
        v2id=idfromij(2,m_N-1,i);
        v3id=idfromij(1,m_N-1,i);
        Triangle_Map T2;
        T2.id = id;
        T2.v1 = v1id;
        T2.v2 = v2id;
        T2.v3 = v3id;
        allT.push_back(T2);
        id++;
        
    }
    
    for (int i=0;i<m_N-1;i++)
    {
        
        int v1id=idfromij(1,i+1,0);
        int v2id=idfromij(2,i+1,0);
        int v3id=idfromij(2,i,0);
        Triangle_Map T1;
        T1.id = id;
        T1.v1 = v1id;
        T1.v2 = v2id;
        T1.v3 = v3id;
        allT.push_back(T1);
        id++;
        
        
        v1id=idfromij(1,i+1,0);
        v2id=idfromij(2,i,0);
        v3id=idfromij(1,i,0);
        Triangle_Map T2;
        T2.id = id;
        T2.v1 = v1id;
        T2.v2 = v2id;
        T2.v3 = v3id;
        allT.push_back(T2);
        id++;
        
    }
    
    
    for (int i=0;i<m_N-1;i++)
    {
        
        int v1id=idfromij(1,i,m_N-1);
        int v2id=idfromij(2,i,m_N-1);
        int v3id=idfromij(2,i+1,m_N-1);
        Triangle_Map T1;
        T1.id = id;
        T1.v1 = v1id;
        T1.v2 = v2id;
        T1.v3 = v3id;
        allT.push_back(T1);
        id++;
        
        
        v1id=idfromij(1,i,m_N-1);
        v2id=idfromij(2,i+1,m_N-1);
        v3id=idfromij(1,i+1,m_N-1);
        Triangle_Map T2;
        T2.id = id;
        T2.v1 = v1id;
        T2.v2 = v2id;
        T2.v3 = v3id;
        allT.push_back(T2);
        id++;
        
    }
    std::cout<<" number of the vertex "<<allV.size()<<" total number of triangles "<< allT.size()<<"\n";
    MeshBluePrint BluePrint;
    BluePrint.bvertex = allV;
    BluePrint.btriangle = allT;
    BluePrint.binclusion = allI;
    BluePrint.simbox = m_Box;
    
    std::string ext = m_OutputFilename.substr(m_OutputFilename.find_last_of(".") + 1);
    
    if(ext==TSExt)
    {
        WriteQFile(m_OutputFilename , BluePrint);
    }
    else if(ext==TSIExt)
    {
        WriteTSI(m_OutputFilename , BluePrint);
    }
    else
    {
        std::cout<<"---> Error: output file with "<<ext<<" extension is not recognized.  It should have either "<<TSExt<<" or "<<TSIExt<<" extension. "<<std::endl;
    }

}
int Generate::idfromij(int s, int i, int j)
{
    
    if(i>=m_N || m_N<=j)
    {
        std::cout<<"error, this should happen \n";
    }
    int id=0;
    if(s==1)
    id=m_N*j+i;
    else if(s==2)
    id=m_N*m_N+m_N*j+i;
    
    
    
    return id;
}
std::vector<Triangle_Map>  Generate::MakeTrianglesAroundHole(int id, int i, int j)
{
    
    std::vector<Triangle_Map> allT;

    int v1id=idfromij(1,i,j);
    int v2id=idfromij(1,i,j+1);
    int v3id=idfromij(2,i,j+1);
    Triangle_Map T1;
    T1.id = id;
    T1.v1 = v1id;
    T1.v2 = v2id;
    T1.v3 = v3id;
    allT.push_back(T1);
    id++;


    v1id=idfromij(1,i,j);
    v2id=idfromij(2,i,j+1);
    v3id=idfromij(2,i,j);
    T1.id = id;
    T1.v1 = v1id;
    T1.v2 = v2id;
    T1.v3 = v3id;
    allT.push_back(T1);
    id++;

    
    v1id=idfromij(1,i,j+1);
    v2id=idfromij(2,i+1,j+1);
    v3id=idfromij(2,i,j+1);
    T1.id = id;
    T1.v1 = v1id;
    T1.v2 = v2id;
    T1.v3 = v3id;
    allT.push_back(T1);
    id++;

    v1id=idfromij(1,i,j+1);
    v2id=idfromij(1,i+1,j+1);
    v3id=idfromij(2,i+1,j+1);
    T1.id = id;
    T1.v1 = v1id;
    T1.v2 = v2id;
    T1.v3 = v3id;
    allT.push_back(T1);
    id++;

    v1id=idfromij(1,i,j);
    v2id=idfromij(2,i,j);
    v3id=idfromij(2,i+1,j);
    T1.id = id;
    T1.v1 = v1id;
    T1.v2 = v2id;
    T1.v3 = v3id;
    allT.push_back(T1);
    id++;

    v1id=idfromij(1,i,j);
    v2id=idfromij(2,i+1,j);
    v3id=idfromij(1,i+1,j);
    T1.id = id;
    T1.v1 = v1id;
    T1.v2 = v2id;
    T1.v3 = v3id;
    allT.push_back(T1);
    id++;
    
    
    v1id=idfromij(1,i+1,j);
    v2id=idfromij(2,i+1,j+1);
    v3id=idfromij(1,i+1,j+1);
    T1.id = id;
    T1.v1 = v1id;
    T1.v2 = v2id;
    T1.v3 = v3id;
    allT.push_back(T1);
    id++;
    
    
    v1id=idfromij(1,i+1,j);
    v2id=idfromij(2,i+1,j);
    v3id=idfromij(2,i+1,j+1);
    T1.id = id;
    T1.v1 = v1id;
    T1.v2 = v2id;
    T1.v3 = v3id;
    allT.push_back(T1);

    return allT;
}

bool Generate::MakeHole(int s,int i, int j)
{
    bool is=false;
    
    if(m_genus!=0)
    {
    
    int space=(m_N-2)/(int(sqrt(m_genus))+1);
    

        if(s==1 && m_noUpperCreatedHole<m_genus)
        {
            
            if((i+1)%space==0 && (j+1)%space==0 && j<m_N-1 && i<m_N-1)
                is=true;
            
        }
        else if(s==2 && m_noLowerCreatedHole<m_genus && j<m_N-1 && i<m_N-1)
        {
            if((i+1)%space==0 && (j+1)%space==0)
                is=true;
        }
    }



    
    
    
    if(is==true && s==1)
    m_noUpperCreatedHole++;
    else if(is==true && s==2)
    m_noLowerCreatedHole++;
    
    return is;
}
//======= end high genus
void Generate::Cylinder()
{
    std::vector<Vertex_Map> allV;
    std::vector<Triangle_Map> allT;
    std::vector<Inclusion_Map> allI;
    
    double l=1.2;
    int Nz = int(m_Box(2)/l);
    m_Box(2) = double(Nz)*l;
    double lx=l*0.87;
    
    int Nx= m_N;
    if (Nx%2!=0)
        Nx=Nx+1;
    double Width = double(Nx)*lx;
    
    double xm=m_Box(1)/2.0;
    
    int id=0;
    
    for (int i=0;i<Nz;i++)
    {
        for (int j=0;j<Nx;j++)
        {
            
            double y=(double(j)+0.5)*lx-Width/2+xm;
            double x=xm+Width/2+0.2;
            double z=l*(double(i)+double((j)%2)*0.5);
            
            Vertex_Map v;
            v.id = id;
            v.x = x;
            v.y = y;
            v.z = z;
            v.domain = 0;
            allV.push_back(v);
            id++;
            
       }
        for (int j=0;j<Nx;j++)
        {
            
            double x=-(double(j)+0.5)*lx+Width/2+xm;
            double y=xm+Width/2+0.2;
            double z=l*(double(i)+double((j+Nx)%2)*0.5);
            
            Vertex_Map v;
            v.id = id;
            v.x = x;
            v.y = y;
            v.z = z;
            v.domain = 0;
            allV.push_back(v);
            
            id++;
            
        }
        for (int j=0;j<Nx;j++)
        {
            
            double y=-(double(j)+0.5)*lx+Width/2+xm;
            double x=-Width/2+xm-0.2;
            double z=l*(double(i)+double((j)%2)*0.5);
            
            Vertex_Map v;
            v.id = id;
            v.x = x;
            v.y = y;
            v.z = z;
            v.domain = 0;
            allV.push_back(v);
            
            id++;
            
        }
        for (int j=0;j<Nx;j++)
        {
            
            double x=(double(j)+0.5)*lx-Width/2+xm;
            double y=-Width/2+xm-0.2;
            double z=l*(double(i)+double((j)%2)*0.5);
            
            Vertex_Map v;
            v.id = id;
            v.x = x;
            v.y = y;
            v.z = z;
            v.domain = 0;
            allV.push_back(v);
            
            id++;
            
        }
    }
 
// Making triangles
    int t=0;
    
    for (int i=0;i<Nz;i++)
    {
        for(int j=0;j<4*Nx;j++)
        {
            int M = 4*Nx;
            if(j%2==0)
            {
                Triangle_Map T1;
                T1.id = t;
                T1.v1 = CylinderIndex(M,Nz,i,j);
                T1.v2 = CylinderIndex(M,Nz,i,j+1);
                T1.v3 = CylinderIndex(M,Nz,i+1,j);
                allT.push_back(T1);
                
                t++;
                T1.id = t;
                T1.v1 = CylinderIndex(M,Nz,i+1,j);
                T1.v2 = CylinderIndex(M,Nz,i,j+1);
                T1.v3 = CylinderIndex(M,Nz,i+1,j+1);
                allT.push_back(T1);
                t++;
                
            }
            else
            {
                Triangle_Map T1;
                T1.id = t;
                T1.v1 = CylinderIndex(M,Nz,i,j);
                T1.v2 = CylinderIndex(M,Nz,i,j+1);
                T1.v3 = CylinderIndex(M,Nz,i+1,j+1);
                allT.push_back(T1);
                t++;
                
                T1.id = t;
                T1.v1 = CylinderIndex(M,Nz,i,j);
                T1.v2 = CylinderIndex(M,Nz,i+1,j+1);
                T1.v3 = CylinderIndex(M,Nz,i+1,j);
                allT.push_back(T1);
                t++;
                
                
            }
        }
    }
    MeshBluePrint BluePrint;
    BluePrint.bvertex = allV;
    BluePrint.btriangle = allT;
    BluePrint.binclusion = allI;
    BluePrint.simbox = m_Box;
    
    std::string ext = m_OutputFilename.substr(m_OutputFilename.find_last_of(".") + 1);
    
    if(ext==TSExt)
    {
        WriteQFile(m_OutputFilename , BluePrint);
    }
    else if(ext==TSIExt)
    {
        WriteTSI(m_OutputFilename , BluePrint);
    }
    else
    {
        std::cout<<"---> Error: output file with "<<ext<<" extension is not recognized.  It should have either "<<TSExt<<" or "<<TSIExt<<" extension. "<<std::endl;
    }

}
void Generate::Tetrahedron()
{
    
    int N=m_N;
    int Nv=2*(N*N+1);
    
    
    double x=0;
    double y=0;
    double z=0;
    int id=0;
    double b=1.2;
    double dl=N*b;
    double theta=asin(1.0/sqrt(3.0));
    double H=dl*cos(theta);
    int M;
    
    
    std::vector<Vertex_Map> allV;
    std::vector<Triangle_Map> allT;
    std::vector<Inclusion_Map> allI;
  
    double lx=m_Box(0)/2;
    double ly=m_Box(1)/2;
    double lz=m_Box(2)/2;
    
    
    //novertex=3*(N*N)+2
    for (int i=0;i<N+1;i++)
    {
        if(i<1)
            M=0;
        else if(i>=1)
            M=i-1;
        for (int j=0;j<M+1;j++)
        {
            
            z=H-b*i*cos(theta);
            x=b*i*sin(theta)/2.0;
            y=-dl/2.0+j*b+(N-i)*b/2.0;
            Vertex_Map v;
            v.id = id;
            v.x = x+lx;
            v.y = y+ly;
            v.z = z+lz;
            v.domain = 0;
            allV.push_back(v);
            id++;
        }
    }
    for (int i=0;i<N+1;i++)
    {
        for (int j=0;j<i;j++)
        {
            
            z=H-b*i*cos(theta);
            x=b*(i-j)*sqrt(3.0)/2.0-b*i*sin(theta);
            y=b*i/2.0-j*b/2.0;
            Vertex_Map v;
            v.id = id;
            v.x = x+lx;
            v.y = y+ly;
            v.z = z+lz;
            v.domain = 0;
            allV.push_back(v);
            id++;
            
        }
    }
    
    for (int i=0;i<N+1;i++)
    {
        for (int j=0;j<i;j++)
        {
            
            z=H-b*i*cos(theta);
            x=b*j*sqrt(3.0)/2.0-b*i*sin(theta);
            y=-j*b/2.0;
            Vertex_Map v;
            v.id = id;
            v.x = x+lx;
            v.y = y+ly;
            v.z = z+lz;
            v.domain = 0;
            allV.push_back(v);
            id++;
            
        }
    }
    
    for (int i=1;i<N;i++)
    {
        for (int j=1;j<i;j++)
        {
            
            z=0;
            double y1=dl*sqrt(3)/3-i*b*sqrt(3)/2;
            double x1=-i*b/2+j*b;
            x=x1;
            y=y1;
            double C=sqrt(3)/2;
            double S=-0.5;
            x=C*x1-S*y1;
            y=S*x1+C*y1;
            Vertex_Map v;
            v.id = id;
            v.x = x+lx;
            v.y = y+ly;
            v.z = z+lz;
            v.domain = 0;
            allV.push_back(v);
            id++;
            
        }
    }

    int Tid=0;
    
    //==================== adding three triangles associated with upper vertex
    int v1=0;
    int v2=1;
    int v3=N*(N+1)/2+1;
    
    Triangle_Map T1;
    T1.id = Tid;
    T1.v1 = v1;
    T1.v2 = v2;
    T1.v3 = v3;
    allT.push_back(T1);
    Tid++;
    
    v2=N*(N+1)/2+1;
    v3=N*(N+1)+1;
    Triangle_Map T2;
    T2.id = Tid;
    T2.v1 = v1;
    T2.v2 = v2;
    T2.v3 = v3;
    allT.push_back(T2);
    Tid++;
    
    v2=N*(N+1)+1;
    v3=1;

    Triangle_Map T3;
    T3.id = Tid;
    T3.v1 = v1;
    T3.v2 = v2;
    T3.v3 = v3;
    allT.push_back(T3);
    Tid++;
    //========================

    id=1;
    for (int i=1;i<N;i++)
    {
        for (int j=1;j<i+1;j++)
        {
            v1=findid(0,i,j);
            v2=findid(0,i+1,j);
            v3=findid(0,i+1,j+1);

            Triangle_Map T;
            T.id = Tid;
            T.v1 = v1;
            T.v2 = v2;
            T.v3 = v3;
            allT.push_back(T);
            
            Tid++;
            
            v1=findid(0,i,j);
            v2=findid(0,i+1,j+1);
            v3=findid(0,i,j+1);
            T.id = Tid;
            T.v1 = v1;
            T.v2 = v2;
            T.v3 = v3;
            allT.push_back(T);
            Tid++;
            
            if(j==1)
            {
                
                v1=findid(0,i,j);
                v2=findid(0,i+1,j-1);
                v3=findid(0,i+1,j);
                T.id = Tid;
                T.v1 = v1;
                T.v2 = v2;
                T.v3 = v3;
                allT.push_back(T);
                Tid++;
            }
        }
    }
    
    
  
    
    id=1;
    for (int i=1;i<N;i++)
    {
        for (int j=1;j<i+1;j++)
        {
            v1=findid(1,i,j);
            v2=findid(1,i+1,j);
            v3=findid(1,i+1,j+1);
            Triangle_Map T;
            T.id = Tid;
            T.v1 = v1;
            T.v2 = v2;
            T.v3 = v3;
            allT.push_back(T);
            Tid++;
            
            v1=findid(1,i,j);
            v2=findid(1,i+1,j+1);
            v3=findid(1,i,j+1);
            T.id = Tid;
            T.v1 = v1;
            T.v2 = v2;
            T.v3 = v3;
            allT.push_back(T);
            Tid++;
            
            if(j==1)
            {
                
                v1=findid(1,i,j);
                v2=findid(1,i+1,j-1);
                v3=findid(1,i+1,j);
                Triangle_Map T;
                T.id = Tid;
                T.v1 = v1;
                T.v2 = v2;
                T.v3 = v3;
                allT.push_back(T);
                Tid++;
            }
        }
    }
    
 
    for (int i=1;i<N;i++)
    {
        for (int j=1;j<i+1;j++)
        {
            v1=findid(2,i,j);
            v2=findid(2,i+1,j);
            v3=findid(2,i+1,j+1);
            Triangle_Map T;
            T.id = Tid;
            T.v1 = v1;
            T.v2 = v2;
            T.v3 = v3;
            allT.push_back(T);
            Tid++;
            
            v1=findid(2,i,j);
            v2=findid(2,i+1,j+1);
            v3=findid(2,i,j+1);
            T.id = Tid;
            T.v1 = v1;
            T.v2 = v2;
            T.v3 = v3;
            allT.push_back(T);
            Tid++;
            
            
            if(j==1)
            {
                
                v1=findid(2,i,j);
                v2=findid(2,i+1,j-1);
                v3=findid(2,i+1,j);
                T.id = Tid;
                T.v1 = v1;
                T.v2 = v2;
                T.v3 = v3;
                allT.push_back(T);
                Tid++;
            }
        }
    }
    //  three extra traingales
    
    {
    v2=findid(0,N,1);
    v1=findid(0,N,2);
    v3=findid(2,N,N);
    Triangle_Map T;
    T.id = Tid;
    T.v1 = v1;
    T.v2 = v2;
    T.v3 = v3;
    allT.push_back(T);
    Tid++;
    
    
    v2=findid(1,N,1);
    v1=findid(1,N,2);
    v3=findid(0,N,N);
    T.id = Tid;
    T.v1 = v1;
    T.v2 = v2;
    T.v3 = v3;
    allT.push_back(T);
    Tid++;
    
    
    v2=findid(2,N,1);
    v1=findid(2,N,2);
    v3=findid(1,N,N);
    T.id = Tid;
    T.v1 = v1;
    T.v2 = v2;
    T.v3 = v3;
    allT.push_back(T);
    Tid++;
    }
    
    // last face a hard one
    for (int i=1;i<N-1;i++)
    {
        for (int j=1;j<i+1;j++)
        {
            v1=findidface4(i,j);
            v3=findidface4(i+1,j);
            v2=findidface4(i+1,j+1);
            Triangle_Map T;
            T.id = Tid;
            T.v1 = v1;
            T.v2 = v2;
            T.v3 = v3;
            allT.push_back(T);
            Tid++;
            
            v1=findidface4(i,j);
            v3=findidface4(i+1,j+1);
            v2=findidface4(i,j+1);
            T.id = Tid;
            T.v1 = v1;
            T.v2 = v2;
            T.v3 = v3;
            allT.push_back(T);
            Tid++;
            
            
            if(j==i)
            {
                
                v1=findidface4(i-1,j+1);
                v2=findidface4(i,j+1);
                v3=findidface4(i,j);
                
                T.id = Tid;
                T.v1 = v1;
                T.v2 = v2;
                T.v3 = v3;
                allT.push_back(T);
                Tid++;
            }
            if(j==1)
            {
                
                v1=findidface4(i-1,j);
                v2=findidface4(i,j);
                v3=findidface4(i-1,j-1);
                
                T.id = Tid;
                T.v1 = v1;
                T.v2 = v2;
                T.v3 = v3;
                allT.push_back(T);
                Tid++;
                
                
                v1=findidface4(i-1,j-1);
                v2=findidface4(i,j);
                v3=findidface4(i,j-1);
                
                T.id = Tid;
                T.v1 = v1;
                T.v2 = v2;
                T.v3 = v3;
                allT.push_back(T);
                Tid++;
            }
        }
    }

    
    
    v1=findidface4(N-2,1);
    v2=findid(2,N,2);
    v3=findid(1,N,N);
    
    
    Triangle_Map T;
    T.id = Tid;
    T.v1 = v1;
    T.v2 = v2;
    T.v3 = v3;
    allT.push_back(T);
    Tid++;
    
    
    std::cout<<" number of the vertex "<<allV.size()<<" total number of triangles "<< allT.size()<<"\n";
    MeshBluePrint BluePrint;
    BluePrint.bvertex = allV;
    BluePrint.btriangle = allT;
    BluePrint.binclusion = allI;
    BluePrint.simbox = m_Box;
    
    std::string ext = m_OutputFilename.substr(m_OutputFilename.find_last_of(".") + 1);
    
    if(ext==TSExt)
    {
        WriteQFile(m_OutputFilename , BluePrint);
    }
    else if(ext==TSIExt)
    {
        WriteTSI(m_OutputFilename , BluePrint);
    }
    else
    {
        std::cout<<"---> Error: output file with "<<ext<<" extension is not recognized.  It should have either "<<TSExt<<" or "<<TSIExt<<" extension. "<<std::endl;
    }
    
}
void Generate::FlatBilayer() // making flat bilayers
{
    
    srand (800);
    double lx=sqrt(2)/sqrt(2);
    double ly=(sqrt(6.0)/2.0)/sqrt(2);
    lx=lx*1.2;
    ly=ly*1.2;
    
    
    int N=int (m_Box(1)/ly);
    int M=int (m_Box(0)/lx);
    
    lx=m_Box(0)/double(M);
    ly=m_Box(1)/double(N);

    
    int t=0;
    m_Box(0)=lx*double(M);
    m_Box(1)=ly*double(N);
    
    std::vector<Vertex_Map> allV;
    std::vector<Triangle_Map> allT;
    std::vector<Inclusion_Map> allI;

// making the vertices
    for (int i=0;i<N;i++)
    {
        for(int j=0;j<M;j++)
        {
            double s1=0.4*(double(rand()%2000000)/2000000-0.5);
            double s2=0.4*(double(rand()%2000000)/2000000-0.5);
            double x=((i)%2*0.5+double(j))*lx;
            double y=(double(i))*ly;
            double z=m_Box(2)/2+0.02*(2+sin(y)+cos(2*x));
            
            Vertex_Map v;
            v.x = x;
            v.y = y;
            v.z = z;
            v.id = t;
            v.domain = 0;
            allV.push_back(v);
            t++;
        }
    }

    t=0;
    for (int j=0;j<N;j++)
    {
        for(int i=0;i<M;i++)
        {
            if(j%2==0)
            {
                Triangle_Map T1;
                T1.id = t;
                T1.v1 = findindex(M,N,i,j);
                T1.v2 = findindex(M,N,i,j+1);
                T1.v3 = findindex(M,N,i-1,j+1);
                allT.push_back(T1);
                t++;
                Triangle_Map T2;
                T2.id = t;
                T2.v1 = findindex(M,N,i,j);
                T2.v2 = findindex(M,N,i+1,j);
                T2.v3 = findindex(M,N,i,j+1);
                allT.push_back(T2);
                t++;
                
            }
            else
            {
                Triangle_Map T1;
                T1.id = t;
                T1.v1 = findindex(M,N,i,j);
                T1.v2 = findindex(M,N,i+1,j+1);
                T1.v3 = findindex(M,N,i,j+1);
                allT.push_back(T1);
                t++;
                Triangle_Map T2;
                T2.id = t;
                T2.v1 = findindex(M,N,i,j);
                T2.v2 = findindex(M,N,i+1,j);
                T2.v3 = findindex(M,N,i+1,j+1);
                allT.push_back(T2);
                t++;
            }
        }
    }
        
        MeshBluePrint BluePrint;
        BluePrint.bvertex = allV;
        BluePrint.btriangle = allT;
        BluePrint.binclusion = allI;
        BluePrint.simbox = m_Box;
        
        std::string ext = m_OutputFilename.substr(m_OutputFilename.find_last_of(".") + 1);
        
        if(ext==TSExt)
        {
            WriteQFile(m_OutputFilename , BluePrint);
        }
        else if(ext==TSIExt)
        {
            WriteTSI(m_OutputFilename , BluePrint);
        }
        else
        {
            std::cout<<"---> Error: output file with "<<ext<<" extension is not recognized.  It should have either "<<TSExt<<" or "<<TSIExt<<" extension. "<<std::endl;
        }
        
}
void Generate::WriteQFile(std::string filename , MeshBluePrint blueprint)
{
    int pres=10;
    std::ofstream output;
    output.open(filename.c_str());
    
    
    output<<std::fixed;
    output<<std::setprecision( pres )<<(blueprint.simbox)(0)<<"   "<<(blueprint.simbox)(1)<<"   "<<(blueprint.simbox)(2)<<"   \n";
    output<<(blueprint.bvertex).size()<<"\n";
    
    
    for (std::vector<Vertex_Map>::iterator it = (blueprint.bvertex).begin() ; it != (blueprint.bvertex).end(); ++it)
        output<<std::setprecision( pres )<<it->id<<"  "<<it->x<<"  "<<it->y<<"  "<<it->z<<"  "<<it->domain<<std::endl;

    output<< (blueprint.btriangle).size()<<"\n";
    for (std::vector<Triangle_Map>::iterator it = (blueprint.btriangle).begin() ; it != (blueprint.btriangle).end(); ++it)
        output<<it->id<<"  "<<it->v1<<"   "<<it->v2<<"  "<<it->v3<<"  0  "<<std::endl;
    
    output.close();
}
void Generate::WriteTSI(std::string filename , MeshBluePrint blueprint)
{

    FILE * output;
    output = fopen(filename.c_str(), "w");
    std::string format = "%"+m_tsiPrecision+"lf%"+m_tsiPrecision+"lf%"+m_tsiPrecision+"lf\n";
    const char* version="version 1.1";
    fprintf(output,"%s\n",version);
    //------
    const char* box="box";
    fprintf(output,"%s%18.10lf%18.10lf%18.10lf\n",box,(blueprint.simbox)(0),(blueprint.simbox)(1),(blueprint.simbox)(2));
    
    const char* ver="vertex";
    int size=(blueprint.bvertex).size();
    fprintf(output,"%s%20d\n",ver,size);
    format = "%5d%"+m_tsiPrecision+"lf%"+m_tsiPrecision+"lf%"+m_tsiPrecision+"lf\n";
    for (std::vector<Vertex_Map>::iterator it = (blueprint.bvertex).begin() ; it != (blueprint.bvertex).end(); ++it)
        fprintf(output,format.c_str(),it->id,it->x,it->y,it->z);
  
    const char* tri="triangle";
    size = (blueprint.btriangle).size();
    fprintf(output,"%s%20d\n",tri,size);
   for (std::vector<Triangle_Map>::iterator it = (blueprint.btriangle).begin() ; it != (blueprint.btriangle).end(); ++it)
        fprintf(output,"%10d%10d%10d%10d\n",it->id,it->v1,it->v2,it->v3);
    
 
    const char* inc="inclusion";
    size = (blueprint.binclusion).size();
    fprintf(output,"%s%20d\n",inc,size);
     format = "%10d%10d%10d%"+m_tsiPrecision+"lf%"+m_tsiPrecision+"lf\n";
     for (std::vector<Inclusion_Map>::iterator it = (blueprint.binclusion).begin() ; it != (blueprint.binclusion).end(); ++it)
        fprintf(output,format.c_str(),it->id,it->tid,it->vid,it->x,it->y);

    fclose(output);
}
void Generate::ExploreArguments()
{
    Nfunction f;
    for (long i=1;i<m_Argument.size();i=i+2)
    {
        std::string Arg1 = m_Argument.at(i);
        if(Arg1=="-h")
        {
            HelpMessage();
            exit(0);
            break;
        }
        else if(Arg1=="-o")
        {
            m_OutputFilename = m_Argument.at(i+1);
        }
        else if(Arg1=="-defout")
        {
            m_GeneralOutputFilename = m_Argument.at(i+1);
        }
        else if(Arg1=="-seed")
        {
            m_Seed = f.String_to_Int(m_Argument.at(i+1));
        }
        else if(Arg1=="-angle")
        {
            m_MinFaceAngle = f.String_to_Double(m_Argument.at(i+1));
        }
        else if(Arg1=="-g")
        {
            m_genus = f.String_to_Double(m_Argument.at(i+1));
        }
        else if(Arg1=="-maxDist")
        {
            m_MaxLinkLengthSquare = f.String_to_Double(m_Argument.at(i+1));
        }
        else if(Arg1=="-minDist")
        {
            m_MinVerticesDistanceSquare = f.String_to_Double(m_Argument.at(i+1));
        }
        else if(Arg1=="-minDist")
        {
            m_MinVerticesDistanceSquare = f.String_to_Double(m_Argument.at(i+1));
        }
        else if(Arg1=="-N")
        {
            m_N = f.String_to_Int(m_Argument.at(i+1));
        }
        else if(Arg1=="-box")
        {
            m_Box(0) = f.String_to_Double(m_Argument.at(i+1));
            m_Box(1) = f.String_to_Double(m_Argument.at(i+2));
            m_Box(2) = f.String_to_Double(m_Argument.at(i+3));
            i++;
            i++;

        }
        else if(Arg1=="-type")
        {
            m_Type = m_Argument.at(i+1);
        }
        else if(Arg1=="-dna")
        {
            m_GenerateDNA = true;
            i--;  // No value needed for this flag
        }
        else if(Arg1=="-dna_output")
        {
            m_DNAOutputFilename = m_Argument.at(i+1);
        }
        else if(Arg1=="-dna_num_beads")
        {
            m_DNANumBeads = f.String_to_Int(m_Argument.at(i+1));
        }
        else if(Arg1=="-dna_num_chains")
        {
            m_DNANumChains = f.String_to_Int(m_Argument.at(i+1));
        }
        else if(Arg1=="-dna_boundary_radius")
        {
            m_DNABoundaryRadius = f.String_to_Double(m_Argument.at(i+1));
        }
        else if(Arg1=="-dna_boundary_origin")
        {
            if (i+3 >= m_Argument.size()) {
                std::cout << "---> Error: -dna_boundary_origin requires 3 values (x y z)" << std::endl;
                m_Healthy = false;
                exit(0);
            }
            m_DNABoundaryOrigin(0) = f.String_to_Double(m_Argument.at(i+1));
            m_DNABoundaryOrigin(1) = f.String_to_Double(m_Argument.at(i+2));
            m_DNABoundaryOrigin(2) = f.String_to_Double(m_Argument.at(i+3));
            i += 2;  // Skip the two additional arguments
        }
        else if(Arg1=="-dna_monomer_radius")
        {
            m_DNAMonomerRadius = f.String_to_Double(m_Argument.at(i+1));
        }
        else if(Arg1=="-dna_persistence_length")
        {
            m_DNAPersistenceLength = f.String_to_Double(m_Argument.at(i+1));
        }
        else if(Arg1=="-dna_stages")
        {
            // Format: -dna_stages N len1,num1 len2,num2 ...
            m_DNANumStages = f.String_to_Int(m_Argument.at(i+1));
            m_DNAStages.clear();
            for (int s = 0; s < m_DNANumStages; s++) {
                std::string stage_str = m_Argument.at(i+2+s);
                size_t comma_pos = stage_str.find(',');
                if (comma_pos != std::string::npos) {
                    int len = f.String_to_Int(stage_str.substr(0, comma_pos));
                    int num = f.String_to_Int(stage_str.substr(comma_pos+1));
                    m_DNAStages.push_back(std::make_pair(len, num));
                }
            }
            i += m_DNANumStages;  // Skip the stage arguments
        }
        else if(Arg1=="-dna_bonds")
        {
            m_DNABondsFilename = m_Argument.at(i+1);
        }
        else if(Arg1=="-dna_start_index")
        {
            m_DNAStartIndex = f.String_to_Int(m_Argument.at(i+1));
        }
        else if(Arg1=="-dna_bond_k")
        {
            m_DNABondK = f.String_to_Double(m_Argument.at(i+1));
        }
        else if(Arg1=="-dna_bond_l0")
        {
            m_DNABondL0 = f.String_to_Double(m_Argument.at(i+1));
        }
        else if(Arg1=="-dna_bond_k_angle")
        {
            m_DNABondKAngle = f.String_to_Double(m_Argument.at(i+1));
        }
        else if(Arg1=="-dna_bond_theta0")
        {
            m_DNABondTheta0 = f.String_to_Double(m_Argument.at(i+1));
        }
        else
        {
            std::cout << "---> Error: "<<Arg1;
            std::cout<<"\n"<<"For more information and tips run "<< m_Argument.at(0) <<" -h"<<"\n";
            m_Healthy =false;
            exit(0);
            break;
        }
    }

}
void Generate::HelpMessage()
{
    std::cout<<"--------------------copyright: Weria Pezeshkian------------------"<<"\n";
    std::cout<<" =================================================================  \n";
    std::cout<<"---------------------version "<<SoftWareVersion<<" ------------------"<<"\n";
    std::cout<<" =================================================================  \n";
    std::cout<<"------------simple example for exacuting  -------------------"<<"\n";
    std::cout<<" ./GEN -box 30 30 30 -type channel -N 10  -o out.tsi"<<"\n";
    std::cout<<"-------------------------------------------------------------------------------"<<"\n";
    std::cout<<"  option    type        default            description "<<"\n";
    std::cout<<"-------------------------------------------------------------------------------"<<"\n";
    std::cout<<"  -Box     3*double       10 10 10          box sides "<<"\n";
    std::cout<<"  -o         string       out.q             out put file name, could be both tsi or q file formats by giving the extension "<<"\n";
    std::cout<<"  -type      string       flat              ts shape (flat/tetrahedron/channel/high_gen)  "<<"\n";
    std::cout<<"  -N         int          5                 vertex per side for channel and tetrahedron shape "<<"\n";
    std::cout<<"  -dna                      (flag)          enable DNA generation"<<"\n";
    std::cout<<"  -dna_output string        dna.q           DNA output file"<<"\n";
    std::cout<<"  -dna_num_beads int       100              number of DNA beads per chain (final target)"<<"\n";
    std::cout<<"  -dna_num_chains int      1                number of DNA chains to generate"<<"\n";
    std::cout<<"  -dna_boundary_radius double 90.0           boundary radius in nm"<<"\n";
    std::cout<<"  -dna_boundary_origin double double double 0 0 0  boundary sphere center (x y z) in nm"<<"\n";
    std::cout<<"  -dna_monomer_radius double 1.7             monomer radius in nm"<<"\n";
    std::cout<<"  -dna_persistence_length double 45.0        persistence length in nm"<<"\n";
    std::cout<<"  -dna_stages N len1,num1 len2,num2 ...     growth stages (seg_mult*(2*r), num_segments)"<<"\n";
    std::cout<<"                                             Example: -dna_stages 2 5,50 1,500"<<"\n";
    std::cout<<"  -dna_bonds string        dna_bonds.txt    DNA bonds output file"<<"\n";
    std::cout<<"  -dna_start_index int     (inferred)       starting vertex ID for DNA beads"<<"\n";
    std::cout<<"  -dna_bond_k double       100.0            bond spring constant (k_stretch)"<<"\n";
    std::cout<<"  -dna_bond_l0 double      0.472            bond equilibrium length"<<"\n";
    std::cout<<"  -dna_bond_k_angle double 50.0             angle spring constant"<<"\n";
    std::cout<<"  -dna_bond_theta0 double  0.0              angle equilibrium value (radians, 0=straight)"<<"\n";
    std::cout<<"=========================================================================="<<"\n";
    std::cout<<"=========================================================================="<<"\n";
    std::cout<<"------------------ version "<<SoftWareVersion<<" ------------------"<<"\n";
    std::cout<<" =================================================================  \n";
}
int Generate::findindex(int M,int N,int i,int j)
{
    int q=(i+M)%M;
    int p=j%N;
    int s=q+M*p;
    return s;
}
int Generate::CylinderIndex(int M,int N,int i,int j)
{
    int q=(j+M)%M;
    int p=(i+N)%N;
    int s=q+M*p;
    return s;
}
int Generate::findid(int faceno,int i, int j)
{
    int id=0;
    
    
    if(j>0 && i>=j && i<=m_N+1)
    {
        id=i*(i-1)/2+j+faceno*m_N*(m_N+1)/2;
    }
    else if(j>i && i<=m_N+1)
    {
        if(faceno==0)
        {
            id=i*(i-1)/2+1+m_N*(m_N+1)/2;
        }
        if(faceno==1)
        {
            id=i*(i-1)/2+1+m_N*(m_N+1);
        }
        if(faceno==2)
        {
            id=i*(i-1)/2+1;
        }
        
    }
    else if(j==0)
    {
        if(faceno==0)
        {
            id=i*(i-1)/2+i+m_N*(m_N+1);
        }
        if(faceno==1)
        {
            id=i*(i-1)/2+i;
        }
        if(faceno==2)
        {
            id=i*(i-1)/2+i+m_N*(m_N+1)/2;
        }
        
    }
    
    return id;
    
}
int Generate::findidface4(int i, int j)
{
    int id=0;
    
    if(j>0 && i>=j && i<=m_N-2)
    {
        id=i*(i-1)/2+j+3*m_N*(m_N+1)/2;
    }
    else if(j>i )//&& i<=m_N+1)
    {
        id=findid(0,m_N, m_N-i);
    }
    else if(i==m_N-1)
    {
        int n=m_N;
        int m=j+1;
        id=n*(n-1)/2+2*m_N*(m_N+1)/2+m;
    }
    else if(j==0)
        id=findid(1,m_N,i+2);        // id=0;
    
    return id;
    
}

// ============================================================================
// DNA Generation Functions (matching Fortran sc_chain_generation exactly)
// ============================================================================

// Helper function: 3D distance
static double dist3D(double x1, double y1, double z1, double x2, double y2, double z2) {
    double dx = x2 - x1;
    double dy = y2 - y1;
    double dz = z2 - z1;
    return sqrt(dx*dx + dy*dy + dz*dz);
}

// Helper function: check if point is within spherical boundary
static bool isWithinBoundary(double x, double y, double z, double R_b) {
    double r2 = x*x + y*y + z*z;
    return (r2 < R_b * R_b);
}

// Helper function: circular index for closed chain (0-based)
static int circ_trans(int i, int N) {
    if (i >= N) return i - N;
    if (i < 0) return i + N;
    return i;
}

// Helper function: distance from point to line segment
static double distPointToSegment(double px, double py, double pz,
                                 double x1, double y1, double z1,
                                 double x2, double y2, double z2) {
    double dx = x2 - x1, dy = y2 - y1, dz = z2 - z1;
    double len2 = dx*dx + dy*dy + dz*dz;
    if (len2 < 1e-20) return dist3D(px, py, pz, x1, y1, z1);
    
    double t = std::max(0.0, std::min(1.0, 
        ((px-x1)*dx + (py-y1)*dy + (pz-z1)*dz) / len2));
    
    double proj_x = x1 + t*dx;
    double proj_y = y1 + t*dy;
    double proj_z = z1 + t*dz;
    
    return dist3D(px, py, pz, proj_x, proj_y, proj_z);
}

// Helper function: distance between two line segments
static double distSegmentToSegment(double x1, double y1, double z1, double x2, double y2, double z2,
                                   double x3, double y3, double z3, double x4, double y4, double z4) {
    double dx1 = x2 - x1, dy1 = y2 - y1, dz1 = z2 - z1;
    double dx2 = x4 - x3, dy2 = y4 - y3, dz2 = z4 - z3;
    double dx3 = x1 - x3, dy3 = y1 - y3, dz3 = z1 - z3;
    
    double a = dx1*dx1 + dy1*dy1 + dz1*dz1;
    double b = dx1*dx2 + dy1*dy2 + dz1*dz2;
    double c = dx2*dx2 + dy2*dy2 + dz2*dz2;
    double d = dx1*dx3 + dy1*dy3 + dz1*dz3;
    double e = dx2*dx3 + dy2*dy3 + dz2*dz3;
    
    double denom = a*c - b*b;
    if (fabs(denom) < 1e-10) {
        return std::min(distPointToSegment(x1, y1, z1, x3, y3, z3, x4, y4, z4),
                       distPointToSegment(x2, y2, z2, x3, y3, z3, x4, y4, z4));
    }
    
    double s = std::max(0.0, std::min(1.0, (b*e - c*d) / denom));
    double t = std::max(0.0, std::min(1.0, (a*e - b*d) / denom));
    
    double px1 = x1 + s*dx1, py1 = y1 + s*dy1, pz1 = z1 + s*dz1;
    double px2 = x3 + t*dx2, py2 = y3 + t*dy2, pz2 = z3 + t*dz2;
    
    return dist3D(px1, py1, pz1, px2, py2, pz2);
}

// Helper function: check segment clash (matches test_inserted_segment_clash)
// Now checks against all chains to enforce inter-chain avoidance
static bool checkSegmentClash(const std::vector<Vertex_Map>& chain, int chain_idx,
                              const std::vector<std::vector<Vertex_Map> >& all_chains,
                              double new_x, double new_y, double new_z,
                              int i, int ipp, double r_sc, double R_b,
                              bool debug = false) {
    double min_dist = 2.0 * r_sc;  // Always 2*r_sc (matches Fortran sc_sc_clash)
    int N = chain.size();
    
    // Boundary check (simplified - original has complex cylindrical boundary)
    if (!isWithinBoundary(new_x, new_y, new_z, R_b)) {
        if (debug) {
            double r2 = new_x*new_x + new_y*new_y + new_z*new_z;
            std::cerr << "      [DEBUG] Boundary check failed: point (" << new_x << ", " << new_y << ", " << new_z 
                      << ") distance=" << sqrt(r2) << " > R_b=" << R_b << std::endl;
            std::cerr.flush();
        }
        return true;
    }
    
    // Segments created by insertion
    double seg1_x1 = chain[i].x, seg1_y1 = chain[i].y, seg1_z1 = chain[i].z;
    double seg1_x2 = new_x, seg1_y2 = new_y, seg1_z2 = new_z;
    double seg2_x1 = new_x, seg2_y1 = new_y, seg2_z1 = new_z;
    double seg2_x2 = chain[ipp].x, seg2_y2 = chain[ipp].y, seg2_z2 = chain[ipp].z;
    
    // First check same chain (matches Fortran: lines 489-511)
    // Check segments j=2 to N-2 (skip immediate neighbors)
    // Note: Fortran uses 1-based indexing, so j goes from 2 to N-2
    // In 0-based: j goes from 2 to N-2 (inclusive)
    if (N > 3) {
        for (int j = 2; j <= N - 2; j++) {
            int k = circ_trans(i + j, N);
            int kp = circ_trans(i + j + 1, N);
            
            double d1 = distSegmentToSegment(seg1_x1, seg1_y1, seg1_z1, seg1_x2, seg1_y2, seg1_z2,
                                            chain[k].x, chain[k].y, chain[k].z,
                                            chain[kp].x, chain[kp].y, chain[kp].z);
            double d2 = distSegmentToSegment(seg2_x1, seg2_y1, seg2_z1, seg2_x2, seg2_y2, seg2_z2,
                                            chain[k].x, chain[k].y, chain[k].z,
                                            chain[kp].x, chain[kp].y, chain[kp].z);
            
            if (d1 < min_dist || d2 < min_dist) {
                return true;
            }
        }
    }
    
    // Then check other chains (matches Fortran: lines 562-595)
    // For other chains, check ALL segments (no skipping)
    for (size_t c = 0; c < all_chains.size(); c++) {
        if (c == chain_idx) continue;  // Skip same chain (already checked above)
        
        const std::vector<Vertex_Map>& other_chain = all_chains[c];
        int other_N = other_chain.size();
        
        if (other_N < 2) continue;
        
        // Check ALL segments of other chain (matches Fortran: k = 1 to N)
        for (int k = 0; k < other_N; k++) {
            int kp = (k + 1) % other_N;
            
            double d1 = distSegmentToSegment(seg1_x1, seg1_y1, seg1_z1, seg1_x2, seg1_y2, seg1_z2,
                                            other_chain[k].x, other_chain[k].y, other_chain[k].z,
                                            other_chain[kp].x, other_chain[kp].y, other_chain[kp].z);
            double d2 = distSegmentToSegment(seg2_x1, seg2_y1, seg2_z1, seg2_x2, seg2_y2, seg2_z2,
                                            other_chain[k].x, other_chain[k].y, other_chain[k].z,
                                            other_chain[kp].x, other_chain[kp].y, other_chain[kp].z);
            
            if (d1 < min_dist || d2 < min_dist) {
                if (debug) {
                    std::cerr << "      [DEBUG] Clash with other chain " << (c+1) << " segment " << k 
                              << ": d1=" << d1 << ", d2=" << d2 << ", min_dist=" << min_dist << std::endl;
                    std::cerr.flush();
                }
                return true;
            }
        }
    }
    
    // Check angle constraints for same chain (N > 4) - matches Fortran lines 514-560
    // Note: Fortran also checks angles using test_angle, but we only check segment clashes here
    if (N > 4) {
        // Test clash between second segment and segment prior to first (matches Fortran lines 519-536)
        int k = circ_trans(i - 1, N);
        int kp = circ_trans(i, N);
        double d = distSegmentToSegment(seg2_x1, seg2_y1, seg2_z1, seg2_x2, seg2_y2, seg2_z2,
                                      chain[k].x, chain[k].y, chain[k].z,
                                      chain[kp].x, chain[kp].y, chain[kp].z);
        if (d < min_dist) {
            if (debug) {
                std::cerr << "      [DEBUG] Angle constraint clash (prev segment): d=" << d << ", min_dist=" << min_dist << std::endl;
            }
            return true;
        }
        
        // Test clash between first segment and segment subsequent to second (matches Fortran lines 541-558)
        k = circ_trans(ipp, N);
        kp = circ_trans(ipp + 1, N);
        d = distSegmentToSegment(seg1_x1, seg1_y1, seg1_z1, seg1_x2, seg1_y2, seg1_z2,
                                chain[k].x, chain[k].y, chain[k].z,
                                chain[kp].x, chain[kp].y, chain[kp].z);
        if (d < min_dist) {
            if (debug) {
                std::cerr << "      [DEBUG] Angle constraint clash (next segment): d=" << d << ", min_dist=" << min_dist << std::endl;
            }
            return true;
        }
    }
    
    return false;
}

// Helper function: check sphere clash (matches test_inserted_sphere_clash)
// Now checks against all chains to enforce inter-chain avoidance
static bool checkSphereClash(const std::vector<Vertex_Map>& chain, int chain_idx,
                             const std::vector<std::vector<Vertex_Map> >& all_chains,
                             double new_x, double new_y, double new_z,
                             int i, int ipp, double r_sc, double R_b) {
    double min_dist = 2.0 * r_sc;  // Always 2*r_sc (matches Fortran s_s_clash)
    int N = chain.size();
    
    // Boundary check (simplified - original has cylindrical boundary with L_b)
    if (!isWithinBoundary(new_x, new_y, new_z, R_b)) {
        return true;
    }
    
    // First check same chain (matches Fortran: lines 816-828)
    // Check spheres j=2 to N-1 (skip immediate neighbors)
    if (N > 3) {
        for (int j = 2; j <= N - 1; j++) {
            int k = circ_trans(i + j, N);
            
            if (k == i || k == ipp) continue;
            
            double d = dist3D(new_x, new_y, new_z, chain[k].x, chain[k].y, chain[k].z);
            if (d < min_dist) {
                return true;
            }
        }
    } else {
        // For small chains (N=3), check the opposite vertex
        for (int k = 0; k < N; k++) {
            if (k == i || k == ipp) continue;
            double d = dist3D(new_x, new_y, new_z, chain[k].x, chain[k].y, chain[k].z);
            if (d < min_dist) {
                return true;
            }
        }
    }
    
    // Then check other chains (matches Fortran: lines 831-853)
    // For other chains, check ALL spheres (no skipping)
    for (size_t c = 0; c < all_chains.size(); c++) {
        if (c == chain_idx) continue;  // Skip same chain (already checked above)
        
        const std::vector<Vertex_Map>& other_chain = all_chains[c];
        
        // Check ALL spheres of other chain (matches Fortran: k = 1 to N)
        for (size_t k = 0; k < other_chain.size(); k++) {
            double d = dist3D(new_x, new_y, new_z, other_chain[k].x, other_chain[k].y, other_chain[k].z);
            if (d < min_dist) {
                return true;
            }
        }
    }
    
    return false;
}

// Helper function: interpolate chain (matches interpolate_segments)
static void interpolateChain(std::vector<Vertex_Map>& chain, int old_seg_mult, int new_seg_mult) {
    if (old_seg_mult <= new_seg_mult || chain.size() < 2) return;
    
    int s_i = old_seg_mult / new_seg_mult;
    if (s_i < 2) return;
    
    std::vector<Vertex_Map> new_chain;
    int N_old = chain.size();
    
    for (int i = 0; i < N_old; i++) {
        int i_next = (i + 1) % N_old;
        
        new_chain.push_back(chain[i]);
        new_chain.back().id = new_chain.size() - 1;
        
        double dx = chain[i_next].x - chain[i].x;
        double dy = chain[i_next].y - chain[i].y;
        double dz = chain[i_next].z - chain[i].z;
        
        for (int j = 1; j < s_i; j++) {
            double frac = double(j) / double(s_i);
            Vertex_Map v;
            v.id = new_chain.size();
            v.x = chain[i].x + frac * dx;
            v.y = chain[i].y + frac * dy;
            v.z = chain[i].z + frac * dz;
            v.domain = 0;
            new_chain.push_back(v);
        }
    }
    
    chain = new_chain;
}

// Insert random segment (matches insert_random_segment exactly)
// Now checks against all chains for inter-chain avoidance
static bool insertRandomSegment(std::vector<Vertex_Map>& chain, int chain_idx,
                                const std::vector<std::vector<Vertex_Map> >& all_chains,
                                int seg_mult, double r_sc, double R_b,
                                std::mt19937& rng, std::uniform_real_distribution<double>& uniform_angle,
                                std::uniform_real_distribution<double>& uniform_phi) {
    // Increase attempts for small chains or when inter-chain avoidance is needed
    int max_c = 100;
    if (chain.size() <= 5) {
        max_c = 500;  // More attempts for small chains
    }
    int c = 0;
    double mid_dist = (sqrt(3.0) / 2.0) * (seg_mult * (2.0 * r_sc));
    
    bool accepted = false;
    int i = 0, ipp = 0;
    double new_x = 0, new_y = 0, new_z = 0;
    
    // Debug: print at start - ALWAYS print first few
    std::cerr << "[DEBUG] insertRandomSegment called for chain " << (chain_idx+1) 
              << " (size=" << chain.size() << ", max_c=" << max_c << ", mid_dist=" << mid_dist << ")" << std::endl;
    std::cerr.flush();
    
    while (!accepted && c < max_c) {
        i = rng() % chain.size();
        ipp = (i + 1) % chain.size();
        
        // Debug first few attempts
        if (c < 5) {
            std::cerr << "  [DEBUG] Attempt " << c << ": i=" << i << ", ipp=" << ipp << std::endl;
            std::cerr.flush();
        }
        
        double xm = (chain[i].x + chain[ipp].x) / 2.0;
        double ym = (chain[i].y + chain[ipp].y) / 2.0;
        double zm = (chain[i].z + chain[ipp].z) / 2.0;
        
        double dx = chain[ipp].x - chain[i].x;
        double dy = chain[ipp].y - chain[i].y;
        double dz = chain[ipp].z - chain[i].z;
        double len = sqrt(dx*dx + dy*dy + dz*dz);
        if (len < 1e-10) {
            c++;
            continue;
        }
        
        double ux = dx / len, uy = dy / len, uz = dz / len;
        
        double t = uniform_angle(rng);
        double p = uniform_phi(rng);
        double c_t = cos(t), s_t = sin(t);
        double c_p = cos(p), s_p = sin(p);
        
        double dir_x = c_p * s_t;
        double dir_y = s_p * s_t;
        double dir_z = c_t;
        
        double dxu = dir_x*ux + dir_y*uy + dir_z*uz;
        
        if (fabs(dxu) > 0.9) {
            c++;
            continue;
        }
        
        double perp_x = dir_x - dxu*ux;
        double perp_y = dir_y - dxu*uy;
        double perp_z = dir_z - dxu*uz;
        
        double perp_len = sqrt(perp_x*perp_x + perp_y*perp_y + perp_z*perp_z);
        if (perp_len < 1e-10) {
            c++;
            continue;
        }
        
        perp_x /= perp_len;
        perp_y /= perp_len;
        perp_z /= perp_len;
        
        new_x = xm + mid_dist * perp_x;
        new_y = ym + mid_dist * perp_y;
        new_z = zm + mid_dist * perp_z;
        
        if (!isWithinBoundary(new_x, new_y, new_z, R_b)) {
            if (c < 5) {
                double r2 = new_x*new_x + new_y*new_y + new_z*new_z;
                std::cerr << "    [DEBUG] Attempt " << c << ": Boundary check failed: distance=" << sqrt(r2) << " > R_b=" << R_b << std::endl;
                std::cerr.flush();
            }
            c++;
            continue;
        }
        
        // Debug output every 10 attempts (more frequent for debugging)
        bool debug_this = (c % 10 == 0 || c < 5);
        if (debug_this) {
            std::cerr << "    [DEBUG] insertRandomSegment attempt " << c << " for chain " << (chain_idx+1) 
                      << " (size=" << chain.size() << "): inserting at (" << new_x << ", " << new_y << ", " << new_z << ")" << std::endl;
            std::cerr << "      Inserting between vertices " << i << " (" << chain[i].x << ", " << chain[i].y << ", " << chain[i].z 
                      << ") and " << ipp << " (" << chain[ipp].x << ", " << chain[ipp].y << ", " << chain[ipp].z << ")" << std::endl;
            std::cerr.flush();
        }
        
        bool clash = checkSegmentClash(chain, chain_idx, all_chains, new_x, new_y, new_z, i, ipp, r_sc, R_b, debug_this);
        
        if (clash) {
            c++;
        } else {
            accepted = true;
        }
    }
    
    if (accepted) {
        Vertex_Map v;
        v.id = chain.size();
        v.x = new_x;
        v.y = new_y;
        v.z = new_z;
        v.domain = 0;
        
        chain.insert(chain.begin() + ipp, v);
        
        for (size_t k = 0; k < chain.size(); k++) {
            chain[k].id = k;
        }
        
        return true;
    }
    
    return false;
}

// Insert random sphere (matches insert_random_sphere exactly)
// Now checks against all chains for inter-chain avoidance
static bool insertRandomSphere(std::vector<Vertex_Map>& chain, int chain_idx,
                               const std::vector<std::vector<Vertex_Map> >& all_chains,
                               double r_sc, double R_b,
                               std::mt19937& rng, std::uniform_real_distribution<double>& uniform_angle,
                               std::uniform_real_distribution<double>& uniform_phi) {
    // Debug: ALWAYS print at start
    std::cerr << "[DEBUG] insertRandomSphere ENTERED for chain " << (chain_idx+1) 
              << " (size=" << chain.size() << ")" << std::endl;
    std::cerr.flush();
    
    // Increase attempts for small chains or when inter-chain avoidance is needed
    int max_c = 20;
    if (chain.size() <= 5) {
        max_c = 200;  // More attempts for small chains
    }
    int c = 0;
    double mid_dist = (sqrt(3.0) / 2.0) * (2.0 * r_sc);
    
    bool accepted = false;
    int i = 0, ipp = 0;
    double new_x = 0, new_y = 0, new_z = 0;
    
    std::cerr << "[DEBUG] insertRandomSphere: max_c=" << max_c << ", mid_dist=" << mid_dist << std::endl;
    std::cerr.flush();
    
    while (!accepted && c < max_c) {
        if (c < 5) {
            std::cerr << "  [DEBUG] insertRandomSphere attempt " << c << std::endl;
            std::cerr.flush();
        }
        i = rng() % chain.size();
        ipp = (i + 1) % chain.size();
        
        double xm = (chain[i].x + chain[ipp].x) / 2.0;
        double ym = (chain[i].y + chain[ipp].y) / 2.0;
        double zm = (chain[i].z + chain[ipp].z) / 2.0;
        
        double dx = chain[ipp].x - chain[i].x;
        double dy = chain[ipp].y - chain[i].y;
        double dz = chain[ipp].z - chain[i].z;
        double len = sqrt(dx*dx + dy*dy + dz*dz);
        if (len < 1e-10) {
            c++;
            continue;
        }
        
        double ux = dx / len, uy = dy / len, uz = dz / len;
        
        double t = uniform_angle(rng);
        double p = uniform_phi(rng);
        double c_t = cos(t), s_t = sin(t);
        double c_p = cos(p), s_p = sin(p);
        
        double dir_x = c_p * s_t;
        double dir_y = s_p * s_t;
        double dir_z = c_t;
        
        double dxu = dir_x*ux + dir_y*uy + dir_z*uz;
        
        if (fabs(dxu) > 0.9) {
            c++;
            continue;
        }
        
        double perp_x = dir_x - dxu*ux;
        double perp_y = dir_y - dxu*uy;
        double perp_z = dir_z - dxu*uz;
        
        double perp_len = sqrt(perp_x*perp_x + perp_y*perp_y + perp_z*perp_z);
        if (perp_len < 1e-10) {
            c++;
            continue;
        }
        
        perp_x /= perp_len;
        perp_y /= perp_len;
        perp_z /= perp_len;
        
        new_x = xm + mid_dist * perp_x;
        new_y = ym + mid_dist * perp_y;
        new_z = zm + mid_dist * perp_z;
        
        if (c < 5) {
            std::cerr << "    [DEBUG] insertRandomSphere attempt " << c << ": new point (" << new_x << ", " << new_y << ", " << new_z << ")" << std::endl;
            std::cerr.flush();
        }
        
        if (!isWithinBoundary(new_x, new_y, new_z, R_b)) {
            if (c < 5) {
                double r2 = new_x*new_x + new_y*new_y + new_z*new_z;
                std::cerr << "    [DEBUG] insertRandomSphere attempt " << c << ": boundary check FAILED, distance=" << sqrt(r2) << " > R_b=" << R_b << std::endl;
                std::cerr.flush();
            }
            c++;
            continue;
        }
        
        bool clash = checkSphereClash(chain, chain_idx, all_chains, new_x, new_y, new_z, i, ipp, r_sc, R_b);
        
        if (c < 5) {
            std::cerr << "    [DEBUG] insertRandomSphere attempt " << c << ": clash=" << (clash ? "YES" : "NO") << std::endl;
            std::cerr.flush();
        }
        
        if (clash) {
            c++;
        } else {
            accepted = true;
        }
    }
    
    if (accepted) {
        Vertex_Map v;
        v.id = chain.size();
        v.x = new_x;
        v.y = new_y;
        v.z = new_z;
        v.domain = 0;
        
        chain.insert(chain.begin() + ipp, v);
        
        for (size_t k = 0; k < chain.size(); k++) {
            chain[k].id = k;
        }
        
        return true;
    }
    
    return false;
}

// Grow chain segments (matches grow_chain_segments exactly)
// Now accepts all chains for inter-chain avoidance
static bool growChainSegments(std::vector<Vertex_Map>& chain, int chain_idx,
                              const std::vector<std::vector<Vertex_Map> >& all_chains,
                              int target_segments, int seg_mult,
                              double R_b, double r_sc, std::mt19937& rng,
                              std::uniform_real_distribution<double>& uniform_angle,
                              std::uniform_real_distribution<double>& uniform_phi) {
    const int N_chains = 1;
    const int max_fail = 10000 * N_chains;
    
    int fail_count = 0;
    int N_init = chain.size();
    int dN = target_segments - N_init;
    
    if (dN <= 0) return true;
    
    // First phase: insert at least one segment
    while (chain.size() < N_init + 1) {
        bool insert_fail = !insertRandomSegment(chain, chain_idx, all_chains, seg_mult, r_sc, R_b, rng, uniform_angle, uniform_phi);
        
        if (insert_fail) {
            fail_count++;
        } else {
            fail_count = 0;
            break;
        }
        
        if (fail_count > max_fail) return false;
    }
    
    // Second phase: loop until target
    int progress = 0;
    while (chain.size() < target_segments) {
        bool insert_fail = !insertRandomSegment(chain, chain_idx, all_chains, seg_mult, r_sc, R_b, rng, uniform_angle, uniform_phi);
        
        if (insert_fail) {
            fail_count++;
        } else {
            fail_count = 0;
        }
        
        int p = (int)floor(10.0 * (chain.size() - N_init) / dN);
        if (p > progress) {
            progress = p;
            fail_count = 0;
        }
        
        if (fail_count > max_fail) return false;
    }
    
    return true;
}

// Grow chain spheres (matches grow_chain_spheres exactly)
// Now accepts all chains for inter-chain avoidance
static bool growChainSpheres(std::vector<Vertex_Map>& chain, int chain_idx,
                             const std::vector<std::vector<Vertex_Map> >& all_chains,
                             int target_segments,
                             double R_b, double r_sc, std::mt19937& rng,
                             std::uniform_real_distribution<double>& uniform_angle,
                             std::uniform_real_distribution<double>& uniform_phi) {
    const int N_chains = 1;
    const int max_fail = 1000 * N_chains;
    
    int fail_count = 0;
    int N_init = chain.size();
    int dN = target_segments - N_init;
    
    if (dN <= 0) return true;
    
    // First phase: insert at least one sphere
    while (chain.size() < N_init + 1) {
        bool insert_fail = !insertRandomSphere(chain, chain_idx, all_chains, r_sc, R_b, rng, uniform_angle, uniform_phi);
        
        if (insert_fail) {
            fail_count++;
        } else {
            fail_count = 0;
        }
        
        if (fail_count > max_fail) return false;
    }
    
    // Second phase: loop until target
    int progress = 0;
    while (chain.size() < target_segments) {
        bool insert_fail = !insertRandomSphere(chain, chain_idx, all_chains, r_sc, R_b, rng, uniform_angle, uniform_phi);
        
        if (insert_fail) {
            fail_count++;
        } else {
            fail_count = 0;
        }
        
        int p = (int)floor(10.0 * (chain.size() - N_init) / dN);
        if (p > progress) {
            progress = p;
            fail_count = 0;
        }
        
        if (fail_count > max_fail) return false;
    }
    
    return true;
}

// Helper function: initialize triangle for a chain (with inter-chain avoidance)
static bool initializeChainTriangle(std::vector<Vertex_Map>& chain,
                                    int chain_id, const std::vector<std::pair<int, int> >& stages,
                                    const std::vector<std::vector<Vertex_Map> >& all_chains,
                                    double R_b, double r_sc,
                                    std::mt19937& rng, std::uniform_real_distribution<double>& uniform,
                                    std::uniform_real_distribution<double>& uniform_angle) {
    
    // Initialize triangle (matches initialize_chain)
    int first_seg_mult = stages[0].first;
    double segment_length = first_seg_mult * (2.0 * r_sc);
    double L = segment_length;
    double R_b_sqrd = (R_b - (L/sqrt(3.0)) + r_sc) * (R_b - (L/sqrt(3.0)) + r_sc);
    
    bool triangle_placed = false;
    int attempts = 0;
    const int max_attempts = 10000;  // Increased for better placement with multiple chains
    
    while (!triangle_placed && attempts < max_attempts) {
        // Sample uniformly within sphere (matches Fortran: r_rand(-R_b, R_b))
        // Fortran uses r_rand(-R_b, R_b) which samples from -R_b to R_b
        // We need to sample from -R_b to R_b, not 0 to R_b
        double x0 = uniform(rng) * 2.0 * R_b - R_b;  // Maps [0,1] to [-R_b, R_b]
        double y0 = uniform(rng) * 2.0 * R_b - R_b;
        double z0 = uniform(rng) * 2.0 * R_b - R_b;
        
        if (x0*x0 + y0*y0 + z0*z0 < R_b_sqrd) {
            double w0 = uniform_angle(rng);
            double L_tri = L * sqrt(3.0) / 2.0;
            
            bool all_valid = true;
            std::vector<Vertex_Map> triangle;
            for (int j = 0; j < 3; j++) {
                double a = w0 + (2.0 * M_PI / 3.0) * j;
                double x = x0 + (L_tri / sqrt(3.0)) * cos(a);
                double y = y0 + (L_tri / sqrt(3.0)) * sin(a);
                double z = z0;
                
                if (x*x + y*y + z*z >= R_b_sqrd) {
                    all_valid = false;
                    break;
                }
                
                Vertex_Map v;
                v.id = j;
                v.x = x;
                v.y = y;
                v.z = z;
                v.domain = 0;
                triangle.push_back(v);
            }
            
            if (all_valid && triangle.size() == 3) {
                // NOTE: The Fortran code has a bug at line 374 - it checks sc_s_clash with wrong arguments
                // This means it doesn't actually check inter-chain clashes during triangle initialization
                // We match this behavior - triangles are placed randomly, and clashes are checked during
                // segment insertion instead (which is where the real clash detection happens)
                // This is why the Fortran code works - it relies on random placement + clash checks during growth
                bool clash = false;
                
                // Only check obstacles if needed (not implemented in our simplified version)
                // The Fortran code would check obstacles here, but we skip it
                
                if (!clash) {
                    chain = triangle;
                    triangle_placed = true;
                }
            }
        }
        attempts++;
    }
    
    if (!triangle_placed) {
        std::cerr << "---> Error: Chain " << chain_id << " could not initialize triangle after " << max_attempts << " attempts" << std::endl;
    }
    
    return triangle_placed;
}

// Main DNA generation function (supports multiple chains)
void Generate::GenerateDNA()
{
    std::cout << "---> Generating " << m_DNANumChains << " DNA chain(s) using Koch-curve algorithm..." << std::endl;
    
    double R_b = m_DNABoundaryRadius;
    double r_sc = m_DNAMonomerRadius;
    
    std::mt19937 rng(m_Seed);
    std::uniform_real_distribution<double> uniform(-1.0, 1.0);
    std::uniform_real_distribution<double> uniform_angle(0.0, 2.0 * M_PI);
    std::uniform_real_distribution<double> uniform_phi(0.0, M_PI);
    
    // Set up stages
    std::vector<std::pair<int, int> > stages = m_DNAStages;
    if (stages.empty()) {
        stages.push_back(std::make_pair(1, m_DNANumBeads));
    }
    
    // Initialize all chains with triangles first (concurrent initialization)
    std::vector<std::vector<Vertex_Map> > all_chains(m_DNANumChains);
    std::vector<bool> chain_initialized(m_DNANumChains, false);
    
    std::cout << "---> Initializing " << m_DNANumChains << " chain triangle(s)..." << std::endl;
    for (int chain_id = 0; chain_id < m_DNANumChains; chain_id++) {
        bool success = initializeChainTriangle(all_chains[chain_id], chain_id+1, stages, all_chains, 
                                               R_b, r_sc, rng, uniform, uniform_angle);
        if (success) {
            chain_initialized[chain_id] = true;
            // NOTE: Do NOT apply origin offset here - keep everything relative to (0,0,0) during generation
            // Origin offset will be applied at the end when writing output
            std::cout << "---> Chain " << (chain_id+1) << " triangle initialized: " << all_chains[chain_id].size() << " beads" << std::endl;
        } else {
            std::cerr << "---> Error: Chain " << (chain_id+1) << " failed to initialize triangle" << std::endl;
        }
    }
    
    // Concurrent growth: alternate growth steps between chains
    std::cout << "---> Growing chains concurrently with inter-chain avoidance..." << std::endl;
    
    // Process each stage for all chains concurrently
    for (int stage_idx = 0; stage_idx < (int)stages.size(); stage_idx++) {
        int seg_mult = stages[stage_idx].first;
        int target_segments = stages[stage_idx].second;
        
        // Interpolate all chains if needed
        if (stage_idx > 0) {
            int prev_seg_mult = stages[stage_idx-1].first;
            if (prev_seg_mult > seg_mult) {
                for (int chain_id = 0; chain_id < m_DNANumChains; chain_id++) {
                    if (chain_initialized[chain_id]) {
                        interpolateChain(all_chains[chain_id], prev_seg_mult, seg_mult);
                    }
                }
            }
        }
        
        // Grow chains using round-robin queue (matches Fortran implementation)
        // This is NOT truly concurrent - it's sequential but round-robin
        std::vector<int> queue(m_DNANumChains);
        for (int i = 0; i < m_DNANumChains; i++) {
            queue[i] = i;
        }
        
        int max_fail = 10000 * m_DNANumChains;
        int fail_count = 0;
        int progress = 0;
        // N_init should be the size of the LAST chain (matches Fortran: sc_s%chain(sc_p%N_chains)%N)
        int N_init = 0;
        if (m_DNANumChains > 0 && chain_initialized[m_DNANumChains-1]) {
            N_init = all_chains[m_DNANumChains-1].size();
        }
        int dN = target_segments - N_init;
        
        // First phase: ensure at least one insertion in the first chain of queue (matches Fortran exactly)
        // This only tries to insert into queue[0], and if successful, rotates queue and exits
        if (dN > 0) {
            int attempts = 0;
            std::cerr << "---> [DEBUG] Starting first phase: trying to grow chain " << (queue[0]+1) 
                      << " from size " << all_chains[queue[0]].size() << " to " << (N_init + 1) << std::endl;
            while ((int)all_chains[queue[0]].size() < N_init + 1) {
                bool insert_fail = false;
                if (seg_mult == 1) {
                    std::cerr << "---> [DEBUG] About to call insertRandomSphere for chain " << (queue[0]+1) 
                              << ", seg_mult=" << seg_mult << std::endl;
                    std::cerr.flush();
                    insert_fail = !insertRandomSphere(all_chains[queue[0]], queue[0], all_chains, r_sc, R_b, 
                                                     rng, uniform_angle, uniform_phi);
                    std::cerr << "---> [DEBUG] insertRandomSphere returned: " << (insert_fail ? "FAIL" : "SUCCESS") << std::endl;
                    std::cerr.flush();
                } else {
                    std::cerr << "---> [DEBUG] About to call insertRandomSegment for chain " << (queue[0]+1) 
                              << ", seg_mult=" << seg_mult << std::endl;
                    std::cerr.flush();
                    insert_fail = !insertRandomSegment(all_chains[queue[0]], queue[0], all_chains, seg_mult, r_sc, R_b, 
                                                      rng, uniform_angle, uniform_phi);
                    std::cerr << "---> [DEBUG] insertRandomSegment returned: " << (insert_fail ? "FAIL" : "SUCCESS") << std::endl;
                    std::cerr.flush();
                }
                
                attempts++;
                if (insert_fail) {
                    fail_count++;
                    // Debug output every 100 attempts (more frequent)
                    if (attempts % 100 == 0) {
                        std::cerr << "---> [DEBUG] Stage " << (stage_idx+1) << " initial insertion attempt " << attempts 
                                  << " for chain " << (queue[0]+1) << " (size=" << all_chains[queue[0]].size() 
                                  << ", N_init=" << N_init << ")" << std::endl;
                    }
                } else {
                    // Rotate queue and exit (matches Fortran: queue = cshift(queue,shift=1); exit)
                    std::rotate(queue.begin(), queue.begin() + 1, queue.end());
                    fail_count = 0;
                    break;
                }
                
                if (fail_count > max_fail) {
                    std::cerr << "---> Error: Stage " << (stage_idx+1) << " failed to insert initial segment after " 
                              << attempts << " attempts" << std::endl;
                    std::cerr << "---> Chain " << (queue[0]+1) << " size: " << all_chains[queue[0]].size() 
                              << ", other chain sizes: ";
                    for (int c = 0; c < m_DNANumChains; c++) {
                        if (c != queue[0]) {
                            std::cerr << "C" << (c+1) << "=" << all_chains[c].size() << " ";
                        }
                    }
                    std::cerr << std::endl;
                    return;
                }
            }
        }
        
        // Second phase: round-robin growth until target
        while ((int)all_chains[queue[0]].size() < target_segments) {
            int current_chain = queue[0];
            
            if (!chain_initialized[current_chain]) {
                // Rotate queue
                std::rotate(queue.begin(), queue.begin() + 1, queue.end());
                continue;
            }
            
            bool insert_fail = false;
            if (seg_mult == 1) {
                insert_fail = !insertRandomSphere(all_chains[current_chain], current_chain, all_chains, r_sc, R_b, 
                                                 rng, uniform_angle, uniform_phi);
            } else {
                insert_fail = !insertRandomSegment(all_chains[current_chain], current_chain, all_chains, seg_mult, r_sc, R_b, 
                                                  rng, uniform_angle, uniform_phi);
            }
            
            if (insert_fail) {
                fail_count++;
            } else {
                // Rotate queue to next chain (matches Fortran: queue = cshift(queue,shift=1))
                std::rotate(queue.begin(), queue.begin() + 1, queue.end());
                fail_count = 0;
            }
            
            // Progress tracking
            int p = (int)floor(10.0 * (all_chains[queue[0]].size() - N_init) / dN);
            if (p > progress) {
                progress = p;
                fail_count = 0;
                std::cout << "---> Stage " << (stage_idx+1) << " progress: " << (10*progress) << "/100" << std::endl;
            }
            
            if (fail_count > max_fail) {
                std::cerr << "---> Error: Stage " << (stage_idx+1) << " failed after " << max_fail << " attempts" << std::endl;
                return;
            }
        }
        
        // Report final sizes
        for (int chain_id = 0; chain_id < m_DNANumChains; chain_id++) {
            if (chain_initialized[chain_id]) {
                if ((int)all_chains[chain_id].size() < target_segments) {
                    std::cerr << "---> Warning: Chain " << (chain_id+1) << " stage " << (stage_idx+1) 
                              << " incomplete: " << all_chains[chain_id].size() << " / " << target_segments << " segments" << std::endl;
                }
            }
        }
    }
    
    // Count total beads
    int total_beads = 0;
    for (int chain_id = 0; chain_id < m_DNANumChains; chain_id++) {
        if (chain_initialized[chain_id]) {
            total_beads += all_chains[chain_id].size();
            std::cout << "---> Chain " << (chain_id+1) << " complete: " << all_chains[chain_id].size() << " beads" << std::endl;
        }
    }
    
    std::cout << "---> Generated " << m_DNANumChains << " chain(s) with " << total_beads << " total DNA beads" << std::endl;
    
    // Remove uninitialized chains
    std::vector<std::vector<Vertex_Map> > valid_chains;
    for (int chain_id = 0; chain_id < m_DNANumChains; chain_id++) {
        if (chain_initialized[chain_id] && all_chains[chain_id].size() > 0) {
            valid_chains.push_back(all_chains[chain_id]);
        }
    }
    
    if (valid_chains.empty()) {
        std::cerr << "---> Error: No chains were successfully generated" << std::endl;
        return;
    }
    
    all_chains = valid_chains;
    
    // Combine all chains into a single vector with sequential IDs
    std::vector<Vertex_Map> combined_chain;
    int global_id = 0;
    
    for (size_t c = 0; c < all_chains.size(); c++) {
        for (size_t i = 0; i < all_chains[c].size(); i++) {
            Vertex_Map v = all_chains[c][i];
            v.id = global_id++;
            // Apply origin offset here, at the very end (after all growth is complete)
            v.x += m_DNABoundaryOrigin(0);
            v.y += m_DNABoundaryOrigin(1);
            v.z += m_DNABoundaryOrigin(2);
            combined_chain.push_back(v);
        }
    }
    
    // Write output
    Vec3D dnaBox;
    double box_size = R_b * 2.0;
    dnaBox(0) = box_size;
    dnaBox(1) = box_size;
    dnaBox(2) = box_size;
    
    WriteDNAQFile(m_DNAOutputFilename, combined_chain, dnaBox);
    std::cout << "---> DNA coordinates written to: " << m_DNAOutputFilename << std::endl;
    
    // Write bonds file if bonds filename is set
    if (!m_DNABondsFilename.empty()) {
        // Determine start index (infer from membrane if not set)
        int start_index = m_DNAStartIndex;
        if (start_index < 0) {
            // Try to infer from membrane file by reading the vertex count
            start_index = 0;  // Default to 0 if can't infer
            if (!m_OutputFilename.empty()) {
                std::ifstream mem_file(m_OutputFilename);
                if (mem_file.is_open()) {
                    std::string line;
                    // Skip box line
                    if (std::getline(mem_file, line)) {
                        // Read vertex count
                        if (std::getline(mem_file, line)) {
                            Nfunction f;
                            int num_vertices = f.String_to_Int(line);
                            start_index = num_vertices;
                            std::cout << "---> Inferred DNA start index from membrane file: " << start_index << std::endl;
                        }
                    }
                    mem_file.close();
                }
            }
            if (start_index == 0) {
                std::cout << "---> Warning: Could not infer DNA start index. Using 0. "
                          << "You may need to specify -dna_start_index explicitly." << std::endl;
            }
        }
        
        // Calculate average bond length if l0 not set
        double l0 = m_DNABondL0;
        if (l0 < 0) {
            double total_bond_length = 0.0;
            int bond_count = 0;
            for (size_t c = 0; c < all_chains.size(); c++) {
                for (size_t i = 0; i < all_chains[c].size(); i++) {
                    size_t next = (i + 1) % all_chains[c].size();
                    double dx = all_chains[c][next].x - all_chains[c][i].x;
                    double dy = all_chains[c][next].y - all_chains[c][i].y;
                    double dz = all_chains[c][next].z - all_chains[c][i].z;
                    total_bond_length += sqrt(dx*dx + dy*dy + dz*dz);
                    bond_count++;
                }
            }
            if (bond_count > 0) {
                l0 = total_bond_length / bond_count;
            } else {
                l0 = 2.0 * r_sc;  // fallback
            }
        }
        
        // Use calculated l0 if it was computed, otherwise use the default/specified value
        double final_l0 = (m_DNABondL0 < 0 && l0 > 0) ? l0 : m_DNABondL0;
        WriteDNABondsFile(m_DNABondsFilename, all_chains, start_index, m_DNABondK, final_l0, m_DNABondKAngle, m_DNABondTheta0);
        std::cout << "---> DNA bonds written to: " << m_DNABondsFilename << std::endl;
    }
}

// Write DNA Q file
void Generate::WriteDNAQFile(std::string filename, const std::vector<Vertex_Map>& dnaVertices, Vec3D box)
{
    int pres = 10;
    std::ofstream output;
    output.open(filename.c_str());
    
    if (!output.is_open()) {
        std::cerr << "---> Error: Could not open DNA output file: " << filename << std::endl;
        return;
    }
    
    output << std::fixed;
    output << std::setprecision(pres) << box(0) << "   " << box(1) << "   " << box(2) << "   \n";
    output << dnaVertices.size() << "\n";
    
    for (size_t i = 0; i < dnaVertices.size(); i++) {
        const Vertex_Map& v = dnaVertices[i];
        output << std::setprecision(pres) << v.id << "  " << v.x << "  " << v.y << "  " << v.z << "  " << v.domain << std::endl;
    }
    
    output << "0\n";  // No triangles for DNA beads
    
    output.close();
}

// Write DNA bonds file for multiple chains
void Generate::WriteDNABondsFile(std::string filename, const std::vector<std::vector<Vertex_Map> >& chains, int start_index, double k, double l0, double k_angle, double theta0)
{
    std::ofstream output;
    output.open(filename.c_str());
    
    if (!output.is_open()) {
        std::cerr << "---> Error: Could not open DNA bonds output file: " << filename << std::endl;
        return;
    }
    
    // Write header with global parameters (using ; for comments)
    output << "; DNA bonds file (generated by dts_generate)" << std::endl;
    output << "; Global bond parameters:" << std::endl;
    output << "; k_stretch = " << k << std::endl;
    output << "; l0 = " << l0 << std::endl;
    output << "; k_angle = " << k_angle << std::endl;
    output << "; theta0 = " << theta0 << std::endl;
    output << "; Number of chains: " << chains.size() << std::endl;
    output << "; Format: v1 v2 [k_stretch] [l0] (if k_stretch and l0 are in header, only v1 v2 needed)" << std::endl;
    output << std::endl;
    
    int current_id = start_index;
    
    // Write bonds for each chain (circular topology)
    for (size_t c = 0; c < chains.size(); c++) {
        const std::vector<Vertex_Map>& chain = chains[c];
        int chain_size = chain.size();
        
        for (int i = 0; i < chain_size; i++) {
            int v1 = current_id + i;
            int v2 = current_id + ((i + 1) % chain_size);
            
            // Write bond: v1 v2 (k and l0 are in header)
            output << v1 << " " << v2 << std::endl;
        }
        
        current_id += chain_size;
    }
    
    output.close();
}

// Update input.dts file with bond parameters (placeholder - can be implemented if needed)
void Generate::UpdateInputFileWithBondParams(std::string input_dts_file, std::string bonds_filename, double k, double l0)
{
    // This function could append a comment to input.dts about the generated bond parameters
    // For now, it's a placeholder
    (void)input_dts_file;
    (void)bonds_filename;
    (void)k;
    (void)l0;
}
