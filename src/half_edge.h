
#include "utils.h"
#include <iostream>
#include <math.h>
#include <map>
#include<bits/stdc++.h> 
#include<algorithm>
using namespace std;
using namespace glm;
using namespace p2t;

extern map<pair<int ,int>,struct halfedge *> hashedge;
extern vector<Point*> points;


// Half Edge Struct
typedef struct halfedge{
    struct halfvertex *v;
    struct halfface *f;
    struct halfedge *prev_prev;
    struct halfedge *next;
    struct halfedge *nex_next;
    struct halfedge *opp;
    struct halfedge *prev;
    
    
} edge;

// Half Face Struct
typedef struct halfface{
    edge *e;
    int triangleType; //2 for terminal, 1 for sleeve and 0 for junction
    int visit;
} face;

// Half Vertex Struct
typedef struct halfvertex{
    float x;
    float y;
    float z;
    struct halfedge *e;
    int vNum;
    int alpha;
    int samplePoints;
    bool boundary;
    
} vertex;



static void deleteFaces(vector<face*> &faces) {

    int totalFaces = faces.size();

    for(auto f:faces){
        edge *del_edge = new edge;
        del_edge = f->e->next->next;
        delete del_edge;
        del_edge = f->e->next;
        delete del_edge;
        del_edge = f->e;
        delete f->e;
        delete f;
    }
    faces.clear();
}
  

// Performing the Cross product using the 2  directions 
static GLfloat crossp(const vector<vertex *> &vertices,int aa, int bb,int  cc){

    float p1x = vertices.at(cc)->x-vertices.at(aa)->x;
    float p1y = vertices.at(cc)->y-vertices.at(aa)->y;
    vec3 direction_2(p1x, p1y,0.0);
    float p2x = vertices.at(bb)->x-vertices.at(aa)->x;
    float p2y = vertices.at(bb)->y-vertices.at(aa)->y;
    vec3 direction_1(p2x, p2y, 0.0);
    GLfloat cross_prod_val = glm::cross(direction_1, direction_2).z;
    return cross_prod_val;
}


// Creating the Half Edge Buffer of face
static void half_edgeface(int i1,int i2,int i3, const vector<vertex *> &vertices, vector<face *> &faces){

    face *f = new face; 
    edge *e1 = new edge;
    edge *e2 = new edge;
    edge *e3 = new edge;
    f->visit = 0;
    f->triangleType = -1;
    f->e = e1;

    e1->f = f;
    e1->opp = NULL;
    e1->next = e2;
    e1->prev = e3;

    e2->f = f;
    e2->opp = NULL;
    e2->next = e3;
    e2->prev = e3;

    e3->f = f;
    e3->opp = NULL;
    e3->next = e1;
    e3->prev = e2;

    edge *edg[3];

    edg[0] = e1;
    edg[1] = e2;
    edg[2] = e3;

    for (int i =0;i<3;i++){
        edg[i]->f = f;
        edg[i]->opp = NULL;
        if(i==2){ edg[i]->prev = edg[i-1];}
        else{ edg[i]->prev = edg[2];}
        edg[i]->next = edg[(i+1)%3];
    }
   
    if (crossp(vertices,i1, i2, i3)<0.0){

        int tt = i1;
        i1 = i3;
        i3 = tt;
    }

    int indices[] = {i1,i2,i3};
    pair<int,int> pair;

    e1->v = vertices[i1];
    vertices[i1]->e = e1;
    pair.first = std::min(i1,i2);
    pair.second = std::max(i1,i2);
    if(hashedge.find(pair)==hashedge.end()) hashedge[pair] = e1;
    else{
        hashedge[pair]->opp = e1;
        e1->opp = hashedge[pair];
    } 

    e2->v = vertices[i2];
    vertices[i2]->e = e2;
    pair.first =std:: min(i2,i3);
    pair.second =std:: max(i2,i3);
    if(hashedge.find(pair)==hashedge.end()) hashedge[pair] = e2;
    else{
        hashedge[pair]->opp = e2;
        e2->opp = hashedge[pair];
    } 

    e3->v = vertices[i3];
    vertices[i3]->e = e3;
    pair.first = std::min(i1,i3);
    pair.second = std::max(i1,i3);
    if(hashedge.find(pair)==hashedge.end()) hashedge[pair] = e3;
    else{
        hashedge[pair]->opp = e3;
        e3->opp = hashedge[pair];
    } 

    faces.push_back(f);

}



static float cal ( float x, float y, Point *points ){
    double x_pts =  pow( x*x - points->x*points->x , 2 );
    double y_pts = pow( y*y - points->y*points->y , 2 );
    float ans = sqrt(x_pts+y_pts);
    return ans;
}



static float cal2 ( float x, float y, vertex *points ){
    double x_pts =  pow( x*x - points->x*points->x , 2 );
    double y_pts = pow( y*y - points->y*points->y , 2 );
    float ans = sqrt(x_pts+y_pts);
    return ans;
        
}


static void half_edgebuff(vector<Point*> &points, vector<Triangle*> &triangles,vector<vertex *> &vertices, vector<face *> &faces){

    vector<int> indices;
    
    float d_val;
    int f_idx;
    int idx;
    float x,y,z;
    for(int i=0;i<points.size();i++){
        vertex *vv = new vertex;
        vv->e=NULL;
        vv->vNum = vertices.size();
        vv->boundary = true;
        vv->x = points[i]->x;
        vv->y = points[i]->y;
        vv->z = 0;
        vertices.push_back(vv);
    }


    for (int j =0; j<triangles.size();j++){
        for (int i =0;i<3;i++){
            f_idx = 0;
            idx= 0;
            x = triangles[j]->GetPoint(i)->x;
            y = triangles[j]->GetPoint(i)->y;

            for (int k=0; k<points.size(); k++){
               d_val= cal(x,y,points[k]);
                if (d_val<=0.00001){
                    f_idx = 1;
                    break;
                }
                idx++;
            }
            indices.push_back(idx);
        }
        half_edgeface(indices[0],indices[1],indices[2],vertices,faces);
        indices.clear();
    }

}
static void marking(vector<face *>& faces){
    
    int c = -1;
    edge *e;
    for(int i=0;i<faces.size();i++){
        if(c==-1) {c += 1;}
        else {c = 0;}
        if(faces[i]->e->opp==NULL) c++;
        e = faces[i]->e->next;
        while(e!=faces[i]->e){
            if(e->opp==NULL) c++;
            e = e->next;
        }
        faces[i]->triangleType=c;
    }
}

// This is a circle checking for the point lying outside the circle
static bool checkout(vertex *v3,float cntr_x,float cntr_y, float rad){
    double center_x =  pow(v3->x-cntr_x,2);
    double center_y = pow(v3->y-cntr_y,2);
    double radius = pow(rad,2);
    
    return center_x+center_y>radius;
}

// This is getting the vertex of each and every index that is taken from the mouseinput
static int vert_idx(float x, float y, vector<vertex *> *vertices){

    int idx=0;
    int f_idx=0;
    int centre;
    float d_val;
    vertex *my_v = new vertex;
    for (auto point: *vertices){
        d_val = sqrt(pow( pow(x,2)  - pow(point->x,2) , 2 )+pow( pow(y,2)  - pow(point->y,2) , 2 ));
        if (d_val<=0.00001){
            f_idx = true;
            break;
        }
        idx++;
    }
    if (!f_idx){
        centre = vertices->size();
            
            my_v->e=NULL;
            my_v->vNum = centre;
            my_v->boundary = false;
            my_v->x = x;
            my_v->y = y;
            my_v->z = 0.0;
            vertices->push_back(my_v);

    }
    else{
        centre = idx;
    }

    return centre;

}
static int centerr(float a, float b, vector<vertex *> *vertt){
    return vert_idx(a,b,vertt);
}


static void addingFaces(int i1,int i2,int i3, const vector<vertex *> &vertices, vector<face *> &prunedFaces){
    
    face *f = new face; 
    edge *e1 = new edge;
    edge *e2 = new edge;
    edge *e3 = new edge;
    f->visit = 0;
    f->triangleType = -1;
    f->e = e1;

   
    edge *edg[3];
    edg[0] = e1;
    edg[1] = e2;
    edg[2] = e3;

    for (int i =0;i<3;i++){
        edg[i]->f = f;
        edg[i]->opp = NULL;

        if(i==2){ edg[i]->prev = edg[i-1];}
        else{ edg[i]->prev = edg[2];}
        edg[i]->next = edg[(i+1)%3];
        
    }
   

    if (crossp(vertices,i1,i2,i3)){
        int t = i1;
        i1 = i3;
        i3 = t;
        
    }
    int indices[] = {i1,i2,i3};

    for (int i=0;i<3;i++){
        edg[i]->v = vertices[indices[i]];
    }
    prunedFaces.push_back(f);
}


static vector<face *> prunningg(vector<vertex *> &vertices, vector<face *> &faces){
    edge *prevE;
    edge *e;
    edge *opp;
    int idx;

    float fanp1;
    float fanp2;
    float fanp3;

    float cntr_x;
    float cntr_y;
    float rad;

    float d;
    int centre;
    bool outside=false;
    int f_idx =0;
    int cnt;
    
    vertex *vert[3];
    vector<face *> pruned_faces; 
    deque<vertex *> uneccVertices;
    for (int i =0; i<faces.size();i++){
        if (faces[i]->triangleType==2){ 
            prevE = faces[i]->e; 
            e = faces[i]->e;
            uneccVertices.clear();
            outside=false;
            while (!outside){

                e = e->next;
                cnt =0;
                if(e->opp==NULL){
                    while (e->opp==NULL && cnt<=3){
                    e = e->next;
                    cnt++;
                    }
                }
                

                faces[i]->visit=1;
                e->f->visit=1;

                vert[0]= e->v;
                vert[1]= e->next->v;
                vert[2]= e->next->next->v;
                if (e->f->triangleType ==0){
            
                    fanp1 = (vert[0]->x+vert[1]->x + vert[2]->x)/3;
                    fanp2 = (vert[0]->y+vert[1]->y + vert[2]->y)/3;
                    fanp3 = 0.0;
                    centre = vert_idx(fanp1,fanp2,&vertices);
                    break;
                }

                int c0=0;
                for ( int i =0; i<uneccVertices.size(); i++){
                    if ( uneccVertices[i]==vert[2]){
                        c0++;
                    }
                }
                if( c0==0){
                    uneccVertices.push_front(vert[2]);
                }
            
                cntr_x = (vert[0]->x+vert[1]->x)/2.0;
                cntr_y = (vert[0]->y+vert[1]->y)/2.0;
                rad = sqrt(pow(vert[0]->x-vert[1]->x,2) + pow(vert[0]->y-vert[1]->y,2))/2.0;
                for (int i=1;i<=uneccVertices.size()-1;i++){
                    vertex * v = uneccVertices[i];
                    if (checkout(v,cntr_x,cntr_y,rad)){
                        centre = vertices.size();
                        // vertices.push_back(makeHalfEdgeVertex(cntr_x,cntr_y,0.0,vertices.size(),false));
                        vertex *v = new vertex;
                        v->x = cntr_x; 
                        v->y = cntr_y; 
                        v->z = 0.0;
                        v->e = NULL;  
                        v->vNum = vertices.size();
                        v->boundary=false;
                        // static vertex *v = new vertex;
                        vertices.push_back(v);
                        outside = true;
                        break;
                    }
                }
                
//  for vert[1]
                int c1=0;
                for ( int i =0; i<uneccVertices.size(); i++){
                    if ( uneccVertices[i]==vert[1]){
                        c1++;
                    }
                }
                if( c1==0){
                    uneccVertices.push_front(vert[1]);
                }
//  for vert[0]
                int c2=0;
                for ( int i =0; i<uneccVertices.size(); i++){
                    if ( uneccVertices[i]==vert[0]){
                        c2++;
                    }
                }
                if( c2==0){
                    uneccVertices.push_back(vert[0]);
                }
                e = e->opp;  
            }

            if (uneccVertices.size() !=0){
                for (int i=0;i<uneccVertices.size()-1;i++){
                    int num1=uneccVertices[i]->vNum;
                    int num2=uneccVertices[i+1]->vNum;
                    half_edgeface(num1,num2,centre,vertices, pruned_faces);
                }
            }
        }   
    }
    int j =0;
    while(j<faces.size()){
        if (faces[j]->triangleType==0){
            int centroid;            
            e = faces[j]->e;
            for (int i =0;i<3;i++){
                fanp1 = (e->v->x + e->next->v->x)/2.0;
                fanp2 = (e->v->y + e->next->v->y)/2.0;
                fanp3 = 0.0;
                f_idx = 0;
                vert[0] = e->v;
                vert[1] = e->next->v;
                vert[2] = e->next->next->v;
                for (int i=0;i<vertices.size();i++){
                    double x_pts =  pow( fanp1*fanp1  - vertices[i]->x*vertices[i]->x , 2 );
                    double y_pts = pow( fanp2*fanp2  - vertices[i]->y*vertices[i]->y , 2 );
                    d = sqrt(x_pts+y_pts);
                    if (d<=0.00001){
                        f_idx = 1;
                        break;
                    }
                }
                if(e->opp!=NULL){
                    if((e->opp->f->triangleType==0 || e->opp->f->visit==0 || f_idx)){
                        fanp1 = (vert[0]->x + vert[1]->x + vert[2]->x)/3.0;
                        fanp2 = (vert[0]->y + vert[1]->y + vert[2]->y)/3.0;
                        fanp3 = 0.0;

                        float aa=(e->v->x+e->next->v->x)/2.0;
                        float bb=(e->v->y+e->next->v->y)/2.0;

                        half_edgeface(e->v->vNum,centerr(aa,bb,&vertices),centerr(fanp1,fanp2,&vertices),vertices,pruned_faces);
                        half_edgeface(centerr(aa,bb,&vertices),e->next->v->vNum,centerr(fanp1,fanp2,&vertices),vertices,pruned_faces);
                    }
                    
                }
                e = e->next;
            }
        }

        else if (faces[j]->visit==0){
            e = faces[j]->e;
            cnt =0;
            e = e->next;
            
            
            while (e->opp!=NULL && cnt<=5 ){
            e = e->next;
            cnt +=1;
            }
            

            faces[j]->visit=1;
            e->f->visit=1;
            vert[0] = e->v;
            vert[1] = e->next->v;
            vert[2] = e->next->next->v;
            
            float aa=(vert[1]->x+vert[2]->x)/2.0;
            float bb=(vert[1]->y+vert[2]->y)/2.0;
            float cc=(vert[0]->x+vert[2]->x)/2.0;
            float dd=(vert[0]->y+vert[2]->y)/2.0;
        
            half_edgeface(vert[0]->vNum,vert[1]->vNum,centerr(aa,bb,&vertices),vertices, pruned_faces);
            half_edgeface(vert[0]->vNum,centerr(aa,bb,&vertices),centerr(cc,dd,&vertices),vertices, pruned_faces);
            half_edgeface(centerr(aa,bb,&vertices),vert[2]->vNum,centerr(cc,dd,&vertices), vertices, pruned_faces);
        }
        j++;
    }

    for(int i =0; i< faces.size(); i++){
        edge *del_edge = new edge;
        del_edge = faces[i]->e->next->next;
        delete del_edge;
        del_edge = faces[i]->e->next;
        delete del_edge;
        del_edge = faces[i]->e;
        delete faces[i]->e;
        delete faces[i];
    }
    faces.clear();

    return pruned_faces;
}


static void l_elevate(edge **hlfedg, float &length, int &n){
    // float dir[2];
    float direction1;
    float direction2;

    bool condition = (*hlfedg)->next->v->boundary;

    if(condition==true){
        // Reversing it check it please
        direction1=(*hlfedg)->v->x -(*hlfedg)->next->v->x ;
        direction1 = pow(direction1,2);
        direction2=  (*hlfedg)->v->y - (*hlfedg)->next->v->y ;
        direction2 = pow(direction2,2);
        length += sqrt(direction1+direction2);
        // len+=sqrt(dir[0]*dir[0]+dir[1]*dir[1]);
        n+=1;
    }
    (*hlfedg)=(*hlfedg)->opp->next;
}



// static void ellipsoid_vert(bool hEdge->v->boundary,bool hEdge->next->v->boundary, vertex* vert[] , int samplePoints, vector<vertex *> n_points, map<pair<int, int>, vector<int>> samplePtsPerEdge, edge* hEdge, vector<int> sampleIndices){
//     if( (hEdge->v->boundary && (!hEdge->next->v->boundary)) || ((!hEdge->v->boundary) && hEdge->next->v->boundary)){
//         if(!hEdge->v->boundary){
//             vert[0]=(hEdge)->v;
//             vert[1]=(hEdge)->next->v;
//         }
//         else{
//             vert[1]=(hEdge)->v;
//             vert[0]=(hEdge)->next->v;
//         }

//         float x_cord = vert[1]->x-vert[0]->x;
//         float y_cord = vert[1]->y-vert[0]->y;
//         vec3 ellipsoid_axis( x_cord , y_cord , 0.0 );
//         // float a= ellipsoid_axis.length();


//         // float b= vert[0]->z;
//         float t = 1.0;
//         for(int i=1;i<samplePoints+1;i++){
//             t=((float)i/(float)(samplePoints+1));
//             if(t<1){
//                 vertex *v = new vertex;
//                 v->x = vert[0]->x + t*ellipsoid_axis[0];
//                 v->y = vert[0]->y + t*ellipsoid_axis[1];
//                 v->z = vert[0]->z*sqrt(1-pow(t,2));
//                 v->vNum = n_points.size();
//                 n_points.push_back(v);
//                 sampleIndices.push_back(n_points.size()-1);   
//             }
             
//         }

//         pair<int,int> p;
//         p.first = std::min(vert[0]->vNum, vert[1]->vNum);
//         p.second = std::max(vert[0]->vNum, vert[1]->vNum);
//         samplePtsPerEdge[p]=sampleIndices;
//         }
//         hEdge=hEdge->next;
// }

// static void update_n_points()

static void addingg(vector<int> e1SamplePts,vector<int> e2SamplePts, int pointts,int pointts2, std::vector<vertex *> n_points, std::vector<face *> elevated_f){
    addingFaces(e1SamplePts[pointts+1], e2SamplePts[pointts+1], e1SamplePts[pointts+2], n_points, elevated_f);
    addingFaces(e2SamplePts[pointts+1], e2SamplePts[pointts+2], e1SamplePts[pointts+2], n_points, elevated_f);
    pointts2++;
}

static vector<face *> elevation_sewing(vector<vertex *> &vertices, vector<face *> &faces){
    
    int samplePoints=0;

    vector<face *> elevated_f; 
    edge *hEdge;
    
    
    for(int i =0 ; i < vertices.size(); i++){
        float length=0.0;
        int count = 0;
        int n=0;

        if(!vertices[i]->boundary){
            if(vertices[i]->e->opp!=NULL){
               
                hEdge=vertices[i]->e;
                // hellliuuuu
                bool x = hEdge==vertices[i]->e;

                l_elevate(&hEdge, length, n);
                while(hEdge->opp!=NULL && hEdge!=vertices[i]->e){
                    l_elevate(&hEdge, length, n);
                    count ++;
                }
                vertices[i]->z=length/n;
                samplePoints=int(length/n/0.15);
            }
            
        }
    }

    vector<vertex *> n_points(vertices);

    std::map<std::pair<int, int>, std::vector<int>> samplePtsPerEdge;
    // bool hEdge->v->boundary = hEdge->v->boundary;
    // bool hEdge->next->v->boundary = hEdge->next->v->boundary; 

    for(int i =0 ; i< faces.size(); i++){
        hEdge=faces[i]->e;
        vertex *vert[2];
        int count=0;
        vector<int> sampleIndices;
        
        // ellipsoid_vert(hEdge->v->boundary, hEdge->next->v->boundary, vert, samplePoints, n_points,samplePtsPerEdge, hEdge, sampleIndices);
        if( (hEdge->v->boundary && (!hEdge->next->v->boundary)) || ((!hEdge->v->boundary) && hEdge->next->v->boundary)){
        if(!hEdge->v->boundary){
            vert[0]=(hEdge)->v;
            vert[1]=(hEdge)->next->v;
        }
        else{
            vert[1]=(hEdge)->v;
            vert[0]=(hEdge)->next->v;
        }

        float x_cord = vert[1]->x-vert[0]->x;
        float y_cord = vert[1]->y-vert[0]->y;
        vec3 ellipsoid_axis( x_cord , y_cord , 0.0 );
        // float a= ellipsoid_axis.length();


        // float b= vert[0]->z;
        float t = 1.0;
        float tt=0.0;
        for(int i=1;i<samplePoints+1;i++){
            float samp_pt= samplePoints+1;
            t=((float)i/(float)samp_pt);
            if(t<1){
                vertex *v = new vertex;
                v->x = vert[0]->x + t*ellipsoid_axis[0];
                v->y = vert[0]->y + t*ellipsoid_axis[1];
                v->z = vert[0]->z*sqrt(1-pow(t,2));
                v->vNum = n_points.size();
                n_points.push_back(v);
                sampleIndices.push_back(n_points.size()-1);   
            }
             
        }

        pair<int,int> p;
        p.first = std::min(vert[0]->vNum, vert[1]->vNum);
        p.second = std::max(vert[0]->vNum, vert[1]->vNum);
        samplePtsPerEdge[p]=sampleIndices;
        }
        hEdge=hEdge->next;

        while(hEdge!=faces[i]->e && count++<3){
            sampleIndices.clear();
            // ellipsoid_vert(hEdge->v->boundary, hEdge->next->v->boundary, vert, samplePoints, n_points,samplePtsPerEdge, hEdge, sampleIndices);
            if( (hEdge->v->boundary && (!hEdge->next->v->boundary)) || ((!hEdge->v->boundary) && hEdge->next->v->boundary)){
        if(!hEdge->v->boundary){
            vert[0]=(hEdge)->v;
            vert[1]=(hEdge)->next->v;
        }
        else{
            vert[1]=(hEdge)->v;
            vert[0]=(hEdge)->next->v;
        }

        float x_cord = vert[1]->x-vert[0]->x;
        float y_cord = vert[1]->y-vert[0]->y;
        vec3 ellipsoid_axis( x_cord , y_cord , 0.0 );
        // float a= ellipsoid_axis.length();


        // float b= vert[0]->z;
        float t = 1.0;
        float bb=0.0;
        for(int i=1;i<samplePoints+1;i++){
            float samp_pt= samplePoints+1;
            t=((float)i/(float)samp_pt);
            if(t<1){
                vertex *v = new vertex;
                v->x = vert[0]->x + t*ellipsoid_axis[0];
                v->y = vert[0]->y + t*ellipsoid_axis[1];
                v->z = vert[0]->z*sqrt(1-pow(t,2));
                v->vNum = n_points.size();
                n_points.push_back(v);
                sampleIndices.push_back(n_points.size()-1);   
            }
             
        }

        pair<int,int> p;
        p.first = std::min(vert[0]->vNum, vert[1]->vNum);
        p.second = std::max(vert[0]->vNum, vert[1]->vNum);
        samplePtsPerEdge[p]=sampleIndices;
        }
        hEdge=hEdge->next;
        }
    }

    int cnt=0;
    vector<int> edge1SamplePts;
    vector<int> edge2SamplePts;

    for(int i=0; i<faces.size(); i++){
        hEdge=faces[i]->e;
        
        cnt=0;
        hEdge=hEdge->next;

        int counter = 0;
      
        while (hEdge != faces[i]->e && counter < 3) {
            bool condd=hEdge->v->boundary;
            bool condd2=hEdge->next->v->boundary;
            if ((condd && condd2) ||
                (!condd && !condd2)) {
                break;
            }
            hEdge = hEdge->next;
            counter++;
        }
        std::pair<int,int> p;
        if(hEdge->next->v->vNum < hEdge->next->next->v->vNum){
            p.first =hEdge->next->v->vNum;
            p.second =hEdge->next->next->v->vNum;
        }
        else if (hEdge->next->v->vNum > hEdge->next->next->v->vNum){
            p.second =hEdge->next->v->vNum;
            p.first =hEdge->next->next->v->vNum;
        }
        // std::vector<int> edge_at_p=

        vector<int> e1SamplePts(samplePtsPerEdge.at(p));
        
        // std::vector<int> edge_at_p=samplePtsPerEdge.at(p);
        // std::vector<int> e1SamplePts(edge_at_p);
        pair<int,int> p1;

        if(hEdge->next->next->v->vNum  < hEdge->next->next->next->v->vNum){
            p1.first =hEdge->next->next->v->vNum ;
            p1.second =hEdge->next->next->next->v->vNum;
        }
        else if (hEdge->next->next->v->vNum > hEdge->next->next->next->v->vNum){
            p1.second =hEdge->next->next->v->vNum;
            p1.first =hEdge->next->next->next->v->vNum;
        }

        vector<int> e2SamplePts(samplePtsPerEdge.at(p1));

        bool cond1=hEdge->v->boundary;
        bool cond2= hEdge->next->v->boundary;
        if(cond1 && cond2){

            int j;
            for ( int i = 0,j = e1SamplePts.size() - 1; i < e1SamplePts.size()/2; i++, j--)  
            {     
                int temp = e1SamplePts[i];  
                e1SamplePts[i] = e1SamplePts[j];  
                e1SamplePts[j] = temp;  
            }
            for ( int i = 0,j = e2SamplePts.size() - 1; i < e2SamplePts.size()/2; i++, j--)  
            {     
                int temp = e2SamplePts[i];  
                e2SamplePts[i] = e2SamplePts[j];  
                e2SamplePts[j] = temp;  
            }
          
        }
       
        addingFaces(hEdge->v->vNum, hEdge->next->v->vNum, e1SamplePts[0], n_points, elevated_f);
        int pointts=0;
        int pointts2=0;

        addingFaces(e2SamplePts[0], e1SamplePts[0], hEdge->next->v->vNum,n_points, elevated_f);
        

     while(pointts<samplePoints-1)
        addingg(e1SamplePts,e2SamplePts,pointts, pointts2, n_points, elevated_f);
        addingFaces(e1SamplePts[pointts], hEdge->next->next->v->vNum, e2SamplePts[pointts], n_points, elevated_f);
    }
   
    deleteFaces(faces);
    // sleep(2);
    return elevated_f;

}
