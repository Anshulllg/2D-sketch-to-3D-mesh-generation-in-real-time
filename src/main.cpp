#include "utils.h"
#include "math.h"
#include "poly2tri.h"
#include <set>
#include <iostream>
#include <math.h>
#include "half_edge.h"
#include <map>
#include<bits/stdc++.h> 

using namespace std;
using namespace glm;
using namespace p2t;

// GLobal variables

vector<float> c_points;
vector<float> verticesToDraw;
int width = 1000, height = 1000; 
bool c_pointsUpdated = false;
vector<Point*> points;
vector<Triangle*> triangles;
vector<float> triarry;
int tarr[10000];
vector<vertex *> vertices;
vector<face *> faces;
map<pair<int ,int>,struct halfedge *> hashedge;
CDT* cdt;
bool flag=false;
int checkk=0;
unsigned int cp_VBO;
unsigned int p2d_VBO;
unsigned int triangle_VBO;

unsigned int cp_VAO;
unsigned int p2d_VAO;
unsigned int traingle_VAO;

int board_w, board_h;

int trisize= triarry.size();
float translatee[] = {0.0,0.1};

void func1(Point *t) {
    triarry.push_back(t->x);
    triarry.push_back(t->y);
    triarry.push_back(0);
    
}

void tribuff_add(){
    for (auto triangle: triangles){
        for (int i =0;i<3;i++){
            func1(triangle->GetPoint(i));
        }
        for(int k=0; k<3; k++){
            func1(triangle->GetPoint(k));
        }
    }
}

void func(edge *e,int flagg) {
    triarry.push_back(e->v->x);
    triarry.push_back(e->v->y);
    if (flagg==0){
         triarry.push_back(e->v->z);
    }else{
        triarry.push_back(-1*e->v->z);
    } 
}

void facebuff(){
    triarry.clear();
    
    for (int i=0; i<faces.size();i++){

        func(faces[i]->e,0);
        func(faces[i]->e->next,0);
        func(faces[i]->e->next->next,0);

        func(faces[i]->e,0);
        func(faces[i]->e->next->next,0);
        func(faces[i]->e->next,0);
    }

    for (int i=0; i<faces.size();i++){
        func(faces[i]->e,1);
        func(faces[i]->e->next,1);
        func(faces[i]->e->next->next,1);

        func(faces[i]->e,1);
        func(faces[i]->e->next->next,1);
        func(faces[i]->e->next,1);
    }

}


void erase_lines(vector<float> &points, vector<Point*>& p2tPoints,vector<Triangle*>& triangles,vector<float>& triarry,bool &mouseDowned){

    c_points.clear();
    points.clear();
    triangles.clear();
    triarry.clear();
    c_pointsUpdated = false;
    mouseDowned = false;
    
}

int main(int, char* argv[])
{
    int mouseDowned = 0;
    bool artboard=true;
    bool check =false;
    int check2= 0;
    GLFWwindow* window = setupWindow(width, height);
    ImGuiIO& io = ImGui::GetIO(); // Create IO object

    // Create VBOs, VAOs  
    glGenBuffers(1, &cp_VBO);
    glGenVertexArrays(1, &cp_VAO);
    glGenBuffers(1, &p2d_VBO);
    glGenVertexArrays(1, &p2d_VAO);
    glGenBuffers(1, &triangle_VBO);
    glGenVertexArrays(1, &traingle_VAO);

    mat4 view = lookAt(vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 0.0), vec3(0.0, 1.0, 0.0));
    unsigned int shaderProgram = createProgram("./shaders/vshader.vs", "./shaders/fshader.fs");
    glUseProgram(shaderProgram);


    if(glGetUniformLocation(shaderProgram, "vView") == -1){
        fprintf(stderr, "Could not bind location: vView\n");
        exit(0);
    }
    glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "vView"), 1, GL_FALSE, glm::value_ptr(view));    



    while (!glfwWindowShouldClose(window))
    {
        glfwPollEvents();
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        ImGui::SetNextWindowPos(ImVec2(0, 0));
            ImVec4 clear_color = ImVec4(1.0f, 1.0f, 1.0f,1.0f);
            ImGui::Begin("Toolbar", NULL, ImGuiWindowFlags_AlwaysAutoResize);

            if (ImGui::SliderFloat2("translation", translatee, -1.0, 1.0)){
                std::cout<<"translation:"<<translatee[0]<<std::endl;

            }
            
            if(ImGui::Button("Erase")){
                erase_lines(c_points,points,triangles,triarry,artboard);
        
            }
            // ImVec4 clear_color = ImColor(114, 144, 154);
            ImGui::ColorEdit3("clear color", (float*)&clear_color);
            ImGui::End();

        // Rendering
        trans(translatee, check2);

       
        view = lookAt(normalize(vec3(0.0, translatee[0],translatee[1])), vec3(0.0, 0.0, 0.0), vec3(0.0, 1.0, 0.0));
        glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "vView"), 1, GL_FALSE, value_ptr(view));

        ImGui::Render();
        // Add a new point on mouse click
        float x,y ;
        ImDrawData * draw= ImGui::GetDrawData();
        if (!artboard){
            artboard=true;
            mouseDowned=0;
            checkk=1;
            ImGui_ImplOpenGL3_RenderDrawData(draw);
            glfwSwapBuffers(window);
        }
        
        
        glfwGetFramebufferSize(window, &board_w, &board_h);
        glViewport(0, 0, board_w, board_h);
        glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
        glClear(GL_COLOR_BUFFER_BIT);

        if (io.MouseDown[0]){
            x = io.MousePos.x;
            y = io.MousePos.y;
            
            c_points.push_back(scaling_X(x,y,height,width));
            c_points.push_back(scaling_Y(x,y,height,width));
            c_points.push_back(0.0);
            if (mouseDowned){
                
                c_points.push_back(scaling_X(x,y,height,width));
                c_points.push_back(scaling_Y(x,y,height,width));
                c_points.push_back(0.0);
            }
            ptpush(x,y,height,width);
            mouseDowned = 1;
            c_pointsUpdated = true;
        }

        if (io.MouseReleased[0]   ){
            x = io.MousePos.x;
            y = io.MousePos.y;
            
            c_points.push_back(scaling_X(x,y,height,width));
            c_points.push_back(scaling_Y(x,y,height,width));
            c_points.push_back(0.0);

            ptpush(x,y,height,width);
            
            cdt = new CDT(points);
            cdt->Triangulate();
    		triangles = cdt->GetTriangles();
            // adding too triangle buffer
            triarry.clear();    //helloo1
            tribuff_add();
            // splitting faces in half  at the codial axis
            half_edgebuff(points,triangles,vertices,faces);
            marking(faces);
            faces= prunningg(vertices, faces);
            c_pointsUpdated = true;
            mouseDowned = 0;
        }
        ///hellooooooooo
        else if ((io.KeyCtrl || io.MouseReleased[1]) ){
            if(io.MouseReleased[1])
                faces= elevation_sewing(vertices, faces);
            facebuff();
            mouseDowned = 0;
            c_pointsUpdated = true;
            
        }

        else if (io.KeyAlt){
            // verticesToDraw.clear();
            c_pointsUpdated=true;
        }

        if(c_pointsUpdated) {
            flag=true;
 // drawww
            if (io.MouseDown[0]){
                glBindVertexArray(cp_VAO);
                glBindBuffer(GL_ARRAY_BUFFER, cp_VBO);
          
                glBufferData(GL_ARRAY_BUFFER, c_points.size()*sizeof(GLfloat), &c_points[0], GL_DYNAMIC_DRAW);
                glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float),0);
                glEnableVertexAttribArray(0); 
                glBindVertexArray(cp_VAO);
                glDrawArrays(GL_LINES, 0, c_points.size()/3);
                glUseProgram(0);

                ImGui_ImplOpenGL3_RenderDrawData(draw);
                glfwSwapBuffers(window);


            }
 // triangulation   1st click, 2nd swee, 3d elevate
            if (( io.MouseReleased[1] ||io.KeyAlt || io.MouseReleased[0] || io.KeyCtrl ) ){
                if(io.KeyAlt)
                    check2=1;

                glBindBuffer(GL_ARRAY_BUFFER, traingle_VAO);
                glBindVertexArray(traingle_VAO);
                glBufferData(GL_ARRAY_BUFFER, triarry.size()*sizeof(GLfloat), &triarry[0], GL_DYNAMIC_DRAW);
                glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3*sizeof(float), 0);
                glEnableVertexAttribArray(0);
                glBindVertexArray(traingle_VAO);
                glDrawArrays(GL_LINES, 0, triarry.size()/3);
                glUseProgram(0);

                points.clear();
                ImGui_ImplOpenGL3_RenderDrawData(draw);
                glfwSwapBuffers(window);
                c_pointsUpdated = false;
            }
        }
        if(flag==false){
            ImGui_ImplOpenGL3_RenderDrawData(draw);
            glfwSwapBuffers(window);
        }
        glUseProgram(shaderProgram);
    }
    
    glDeleteBuffers(1, &triangle_VBO);
    glDeleteBuffers(1,&cp_VBO);

    return 0;

    
}
