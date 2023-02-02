// This is the starter code for the CS 6491 Ray Tracing project.
//
// The most important part of this code is the interpreter, which will
// help you parse the scene description (.cli) files.
import java.util.*;
public class Color{
  float r;
  float g;
  float b;
}
public class Surface{
  float dr;
  float dg;
  float db;
}
public class Light{
  float x;
  float y;
  float z;
  float r;
  float g;
  float b;
}
public class Triangle{
  PVector vec1;
  PVector vec2;
  PVector vec3;
  Surface sur;
  Triangle(PVector a, PVector b, PVector c){
    vec1 = a;
    vec2 = b;
    vec3 = c;
  }
}

public class Ray{
  PVector origin;
  PVector slope;
}

public class matrix{
  public float[][] data;
  public matrix() {
    data = new float[4][4];
  }
  
  public matrix plus(matrix B){
    matrix A = this;    
    matrix C = new matrix();
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            C.data[i][j] = A.data[i][j] + B.data[i][j];
    return C;
  }

  public matrix times(matrix B){
    matrix A = this;      
    matrix C = new matrix();
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            for (int k = 0; k < 4; k++)
                C.data[i][j] += (A.data[i][k] * B.data[k][j]);
    return C;
  }
  public float[][] timesv(float B[][]){
    matrix A = this;      
    float C[][]={{0},{0},{0},{0}};
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 1; j++)
            for (int k = 0; k < 4; k++)
                C[i][j] += (A.data[i][k] * B[k][j]);
    return C;
  }

  

  public void show(){
    for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++)
        println(data[i][j]);
    println();
    }
  }
}
matrix identity(){
    matrix iden= new matrix();
    for (int i = 0; i < 4; i++)
            iden.data[i][i] = 1;
        return iden;
  }
matrix translation_mat(float tx, float ty, float tz){
    matrix mat=identity();
    mat.data[0][3]+=tx;
    mat.data[1][3]+=ty;
    mat.data[2][3]+=tz;
    
    return mat;

  }

matrix rotatex_mat(float angle){
    matrix mat=identity();
    mat.data[1][1]=cos(radians(angle));
    mat.data[1][2]=-sin(radians(angle));
    mat.data[2][1]=sin(radians(angle));
    mat.data[2][2]=cos(radians(angle));
    return mat;

  }

matrix rotatey_mat(float angle){
    matrix mat=identity();
    mat.data[0][0]=cos(radians(angle));
    mat.data[0][2]=sin(radians(angle));
    mat.data[2][0]=-sin(radians(angle));
    mat.data[2][2]=cos(radians(angle));
    return mat;

  }

matrix rotatez_mat(float angle){
    matrix mat=identity();
    mat.data[0][0]=cos(radians(angle));
    mat.data[0][1]=-sin(radians(angle));
    mat.data[1][0]=sin(radians(angle));
    mat.data[1][1]=cos(radians(angle));
    return mat;

  }

matrix scale_mat(float sx, float sy, float sz){
    matrix mat=identity();
    mat.data[0][0]*=sx;
    mat.data[1][1]*=sy;
    mat.data[2][2]*=sz;
    return mat;

  }
public class Hit{
  PVector hitCoord;
  Triangle hitObj;
  Surface sur;
  PVector snormal;
  float t;
}


Hit gettriangleintersect(Triangle tri, Ray r){
  Hit answer=new Hit();
  PVector tri_normal=(PVector.sub(tri.vec2,tri.vec1)).cross(PVector.sub(tri.vec3,tri.vec1)).normalize();
  float t_denom=PVector.dot(tri_normal,r.slope);
  
  if (t_denom!=0)
  {
    float t=PVector.dot(tri_normal,PVector.sub(tri.vec1,r.origin))/t_denom;
    if(t>=0)
    {
      PVector p=PVector.add(r.origin,PVector.mult(r.slope,t));
      float triple1=PVector.dot((PVector.sub(p,tri.vec1)).cross(PVector.sub(tri.vec2,tri.vec1)),tri_normal);
      float triple2=PVector.dot((PVector.sub(p,tri.vec2)).cross(PVector.sub(tri.vec3,tri.vec2)),tri_normal);
      float triple3=PVector.dot((PVector.sub(p,tri.vec3)).cross(PVector.sub(tri.vec1,tri.vec3)),tri_normal);

      if((triple1>=0 && triple2>=0 && triple3>=0)||(triple1<=0 && triple2<=0 && triple3<=0))
      {
        answer.hitCoord=p;
        answer.hitObj=tri;
        answer.sur=tri.sur;
        answer.t=t;

        if(PVector.dot(tri_normal,r.slope)<0){
          answer.snormal=tri_normal;
        }
        else{
          answer.snormal=PVector.mult(tri_normal,-1);
        }
      }

    
    }
  }
  return answer;

}

boolean debug_flag = false;
float fov;
Color bg;
PVector first;
PVector second;
ArrayList<Light> lightlist=new ArrayList<Light>();
ArrayList<Triangle> trianglelist=new ArrayList<Triangle>();
Stack<matrix> matrix_stack= new Stack<matrix>();
Surface curr_surf = new Surface();

void setup() {

  reset_scene();
  size (300, 300);  
  noStroke();
  colorMode(RGB, 1.0);
  background (0, 0, 0);
}

void keyPressed() {
  reset_scene();
  matrix_stack.push(identity());
  switch(key) {
    case '1': interpreter("s1.cli"); break;
    case '2': interpreter("s2.cli"); break;
    case '3': interpreter("s3.cli"); break;
    case '4': interpreter("s4.cli"); break;
    case '5': interpreter("s5.cli"); break;
    case '6': interpreter("s6.cli"); break;
  }
}

// this routine helps parse the text in a scene description file
void interpreter(String file) {

  println("Parsing '" + file + "'");
  String str[] = loadStrings (file);
  if (str == null) println ("Error! Failed to read the file.");
  
  for (int i = 0; i < str.length; i++) {
    
    String[] token = splitTokens (str[i], " ");   // get a line and separate the tokens
    if (token.length == 0) continue;              // skip blank lines

    if (token[0].equals("fov")) {
        fov=float(token[1]);
        println("fov "+fov);
    }
    else if (token[0].equals("background")) {
      bg= new Color();
      bg.r = float(token[1]);  // this is how to get a float value from a line in the scene description file
      bg.g = float(token[2]);
      bg.b = float(token[3]);
      println ("background = " + bg.r + " " + bg.g + " " + bg.b);
    }
    else if (token[0].equals("light")) {
      Light li = new Light();
      li.x =float(token[1]);
      li.y =float(token[2]);
      li.z =float(token[3]);
      li.r =float(token[4]);
      li.g =float(token[5]);
      li.b =float(token[6]);

      lightlist.add(li);
      println ("light = " + li.x + " " + li.y + " " + li.z + " "+li.r+ " "+li.g+ " "+li.b );
    }
    else if (token[0].equals("read")) {

        interpreter (token[1]);

    }
    else if (token[0].equals("surface")) {
      Surface new_surf= new Surface();
      curr_surf=new_surf;
      curr_surf.dr =float(token[1]);
      curr_surf.dg =float(token[2]);
      curr_surf.db =float(token[3]);
      
      println ("surface = " + curr_surf.dr + " " + curr_surf.dg + " " + curr_surf.db);
    }    
    else if (token[0].equals("begin")) {
      first = null;
      second = null;
    }
    else if (token[0].equals("vertex")) {
      matrix temp=new matrix();
      temp=matrix_stack.peek();

      float vert[][]={{float(token[1])},{float(token[2])},{float(token[3])},{1}};

      float trans[][] =temp.timesv(vert);
      // trans.show();

      if(first == null && second == null) 
      {
        first = new PVector(trans[0][0], trans[1][0], trans[2][0]);
      }
      else if (first != null && second == null)
      {
        second = new PVector(trans[0][0], trans[1][0], trans[2][0]);
      } 
      else{
        Triangle t = new Triangle(first, second, new PVector( trans[0][0], trans[1][0], trans[2][0] ));
        t.sur=curr_surf;
        trianglelist.add(t);
        
        println("Triangle "+ t.vec1+" "+t.vec2+" "+t.vec3+" ");
      }
      
    }
    else if (token[0].equals("end")) {
      first = null;
      second = null;
    }

    else if(token[0].equals("translate")){
      matrix mat=translation_mat(float(token[1]),float(token[2]),float(token[3]));
      matrix top=matrix_stack.pop();
      matrix prod=top.times(mat);
      matrix_stack.push(prod);

    }

    else if(token[0].equals("scale")){
      matrix mat=scale_mat(float(token[1]),float(token[2]),float(token[3]));
      matrix top=matrix_stack.pop();
      matrix prod=top.times(mat);
      matrix_stack.push(prod);

    }
    else if(token[0].equals("rotatex")){
      matrix mat=rotatex_mat(float(token[1]));
      matrix top=matrix_stack.pop();
      matrix prod=top.times(mat);
      matrix_stack.push(prod);

    }
    else if(token[0].equals("rotatey")){
      matrix mat=rotatey_mat(float(token[1]));
      matrix top=matrix_stack.pop();
      matrix prod=top.times(mat);
      matrix_stack.push(prod);

    }
    else if(token[0].equals("rotatez")){
      matrix mat=rotatez_mat(float(token[1]));
      matrix top=matrix_stack.pop();
      matrix prod=top.times(mat);
      matrix_stack.push(prod);

    }
    else if(token[0].equals("push")){
      matrix dup=new matrix();
      dup=matrix_stack.peek();
      matrix_stack.push(dup);
    }
    else if(token[0].equals("pop")){
      if(matrix_stack.size()==1)
      {
        println("error msg: only 1 element");
        break;
      }
      matrix_stack.pop();
    }


    else if (token[0].equals("render")) {
      draw_scene();   // this is where you actually perform the scene rendering
      println("Saving image to '" + file+"_res'");
      save(file+"res");
    }
    else if (token[0].equals("#")) {
      // comment (ignore)
    }

    else {
      println ("unknown command: " + token[0]);
    }
  }
  
}

void reset_scene() {
  // reset your scene variables here
  bg=null;
  fov=0;
  first=null;
  second=null;
  lightlist.clear();
  trianglelist.clear();
  matrix_stack.clear();
  
}

void draw_scene() {

  for(int y = 0; y < height; y++) {
    for(int x = 0; x < width; x++) {
      
      float k = tan(radians(fov)/2);
      float x_ = ((x-(width/2))*((2*k)/width));
      float y_ = -((y-(height/2))*((2*k)/height));
      
      Ray eyeray=new Ray();
      eyeray.origin=new PVector(0,0,0);
      eyeray.slope=new PVector(x_,y_,-1);

      ArrayList<Triangle> cand_shapes=new ArrayList<Triangle>();
      ArrayList<Float> cand_dis=new ArrayList<Float>();
      ArrayList<Float> cand_t=new ArrayList<Float>();
      ArrayList<Hit> cand_hit=new ArrayList<Hit>();

      Hit curr_hit=new Hit();
      Hit best_hit=new Hit();

      Triangle best_tri=null;

      for (int i = 0; i < trianglelist.size(); i++)
      {
        Triangle curr_tri=trianglelist.get(i);
        curr_hit=gettriangleintersect(curr_tri,eyeray);
        if(curr_hit.hitCoord!=null){
          float distance=PVector.dist(curr_hit.hitCoord,eyeray.origin);
          cand_shapes.add(curr_tri);
          cand_dis.add(distance);
          cand_hit.add(curr_hit);
          cand_t.add(curr_hit.t);
        }
      }
      if(cand_shapes.size()==0){
        if(bg==null)
        {
          println("bg is null");
        }
        color pixcolor = color(bg.r, bg.g, bg.b);
        set (x, y, pixcolor);
      }
      else{

        int min_t_idx = 0;

        for (int i=1; i < cand_t.size(); ++i) {
          if (cand_t.get(i) < cand_t.get(min_t_idx)) {
            min_t_idx = i;
          }
        }
        best_hit=cand_hit.get(min_t_idx);
        best_tri=cand_shapes.get(min_t_idx);

        PVector temp=getcolor(best_tri,best_hit,eyeray);
        color pixcolor=color(temp.x,temp.y,temp.z);
        set (x, y, pixcolor);
      }

    }
  }
}

PVector getcolor(Triangle tri, Hit hit, Ray eyeray){
  
  if(hit.hitObj == null){
    if (bg != null)  return new PVector(bg.r, bg.g, bg.b);
    else return new PVector(0,0,0);
  }
  PVector hitpoint = hit.hitCoord;  //answer.hitCoord = intersection point
  PVector result = new PVector(0,0,0);
  for (int i = 0; i < lightlist.size(); i++)
  {
    //diffuse lighting
    PVector light = new PVector(lightlist.get(i).x, lightlist.get(i).y, lightlist.get(i).z);
    PVector lightsource = PVector.sub(light, hitpoint);
    
    lightsource.normalize();
    float dotprod = PVector.dot(lightsource, hit.snormal);
    float maxi = max(dotprod, 0);
    //shadows
    Ray shadowray = new Ray();
    shadowray.origin=hitpoint;
    shadowray.slope=lightsource;
    int castshadow=1;
    for(int j=0;j<trianglelist.size();j++)
    {
      Hit checkobj=gettriangleintersect(trianglelist.get(j),shadowray);
      if(checkobj.hitCoord!=null && checkobj.t>0.0001 && checkobj.t< hitpoint.dist(light))
      {
        castshadow=0;
        break;
      }
      
    }
      result.x += (lightlist.get(i).r * maxi * hit.sur.dr * castshadow); 
      result.y += (lightlist.get(i).g * maxi * hit.sur.dg * castshadow) ;
      result.z += (lightlist.get(i).b * maxi * hit.sur.db * castshadow) ;
  }

  return result;
}
// prints mouse location clicks, for help debugging
void mousePressed() {
  println ("You pressed the mouse at " + mouseX + " " + mouseY);
}

// you don't need to add anything here
void draw() {
}
