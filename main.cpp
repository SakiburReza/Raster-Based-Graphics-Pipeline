#include<iostream>
#include<cmath>
#include<cassert>
#include<bits/stdc++.h>
#include "bitmap_image.hpp"

#define pi (2*acos(0.0))

using namespace std;


class Homogenous_Point{
    public:
    double x,y,z,w;
    Homogenous_Point(){}
    Homogenous_Point(double x, double y,double z){
        this->x= x;
        this->y=y;
        this->z = z;
        this->w=1;
    }
    Homogenous_Point(double x, double y,double z, double w){
        this->x= x;
        this->y=y;
        this->z = z;
        this->w=w;
    }
};

class Vector3D{
	public:
	double x,y,z;
	Vector3D(){
		x=0;
		y=0;
		z=0;
	}
	Vector3D(double x, double y, double z){
		this->x = x;
		this->y =y;
		this->z = z;
	}
	double magnitude(){
		return sqrt((this->x * this->x) + (this->y * this->y)+(this->z * this->z));
	}
    void normalize(){
        double val = this->magnitude();
        x = x/val;
        y = y/val;
        z = z/val;
    }

	Vector3D operator+(const Vector3D a){
		return Vector3D(this->x+a.x, this->y+a.y ,this->z+a.z);
	}
	Vector3D operator-(const Vector3D a){
		return Vector3D(this->x-a.x, this->y-a.y ,this->z-a.z);
	}
	double operator*(const Vector3D a){
		return (this->x * a.x+ this->y * a.y + this->z*a.z);
	}
	Vector3D operator*(const double a){
		return Vector3D(this->x * a, this->y * a , this->z*a);
	}
	Vector3D operator^(const Vector3D a){
		return Vector3D((this->y*a.z - this->z*a.y), (this->z*a.x-this->x*a.z) , (this->x*a.y-this->y*a.x));
	}
	void print(){
		printf("x = %lf , y = %lf , z = %lf\n",x,y,z);
	}

};

class Color{
    public:
    double r, g, b;
    Color(){
        r=10;
        g=10;
        b=10;
    }
    Color(double r, double g, double b){
        this->r = r;
        this->g = g;
        this->b = b;
    }
};

class Triangle{
    public:
    Homogenous_Point points[3];
    Color color;
    Triangle(Homogenous_Point pointA,Homogenous_Point pointB, Homogenous_Point pointC,Color color){
        points[0] = pointA;
        points[1] = pointB;
        points[2] = pointC;
        this->color = color;
    }
    Triangle(){}  
    void decendingSortByYValue(){

         for(int i=0;i<2;i++){
            for(int j=i+1; j < 3 ;j++){
                if(points[i].y < points[j].y){
                    Homogenous_Point temp = points[i];
                    points[i]=points[j];
                    points[j]=temp;//swap
                }
            }
         }
     }
     int uniqueIndexForY(){

       if(points[1].y == points[2].y ) return 0;
       if(points[2].y == points[0].y ) return 1;
       if(points[0].y == points[1].y) return 2;
       return 0;

     }

      double xmin(){
        double mx = points[0].x;
        for ( int i = 1  ; i < 3 ; i++){
            if( mx > points[i].x ) mx = points[i].x;
        }
       return mx;
     }

     double xmax(){

        double mx = points[0].x;
        for ( int i = 1  ; i < 3 ; i++){
            if( mx < points[i].x ) mx = points[i].x;
        }
       return mx;
     }   

    
};

class Matrix{
    public:
        double values[4][4];
        int rows, columns;
    
    Matrix(int n){
        rows = n;
        columns = n;
    }
    Matrix( int row, int column){
        rows = row;
        columns = column;
    }
    static Matrix identity_matrix(int n){
        Matrix m(n);

        for(int i=0;i<n;i++){
            for (int j=0;j<n;j++){
                if(i==j) m.values[i][j]=1;
                else m.values[i][j]=0;
            }
        }
        return m;
    }
    Matrix operator+(const Matrix mat){
        assert(rows==mat.rows && columns == mat.columns);

        Matrix m(rows,columns);
        for(int i=0;i<rows;i++){
            for (int j=0;j<columns;j++){
                m.values[i][j] = this->values[i][j] + mat.values[i][j];
            }
        }
        return m;
    }
    Matrix operator-(const Matrix mat){
        assert(rows==mat.rows && columns == mat.columns);

        Matrix m(rows,columns);
        for(int i=0;i<rows;i++){
            for (int j=0;j<columns;j++){
                m.values[i][j] = this->values[i][j] - mat.values[i][j];
            }
        }
        return m;
    }
    Matrix operator*(const Matrix mat){
        assert(columns==mat.rows);

        Matrix m(rows,mat.columns);

        for (int i = 0; i < m.rows; i++) {
            for (int j = 0; j < m.columns; j++) {
                double val = 0;
                for (int k = 0; k < columns; k++) {
                    val += this->values[i][k] * mat.values[k][j];
                }
                m.values[i][j] = val;
            }
        }
        return m;
    }
    Matrix operator*(const double k){
        Matrix m(rows,columns);

        for(int i=0;i<m.rows;i++){
            for(int j = 0; j<m.columns;j++){
                m.values[i][j] = k * values[i][j];
            }
        }
        return m;
    }
    Homogenous_Point operator*(const Homogenous_Point p){
        assert(rows==4 && columns == 4);

        Matrix m(4,1);  //What about grabage?

        m.values[0][0]=p.x;
        m.values[1][0]=p.y;
        m.values[2][0]=p.z;
        m.values[3][0]=p.w;

        Matrix mat = (*this) * m;

        return Homogenous_Point(mat.values[0][0]/mat.values[3][0],mat.values[1][0]/mat.values[3][0],mat.values[2][0]/mat.values[3][0]);
    }
    void print(){
        for(int i = 0; i<rows;i++){
            for (int j = 0; j < columns ; j++)
            {
                cout<<values[i][j]<<"\t";
            }
            cout<<endl;
        }
    }
};

Matrix translation_Matrix(double tx, double ty, double tz)
{
    Matrix m = Matrix::identity_matrix(4);
    m.values[0][3]=tx;
    m.values[1][3]=ty;
    m.values[2][3]=tz;
    //m.values[3][3]=1;
    return m;
}
Matrix scale_Matrix(double sx, double sy, double sz)
{
    Matrix m = Matrix::identity_matrix(4);
    m.values[0][0]=sx;
    m.values[1][1]=sy;
    m.values[2][2]=sz;
    return m;
}
Matrix rotate_Matrix(double angle, double ax, double ay, double az)
{
    Vector3D a(ax,ay,az);
    a.normalize();

    Vector3D i(1,0,0), j(0,1,0), k(0,0,1);

    double c = cos(angle * pi/180);
    double s = sin(angle * pi/180);

    Vector3D c1 = i*c + a * (a*i) * (1-c) + (a^i)*s;
    Vector3D c2 = j*c + a * (a*j) * (1-c) + (a^j)*s;
    Vector3D c3 = k*c + a * (a*k) * (1-c) + (a^k)*s;
   
   Matrix rotationMatrix = Matrix::identity_matrix(4);

   rotationMatrix.values[0][0] = c1.x;
   rotationMatrix.values[1][0] = c1.y;
   rotationMatrix.values[2][0] = c1.z;

   rotationMatrix.values[0][1] = c2.x;
   rotationMatrix.values[1][1] = c2.y;
   rotationMatrix.values[2][1] = c2.z;

   rotationMatrix.values[0][2] = c3.x;
   rotationMatrix.values[1][2] = c3.y;
   rotationMatrix.values[2][2] = c3.z;



    return rotationMatrix;


}
double eyeX, eyeY, eyeZ;
double lookX, lookY, lookZ;
double upX, upY, upZ;
double fovX,fovY, aspectRatio, near, far;
int screenHeight, screenWidth;
double x_left, y_bottom, x_right, y_top, front, rear;
Color background(0,0,0);


Matrix matrix = Matrix::identity_matrix(4);
stack <Matrix> matrix_stack;

void Modeling_Transformations()
{
    freopen("scene.txt","r",stdin);
    freopen("stage1.txt","w",stdout);

    cout<<std::fixed;
    cout<<std::setprecision(7);

    cin>>eyeX>>eyeY>>eyeZ;
    cin>>lookX>>lookY>>lookZ;
    cin>>upX>>upY>>upZ;
    cin>>fovY>>aspectRatio>>near>>far;


    string temp;

    while(1){
        cin>>temp;
       
        if(temp=="triangle"){
            double x, y, z;
            for(int i =0;i<3;i++){
                cin>>x>>y>>z;
                Homogenous_Point p(x,y,z);
                p = matrix*p;
                cout<<p.x<<" "<<p.y<<" "<<p.z<<endl;
            }
            cout<<endl;
        }
        else if(temp=="translate"){
            double tx, ty, tz;

            cin>>tx>>ty>>tz;
            matrix = matrix * translation_Matrix(tx,ty,tz);

        }
        else if(temp=="scale"){
            double sx, sy, sz;
            cin>>sx>> sy>>sz;
            matrix = matrix * scale_Matrix(sx,sy,sz);
        }
        else if(temp=="rotate"){
            double angle, ax, ay, az;
            cin>>angle>>ax>>ay>>az;
            matrix = matrix * rotate_Matrix(angle,ax,ay,az);
        }
        else if(temp=="push"){
            matrix_stack.push(matrix);
        }
        else if(temp=="pop"){
            if(!matrix_stack.empty()){
                matrix = matrix_stack.top();
                matrix_stack.pop();
            }
        }
        else if(temp=="end") break;
    }
    fclose(stdin);
    fclose(stdout);

    
}
void View_Transformation(){
    ifstream stage1;
    ofstream stage2;
    stage1.open ("stage1.txt");
    stage2.open ("stage2.txt");
    stage2 << std::fixed;
    stage2 << std::setprecision(7);

    Vector3D l(lookX-eyeX, lookY-eyeY,lookZ-eyeZ);
    l.normalize();
    Vector3D u(upX,upY,upZ);
    Vector3D r = l^u;
    r.normalize();
    u = r ^ l;

    Matrix T = translation_Matrix(-eyeX,-eyeY,-eyeZ);
    Matrix R = Matrix::identity_matrix(4);

    R.values[0][0] = r.x;
    R.values[0][1] = r.y;
    R.values[0][2] = r.z;

    R.values[1][0] = u.x;
    R.values[1][1] = u.y;
    R.values[1][2] = u.z;

    R.values[2][0] = -l.x;
    R.values[2][1] = -l.y;
    R.values[2][2] = -l.z;

    Matrix V = R * T;

    while(1){
        double x, y,z;
        bool eof_found = false;
        for(int i =0; i<3; i ++){
            stage1>>x;
            if(stage1.eof()) {
                eof_found = true;
                break;
            }
            stage1>>y>>z;
            Homogenous_Point p(x,y,z);
            p = V * p;
            stage2<<p.x<<" "<<p.y<<" "<<p.z<<endl;
        }
        if(eof_found) break;
        stage2<<endl;
    }
    stage1.close();
    stage2.close();
}
void Projection_Transformation()
{
    ifstream stage2;
    ofstream stage3;
    stage2.open ("stage2.txt");
    stage3.open ("stage3.txt");
    stage3 << std::fixed;
    stage3 << std::setprecision(7);

    fovX = fovY * aspectRatio;
    double t = near * tan((pi/180.0) * (fovY/2.0));
    double r = near * tan((pi/180.0) * (fovX/2.0));

    Matrix p = Matrix::identity_matrix(4);

    p.values[0][0] = near/r;
    p.values[1][1] = near/t;
    p.values[2][2] = -(far+near)/(far-near);
    p.values[2][3] = -(2*far*near)/(far-near);
    p.values[3][2] = -1;
    p.values[3][3] = 0;

     while(1){
        double x, y,z;
        bool eof_found = false;
        for(int i =0; i<3; i ++){
            stage2>>x;
            if(stage2.eof()) {
                eof_found = true;
                break;
            }
            stage2>>y>>z;
            Homogenous_Point q(x,y,z);
            q = p * q;
            stage3<<q.x<<" "<<q.y<<" "<<q.z<<endl;
        }
        if(eof_found) break;
        stage3<<endl;
    }
    stage2.close();
    stage3.close();
}

void Clipping_ScanConversion(){
    ifstream config, stage3;
    ofstream out;
    config.open("config.txt");
    stage3.open("stage3.txt");
    out.open("z-buffer.txt");

    cout<<std::fixed;
    cout<<std::setprecision(7);
    
    config>>screenWidth>>screenHeight;
    config>>x_left;
    config>>y_bottom;
    config>>front>>rear;

    x_right = -x_left;
    y_top = - y_bottom;

    queue<Triangle> triangles;
    double x,y,z,dx,dy;
    bool eof_found = false;
    Triangle temp;

    for(int i = 0;; i++){
       

        for ( int j=0;j<3;j++){
            stage3>>x;
            if(stage3.eof()){
                eof_found = true;
                break;
            }
            stage3>>y>>z;
            Homogenous_Point p(x,y,z);
            temp.points[j] = p;
        }
        if(eof_found) break;

            temp.color.r = rand() % 256;
            temp.color.g = rand() % 256;
            temp.color.b = rand() % 256;
            triangles.push(temp);
    }
    

    config.close();
    stage3.close();

    Color **frameBuffer = new Color *[screenHeight];
    double **zBuffer = new double *[screenHeight];

    for(int i=0; i<screenHeight; i++){
        frameBuffer[i] = new Color [screenWidth];
        zBuffer[i] = new double [screenWidth];

        for(int j=0;j<screenWidth;j++){
            frameBuffer[i][j] = background;
            zBuffer[i][j] = rear;
        } 
    }

    dx = (x_right - x_left) / screenWidth;
    dy = (y_top - y_bottom) / screenHeight;
    
    while(triangles.size()!=0){
        temp = triangles.front();
        triangles.pop();

        temp.decendingSortByYValue();

        double xMax = temp.xmax();
        double xMin = temp.xmin();

        int yLine[3];

        for(int i=0; i<3 ;i++){
            double ty = y_top - temp.points[i].y;
            yLine[i] = ty/dy;
        }
        int id0 = temp.uniqueIndexForY();
        int id1 = (id0+1) % 3;
        int id2 = (id0+2) % 3;

        for(int line_y = yLine[0]+1 ; line_y<=yLine[2];line_y++){
            if(line_y <= 0 || (line_y>=screenHeight-1)) continue;

            double ys = y_top - (dy*line_y);

            double xA = (ys - temp.points[id1].y) / (temp.points[id0].y - temp.points[id1].y) * (temp.points[id0].x - temp.points[id1].x) +temp.points[id1].x;
            double zA = (ys - temp.points[id1].y) / (temp.points[id0].y - temp.points[id1].y) * (temp.points[id0].z - temp.points[id1].z) +temp.points[id1].z;
             double xB = (ys - temp.points[id2].y) / (temp.points[id0].y - temp.points[id2].y) * (temp.points[id0].x - temp.points[id2].x) +temp.points[id2].x;
              double zB = (ys - temp.points[id2].y) / (temp.points[id0].y - temp.points[id2].y) * (temp.points[id0].z - temp.points[id2].z) +temp.points[id2].z;

              double x_l, x_r, z_l, z_r ; //l=left, r=right
              if(xA < xB){
                x_l = xA;
                z_l = zA;
                x_r = xB;
                z_r = zB;
              }
              else{
                x_l = xB;
                z_l = zB;
                x_r = xA;
                z_r = zA;
              }
              if(x_l < xMin) {
                x_l = xMin + abs(xMin -  x_l);
              }
              if(x_r > xMax){
                x_r = xMax - abs(xMax - x_r);
              }

              int xLLine = (- x_left + x_l)/dx;
              int xRLine = (-x_left+x_r)/dx;
              
              for(int xLine = xLLine+1; xLine<xRLine; xLine++){
                if(xLine < 0 || xLine>=(screenWidth-1)) continue;

                double xs = xLine * dx + x_left ;
                double zs = ((xs-x_r)/(x_l-x_r)) * (z_l - z_r) + z_r;

                if(zs < rear && zBuffer[xLine][line_y]>zs){
                    zBuffer[xLine][line_y] = zs;
                    frameBuffer[xLine][line_y] = temp.color;
                    out<<zs<<" ";
                }
                    
                
              }

        }
    }

    //Creating bmp images
    bitmap_image image(screenWidth,screenHeight);

    for(int i=0;i<screenHeight;i++){
        for(int j=0;j<screenWidth;j++){
            image.set_pixel(i,j,frameBuffer[i][j].r,frameBuffer[i][j].g,frameBuffer[i][j].b);

        }
    }
    image.save_image("out.bmp");
    out<<endl;
    out.close();

    //Freeing the memory

    for(int i = 0; i<screenHeight;i++){
        delete frameBuffer[i];
        delete zBuffer[i];
    }
    delete frameBuffer;
    delete zBuffer;
}
int main()
{
    Modeling_Transformations();
    View_Transformation();
    Projection_Transformation();

    Clipping_ScanConversion();

    return 0;
}