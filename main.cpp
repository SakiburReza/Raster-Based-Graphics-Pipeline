#include<iostream>
#include<cmath>
#include<cassert>
#include<bits/stdc++.h>

#define pi (2*acos(0.0))

using namespace std;


class Homogenous_Point{
    public:
    double x,y,z,w;

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

        return Homogenous_Point(mat.values[0][0],mat.values[1][0],mat.values[2][0],mat.values[3][0]);
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
double fovY, aspectRatio, near, far;

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
int main()
{
    Modeling_Transformations();
    View_Transformation();
    
}