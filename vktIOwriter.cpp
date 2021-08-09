#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

struct Point3D{
    double x; 
    double y;
    double z;
    
    //Constructor for 3D Point
    Point3D(double _x, double _y, double _z)
    {
        x = _x;
        y = _y;
        z = _z;
    }
    ~Point3D(){}
    
//Function to get X, Y and Z 
    double getX() 
    { 
        return x; 
    }
    double getY() 
    { 
        return y; 
    }
    double getZ() 
    { 
        return z;
    }

};


/*DATASET STRUCTURED_POINTS
DIMENSIONS nx ny nz
ORIGIN xyz
SPACING sx sy sz*/
void WriteStructuredPoints( float nx, float ny, float nz, float x, float y, float z, float sx, float sy, float sz)
{
    float Nx = nx, Ny = ny, Nz = nz;
    float X = x, Y = y, Z = z;
    float Sx = sx, Sy = sy, Sz = sz;
    std::ofstream file;
    file.open("StructuredPiont.vtk");
    file << "# vtk DataFile Version 3.0" <<std::endl;
    file << "Volume example " <<std::endl;
    file << "ASCII" <<std::endl;
    file << "DATASET STRUCTURED_POINTS "<<std::endl;
    file << "DIMENSIONS" << " " <<Nx <<" "<<Ny <<" "<<Nz<<std::endl;
    file << "SPACING" << " " <<Sx <<" "<<Sy <<" "<<Sz<<std::endl;
    file << "ORIGIN" << " " <<x <<" "<<y <<" "<<z<<std::endl;
    file.close();
}

void WriteStructuredGrid (float nx, float ny, float nz, int n, std::vector<Point3D> &pv, std::vector<float> &sv)
{
    float Nx = nx, Ny = ny, Nz = nz;
    std::ofstream file;
    file.open("StructuredGrid.vtk");
    file << "# vtk DataFile Version 3.0" <<std::endl;
    file << "Volume example" <<std::endl;
    file << "ASCII" <<std::endl;
    file << "DATASET STRUCTURED_GRID "<<std::endl;
    file << "DIMENSIONS" << " " <<Nx <<" "<<Ny <<" "<<Nz<<std::endl;
    file << "POINTS" <<" "<< n << " "<< "float" <<std::endl;
    // for (int i = 0; i < pv; i++)
    // {
    //     file << pv.at(i).getX() << " " <<pv.at(i).getY() << " "<<pv.at(i).getZ()<< std::endl;
    // }
    for (auto &x: pv)
        file << x.getX() <<" " << x.getY() <<" "<< x.getZ()<< std::endl;
    file << "POINT_DATA" << " " << n << std::endl;
    file << "SCALARS sample_scalars float 1"<<std::endl;
    file << "LOOKUP_TABLE default"<<std::endl;

    for (auto &x :sv)
        file << x <<std::endl;
        
    file.close();
}

void WriteRectilinearGrid(float nx, float ny, float nz, int n, std::vector<Point3D> &pv, std::vector<float> &sv)
{
    float Nx = nx, Ny = ny, Nz = nz;
    std::ofstream file;
    file.open("RectilinearGrid.vtk");
    file << "# vtk DataFile Version 3.0" <<std::endl;
    file << "Volume example" <<std::endl;
    file << "ASCII" <<std::endl;
    file << "DATASET RECTILINEAR_GRID"<<std::endl;
    file << "DIMENSIONS" << " " <<n <<" "<<n <<" "<<n<<std::endl;
    file << "X_COORDINATES" <<" "<< n << " "<< "float" <<std::endl;
    // for (int i = 0; i < pv; i++)
    // {
    //     file << pv.at(i).getX() << " " <<pv.at(i).getY() << " "<<pv.at(i).getZ()<< std::endl;
    // }
    for (auto &x: pv)
        file << x.getX() <<" " ;
    std::cout << std::endl;
    std::cout << std::endl;

    file << "Y_COORDINATES" <<" "<< n << " "<< "float" <<std::endl;
    for (auto &x: pv)
        file << x.getY() <<" " ; 
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    
    file << "Z_COORDINATES" <<" "<< n << " "<< "float" <<std::endl;
    for (auto &x: pv)
        file << x.getZ() <<" " ; 
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    // file << "POINT_DATA" << " " << n << std::endl;
    // file << "SCALARS sample_scalars float 1"<<std::endl;
    // file << "LOOKUP_TABLE default"<<std::endl;

    // for (auto &x :sv)
    //     file << x <<std::endl;
        
    file.close();
}


float circle (Point3D m, Point3D cm, float r_)
{
    Point3D cm_(cm);
    std::vector<float> sv;
    float r = r_;
    float sq_r = r*r;
    float ref;
    Point3D x(m);
    ref = (pow((x.getX() - cm.getX()),2)) + (pow((x.getY() - cm.getY()),2)) + (pow((x.getZ() - cm.getZ()),2));
    
    if (sq_r > ref )
    {
        return 1;
    }
    else 
    {
        return 0;
    }
}

float square (Point3D m, Point3D cm, float r_)
{
    Point3D cm_(cm);
    std::vector<float> sv;
    float r = r_;
    float sq_r = r*r;
    float ref;
    Point3D x(m);
    ref = (abs((x.getX() - cm.getX()) + (x.getY() - cm.getY()) + (x.getZ() - cm.getZ()))) + (abs((x.getX() - cm.getX()) - (x.getY() - cm.getY()) - (x.getZ() - cm.getZ()))) ; 
    if (sq_r > ref )
    {
        return 1;
    }
    else 
    {
        return 0;
    }
}

float triangle_area (Point3D a, Point3D b, Point3D c)
{
    Point3D a_(a), b_(b), c_(c);
    float dX, dY, dZ;

    dX = 0.5 *(abs(a_.getY()*(b_.getZ() - c_.getZ()) - (b_.getY()*(a_.getZ() - c_.getZ())) + (c_.getY()*(a_.getZ() - b_.getZ()))));
    dY = 0.5 *(abs(a_.getZ()*(b_.getX() - c_.getX()) - (b_.getZ()*(a_.getX() - c_.getX())) + (c_.getZ()*(a_.getX() - b_.getX()))));
    dZ = 0.5 *(abs(a_.getX()*(b_.getY() - c_.getY()) - (b_.getX()*(a_.getY() - c_.getY())) + (c_.getX()*(a_.getY() - b_.getY()))));

    float area = sqrt((pow(dX,2)) + (pow(dY,2)) + (pow(dZ,2)));
    return area;
}

float inside_triangle(Point3D a, Point3D b, Point3D c, Point3D d)
{
    float A = triangle_area(a, b, c);
    float A1 = triangle_area(a, b, d);
    float A2 = triangle_area (a, c, d);
    float A3 = triangle_area(b, c, d);

    float A_sum = A1 + A2 + A3;
    if (A == A_sum)
    {
        return 1.0;
    }
    else
    {
        return 0.0;
    }

}

int main()
{
    std::cout <<"vkt IO Writer " << std::endl;
    /*For 2D Nz = 1*/
    

    std::vector<Point3D> input_points;
    float nx = 100.0f, ny =100.0f, nz = 1.0f; //domain size
    float sx = 0.5f, sy = 0.5f, sz = 1.0f; //Spacing size
    float n = (nx/sx)*(ny/sy)*(nz/sz); //total no of points;
    std :: cout << "Creating cartesian mesh is started......." <<std::endl;
    for (float z = 0.0; z < nz; z += sz)
    {
        for (float y = 0.0; y < ny ; y += sy)
        {
            for (float x = 0.0; x < nx ; x += sx)
            {
                // Point3D temp;
                std :: cout <<x << " " <<y <<" " << z <<std::endl;
                
                Point3D temp(x, y, z);
                input_points.push_back(temp);
            }
        }
    }
    std :: cout << "Creating cartesian mesh is done......." <<std::endl;

    std::vector<float> sv;
    float j;
    Point3D cm_(50,50,0);
    float r_= 15.0;
    Point3D a(20, 20, 0), b(70, 20, 0), c(60, 50, 0);
    for (auto &x: input_points)
    {
        Point3D m(x.getX(), x.getY(), x.getZ());

        j = circle(m, cm_, r_);
        // j =  square(m, cm_, r_ );
        j = inside_triangle(a, b, c, m);
        sv.push_back(j);

    }
    
    // WriteStructuredPoints(10, 5, 1, 0, 0, 0, 1, 1, 1); //DIMENSIONS nx ny nz, ORIGIN xyz SPACING sx sy sz
    WriteStructuredGrid(nx, ny, nz, n, input_points, sv );
    // WriteRectilinearGrid(nx, ny, nz, n, input_points, sv);


    


    // Point3D a(3, 8, 1), b(-4, 2, 1), c(5, 1, 2);

    // float area = triangle_area(a, b, c);
    // std::cout << "Triangle area: " << " " <<area <<std::endl; 
    return 0;
        
}