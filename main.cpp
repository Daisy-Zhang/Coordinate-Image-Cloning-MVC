// std and stb_image library
#include <iostream>
#include <bits/stdc++.h>
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

// CGAL includes for triangulation
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>

// CGAL typedefs and definitions
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> DelaunayTriangulation;

typedef CGAL::Delaunay_mesh_size_criteria_2<DelaunayTriangulation> Criteria;

typedef DelaunayTriangulation::Finite_vertices_iterator FiniteVerticesIterator;
typedef DelaunayTriangulation::Finite_faces_iterator FiniteFacesIterator;
typedef DelaunayTriangulation::Finite_edges_iterator FiniteEdgesIterator;
typedef DelaunayTriangulation::Face_circulator FacesCirculator;
typedef DelaunayTriangulation::Point_iterator PointIterator;
typedef DelaunayTriangulation::Vertex Vertex;
typedef DelaunayTriangulation::Vertex_handle VertexHandle;
typedef DelaunayTriangulation::Edge Edge;
typedef DelaunayTriangulation::Face Face;
typedef DelaunayTriangulation::Face_handle FaceHandle;
typedef DelaunayTriangulation::Triangle Triangle;

typedef CGAL::Cartesian<double> Metric;
typedef DelaunayTriangulation::Point Point;
typedef CGAL::Polygon_2<Metric> Polygon2D;
typedef CGAL::Bbox_2 BBox2D;
typedef CGAL::Delaunay_mesher_2<DelaunayTriangulation, Criteria> Mesher;

std::vector<Point> boundary_point;
std::map<Point, int> triangle_vertex_index;
std::vector<Triangle> triangle_meshes;
std::vector<Point> vertex_list;

#define getpixel(img, i, j, k) img.image[int((i) * (img.w) * (img.c) + (j) * (img.c) + (k))]

struct Image{
    double* image;
	unsigned char* buf;
	std::string filename;
	int w, h, c;
	int w0, w1, h0, h1;

    Image(): image(NULL), w(0), h(0), c(0) {}
	Image(int _h, int _w, int _c)
	{
		h = _h, w = _w, c = _c;
		image = new double[w * h * c];
		buf = new unsigned char[w * h * c];
		memset(image, 0, sizeof(double) * w * h * c);
		memset(buf, 0, sizeof(unsigned char) * w * h * c);
	}
	Image(std::string _filename): filename(_filename) {
		buf = stbi_load(filename.c_str(), &w, &h, &c, 0);
		image = new double[w * h * c];
		for (int i = 0; i < w * h * c; i ++)
			image[i] = buf[i];
	}
	void write(const char* output_filename){
		unsigned char* buf_ = new unsigned char[w * h * c];
		for (int i = 0; i < w * h * c; i ++)
			buf_[i] = image[i] < 0 ? 0 : image[i] > 255 ? 255 : image[i];
		stbi_write_png(output_filename, w, h, c, buf_, 0);
	}
};

void getAdaptiveMesh(const std::vector<Point> &boundary) {
	DelaunayTriangulation* adaptive_mesh = new DelaunayTriangulation();
	adaptive_mesh -> clear();
	adaptive_mesh -> insert(boundary.begin(), boundary.end());
	for (int i = 0; i < boundary.size(); i ++) {
		Point cur_point = boundary[i];
		Point next_point = boundary[(i + 1) % boundary.size()];
		adaptive_mesh -> insert_constraint(cur_point, next_point);
	}

	Mesher mesher(*adaptive_mesh);
	mesher.set_criteria(Criteria(0.125, 0));
	mesher.refine_mesh();
	printf("Done with:\nnumber of vertices: %d\nnumber of triangles: %d \n", adaptive_mesh -> number_of_vertices(), adaptive_mesh -> number_of_faces());
	
	int index = 0;
	for (FiniteFacesIterator iter = adaptive_mesh -> finite_faces_begin(); iter != adaptive_mesh -> finite_faces_end(); iter ++) {
		Triangle triangle = adaptive_mesh -> triangle(iter);

		Point a = triangle.vertex(0);
		Point b = triangle.vertex(1);
		Point c = triangle.vertex(2);

		if (triangle_vertex_index.find(a) == triangle_vertex_index.end()) {
			triangle_vertex_index[a] = index ++;
			vertex_list.push_back(a);
		}
		if (triangle_vertex_index.find(b) == triangle_vertex_index.end()) {
			triangle_vertex_index[b] = index ++;
			vertex_list.push_back(b);
		}
		if (triangle_vertex_index.find(c) == triangle_vertex_index.end()) {
			triangle_vertex_index[c] = index ++;
			vertex_list.push_back(c);
		}

		triangle_meshes.push_back(triangle);
	}
}

double getDistance(const Point &a, const Point &b) {
	return (a.x() - b.x()) * (a.x() - b.x()) + (a.y() - b.y()) * (a.y() - b.y());
}

double getDotProduct(const Point &a, const Point &b, const Point &c) {
	return (b.x() - a.x()) * (c.x() - a.x()) + (b.y() - a.y()) * (c.y() - a.y());
}

double getTanValue(const Point &a_, const Point &b_, const Point &c_) {
	// cal tan(a/2)
	double ab = getDistance(a_, b_);
	double ac = getDistance(a_, c_);
	double dot_product = getDotProduct(a_, b_, c_);

	return sqrt((ab*ac - dot_product) / (ab*ac + dot_product));
}

struct TransformWeight{
	double w[3];
	double getRes(double x, double y) {return x * w[0] + y * w[1] + w[2];}
}t_w[5];

void calTransformWeight(TransformWeight &tw, const Point &v1, const Point &v2, const Point &v3, double c0, double c1, double c2)
{
	double a[15][15];
	int n=3;
	a[1][1] = v1.x(), a[2][1] = v2.x(), a[3][1] = v3.x();
	a[1][2] = v1.y(), a[2][2] = v2.y(), a[3][2] = v3.y();
	a[1][3] = a[2][3] = a[3][3] = 1;
	a[1][4] = c0, a[2][4] = c1, a[3][4] = c2;

	int i, j, k, las;
	double t;
	for(i = 1; i <= n; i ++) {
		for(t = 0, las = j = i; j <= n; j ++)
			if(abs(a[j][i])>t)t=abs(a[j][i]),las=j;
		if(j = las, j != i)
			for(k=1;k<=n+1;k++)t=a[i][k],a[i][k]=a[j][k],a[j][k]=t;
		for(j = i + 1; j <= n; j ++)
			for(t=a[j][i]/a[i][i],k=i;k<=n+1;k++)a[j][k]-=a[i][k]*t;
	}
	for(i = n; i >= 1; i --)
	for(a[i][n+1]/=a[i][i],j=i-1;j;j--)a[j][n+1]-=a[j][i]*a[i][n+1];

	tw.w[0] = a[1][4], tw.w[1] = a[2][4], tw.w[2] = a[3][4];
}

double evaluatePointToLine(double x, double y, double x1, double y1, double x2, double y2)
{
	double a = y2 - y1;
	double b = x1 - x2;
	double c = x2 * y1 - x1 * y2;
 
	return a * x + b * y + c;
}

bool inside(const Point &p, const Point &a, const Point &b, const Point &c)
{
	double d1 = evaluatePointToLine(p.x(), p.y(), a.x(), a.y(), b.x(), b.y());
	double d2 = evaluatePointToLine(p.x(), p.y(), b.x(), b.y(), c.x(), c.y());
	if (d1 * d2 < 0)
		return false;
 
	double d3 = evaluatePointToLine(p.x(), p.y(), c.x(), c.y(), a.x(), a.y());
	if (d2 * d3 < 0)
		return false;
 
	return true;
}

int main(int argc, char *argv[]) {
    // para process
    std::string mask_filename, src_filename, target_filename, output_filename;
    double pos_x, pos_y = 0;
    
    mask_filename = argv[1];
    src_filename = argv[2];
    target_filename = argv[3];
    output_filename = argv[4];
    pos_x = atoi(argv[5]);
    pos_y = atoi(argv[6]);

    // read data process
    Image src_image(src_filename), mask_image(mask_filename), target_image(target_filename);
    if(mask_image.w != src_image.w || mask_image.h != src_image.h) {
        printf("input mask image and src image not matched.\n");
        return -1;
    }
	printf("get input and read data finished.\n");
    
    // get image boundary process
	std::map<Point, int> boundary_index;
	std::string cmd = "python get_boundary.py "+ mask_filename;
	system(cmd.c_str());
	FILE *boundary_file = fopen("boundary_point.txt", "r");
	int x, y, ans = 0;
	while (~fscanf(boundary_file, "%d%d", &x, &y)) {
		boundary_index[Point(x, y)] = ans ++;

		boundary_point.push_back(Point(x, y));
	}
	fclose(boundary_file);
	printf("get boundary finished.\n");

    // cal weight process
	getAdaptiveMesh(boundary_point);
	// adaptive mesh output test
	/*FILE *mesh_file = fopen("triangle_vertex.txt", "w");
	for(auto iter = triangle_meshes.cbegin(); iter != triangle_meshes.cend(); iter ++) {
		fprintf(mesh_file, "%d %d %d %d %d %d\n", posTransfer((iter->vertex(0)).x()), posTransfer((iter->vertex(0)).y()), posTransfer((iter->vertex(1)).x()), posTransfer((iter->vertex(1)).y()), posTransfer((iter->vertex(2)).x()), posTransfer((iter->vertex(2)).y()));
	}
	fclose(mesh_file);*/
	printf("get adaptive mesh finished.\n");

    // get weight for vertex
	double boundary_diff[ans + 10][4], weight[ans + 10], tan_val[ans + 10], vertex_color[vertex_list.size() + 10][4];
	int channel = std::min(src_image.c, target_image.c);
	memset(vertex_color, 0, sizeof vertex_color);
	for (int i = 0; i < vertex_list.size(); i ++) {
		if (boundary_index.find(vertex_list[i]) != boundary_index.end()) { // boundary point
			int src_x = int(.5 + vertex_list[i].x());
			int src_y = int(.5 + vertex_list[i].y());
			int tar_x = src_x + pos_x, tar_y = src_y + pos_y;
			for (int k = 0; k < channel; k ++) {
				vertex_color[triangle_vertex_index[vertex_list[i]]][k] = boundary_diff[boundary_index[vertex_list[i]]][k] = getpixel(target_image, tar_x, tar_y, k) - getpixel(src_image, src_x, src_y, k);
				//printf("%f\n", getpixel(target_image, tar_x, tar_y, k) - getpixel(src_image, src_x, src_y, k));
			}
		}
	}
	int tmp = 0;
	for (int i = 0; i < vertex_list.size(); i ++) {
		if (boundary_index.find(vertex_list[i]) == boundary_index.end()) { // internal point
			tmp ++;
			int index = triangle_vertex_index[vertex_list[i]];
			for (int j = 0; j < ans; j ++) {
				tan_val[j] = getTanValue(vertex_list[i], boundary_point[j], boundary_point[(j + 1) % ans]);
				//printf("%f\n", tan_val[j]);
			}
			double w_sum = 0;
			for (int j = 0; j < ans; j ++) {
				weight[j] = (tan_val[(j+ans-1)%ans] + tan_val[j]) / getDistance(vertex_list[i], boundary_point[j]);
				w_sum += weight[j];
			}
			//printf("%f\n", w_sum);
			for (int j = 0; j < ans; j ++) {
				double lmd = weight[j] / w_sum;
				//printf("%f\n", lmd);
				for(int k = 0; k < channel; k ++) {
					vertex_color[index][k] += lmd * boundary_diff[j][k];
					//printf("%f %f\n", lmd, boundary_diff[j][k]);
				}
			}
		}
	}
	//printf("%d %d %d\n", ans, vertex_list.size(), tmp);

	// for each triangular mesh, calc the inside points' color
	for (int i = 0; i < triangle_meshes.size(); ++i)
	{
		Point v1 = triangle_meshes[i].vertex(0);
		Point v2 = triangle_meshes[i].vertex(1);
		Point v3 = triangle_meshes[i].vertex(2);
		double x_min = std::min(std::min(v1.x(), v2.x()), v3.x());
		double x_max = std::max(std::max(v1.x(), v2.x()), v3.x());
		double y_min = std::min(std::min(v1.y(), v2.y()), v3.y());
		double y_max = std::max(std::max(v1.y(), v2.y()), v3.y());
		// calc transform matrix
		int index1 = triangle_vertex_index[v1];
		int index2 = triangle_vertex_index[v2];
		int index3 = triangle_vertex_index[v3];
		for (int k = 0; k < channel; k ++) {
			calTransformWeight(t_w[k], v1, v2, v3, vertex_color[index1][k], vertex_color[index2][k], vertex_color[index3][k]);
			/*if(vertex_color[index1][k] < 0) {
				printf("index1 %f\n", vertex_color[index1][k]);
			}
			if(vertex_color[index2][k] < 0) {
				printf("index2 %f\n", vertex_color[index2][k]);
			}
			if(vertex_color[index3][k] < 0) {
				printf("index3 %f\n", vertex_color[index3][k]);
			}*/
		}
		for (int x = std::max(0., x_min - 1); x <= x_max + 1 && x < src_image.h; x ++) {
			for (int y = std::max(0., y_min - 1); y <= y_max + 1 && y < src_image.w; y ++) {
				if (inside(Point(x, y), v1, v2, v3)) {
					for (int k = 0; k < channel; k ++)
					{
						double c = t_w[k].getRes(x, y) + getpixel(src_image, x, y, k);
						/*if(c < 0) {
							printf("%f %f\n", getpixel(src_image, x, y, k), t_w[k].getRes(x, y));
						}*/
						getpixel(target_image, x + pos_x, y + pos_y, k) = c;
					}
				}
			}
		}
	}
	// output
	printf("result stored in %s\n", output_filename.c_str());
	target_image.write(output_filename.c_str());

    printf("All process finished!\n");
    return 0;
}