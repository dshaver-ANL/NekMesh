#define max_elem 25000
#define max_mesh 250

typedef struct { double x; double y; double z;} point;
typedef struct { int n; point *p;} connector;
typedef struct { int elid; int esid; double curve; char ccurve;} edge;
typedef struct { int vid[4]; char BC[4][2][4];} quad;
typedef struct { int vid[8];} hex;
typedef struct { connector con[4]; int neiid[4]; char BC[4][2][4];} mesh;

int make_quad_space(int,int,point*,char[4][2][4]);
int make_gquad_space(int,int,double,point*,char[4][2][4]);
int make_g2quad_space(int,int,double,double,point*,char[4][2][4]);
int make_cquad_space(int,int,double,double,double,point*,char[4][2][4]);

int make_tri_space(int,point*,char[4][2][4]);

double distance(point,point);
double growth_ratio(int,double,double);
int get_g_side(point,point,double,connector*);
point line_circle_intercept(point,point,point,double);
point line_line_intercept(point,point,point,point);
point circle_center_2pR(point,point,double);
double get_theta_0(point,point,int,double,double,point);
point rotate_point(point,double,point);
point reflect_point(point,point,point);
point translate_point(point,point);
point midpoint(point,point);
point linpoint(double,point,point);
point centroid(point*,int);
int invert(double*,int);
int invert_connector(connector*);
connector make_lin_side(point,point,int);

int write_rea(char*);
int output_pts(point*,int,char*);

int make_vtri_space(int,int,double*,double*);
int make_ctri_space(int,int,double,double,double*,double*);
