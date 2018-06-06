#define max_elem 25000

typedef struct { double x;double y;} point;
typedef struct { int elid; int esid; double curve; char ccurve;} edge;
typedef struct { int vid[4];} quad;
typedef struct { int vid[8];} hex;

int make_quad_space(quad*,point*,int,int,double*,double*,int*,int*);
int make_gquad_space(quad*,point*,int,int,double,double*,double*,int*,int*);
int make_cquad_space(quad*,point*,edge*,int,int,double,double,double*,double*,int*,int*,int*);
int make_tri_space(quad*,point*,int,double*,double*,int*,int*);
int make_vtri_space(quad*,point*,int,int,double*,double*,int*,int*);
int make_ctri_space(quad*,point*,edge*,int,int,double,double,double*,double*,int*,int*);

double growth_ratio(int,double,double);
int get_g_side(double,double,double,double,int,double,double*, double*);
int line_circle_intercept(double,double,double,double,double,double,double,double*,double*);
int circle_center_2pR(double,double,double,double,double,double*,double*);

int write_rea(quad*,point*,edge*,char*);
