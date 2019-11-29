#if !defined(max)
#define abs(x) ( (x) > 0.0 ? x : -(x) )
#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))
#define sqr(x) (x)*(x)
#endif

#define X(ix,iy) (ix)*iNy+ (iy)

char *replace_fno(char *, char *, char *);
char *fno_to_fstr(int, int);
int count_digits(int);

void GetRegionsAverages(double  *,double  *,int ,int ,double  *,int ,double  *,double  *,double  *,double  *,int );

void ChanVeseSegmentationNew (double *, int , int , int , double *, int , int , double , double , double *, double *, double *, double *, double *);
