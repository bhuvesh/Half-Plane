#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/core/core.hpp"
#include <bits/stdc++.h>


using namespace cv;
using namespace std;


int lf=0,rg=0,INF=10000000;
double error2=0.00000001,err=0.00000001;

class point {
public:
    double px, py;

    bool operator <(const point &p) const {
        return px + err < p.px || ( fabs(px - p.px) < err && py + err < p.py);
    }
    bool operator ==(const point &p) const {
        return fabs(px-p.px)<err && fabs(py-p.py)<err;
    }
};


struct node
{
    point p;
    int pos;
    struct node *left;
    struct node *right;
};


struct hullv
{
    vector<point> top;
    int tops;
    vector<point> bot;
    int bots;
    vector<point> full;
    int s;
    Point *allpts;

};


Mat back;
vector<hullv> oni;


long double cross(const point &O, const point &A, const point &B)
{
    return (A.px - O.px) * (B.py - O.py) - (A.py - O.py) * (B.px - O.px);
}

struct lines{
    double m;
    double c;
};
void reset();

double modu(double a)
{
    if(a<0)
    {
        return -a;
    }
    return a;
}

void printit(point p);

/*
Geometric functions
*/
int intline(vector<point> arr,int st,int n,double ll,double rr,lines specs);
double distline(point p,lines specs);
int miniline(vector<point> arr,int n,lines specs,double r);
void multiint(vector<point> arr,int n,lines specs, double rr,double sig);
void pntadd(vector<point> uhull,int us,vector<point> lhull,int ls,lines specs,double sig);
void findpoints (lines qline, double sign);







struct hullt
{
    node *top;
    node *bot;
};

vector<hullt> layers;

vector<point> arr;

node* fromarr(int start , int last)
{
    //cerr<<"fromarr  start: "<<start<<"  end:  "<<last<<endl;
    if(start ==  last)
    {
        node *top;
        top = (node *) malloc(sizeof(node));
        top->p = arr[start];
        top->left = NULL;
        top->pos=start;
        top->right = NULL;
        return top;
    }
    if (start > last)
        return NULL;

    int mid = (start + last)/2;
    node *top;
    top = (node *) malloc(sizeof(node));
    top->p = arr[mid];
    top->pos = mid;
    top->left = fromarr(start , mid - 1);
    top->right = fromarr(mid + 1, last);
    return top;
}




void make_trees(vector <hullv> oni)
{
    for(int i=0; i<oni.size(); i++)
    {
        hullt h;
        arr = oni[i].top;
        h.top = fromarr(0,oni[i].tops-1);
        arr = oni[i].bot;
        h.bot = fromarr(0, oni[i].bots - 1);

        layers.push_back(h);
    }
}



hullv convex_hull(vector<point> P)
{
    hullv x;
    int n = P.size(), k = 0;
    vector<point> H(2*n);
    for (int i = 0; i < n; ++i) {
        while (k >= 2 && cross(H[k-2], H[k-1], P[i]) + err <= 0) k--;
        H[k++] = P[i];
    }
    x.bot = H;
    x.bot.resize(k);
    int last = k;
    for (int i = n-2, t = k+1; i >= 0; i--) {
        while (k >= t && cross(H[k-2], H[k-1], P[i]) + err <= 0) k--;
        H[k++] = P[i];
    }
    H.resize(k);
    for(int i=k-1; i>=last-1; i--)
    {
        x.top.push_back(H[i]);
    }
    x.bots= x.bot.size();
    x.tops=x.top.size();
    x.full = H;
    x.s = H.size();
    x.allpts = (Point*) malloc (x.s * sizeof(Point));
    for(int i=0;i<x.s;i++)
    {
        double x1 = H[i].px;
        double y1 = H[i].py; 
        x.allpts[i] = Point(x1,y1);
    }
    return x;
}




vector<hullv> onion(vector<point> P)
{
    vector <hullv> oni;
    hullv h;
    while (P.size()>1)
    {
        //cerr<<P.size()<<endl;
        h=convex_hull(P);
        oni.push_back(h);
        vector <point> temp;

        int i,j,k,l,m;

        i=0;
        j=0;

        while(i<h.tops&&j<P.size())
        {
            if (h.top[i] == P[j])
            {
                i++;
                j++;
            }
            else if(h.top[i] < P[j])
            {
                i++;
            }
            else
            {
                temp.push_back(P[j]);
                j++;
            }
        }
        //cerr<<"hello\n";
        if (temp.empty())
        {
            P.clear();
        }
        else P = temp;

        temp.clear();

        //cerr<<"hi\n";

        i=0;
        j=0;
        while(i<h.bots&&j<P.size())
        {
            if (h.bot[i] == P[j])
            {
                i++;
                j++;
            }
            else if(h.bot[i] < P[j])
            {
                i++;
            }
            else
            {
                temp.push_back(P[j]);
                j++;
            }
        }

        //cerr<<"hello again\n";

        if (temp.empty())
        {
            P.clear();
        }
        else P = temp;

        //cerr<<"hi again\n";
        temp.clear();
        if(P.empty()) break;
    }

    if (P.size()==1)
    {
        hullv h;
        h.top.push_back(P[0]);
        h.bot.push_back(P[0]);
        h.full.push_back(P[0]);
        h.bots=1;
        h.tops=1;
        h.s=1;
        oni.push_back(h);
    }

    return oni;

}

void print(node* root)
{
    if (root==NULL)
        return;
    print(root->left);
    cout<<"("<<root->p.px<<","<<root->p.py<<") "<<root->pos<<"  ";
    print(root ->right);

}

void drawPoint(int x, int y)
{
        Point pt =  Point(x, y);
        circle( back, pt, 2, Scalar( 0, 0, 0 ), -1 ,3 );    
}

void drawPoint(int x, int y, Scalar S)
{
        Point pt =  Point(x, y);
        circle( back, pt, 3, S , -1 ,3 );   
}

void drawLine(int x1,int y1, int x2, int y2, Scalar S)
{
    int thickness = 1;
  int lineType = 8;
  Point start = Point(x1,y1);
  Point end = Point(x2,y2);
  line( back, start, end, S, thickness, lineType);
}

void drawLine(int x1,int y1, int x2, int y2)
{
    int thickness = 1;
  int lineType = 8;
  Point start = Point(x1,y1);
  Point end = Point(x2,y2);
  line( back, start, end, Scalar( 0, 0, 0 ), thickness, lineType);
}

vector <point> P;
vector <point> pts;
vector <point> L;
lines qline;
int key, drawn, query=0;
int n;
void CallBackFunc(int event, int x, int y, int flags, void* userdata)
{
    
    if ( event == EVENT_MOUSEMOVE )
     {

        ///cout<<query<<endl;
         if (query==2)
        {

            reset();
            drawPoint(L[0].px,L[0].py, Scalar(255,0,0));
           
            drawPoint(x,y, Scalar(255,0,0));
            //L.push_back(inp);
            if (L[0].px == x)
            {
                qline.m = INF;
            }
            else
                qline.m = ((double)(L[0].py - y))/((double)(L[0].px - x));
            qline.c = L[0].py - qline.m*L[0].px;
            cout<<qline.m <<"  "<<qline.c<<endl;
            int x1,y1,x2,y2;
            if (fabs(qline.m) < 45)
            {

                            x1=-5;
               y1 = int (qline.m*x1 + qline.c);
              x2 = 805;
               y2 = int (qline.m*x2 + qline.c);
            }
            else
            {

                cout<<"here\n";
                y1 = 0; 
                x1 = x;
                y2 = 800;
                x2 = x;             
            }

            //query=3;
           // drawLine(L[0].px, L[0].py, L[1].px ,L[1].py, Scalar (50,50,50));
            drawLine(x1, y1, x2 , y2, Scalar (50,50,50));
        }


     }

    if  ( event == EVENT_LBUTTONDOWN)
    {
        cout<<query<<endl;
        if (key!='d' && drawn==0)
        {
            point inp;
            inp.px=x;
            inp.py=y;
            drawPoint(x,y, Scalar(0,0,255));
            P.push_back(inp);

        }

        if (query==1)
        {   
            L.clear();
            pts.clear();
            
            point inp;
            inp.px=x;
            inp.py=y;
            drawPoint(x,y, Scalar(255,0,0));
            L.push_back(inp);
            query=2;
        }

        else if (query==2)
        {
            reset();
            drawPoint(L[0].px,L[0].py, Scalar(255,0,0));
            point inp;
            inp.px=x;
            inp.py=y;
            drawPoint(x,y, Scalar(255,0,0));
            L.push_back(inp);
            if (L[0].px == L[1].px)
            {
                qline.m = INF;
            }
            else
                qline.m = ((double)(L[0].py - L[1].py))/((double)(L[0].px - L[1].px));
            qline.c = L[0].py - qline.m*L[0].px;

            int x1,y1,x2,y2;
             if (fabs(qline.m) < 45)
            {
                            x1=-5;
               y1 = int (qline.m*x1 + qline.c);
              x2 = 805;
               y2 = int (qline.m*x2 + qline.c);
            }
            else
            {
                y1 =-10 ;
                x1 = x;
                y2 = 800;
                x2 = x;             
            }
            query=3;
           // drawLine(L[0].px, L[0].py, L[1].px ,L[1].py, Scalar (50,50,50));
            drawLine(x1, y1, x2 , y2, Scalar (50,50,50));
        }

        else if (query==3)
        {
            point inp;
            inp.px=x;
            inp.py=y;

            
            
                findpoints(qline , distline(inp, qline));

                for(int i=0;i<(int)pts.size();i++)
                {
                 drawPoint((int)pts[i].px,(int)pts[i].py, Scalar(0,153,0));
                }
           query=0;
        }

    }
}


void reset()
{
    back = Scalar::all(255);
       int k = 125/oni.size();
            for(int i=0; i<oni.size(); i++)
            {
                fillConvexPoly(back, oni[i].allpts, oni[i].s -1 , Scalar (230 - k*i, 230 - k*i , 255) );
            }
    for(int i=0; i<oni.size(); i++)
            {
                for(int j=0; j<oni[i].tops - 1; j++)
                {
                    drawLine((int)oni[i].top[j].px , (int) oni[i].top[j].py , (int) oni[i].top[j+1].px, (int) oni[i].top[j+1].py ,Scalar (246,12,4) );
                }

                for(int j=0; j<oni[i].bots - 1; j++)
                {
                    drawLine((int)oni[i].bot[j].px , (int) oni[i].bot[j].py ,(int) oni[i].bot[j+1].px, (int) oni[i].bot[j+1].py , Scalar (246,12,4) );
                }

            }

            for(int i=0;i<n;i++)
            {
                 drawPoint((int)P[i].px,(int)P[i].py, Scalar(0,0,255));
            }
}

void remove_repeat()
{
    vector <point> tem;
    for(int i=0;i<n;i++)
    {
        while(P[i]==P[i+1])
            i++;
        tem.push_back(P[i]);

    }
    P=tem;
}


int main()
{

    drawn=0;
    key = 'f';
    
 

    cerr<<"onion made "<<oni.size()<<endl;

    back = imread("back.jpg");
    namedWindow("Onion",CV_WINDOW_AUTOSIZE);
    setMouseCallback("Onion", CallBackFunc, NULL);

    while(key!=27)
    {
        imshow("Onion", back);
        key = waitKey(5);

        if (key=='d'&& drawn==0)
        {

            sort(P.begin(), P.end());
            n=P.size();
            remove_repeat();
             n=P.size();
            oni = onion(P);


            for(int i=0; i<oni.size(); i++)
            {
                for(int j=0; j<oni[i].tops; j++)
                    cout<<"("<<oni[i].top[j].px<<","<<oni[i].top[j].py<<")  ";
                cout<<endl;

                for(int j=0; j<oni[i].bots; j++)
                    cout<<"("<<oni[i].bot[j].px<<","<<oni[i].bot[j].py<<")  ";
                cout<<endl;

                for(int j=0; j<oni[i].s;j++)
                    cout<<oni[i].allpts[j].x<<","<<oni[i].allpts[j].y<<"    ";
                cout<<endl;
                cout<<endl;
            }
            int k = 125/oni.size();
            for(int i=0; i<oni.size(); i++)
            {
                fillConvexPoly(back, oni[i].allpts, oni[i].s -1 , Scalar (250, 230 - k*i , 230 - k*i) );
            }
            

            for(int i=0; i<oni.size(); i++)
            {
                for(int j=0; j<oni[i].tops - 1; j++)
                {
                    drawLine((int)oni[i].top[j].px , (int) oni[i].top[j].py , (int) oni[i].top[j+1].px, (int) oni[i].top[j+1].py , Scalar (246,12,4));
                }

                for(int j=0; j<oni[i].bots - 1; j++)
                {
                    drawLine((int)oni[i].bot[j].px , (int) oni[i].bot[j].py ,(int) oni[i].bot[j+1].px, (int) oni[i].bot[j+1].py, Scalar (246,12,4) );
                }

            }

            for(int i=0;i<n;i++)
            {
                 drawPoint((int)P[i].px,(int)P[i].py, Scalar(0,0,255));
            }
            drawn=1;

        }

        if(key=='n' && drawn==1 && query==0)
        {
            reset();
            query=1;
            cout<<key<<endl;
        }

    }
    
    




  
    make_trees(oni);
    for(int i=0; i<layers.size(); i++)
    {
        print(layers[i].top);
        cout<<endl;
        print(layers[i].bot);
        cout<<endl<<endl;
    }






}

void printit(point p)
{
    printf("point: %lf %lf\n",p.px,p.py);
}

/*
Geometric functions
*/

int intline(vector<point> arr,int st,int n,double ll,double rr,lines specs)
{
    int l=st,r=n-1,m,i;
    printf("Initial values are : %lf %lf\n",ll,rr);
    double a,b,c;
    while(l<=r)
    {
        m=((l+r)/2);            
        a=distline(arr[m],specs);
        printf("point: %d %lf %lf %lf\n",m,arr[m].px,arr[m].py,a);
        if(modu(a)<=error2)
        {
            return m;
        }
        if((a*ll)>0)
        {
            l=m+1;  
        }
        else if((a*rr)>0)
        {
            r=m-1;
        }
        else
        {
            return m;
        }
    }
    m=((l+r)/2);
    return m;
}

double distline(point p,lines specs)
{
    double x,y,pp,qq,tm,tc;
    x=p.px;
    y=p.py;
    tm=specs.m;
    tc=specs.c;
    pp=(y-(tm*x)-tc);
    qq=sqrt(1+(tm*tm));
    return (pp/qq);
}

int miniline(vector<point> arr,int n,lines specs,double rr)
{
    int l=0,r=n-1,m;
    double a,b,c;
    while(l<=r)
    {
        m=((l+r)/2);
        a=distline(arr[m],specs);
        a*=rr;
        if((m-1)<0)
        {
            c=INF;
        }
        else
        {
            c=distline(arr[m-1],specs);
            c*=rr;
        }
        
        if((m+1)>=n)
        {
            b=INF;
        }
        else
        {
            b=distline(arr[m+1],specs);
            b*=rr;
        }
        

        printf("point: %lf %lf %lf\n",arr[m].px,arr[m].py,a);
        
        /*if a is minima: return a*/    
        if((a<=c)&&(a<=b))
        {
            return m;
        }
            
        if(a>b)
        {
            l=m+1;
        }
        else if(a<b)
        {
            r=m-1;
        }
        else
        {
            return m;
        }
    }
    return -1;
}
    
void multiint(vector<point> arr,int n,lines specs,double rr,double sig)
{
    int a,b,c,i,j;
    double l,r;
    point tmp;  
    a=miniline(arr,n,specs,rr);
    printf("the values are a: %lf rr:%lf sig:%lf\n",distline(arr[a],specs),rr,sig);
    if((distline(arr[a],specs)*rr)>0)
    {
        if((distline(arr[a],specs)*sig)>0)
        {
            //add arr[0] to arr[n-1];
            for(i=0;i<n;i++)
            {
                pts.push_back(arr[i]);
            }   
        }
    }   
    else
    {
        printf("a is %d\n",a);
        //first point of intersection
        l=distline(arr[0],specs);
        r=distline(arr[a],specs);
        b=intline(arr,0,a+1,l,r,specs);
        //second    
        l=distline(arr[a],specs);
        r=distline(arr[n-1],specs);
        c=intline(arr,a,n,l,r,specs);
        printf("%d %d\n",b,c);
        if((distline(arr[c],specs)*sig)>0)
        {
            //add arr[b+1] to arr[c];
            for(i=b+1;i<=c;i++)
            {
                pts.push_back(arr[i]);
            }       
        }
        else
        {
            //add arr[0] to arr[b];
            for(i=0;i<=b;i++)
            {
                pts.push_back(arr[i]);
            }
            //add arr[a+c+1] to arr[n-1];
            for(i=c+1;i<n;i++)
            {
                pts.push_back(arr[i]);
            }       
        }
    }
    //return b;
}

void pntadd(vector<point> uhull,int us,vector<point> lhull,int ls,lines specs,double sig)
{
    double m,c,l,r;
    int t,i;
    
    l=distline(uhull[0],specs);
    r=distline(uhull[us-1],specs);
    
    printf("%lf %lf\n",l,r);
    
    if((l*r)<0)
    {
        t=intline(uhull,0,us,l,r,specs);
        printf("%d ,Good also\n",t);    
        if((distline(uhull[t],specs)*sig)>0)
        {
            //add arr[0] to arr[t];
            for(i=0;i<=t;i++)
            {
                pts.push_back(uhull[i]);
            }
        }
        else
        {
            //add arr[t+1] to arr[us-1];
            for(i=t+1;i<us;i++)
            {
                pts.push_back(uhull[i]);
            }
        }
    }
    else
    {
        multiint(uhull,us,specs,r,sig);
    }

    l=distline(lhull[0],specs);
    r=distline(lhull[ls-1],specs);
    
    printf("%lf %lf\n",l,r);
    
    if((l*r)<0)
    {
        t=intline(lhull,0,ls,l,r,specs);
        printf("%d ,Good also\n",t);
        if((distline(lhull[t],specs)*sig)>0)
        {
            //add arr[0] to arr[t];
            for(i=0;i<=t;i++)
            {
                pts.push_back(lhull[i]);
            }       
        }
        else
        {
            //add arr[t+1] to arr[ls-1];
            for(i=t+1;i<ls;i++)
            {
                pts.push_back(lhull[i]);
            }
        }   
    }
    else
    {
        multiint(lhull,ls,specs,r,sig);
    }
}





void findpoints(lines specs, double sig)
{
    for(int i=0;i<oni.size();i++)
    {
        //double m,c,l,r,sig;
        pntadd(oni[i].top,oni[i].tops,oni[i].bot,oni[i].bots,specs,sig);        
    }
}