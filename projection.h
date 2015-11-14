/********* compute the signed distance according to the curvature filter projection **************/
//FilterType should be integer
//input image should be float
//output image (d_m) should be float and it is the signed distance

//four types of filters
inline float Scheme_TV(int i, float* p_pre, float* p, float* p_nex);
inline float Scheme_MC(int i, float* p_pre, float* p, float* p_nex);
inline float Scheme_GC(int i, float* p_pre, float* p, float* p_nex);
inline float Scheme_LS(int i, float* p_pre, float* p, float* p_nex);

void projection(int Type, Mat& imgF, Mat& d_m)
{
    float (*Local)(int i, float* p_pre, float* p, float* p_nex);
    float *p, *p_pre, *p_down, *p_dm;
    switch(Type)
    {
        case 0:
        {
            Local = Scheme_TV; cout<<"TV Filter: "; break;
        }
        case 1:
        {
            Local = Scheme_MC; cout<<"MC Filter: "; break;
        }
        case 2:
        {
            Local = Scheme_GC; cout<<"GC Filter: "; break;
        }
        //DC filter is not included
        case 4:
        {
            Local = Scheme_LS; cout<<"Bernstein Filter: "; break;
        }
        default:
        {
            cout<<"The filter type is wrong. Do nothing."<<endl; return;
        }
    }

    int M = imgF.rows; int N = imgF.cols;
    if (M!=d_m.rows || N!=d_m.cols)
    {
        cout<<"the input and output images do NOT have the same size."; return;
    }
    /********* four domain sets ***************/
    //black circle
    for (int i = 1; i < M-1; ++i,++i)
    {
        p = imgF.ptr<float>(i);
        p_pre = imgF.ptr<float>(i-1);
        p_down = imgF.ptr<float>(i+1);
        p_dm = d_m.ptr<float>(i);
        for (int j = 1; j < N-1; ++j, ++j)
        {
            p_dm[j] = (*Local)(j,p_pre,p,p_down);
        }
    }

    //black triangle
    for (int i = 2; i < M-1; ++i,++i)
    {
        p = imgF.ptr<float>(i);
        p_pre = imgF.ptr<float>(i-1);
        p_down = imgF.ptr<float>(i+1);
        for (int j = 2; j < N-1; ++j, ++j)
        {
            p_dm[j] = (*Local)(j,p_pre,p,p_down);
        }
    }

    //white circle
    for (int i = 1; i < M-1; ++i,++i)
    {
        p = imgF.ptr<float>(i);
        p_pre = imgF.ptr<float>(i-1);
        p_down = imgF.ptr<float>(i+1);
        for (int j = 2; j < N-1; ++j, ++j)
        {
            p_dm[j] = (*Local)(j,p_pre,p,p_down);
        }
    }

    //white triangle
    for (int i = 2; i < M-1; ++i,++i)
    {
        p = imgF.ptr<float>(i);
        p_pre = imgF.ptr<float>(i-1);
        p_down = imgF.ptr<float>(i+1);
        for (int j = 1; j < N-1; ++j, ++j)
        {
            p_dm[j] = (*Local)(j,p_pre,p,p_down);
        }
    }
}

/************************ projection operation for each pixel ******************************************/
inline float Scheme_GC(int i, float * __restrict p_pre, float * __restrict p, float * __restrict p_nex)
{
    float dist[4];
    dist[2] = 2*p[i];
    dist[0] = (p_pre[i] + p_nex[i]) - dist[2];
    dist[1] = (p[i-1] + p[i+1]) - dist[2];
    if (fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];
    
    dist[1] = (p_pre[i-1] + p_nex[i+1]) - dist[2];
    if (fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];
    dist[1]  = (p_nex[i-1] + p_pre[i+1]) - dist[2];
    if (fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];
    dist[0] /= 2;

    dist[3] = 3*p[i];
    dist[1] = (p_pre[i] + p_pre[i-1] + p[i-1])- dist[3];
    
    dist[2] = (p_pre[i] + p_pre[i+1] + p[i+1])- dist[3];
    if (fabsf(dist[2])<fabsf(dist[1])) dist[1] = dist[2];
    dist[2] = (p_nex[i] + p_nex[i-1] + p[i-1])- dist[3];
    if (fabsf(dist[2])<fabsf(dist[1])) dist[1] = dist[2];
    dist[2] = (p_nex[i] + p_nex[i+1] + p[i+1])- dist[3];
    if (fabsf(dist[2])<fabsf(dist[1])) dist[1] = dist[2];
    dist[1] /= 3;

    if (fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];
    
    return dist[0];
}

inline float Scheme_MC(int i, float* p_pre, float* p, float* p_nex)
{
    //       a   b
    //       I   e
    //       c   d
    // return (2.5(a+c)+5*e)-b-d)/8.0;
    float dist[2];
    float tmp = 8*p[i];

    dist[0] = 2.5f*(p_pre[i]+p_nex[i]) + 5*p[i+1] - p_pre[i+1] - p_nex[i+1] - tmp;
    dist[1] = 2.5f*(p_pre[i]+p_nex[i]) + 5*p[i-1] - p_pre[i-1] - p_nex[i-1] - tmp;
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

    dist[1] = 2.5f*(p[i-1]+p[i+1]) + 5*p_nex[i] - p_nex[i-1] - p_nex[i+1] - tmp;
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];
    dist[1] = 2.5f*(p[i-1]+p[i+1]) + 5*p_pre[i]- p_pre[i-1] - p_pre[i+1] - tmp;
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

    dist[0] /= 8;
    
    return dist[0];
}

inline float Scheme_LS(int i, float* p_pre, float* p, float* p_nex)
{
    //   f   a   b            0 1/2 0               3/7 1/7 -1/7
    //       I   e               -1 0                    -1  1/7
    //       c   d              1/2 0                    0   3/7
    // 
    float dist[4];
    float tmp = 2*p[i];

    dist[0] = p_pre[i]+p_nex[i] - tmp;
    dist[1] = p[i-1] + p[i+1] - tmp;
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];
    dist[1] = p_pre[i-1] + p_nex[i+1] - tmp;
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];
    dist[1] = p_nex[i-1] + p_pre[i+1] - tmp;
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

    tmp *= 3.5f;
    dist[0] *= 3.5f;

    dist[2] = 3*(p_pre[i+1] + p_nex[i-1]) - tmp;
    dist[3] = 3*(p_pre[i-1] + p_nex[i+1]) - tmp;

    dist[1] = p_pre[i] - p_pre[i-1] +  p[i-1] + dist[2];
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

    dist[1] = p_pre[i] - p_pre[i+1] + p[i+1] + dist[3];
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

    dist[1] = p_nex[i] - p_nex[i-1] + p[i-1] + dist[3];
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

    dist[1] = p_nex[i]  - p_nex[i+1] + p[i+1] + dist[2];
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

    dist[0] /= 7;

    return dist[0];
}


inline float Scheme_TV(int i, float* p_pre, float* p, float* p_nex)
{
    //       a   b
    //       I   e
    //       c   d
    // return (a+b+c+d+e)/5.0;
    float dist[2];
    float tmp = 5*p[i];

    dist[0] = p_pre[i]+ p_pre[i+1]+p[i+1]+p_nex[i] + p_nex[i+1] - tmp;

    dist[1] = p_pre[i]+p[i-1]+p_nex[i] + p_pre[i-1] + p_nex[i-1] - tmp;
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

    dist[1] = p[i-1]+p[i+1]+p_nex[i] + p_nex[i-1] + p_nex[i+1] - tmp;
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];
    dist[1] = p[i-1]+p[i+1]+p_pre[i] + p_pre[i-1] + p_pre[i+1] - tmp;
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

    //diag
    dist[1] = p_pre[i-1]+p_pre[i] + p_pre[i+1] + p[i-1] + p_nex[i-1]- tmp;
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];
    dist[1] = p_pre[i-1]+p_pre[i] + p_pre[i+1] + p[i+1] + p_nex[i+1]- tmp;
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];
    dist[1] = p_nex[i-1]+p_nex[i] + p_nex[i+1] + p[i-1] + p_pre[i-1]- tmp;
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];
    dist[1] = p_nex[i-1]+p_nex[i] + p_nex[i+1] + p[i+1] + p_pre[i+1]- tmp;
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];


    dist[0]/=5;
    return dist[0];
}