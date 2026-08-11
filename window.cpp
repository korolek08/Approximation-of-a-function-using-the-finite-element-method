
#include <QPainter>
#include <stdio.h>

#include "window.h"

#define DEFAULT_A -10
#define DEFAULT_B 10
#define DEFAULT_N 10
#define LEN 5000
#define F(I, J) (f(a + (I)*hx, c + (J)*hy))
#define FF(IS, JS, S) (IA_ij(nx, ny, hx, hy, i, j, (IS), (JS), (S), I, A))


#include<cmath>
double f0(double x, double y) {
    (void) x, (void) y;
	return 1;
}
double f1(double x, double y) {
    (void) y;
	return x;
}
double f2(double x, double y) {
    (void) x;
	return y;
}
double f3(double x, double y) {
	return x + y;
}
double f4(double x, double y) {
	return sqrt(x*x + y*y);
}
double f5(double x, double y) {
	return (x*x + y*y);
}
double f6(double x, double y) {
	return exp(x*x - y*y);
}
double f7(double x, double y) {
	return 1/(25*(x*x + y*y) + 1);
}




//thread


double get_cpu_time() {
	struct rusage buf;  
	getrusage(RUSAGE_THREAD, &buf);
	return buf.ru_utime.tv_sec + buf.ru_utime.tv_usec * 1.e-6;
}

double get_full_time() { 
	struct timeval buf;
	gettimeofday(&buf, 0); 
	return buf.tv_sec + buf.tv_usec * 1.e-6;
}

void reduce_sum(int p, double* a, int n) {
	static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
	static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
	static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
	static int t_in=0;
	static int t_out=0;
	static double *r=nullptr;
	int i=0;
	if (p<1) return; //?
	pthread_mutex_lock(&m);
	if (r==nullptr) r=a;
	else for (i=0; i<n; i++) r[i]+=a[i];
	t_in++;
	if (t_in>=p) {
		t_out=0;
		pthread_cond_broadcast(&c_in);
	} else while (t_in<p) pthread_cond_wait(&c_in, &m);
	
	if (r!=a) for (i=0; i<n; i++) a[i]=r[i];
	t_out++;
	if (t_out>=p) {
		t_in=0;
		r=nullptr;
		pthread_cond_broadcast(&c_out); //?
	} else while (t_out<p) pthread_cond_wait(&c_out, &m);
	pthread_mutex_unlock(&m);
}
void reduce_sum_int(int p, int* a, int n) {
	static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
	static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
	static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
	static int t_in=0;
	static int t_out=0;
	static int *r=nullptr;
	int i=0;
	if (p<1) return; //?
	pthread_mutex_lock(&m);
	if (r==nullptr) r=a;
	else for (i=0; i<n; i++) r[i]+=a[i];
	t_in++;
	if (t_in>=p) {
		t_out=0;
		pthread_cond_broadcast(&c_in);
	} else while (t_in<p) pthread_cond_wait(&c_in, &m);
	
	if (r!=a) for (i=0; i<n; i++) a[i]=r[i];
	t_out++;
	if (t_out>=p) {
		t_in=0;
		r=nullptr;
		pthread_cond_broadcast(&c_out); //?
	} else while (t_out<p) pthread_cond_wait(&c_out, &m);
	pthread_mutex_unlock(&m);
}

void reduce_max(int p, double* a, double* mx) {
	static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
	static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
	static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
	static int t_in=0;
	static int t_out=0;
	if (p<1) return; //?
	pthread_mutex_lock(&m);

	t_in++;
	if (t_in>=p) {
		t_out=0;
		pthread_cond_broadcast(&c_in);
	} else while (t_in<p) pthread_cond_wait(&c_in, &m);
	
	double t=0;
	for (int i=0; i<p; i++) {
	    t=a[i];
	    if (*mx < t) {
	        *mx=t;
	    }
	}

	t_out++;
	if (t_out>=p) {
		t_in=0;
		pthread_cond_broadcast(&c_out); //?
	} else while (t_out<p) pthread_cond_wait(&c_out, &m);
	pthread_mutex_unlock(&m);
}




int process_args(Args* a) {
	int t=0;
	for (int i=0; i<a->p; i++) {
		if (a[i].error_type!=io_status::SUCCESS) {
			//printf("Error in file %s\n", a[i].name);
			t=1;
		}
	}
	if (t==1) return 1;
	return 0;
}


//thread

//msr


void ij2l(int nx, int /* ny */, int i, int j, int &l) {
    l = i + j*(nx + 1);
}

void l2ij(int nx, int /* ny */, int &i, int &j, int l) {
    j = l / (nx + 1);
    i = l - j*(nx + 1);
}

int get_len_msr(int nx, int ny) {
    return 6*(nx - 1)*(ny - 1) + 4*(2*(nx - 1) + 2*(ny - 1)) + 3*2 + 2*2;
}

double Pf (double *f, double x0, double y0, int nx, int ny, double hx, double hy, double a, double c) {
    int i = (x0 - a) / hx;
    int j = (y0 - c) / hy;
    int l;
    ij2l (nx, ny, i, j, l);
    if ((y0 - c)/hy - j > (x0 - a)/hx - i) 
        return -f[l]*((y0 - c)/hy - (j + 1)) + f[l+nx+1]*(- (x0 - a)/hx + i + (y0 - c)/hy - j) + f[l+1+nx+1]*((x0 - a)/hx - i);
    else 
        return -f[l]*((x0 - a)/hx - (i + 1)) + f[l+1]*((x0 - a)/hx - i - (y0 - c)/hy + j) + f[l+1+nx+1]*((y0 - c)/hy - j);
}

void thread_rows(int n, int p, int k, int& i1, int& i2) {
    i1 = n*k; i1/=p;
    i2 = n*(k+1); i2/=p;
}

void matrix_mult_vector_msr(double* A, int* I, int n, double* x, double* y, int p, int k) {
    int i, j, l, J;
    int i1, i2;
    double s;
    thread_rows(n, p, k, i1, i2);
    for (i = i1; i < i2; i++) {
        s = A[i]*x[i];
        l = I[i + 1] - I[i];
        J = I[i];
        for (j = 0; j < l; j++)
            s += A[J + j]*x[I[J + j]];
        y[i] = s;
    }
    reduce_sum(p, 0, 0);
}

void mult_sub_vector(int n, double* a, double* b, double coeff, int p, int k) {
    int i, i1, i2;
    thread_rows(n, p, k, i1, i2);
    for (i = i1; i < i2; i++) {
        a[i] -= coeff*b[i];
    }
    reduce_sum(p, 0, 0);
}




int iter = 0;
void apply_precondition_msr_matrix(int n, double* A, int* /* I */, double* v, double* r, int p, int k) {
    int i, i1, i2;
    thread_rows(n, p, k, i1, i2);
    //Jacobi
    
    for (i = i1; i < i2; i++) {
        //if (iter <= p) printf("%e %e\n", r[i], A[i]);
        v[i] = r[i] / A[i];
    }
    iter++;
    reduce_sum(p, 0, 0);
}

double scalar_product(int n, double* x, double* y, int p, int k) {
    int i, i1, i2;
    double s = 0;
    thread_rows(n, p, k, i1, i2);
    for (i = i1; i < i2; i++) s+=x[i]*y[i];
    reduce_sum(p, &s, 1);
    return s;
}





int IA_ij(int nx, int ny, double hx, double hy, int i, int j, int is, int js, int s, int* I, double* A) {
    int l, ls;
    double Sa = hx*hy;
    ij2l(nx, ny, i, j, l);
    ij2l(nx, ny, is, js, ls);
    if (I) I[s] = ls;
    if (A) {
        
        if (i > 0 && i < nx && j > 0 && j < ny) {
            if (l == ls)
                A[s] = 6 * 1/12. * Sa;
            else 
                A[s] = 2 * 1/24. * Sa;
        }
        if (j==0 && i > 0 && i < nx) {
            if (l == ls)
                A[s] = 3 * 1/12. * Sa;
            else if (s==0 || s==1)
                A[s] = 1 * 1/24. * Sa;
            else if (s==2 || s==3)
                A[s] = 2 * 1/24. * Sa;
            else {
                printf("s %d 1\n", s);
                return 1;
            }
        }
        if (j==ny && i > 0 && i < nx) {
            if (l == ls)
                A[s] = 3 * 1/12. * Sa;
            else if (s==0 || s==3)
                A[s] = 1 * 1/24. * Sa;
            else if (s==1 || s==2)
                A[s] = 2 * 1/24. * Sa;
            else {
                printf("s %d 2\n", s);
                return 1;
            }
        }
        if (i==0 && j > 0 && j < ny) {
            if (l == ls)
                A[s] = 3 * 1/12. * Sa;
            else if (s==0 || s==3)
                A[s] = 2 * 1/24. * Sa;
            else if (s==1 || s==2)
                A[s] = 1 * 1/24. * Sa;
            else {
                printf("s %d 3\n", s);
                return 1;
            }
        }
        if (i==nx && j > 0 && j < ny) {
            if (l == ls)
                A[s] = 3 * 1/12. * Sa;
            else if (s==0 || s==3)
                A[s] = 1 * 1/24. * Sa;
            else if (s==1 || s==2)
                A[s] = 2 * 1/24. * Sa;
            else {
                printf("s %d 4\n", s);
                return 1;
            }
        }
        if (i==0 && j==0) {
             if (l == ls)
                A[s] = 2 * 1/12. * Sa;
            else if (s==0 || s==1)
                A[s] = 1 * 1/24. * Sa;
            else if (s==2)
                A[s] = 2 * 1/24. * Sa;
            else {
                printf("s %d 5\n", s);
                return 1;
            }
        }
        if (i==nx && j==ny) {
             if (l == ls)
                A[s] = 2 * 1/12. * Sa;
            else if (s==0 || s==2)
                A[s] = 1 * 1/24. * Sa;
            else if (s==1)
                A[s] = 2 * 1/24. * Sa;
            else {
                printf("s %d 5\n", s);
                return 1;
            }
        }
        
        if ((i==0 && j==ny) || (i==nx && j==0)) {
             if (l == ls)
                A[s] = 1 * 1/12. * Sa;
            else if (s==0 || s==1)
                A[s] = 1 * 1/24. * Sa;
            else {
                printf("s %d 6\n", s);
                return 1;
            }
        }
        
        
        
    }
    return 0;
}


double F_ij(int nx, int ny, double hx, double hy, double a, double c, int i, int j, double (*f)(double x, double y)) {
    double w = hx*hy/192.;
    if (i > 0 && i < nx && j > 0 && j < ny)
        return w*(36*F(i, j) + 20*(F(i + 0.5, j) + F(i, j - 0.5) + F(i - 0.5, j - 0.5) + F(i - 0.5, j) + F(i, j + 0.5) + F(i + 0.5, j + 0.5)) + 4*(F(i + 0.5, j - 0.5) + F(i - 0.5, j - 1) + F(i - 1,  j - 0.5) + F(i - 0.5, j + 0.5) + F(i + 0.5, j + 1) + F(i + 1, j + 0.5)) + 2*(F(i + 1, j) + F(i, j - 1) + F(i - 1, j - 1) + F(i - 1, j) + F(i, j + 1) + F(i + 1, j + 1)));
    if (i > 0 && i < nx && j == 0)
        return w*(18*F(i, j) + 10*(F(i + 0.5, j) + F(i - 0.5, j)) + 20*(F(i, j + 0.5) + F(i + 0.5, j + 0.5)) + 4*(F(i - 0.5, j + 0.5) + F(i + 0.5, j + 1) + F(i + 1, j + 0.5)) + F(i - 1, j) + F(i + 1, j) + 2*(F(i, j + 1) + F(i + 1, j + 1)));
    if (i > 0 && i < nx && j == ny)
        return w*(18*F(i, j) + 10*(F(i + 0.5, j) + F(i - 0.5, j)) + 20*(F(i - 0.5, j - 0.5) + F(i, j - 0.5)) + 4*(F(i - 1, j - 0.5) + F(i - 0.5, j - 1) + F(i + 0.5, j - 0.5)) + F(i - 1, j) + F(i + 1, j) + 2*(F(i - 1, j - 1) + F(i, j - 1)));
    if (j > 0 && j < ny && i == 0)
        return w*(18*F(i, j) + 10*(F(i, j - 0.5) + F(i, j + 0.5)) + 20*(F(i + 0.5, j) + F(i + 0.5, j + 0.5)) + 4*(F(i + 0.5, j - 0.5) + F(i + 0.5, j + 1) + F(i + 1, j + 0.5)) + F(i, j - 1) + F(i, j + 1) + 2*(F(i + 1, j) + F(i + 1, j + 1)));
    if (j > 0 && j < ny && i == nx)
        return w*(18*F(i, j) + 10*(F(i, j - 0.5) + F(i, j + 0.5)) + 20*(F(i - 0.5, j - 0.5) + F(i - 0.5, j)) + 4*(F(i - 0.5, j - 1) + F(i - 1, j - 0.5) + F(i - 0.5, j + 0.5)) + F(i, j - 1) + F(i, j + 1) + 2*(F(i - 1, j - 1) + F(i - 1, j)));
    //corners
    if (i == 0 && j == 0)
        return w*(12*F(i, j) + 10*(F(i + 0.5, j) + F(i, j + 0.5)) + 20*(F(i + 0.5, j + 0.5)) + 4*(F(i + 0.5, j + 1) + F(i + 1, j + 0.5)) + F(i + 1, j) + F(i, j + 1) + 2*(F(i + 1, j + 1)));
    if (i == nx && j == ny)
        return w*(12*F(i, j) + 10*(F(i, j - 0.5) + F(i - 0.5, j)) + 20*(F(i - 0.5, j - 0.5)) + 4*(F(i - 1, j - 0.5) + F(i - 0.5, j - 1)) + F(i, j - 1) + F(i - 1, j) + 2*(F(i - 1, j - 1)));
    //other corners
    if (i == 0 && j == ny)
        return w*(6*F(i, j) + 10*(F(i + 0.5, j) + F(i, j - 0.5)) + 4*(F(i + 0.5, j - 0.5)) + F(i + 1, j) + F(i, j - 1));
    if (i == nx && j == 0)
        return w*(6*F(i, j) + 10*(F(i - 0.5, j) + F(i, j + 0.5)) + 4*(F(i - 0.5, j + 0.5)) + F(i - 1, j) + F(i, j + 1));
    return 1e308;
}


int get_off_diag(int nx, int ny, double hx, double hy, int i, int j, int* I = nullptr, double* A = nullptr) {
    int s = 0;
    if (i < nx) {
        FF(i + 1, j, s);
        s++;
    }
    if (j > 0) {
        FF(i, j - 1, s);
        s++;
    }
    if (i > 0 && j > 0) {
        FF(i - 1, j - 1, s);
        s++;
    }
    if (i > 0) {
        FF(i - 1, j, s);
        s++;
    }
    if (j < ny) {
        FF(i, j + 1, s);
        s++;
    }
    if (i < nx && j < ny) {
        FF(i + 1, j + 1, s);
        s++;
    }
    return s;
}


int get_len_msr_off_diag(int nx, int ny) {
    double hx = 0, hy = 0;
    int i, j, res = 0;
    for (i = 0; i < nx; i++) 
        for (j = 0; j < ny; j++) 
            res += get_off_diag(nx, ny, hx, hy, i, j);
    return res;
}



void fill_I(int nx, int ny, int* I) {
    int i, j, l, N = (nx + 1)*(ny + 1);
    int r = N + 1;
    double hx = 0, hy = 0;
    for (l = 0; l < N; l++) {
        l2ij(nx, ny, i, j, l);
        int s = get_off_diag(nx, ny, hx, hy, i, j);
        I[l] = r;
        r += s;
    }
    I[l] = r;
}

int get_diag(int nx, int ny, double hx, double hy, int i, int j, int* /* I */, double* A) { //
    int l;
    ij2l (nx, ny, i, j, l);
    return IA_ij(nx, ny, hx, hy, i, j, i, j, 0, nullptr, A);
}



int fill_IA(int N, int nx, int ny, double hx, double hy, int* I, double* A, int p, int k) {
    int i, j, l, l1, l2, r, s, t;
    int err = 0, len = 0;
    thread_rows(N, p, k, l1, l2);
    for (l = l1; l < l2; l++) {
        r = I[l];
        s = I[l + 1] - I[l];
        l2ij(nx, ny, i, j, l);
        if (get_diag(nx, ny, hx, hy, i, j, I, A + l) != 0) {
            err = -1;
            break;
        }
        t = get_off_diag(nx, ny, hx, hy, i, j, I + r, A + r);
        if (t != s) {
            err = -1;
            break;
        }
        len += s;
    }
    reduce_sum_int(p, &err, 1);
    if (err < 0) return -1;
    reduce_sum_int(p, &len, 1);
    if (I[N] != N + 1 + len) return -2;
    return 0;
}

void fill_b(int n, int nx, int ny, double hx, double hy, double* b, double a, double c, int p, int k, double (*f)(double x, double y)) {
    int l, l1, l2;
    thread_rows(n, p, k, l1, l2);
    int i, j;
    for (l = l1; l < l2; l++) {
        l2ij (nx, ny, i, j, l);
        b[l] = F_ij(nx, ny, hx, hy, a, c, i, j, f);
    }
    reduce_sum(p, 0, 0);
}

void print_matrix_old(double* a, int n, int m, int p) {
	int np=(n>p ? p : n);
	int mp=(m>p ? p : m);
	int i, j;
	for (i=0; i<np; i++){
		for (j=0; j<mp; j++){
			printf(" %10.3e", a[i*m+j]);
		}
		printf("\n");
	}
}

int minimal_residual_msr_matrix(int n, double* A, int* I, double* b, double* x, double* r, double* u, double* v, double eps, int maxit, int p, int k) {
    double prec, b_norm2, tau, c1, c2;
    int it;
    b_norm2 = scalar_product(n, b, b, p, k);
    //reduce_sum
    prec = b_norm2*eps * eps;
    matrix_mult_vector_msr(A, I, n, x, r, p, k);
    //reduce_sum
    //???
    mult_sub_vector(n, r, b, 1., p, k);
    //r-=1.*b
    //reduce
    
    //print_matrix_old(r, n, 1, n);
    for (it = 0; it < maxit; it++) {
        apply_precondition_msr_matrix(n, A, I, v, r, p, k);
        //if (it == 0) print_matrix_old(v, n, 1, n);
        //reduce_sum
        matrix_mult_vector_msr(A, I, n, v, u, p, k);
        //reduce_sum
        
        c1 = scalar_product(n, u, r, p, k); //ur
        c2 = scalar_product(n, u, u, p, k); //uu
        //printf("%e %e %e \n", c1, c2, prec);
        if (fabs(c1) < prec || fabs(c2) < prec) break;
        tau = c1/c2;
        // x -=tau v
        mult_sub_vector(n, x, v, tau, p, k);
        // r -=tau u
        
        mult_sub_vector(n, r, u, tau, p, k);
    }
    
    if (it > maxit) return -1;
    return it;
}


int minimal_residual_msr_matrix_full(int n, double* A, int* I, double* b, double* x, double* r, double* u, double* v, double eps, int maxit, int maxstep, int p, int k) {
    int step, ret, its=0;
    for (step = 0; step < maxstep; step++) {
        ret = minimal_residual_msr_matrix(n, A, I, b, x, r, u, v, eps, maxit, p, k);
        if (ret >= 0) {
            its+=ret;
            break;
        }
        its+=maxit;
    }
    if (step >= maxstep) return -1;
    return its;
}

double res_1(int n, int nx, int ny, double hx, double hy, double a, double c, int p, int k, double* values, double* r1, double (*f)(double x, double y)) {
    int i, j;
    int l, l1, l2;
    double f_value;
    double value;
    thread_rows(n, p, k, l1, l2);
    double mx = 0, r;
    for (l = l1; l < l2; l++) {
        l2ij(nx, ny, i, j, l);
        if (i < nx && j < ny) {
            f_value = f(a + i*hx + hx/3, c + j*hy + 2*hy/3);
            value = Pf(values, a + i*hx + hx/3, c + j*hy + 2*hy/3, nx, ny, hx, hy, a, c);
            r = fabs(f_value - value);
            if (r > mx) mx = r;
            
            f_value = f(a + i*hx + 2*hx/3, c + j*hy + hy/3);
            value = Pf(values, a + i*hx + 2*hx/3, c + j*hy + hy/3, nx, ny, hx, hy, a, c);
            r = fabs(f_value - value);
            if (r > mx) mx = r;
        }
    }
    r1[k] = mx;
    mx = 0;
    reduce_max(p, r1, &mx);
    r1[k] = mx;
    //printf("%e \n", mx);
    return mx;
}

double res_2(int n, int nx, int ny, double hx, double hy, double a, double c, int p, int k, double* values, double* r2, double (*f)(double x, double y)) {
    int i, j;
    int l, l1, l2;
    double f_value;
    double value;
    thread_rows(n, p, k, l1, l2);
    double sum = 0, r;
    for (l = l1; l < l2; l++) {
        l2ij(nx, ny, i, j, l);
        if (i < nx && j < ny) {
            f_value = f(a + i*hx + hx/3, c + j*hy + 2*hy/3);
            value = Pf(values, a + i*hx + hx/3, c + j*hy + 2*hy/3, nx, ny, hx, hy, a, c);
            r = fabs(f_value - value);
            sum += r;
            
            f_value = f(a + i*hx + 2*hx/3, c + j*hy + hy/3);
            value = Pf(values, a + i*hx + 2*hx/3, c + j*hy + hy/3, nx, ny, hx, hy, a, c);
            r = fabs(f_value - value);
            sum += r;
        }
    }
    sum*=hx*hy/2;
    reduce_sum(p, &sum, 1);
    r2[k] = sum;
    return sum;
}

double res_3(int n, int nx, int ny, double hx, double hy, double a, double c, int p, int k, double* values, double* r3, double (*f)(double x, double y)) {
    int i, j;
    int l, l1, l2;
    double f_value;
    thread_rows(n, p, k, l1, l2);
    double mx = 0, r;
    for (l = l1; l < l2; l++) {
        l2ij(nx, ny, i, j, l);
        if (i < nx && j < ny) {
            f_value = f(a + i*hx, c + j*hy);
            r = fabs(f_value - values[l]);
            if (r > mx) mx = r;
        }
    }
    r3[k] = mx;
    mx = 0;
    reduce_max(p, r3, &mx);
    r3[k] = mx;
    return mx;
}

double res_4(int n, int nx, int ny, double hx, double hy, double a, double c, int p, int k, double* values, double* r4, double (*f)(double x, double y)) {
    int i, j;
    int l, l1, l2;
    double f_value;
    thread_rows(n, p, k, l1, l2);
    double sum = 0, r;
    for (l = l1; l < l2; l++) {
        l2ij(nx, ny, i, j, l);
        if (i < nx && j < ny) {
            f_value = f(a + i*hx, c + j*hy);
            r = fabs(f_value - values[l]);
            sum += r;
        }
    }
    sum*=hx*hy;
    reduce_sum(p, &sum, 1);
    r4[k] = sum;
    return sum;
}

void theArgs::deleting_mem() {
    if (A) delete [] A;
    if (I) delete [] I;
    if (bb) delete [] bb;
    if (x) delete [] x;
    if (r) delete [] r;
    if (u) delete [] u;
    if (v) delete [] v;
}
int theArgs::initialize() {
    deleting_mem();
    N = (nx + 1)*(ny + 1);
    len = get_len_msr(nx, ny) + N + 1;

    A = new double[len];
    I = new int[len];
    bb = new double[N];
    u = new double[N];
    v = new double[N];
    x = new double[N];
    r = new double[N];
    if (!A || !I || !bb || !u || !v || !x || !r) {
        printf ("ERROR!\n");
        deleting_mem();
        return -1;
    }
    memset (A, 0, len*sizeof (double));
    memset (I, 0, len*sizeof (int));
    return 0;
}

void Window::thread_clean() {
    if (theAr->task != task_status::WAIT) {
        pthread_mutex_lock (&theAr->m);
        theAr->deleting = 1;
        pthread_cond_broadcast (&theAr->wait_condition);
        pthread_mutex_unlock (&theAr->m);

        for(int k = 1; k < p; k++) pthread_join(ar[k].tid, 0);
    }
}

void Window::init_func() {
    switch (k) {
        case 0:
            f_name = "f(x) = 1";
            f = f0;
            break;
        case 1:
            f_name = "f(x) = x";
            f = f1;
            break;
        case 2:
            f_name = "f(x) = y";
            f = f2;
            break;
        case 3:
            f_name = "f(x) = x + y";
            f = f3;
            break;
        case 4:
            f_name = "f(x) = sqrt(x*x + y*y)";
            f = f4;
            break;
        case 5:
            f_name = "f(x) = (x*x + y*y)";
            f = f5;
            break;
        case 6:
            f_name = "f(x) = exp(x*x - y*y)";
            f = f6;
            break;
        case 7:
            f_name = "f(x) = 1/(25*(x*x + y*y) + 1)";
            f = f7;
            break;
    }
}

Window::Window (QWidget *parent) : QWidget (parent) {
  a = DEFAULT_A;
  b = DEFAULT_B;
  c = DEFAULT_A;
  d = DEFAULT_B;
  nx = DEFAULT_N;
  ny = DEFAULT_N;
  mx = DEFAULT_N;
  my = DEFAULT_N;
  N = (nx + 1)*(ny + 1);
  //k = 0;

  init_func ();
}

QSize Window::minimumSizeHint() const {
  return QSize (100, 100);
}

QSize Window::sizeHint() const {
  return QSize (1000, 1000);
}



int Window::commandLine(int argc, char **argv)
{
 //double r1, r2;
  if (!((argc==13) && (sscanf(argv[1], "%lf", &a)==1) && (sscanf(argv[2], "%lf", &b)==1) && (sscanf(argv[3], "%lf", &c)==1) && (sscanf(argv[4], "%lf", &d)==1) && (sscanf(argv[5], "%d", &nx)==1) && (sscanf(argv[6], "%d", &ny)==1) && (sscanf(argv[7], "%d", &mx)==1) && (sscanf(argv[8], "%d", &my)==1) && (sscanf(argv[9], "%d", &k)==1) && (sscanf(argv[10], "%lf", &eps)==1) && (sscanf(argv[11], "%d", &maxit)==1) && (sscanf(argv[12], "%d", &p)==1) && p > 0 && k<=7 && k>=0 && nx >0 && ny >0 && mx>0 && my>0 && eps > 0)){
		printf("Usage: %s a b c d nx ny mx my k eps maxit p\n", argv[0]);
		return -2;
	}
  //init_mem ();
  //Building_constructions();
  bool ok;
  a = QString(argv[1]).toDouble(&ok);
  b = QString(argv[2]).toDouble(&ok);
  c = QString(argv[3]).toDouble(&ok);
  d = QString(argv[4]).toDouble(&ok);
  eps = QString(argv[10]).toDouble(&ok);
  if (!((a < b) && (c < d) && (eps > 0))) return -2;
  //printf("CD commandLine %e %e %e %e %d %d %d %d %e\n", a, b, c, d, nx, ny, mx, my, eps);
  N = (nx + 1)*(ny + 1);
  len = get_len_msr(nx, ny) + N + 1;
  init_func();
  return 0;
}


void Window::add_param(Args* ar) {
    for (int kk = 0; kk < p; kk++) {
        ar[kk].theAr = theAr;
        ar[kk].eps = eps;
        ar[kk].maxit = maxit;
        ar[kk].maxstep = maxstep;
        ar[kk].N = N;
        ar[kk].len = theAr->len;
        ar[kk].r1 = r1;
        ar[kk].r2 = r2;
        ar[kk].r3 = r3;
        ar[kk].r4 = r4;
    }
}

void Window::init_r() {
    r1 = new double[p];
    r2 = new double[p];
    r3 = new double[p];
    r4 = new double[p];
}

void Window::delete_mem() {
    if (theAr) delete theAr;
    if (timer) delete timer;
    if (x) delete [] x;
}


void Window::init_mem() {
    N = (nx + 1)*(ny + 1);
    delete_mem();
    theAr = new theArgs ();
    //printf("INIT_MEM\n");
    timer = new QTimer(this);
    theAr->a = a;
    theAr->b = b;
    theAr->c = c;
    theAr->d = d;
    theAr->nx = nx;
    theAr->ny = ny;
    
    x = new double[N];
    for (int i = 0; i < N; i++) x[i] = 0;
    theAr->initialize();
    theAr->f = f;
    theAr->k = k;
    //printf("CD init_mem %e %e \n", c, d);
    
    //f = f0;
    //theAr->k = k;
    init_func();
}

//CHANGE_HERE





int Window::close_tablo() {
    pthread_mutex_lock(&theAr->m);
    int status = theAr->curr_calcul;
    pthread_mutex_unlock (&theAr->m);
    if (status == calculation::PROCESS) {
        QMessageBox::information (0, "Calculation", "Calculation in process, wait please!");
        return 1;
    }
    pthread_mutex_lock(&theAr->m);
    theAr->task = task_status::STOP;
    pthread_cond_broadcast(&theAr->wait_condition);
    pthread_mutex_unlock (&theAr->m);
    return 0;
}

int close_X(Window *pointer) {
    return pointer->close_tablo();
}

int nothing(Window *) {
    return 0;
}

void Window::closeEvent(QCloseEvent *) {
    action(this);
}




void Window::continuous() {
    int task;
    int curr_calcul, curr_plan;
    pthread_mutex_lock (&theAr->m);
    task = theAr->task;
    curr_plan = theAr->curr_plan;
    curr_calcul = theAr->curr_calcul;
    pthread_mutex_unlock (&theAr->m);
    
    switch(task) {
        case WAIT:
            if (curr_calcul==calculation::END) {
                pthread_mutex_lock(&theAr->m);
                if (curr_plan==0) {
                    nx = theAr->nx;
                    ny = theAr->ny;
                    N = theAr->N;
                    a = theAr->a;
                    b = theAr->b;
                    c = theAr->c;
                    d = theAr->d;
                    f = theAr->f;
                    if (x) delete [] x;
                    x = new double[N];
                    for (int i = 0; i < N; i++) x[i] = theAr->x[i];
                    theAr->curr_calcul = calculation::REST;
                    //printf("CD continious %e %e \n", c, d);
                } else {
                    theAr->task = task_status::STOP;
                    printf ("ERROR!\n");
                    window()->close ();
                }
                action = close_X;
                pthread_mutex_unlock(&theAr->m);
                
                update ();
            }
            break;
        case STOP:
            window()->close();
            break;
        default:
            return;
    }
}

void Window::connection_continuous() {
    connect (timer, SIGNAL(timeout()), this, SLOT(continuous()));
    timer->start(50);
}

void Window::restart() {
    //printf("restart\n");
    pthread_mutex_lock (&theAr->m);
    theAr->task = task_status::DO;
    if (theAr->initialize()) {
        printf ("ERROR!\n");
        window()->close ();
    }
    //printf("BRUH\n");


    pthread_cond_broadcast (&theAr->wait_condition);
    pthread_mutex_unlock (&theAr->m);
}



/// change current function for drawing


void Window::change_func () {
    int curr_calcul;
    pthread_mutex_lock (&theAr->m);
    curr_calcul = theAr->curr_calcul;
    pthread_mutex_unlock (&theAr->m);
    if (curr_calcul == calculation::PROCESS) {
        QMessageBox::information (0, "Calculation", "Calculation in process, wait please!");
        return;
    }
    
    pthread_mutex_lock (&theAr->m);
    action = nothing;
    k = (k + 1) % 8;
    update();
    switch (k) {
        case 0:
            f_name = "f(x) = 1";
            f = f0;
            break;
        case 1:
            f_name = "f(x) = x";
            f = f1;
            break;
        case 2:
            f_name = "f(x) = y";
            f = f2;
            break;
        case 3:
            f_name = "f(x) = x + y";
            f = f3;
            break;
        case 4:
            f_name = "f(x) = sqrt(x*x + y*y)";
            f = f4;
            break;
        case 5:
            f_name = "f(x) = (x*x + y*y)";
            f = f5;
            break;
        case 6:
            f_name = "f(x) = exp(x*x - y*y)";
            f = f6;
            break;
        case 7:
            f_name = "f(x) = 1/(25*(x*x + y*y) + 1)";
            f = f7;
            break;
    }
    

    //rebuild x
    //Building_constructions();
    theAr->f = f;
    theAr->k = k;
    
    pthread_mutex_unlock (&theAr->m);
    restart();
    
}

QPointF Window::l2g (double x_loc, double y_loc)
{
  double x_gl = (x_loc - a) / (b - a) * width();
  double y_gl = (d - y_loc) / (d - c) * height();
  return QPointF (x_gl, y_gl);
}

void Window::change_state() {
    state+=1;
    state%=3;
    //Building_constructions();
    update();
}

void Window::up_scale() {
    zoom++;
    double a_new = a, b_new = b, c_new = c, d_new = d;
    a = a_new + (b_new - a_new)/4;
    b = b_new - (b_new - a_new)/4;
    c = c_new + (d_new - c_new)/4;
    d = d_new - (d_new - c_new)/4;
    //Building_constructions();
    update();
}

void Window::down_scale() {
    zoom--;
    double a_new = a, b_new = b, c_new = c, d_new = d;
    a = a_new - (b_new - a_new)/2;
    b = b_new + (b_new - a_new)/2;
    c = c_new - (d_new - c_new)/2;
    d = d_new + (d_new - c_new)/2;
    //Building_constructions();
    update();
    
}

void Window::up_n() {
    pthread_mutex_lock (&theAr->m);
    int curr_calcul = theAr->curr_calcul;
    pthread_mutex_unlock (&theAr->m);
    if (curr_calcul == calculation::PROCESS) {
        QMessageBox::information(0, "Calculation", "Calculation in process, wait please!");
        return;
    }
    pthread_mutex_lock (&theAr->m);
    action = nothing;
    theAr->nx*=2;
    theAr->ny*=2;
    theAr->N = (theAr->nx + 1)*(theAr->ny + 1);
    pthread_mutex_unlock (&theAr->m);
    
    //Building_constructions();

    restart();
    //update();
    
    
}


void Window::down_n() {
    int r=0;
    pthread_mutex_lock (&theAr->m);
    int curr_calcul = theAr->curr_calcul;
    pthread_mutex_unlock (&theAr->m);
    if (curr_calcul == calculation::PROCESS) {
        QMessageBox::information (0, "Calculation", "Calculation in process, wait please!");
        return;
    }

    pthread_mutex_lock (&theAr->m);
    action = nothing;
    theAr->nx/=2;
    theAr->ny/=2;
    if (theAr->nx < 5) {
        r = 1;
        theAr->nx = 5;
        fadingMessage* message = new fadingMessage("MINIMUM 5 POINTS!", 5000);
        message->show();
        
        //action = close_X;
        
    }
    if (theAr->ny < 5) {
        //r = 1;
        
        theAr->ny = 5;
        if (r==0) {
            fadingMessage* message = new fadingMessage("MINIMUM 5 POINTS!", 5000);
            message->show();
        }
        //action = close_X;
    }
    theAr->N = (theAr->nx + 1)*(theAr->ny + 1);
    pthread_mutex_unlock (&theAr->m);

    restart();
    //update ();
    //if (r==1) action = close_X;
    
    
}
void Window::up_perturbation() {
    perturbation++;
    //Building_constructions();
    update();
}
void Window::down_perturbation() {
    perturbation--;
    //Building_constructions();
    update();
}
void Window::up_visualization() {
    mx*=2;
    my*=2;
    update();
}

void Window::down_visualization() {
    mx = (mx+1)/2;
    my = (my+1)/2;
    int r = 0;
    if (mx < 5) {
        mx = 5;
        fadingMessage* message = new fadingMessage("MINIMUM 5 POINTS!", 5000);
        message->show();
        
        r = 1;
    }
    if (my < 5) {
        my = 5;
        if (r==0) {
            fadingMessage* message = new fadingMessage("MINIMUM 5 POINTS!", 5000);
            message->show();
        }
        
    }
    update();
}

void Window::COLORS (double value, int &color1, int &color2, int &color3) { 
    int ar1 = 199, ar2 = 180, ar3 = 70;
    int br1 = 168, br2 = 216, br3 = 125;
    int cr1 = 193, cr2 = 231, cr3 = 176;
    int dr1 = 231, dr2 = 183, dr3 = 176;
    double xr1=0, xr2=0, xr3=0;
    double yr1=0, yr2=0, yr3=0, X=0, Y=0;
    if (value <= l1) {
        X = limit1;
        Y = l1;
        //printf("limit1 %e l1 %e \n", X, Y);
        xr1 = ar1, xr2 = ar2, xr3 = ar3;
        yr1 = br1, yr2 = br2, yr3 = br3;
    } else if (value <= l2) {
        X = l1;
        Y = l2;
        //printf("l1 %e l2 %e \n", X, Y);
        xr1 = br1, xr2 = br2, xr3 = br3;
        yr1 = cr1, yr2 = cr2, yr3 = cr3;
    } else { // if (value <= limit2) 
        X = l2;
        Y = limit2;
        //printf("l2 %e limit2 %e \n", X, Y);
        xr1 = cr1, xr2 = cr2, xr3 = cr3;
        yr1 = dr1, yr2 = dr2, yr3 = dr3;
    }
    if ( fabs(Y - X) < 1e-15 ) {
        color1 = ar1;
        color2 = ar2;
        color3 = ar3;
    } else {
        color1 = (int) ((double) ((xr1*(Y - value) + yr1*(value - X)) / (Y - X) ))%256;
        color2 = (int) ((double) ((xr2*(Y - value) + yr2*(value - X)) / (Y - X) ))%256;
        color3 = (int) ((double) ((xr3*(Y - value) + yr3*(value - X)) / (Y - X) ))%256;
    }
    
      
    
}

void Window::LimitationsF() {
    double hx = (b - a)/mx, hy = (d - c)/my, value;
    double max_y = f(a, c), mx_y, mn_y, min_y = max_y, res;
    for (int i = 0; i < mx; i++) {
        for (int j = 0; j < my; j++) {
            value = f(a + i*hx + hx/3, c + j*hy + 2*hy/3);
            if (value > max_y) max_y = value;
            if (value < min_y) min_y = value;
            
            value = f(a + i*hx + 2*hx/3, c + j*hy + hy/3);
            if (value > max_y) max_y = value;
            if (value < min_y) min_y = value;
        }

      //printf("VALUE %e %e\n", x1, y1);
    }
    mn_y = fabs(min_y), mx_y = fabs(max_y);
    res = (mn_y < mx_y ? mx_y : mn_y);
    
    res+= 0.1*perturbation*res;
    limit = res;
    limit2 = max_y;
    limit1 = min_y;
    l1 = limit1 + (limit2- limit1)/3;
    l2 = limit1 + 2*(limit2 - limit1)/3;
    //printf("limit1 limit2 l1 l2");
}

void Window::LimitationsPf() {
    double h1 = (b - a)/nx, h2 = (d - c)/ny;
    double hx = (b - a)/mx, hy = (d - c)/my, value;
    double max_y = Pf(x, a, c, nx, ny, h1, h2, a, c), mx_y, mn_y, min_y = max_y, res;
    for (int i = 0; i < mx; i++) {
        for (int j = 0; j < my; j++) {
            
            value = Pf(x, a + i*hx + hx/3, c + j*hy + 2*hy/3, nx, ny, h1, h2, a, c);
            if (value > max_y) max_y = value;
            if (value < min_y) min_y = value;
            
            value = Pf(x, a + i*hx + 2*hx/3, c + j*hy + hy/3, nx, ny, h1, h2, a, c);
            if (value > max_y) max_y = value;
            if (value < min_y) min_y = value;
        }

      //printf("VALUE %e %e\n", x1, y1);
    }
    mn_y = fabs(min_y), mx_y = fabs(max_y);
    res = (mn_y < mx_y ? mx_y : mn_y);
    
    if (perturbation) res+=0.1*perturbation*res;
    limit = res;
    limit2 = max_y;
    limit1 = min_y;
    l1 = limit1 + (limit2- limit1)/3;
    l2 = limit1 + 2*(limit2 - limit1)/3;
}



void Window::LimitationsError() {
    double h1 = (b - a)/nx, h2 = (d - c)/ny, mx_y, mn_y;
    double hx = (b - a)/mx, hy = (d - c)/my, value;
    double max_y = f(a, c) - Pf(x, a, c, nx, ny, h1, h2, a, c), min_y = max_y, res;
    for (int i = 0; i < mx; i++) {
        for (int j = 0; j < my; j++) {   
            value = f(a + i*hx + hx/3, c + j*hy + 2*hy/3) - Pf(x, a + i*hx + hx/3, c + j*hy + 2*hy/3, nx, ny, h1, h2, a, c);
            if (value < min_y) min_y = value;
            if (value > max_y) max_y = value;
            
            
            value = f(a + i*hx + 2*hx/3, c + j*hy + hy/3) - Pf(x, a + i*hx + 2*hx/3, c + j*hy + hy/3, nx, ny, h1, h2, a, c);
            if (value < min_y) min_y = value;
            if (value > max_y) max_y = value;
        }

      //printf("VALUE %e %e\n", x1, y1);
    }
    mn_y = fabs(min_y), mx_y = fabs(max_y);
    res = (mn_y < mx_y ? mx_y : mn_y);
    
    res+= 0.1*perturbation*res;
    limit = res;
    limit2 = max_y;
    limit1 = min_y;
    l1 = limit1 + (limit2- limit1)/3;
    l2 = limit1 + 2*(limit2 - limit1)/3;
}



void Window::minicolor(QPainter *painter, QPen pen, QBrush brush, double value, QPointF* vert) {
    int color1, color2, color3;
    COLORS(value, color1, color2, color3);
    pen.setColor(QColor (color1, color2, color3));
    brush.setColor(QColor (color1, color2, color3));
    painter->setPen (pen);
    painter->setBrush (brush);
    painter->drawPolygon(vert, 3);

}


void Window::illustrateF (QPainter *painter) {
  QPen pen(Qt::black, 0, Qt::SolidLine);
  QBrush brush(Qt::black, Qt::SolidPattern);
  int color1, color2, color3;
  double hx = (b - a)/mx, hy = (d - c)/my, value;
  for (int i = 0; i < mx; i++) {
      for (int j = 0; j < my; j++) {
          QPointF vert[3] = {l2g(a + i*hx, c + j*hy), l2g(a + i*hx, c + (j + 1)*hy), l2g (a + (i + 1)*hx, c + (j + 1)*hy)};
          value = f(a + i*hx + hx/3, c + j*hy + 2*hy/3);
          if (i==(mx + 1)/2 && j==(my + 1)/2) {
            value += 0.1*perturbation*limit;          
          }
          COLORS(value, color1, color2, color3);
            pen.setColor(QColor (color1, color2, color3));
            brush.setColor(QColor (color1, color2, color3));
            painter->setPen (pen);
            painter->setBrush (brush);
            painter->drawPolygon(vert, 3);
          
          vert[1] = l2g(a + (i + 1)*hx, c + j*hy);
          value = f(a + i*hx + 2*hx/3, c + j*hy + hy/3);
          if (i==(mx + 1)/2 && j==(my + 1)/2) {
            value += 0.1*perturbation*limit;          
          }
          COLORS(value, color1, color2, color3);
            pen.setColor(QColor (color1, color2, color3));
            brush.setColor(QColor (color1, color2, color3));
            painter->setPen (pen);
            painter->setBrush (brush);
            painter->drawPolygon(vert, 3);
        }
    }
}

void Window::illustratePf(QPainter *painter) {
  double hx = (b - a)/mx, hy = (d - c)/my, value;
  double h1 = (b - a)/nx, h2 = (d - c)/ny;
  QPen pen (Qt::black, 0, Qt::SolidLine);
  QBrush brush (Qt::black, Qt::SolidPattern);
  for (int i = 0; i < mx; i++) {
      for (int j = 0; j < my; j++) {
          QPointF vert[3] = {l2g(a + i*hx, c + j*hy), l2g(a + (i + 1)*hx, c + j*hy), l2g (a + (i + 1)*hx, c + (j + 1)*hy)};
          value = Pf(x, a + i*hx + 2*hx/3, c + j*hy + hy/3, nx, ny, h1, h2, a, c);
          if (i==(mx + 1)/2 && j==(my + 1)/2) {
            value += 0.1*perturbation*limit;          
          }
          minicolor(painter, pen, brush, value, vert);
          
          vert[1] = l2g(a + i*hx, c + (j + 1)*hy);
          value = Pf(x, a + i*hx + hx/3, c + j*hy + 2*hy/3, nx, ny, h1, h2, a, c);
          //printf("PF2 %e \n", value);
          if (i==(mx + 1)/2 && j==(my + 1)/2) {
            value += 0.1*perturbation*limit;          
          }
          minicolor(painter, pen, brush, value, vert);
        }
    }
}

//STOPPED HERE

void Window::illustrateError(QPainter *painter) {
  double hx = (b - a)/mx, hy = (d - c)/my, value;
  double h1 = (b - a)/nx, h2 = (d - c)/ny;
  QPen pen(Qt::black, 0, Qt::SolidLine);
  QBrush brush(Qt::black, Qt::SolidPattern);
  for (int i = 0; i < mx; i++) {
      for (int j = 0; j < my; j++) {
          QPointF vert[3] = {l2g(a + i*hx, c + j*hy), l2g(a + (i + 1)*hx, c + j*hy), l2g (a + (i + 1)*hx, c + (j + 1)*hy)};
          value = f(a + i*hx + 2*hx/3, c + j*hy + hy/3) - Pf(x, a + i*hx + 2*hx/3, c + j*hy + hy/3, nx, ny, h1, h2, a, c);
          if (i==(mx + 1)/2 && j==(my + 1)/2) {
            value += 0.1*perturbation*limit;          
          }
          minicolor(painter, pen, brush, value, vert);
          
          vert[1] = l2g(a + i*hx, c + (j + 1)*hy);
          value = f(a + i*hx + hx/3, c + j*hy + 2*hy/3) - Pf(x, a + i*hx + hx/3, c + j*hy + 2*hy/3, nx, ny, h1, h2, a, c);
          if (i == (mx + 1)/2 && j == (my + 1)/2) {
            value += 0.1*perturbation*limit;          
          }
          minicolor(painter, pen, brush, value, vert);
        }
    }
}

void Window::paintEvent(QPaintEvent * /* event */) {  
  QPainter painter (this);
  char print_state[LEN];
  char print_F_max[LEN];
  char print_scale[LEN];
  char print_n_p[LEN];
  char print_scale_2[LEN];
  char print_n_p_2[LEN];
  
  printf("%s\n", f_name);
  sprintf(print_state, "State = %d", state+1);
  QPen pen_black(Qt::black, 1, Qt::SolidLine); 
  QPen pen_red(Qt::red, 0.5, Qt::SolidLine); 
  QPen pen_blue(Qt::blue, 2, Qt::SolidLine); 
  QPen pen_green(Qt::green, 2, Qt::SolidLine); 
  QPen pen_magenta(Qt::magenta, 0.5, Qt::SolidLine); 

  painter.setPen (pen_black);
  
  if (state == 0) {
    LimitationsF();
    illustrateF(&painter);
  } else if (state == 1) {
    LimitationsPf();
    illustratePf(&painter);
  } else {
    LimitationsError();
    illustrateError(&painter);
  }
  
  //printf("VALUE 1 %e %e\n", min_y, max_y);
  
  sprintf(print_scale, "Current scale: %d. a = %e; b = %e;\n", zoom, a, b);
  sprintf(print_scale_2, "                           c = %e; d = %e\n", c, d);
  sprintf(print_n_p, "nx = %d ny = %d mx = %d my = %d \n", nx, ny, mx, my);
  sprintf(print_n_p_2, "perturbation = %d\n", perturbation);
  
  printf("Current scale: %d. a = %e; b = %e;\n", zoom, a, b);
  printf("                   c = %e; d = %e\n", c, d);
  printf("nx = %d ny = %d mx = %d my = %d \n", nx, ny, mx, my);
  printf("perturbation = %d\n", perturbation);
  

  
  sprintf(print_F_max, "|F| = %e\n", limit);
  printf("|F| = %e\n \n", limit);
  //delta_y = 0.01 * (mxx - mnn + 1e-17);
  //mnn -= delta_y;
  //mxx += delta_y;
  //printf("MXX %e\n", mxx);

  painter.setFont(QFont("Verdana", 14));
  // draw axis
  painter.setPen (pen_blue);
  painter.drawLine (l2g(0, c), l2g(0, d));
  painter.drawLine (l2g(a, 0), l2g(b, 0));

  // render function name
  painter.setPen ("black");
  painter.drawText (0, 20, f_name);
  //painter.setPen("black");
  painter.drawText (0, 40, print_state);
  //painter.setPen("blue");
  painter.drawText (0, 60, print_F_max);
  //painter.setPen("blue");
  painter.drawText (0, 80, print_scale);
  painter.drawText (0, 100, print_scale_2);
  //painter.setPen("blue");
  painter.drawText (0, 120, print_n_p);
  painter.drawText (0, 140, print_n_p_2);
  /*
  painter.setFont(QFont("Verdana", 10));
  painter.setPen ("black");
  painter.drawText (0, 120, color1);
  painter.setPen ("blue");
  painter.drawText (0, 140, color2);
  painter.setPen ("green");
  painter.drawText (0, 160, color3);
  painter.setPen ("red");
  painter.drawText (0, 180, color4);
  painter.setPen ("magenta");
  painter.drawText (0, 200, color5);
  */
  
  

}


void* thread_func(void* arg) { //write here counterRes
    Args* ar=(Args *)arg;
    theArgs* theAr = ar->theAr;
    pthread_mutex_lock (&theAr->m);
    theAr->curr_calcul = calculation::PROCESS;
    pthread_mutex_unlock (&theAr->m);
    double a = theAr->a;
    double b = theAr->b;
    double c = theAr->c;
    double d = theAr->d;
    //printf("CD thread_func %e %e", c, d);
    int nx = theAr->nx;
    int ny = theAr->ny;
    //int mx = ar->mx;
    //int my = ar->my;
    double eps = ar->eps;
    int maxit = ar->maxit;
    int maxstep = ar->maxstep;
    int p = ar->p;
    int N = theAr->N;
    double* A  = theAr->A;
    int* I  = theAr->I;
    double* bb = theAr->bb;
    double* u = theAr->u;
    double* v = theAr->v;
    double* x = theAr->x;
    double* r = theAr->r;
    double* r1 = ar->r1;
    double* r2 = ar->r2;
    double* r3 = ar->r3;
    double* r4 = ar->r4;
    double (*f) (double,double);
    f = theAr->f;
    int n_thread = ar->n_thread;
    //printf("n_thread %d\n", n_thread);
    //printf("WORKS! %d\n", n_thread);

    cpu_set_t cpu;
    CPU_ZERO(&cpu);
    int n_cpus=get_nprocs();
    int cpu_id=n_cpus-1-(n_thread%n_cpus);
    CPU_SET(cpu_id, &cpu);
    pthread_t tid=pthread_self();
    pthread_setaffinity_np(tid, sizeof(cpu), &cpu);
    int t1 = N*n_thread/p;
    int t2 = N*(n_thread + 1)/p;
    memset (bb + t1, 0, (t2 - t1)*sizeof(double));
    memset (u + t1, 0, (t2 - t1)*sizeof(double));
    memset (v + t1, 0, (t2 - t1)*sizeof(double));
    memset (x + t1, 0, (t2 - t1)*sizeof(double));
    memset (r + t1, 0, (t2 - t1)*sizeof(double));

    reduce_sum(p, 0, 0);
    //printf("reduce_sum %d\n", p);
    //printf("WORKS! 1 %d\n", n_thread);
    double hx = (b - a)/nx;
    double hy = (d - c)/ny;
    //printf("fill_I \n");
    double T1 = get_full_time();
    if (n_thread==0) fill_I(nx, ny, I);
    //printf("fill_I 2\n");
    reduce_sum(p, 0, 0);
    //printf("WORKS! 2\n");
    fill_IA (N, nx, ny, hx, hy, I, A, p, n_thread);
    reduce_sum(p, 0, 0);
    //printf("fill_IA \n");
    fill_b (N, nx, ny, hx, hy, bb, a, c, p, n_thread, f);
    reduce_sum(p, 0, 0);
    //printf("fill_b \n");
    
    //printf("WORKS! 3\n");
    ar->it = minimal_residual_msr_matrix_full(N, A, I, bb, x, r, u, v, eps, maxit, maxstep, p, n_thread); 
    T1 = get_full_time() - T1;
    reduce_sum_int(ar->p, &ar->it, 1);
    
    //printf("WORKS! 4\n");
    double T2 = get_full_time();
    res_1(N, nx, ny, hx, hy, a, c, p, n_thread, x, r1, f);
    res_2(N, nx, ny, hx, hy, a, c, p, n_thread, x, r2, f);
    res_3(N, nx, ny, hx, hy, a, c, p, n_thread, x, r3, f);
    res_4(N, nx, ny, hx, hy, a, c, p, n_thread, x, r4, f);
    T2 = get_full_time() - T2;
    
    

    if (n_thread == 0) {
        printf ( "%s : Task = %d R1 = %e R2 = %e R3 = %e R4 = %e T1 = %.2f T2 = %.2f\\It = %d E = %e K = %d Nx = %d Ny = %d P = %d\n",
        ar->arguments, 1, r1[0], r2[0], r3[0], r4[0], T1, T2, ar->it, eps, theAr->k, nx, ny, p);
    }
    pthread_mutex_lock (&theAr->m);
    if (theAr->deleting==1) theAr->task = task_status::STOP;
    else theAr->task = task_status::WAIT;
    theAr->curr_plan = 0;
    theAr->curr_calcul = calculation::END;
    pthread_mutex_unlock (&theAr->m);
    return 0;
}

void* continue_thread (void* arg) {
    Args* ar = (Args *) arg;
    theArgs* theAr = ar->theAr;
    reduce_sum(ar->p, 0, 0);
    while (true) {
        pthread_mutex_lock (&theAr->m);
        while (theAr->task==task_status::WAIT) pthread_cond_wait(&theAr->wait_condition, &theAr->m);
        pthread_mutex_unlock (&theAr->m);
        if (theAr->task==task_status::DO) thread_func(ar);
        else {
            if (theAr->task==task_status::STOP) reduce_sum(ar->p, 0, 0);
            break;
        }
    }
    return 0;
}

