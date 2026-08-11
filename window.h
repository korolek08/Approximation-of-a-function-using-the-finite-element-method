#ifndef WINDOW_H
#define WINDOW_H

#include <QtWidgets/QtWidgets>
#define STEP 100

#include <stdio.h>
#include <sys/sysinfo.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>
#include "io_status.h"



class fadingMessage : public QWidget {
    Q_OBJECT

    public:
        fadingMessage(const QString& message, int duration = 5000, QWidget* parent = nullptr) : QWidget(parent) {
            auto screen = QGuiApplication::primaryScreen(); 
            QRect rect = screen->geometry(); 
            QPoint center = rect.center(); 
            
            
            label = new QLabel(message, this);
            label->setAlignment(Qt::AlignCenter);
            QVBoxLayout* layout = new QVBoxLayout(this);
            layout->addWidget(label);
            setLayout(layout);
            
            setGeometry(0, 0, 200, 50); 
            center.setX(center.x() - 100);
            center.setY(center.y() - 25);

            
            move(center); 
            timer = new QTimer(this);
            connect(timer, &QTimer::timeout, this, &fadingMessage::hideMessage);
            timer->start(duration);
        }

public slots:
    void hideMessage() {
        close();
    }

private:
    QLabel* label;
    QTimer* timer;
};

class Window;
int nothing(Window *);
int close_X(Window *);
int get_len_msr(int nx, int ny);


class theArgs {
    public:
        int k;
        int deleting = 0;
        int curr_plan = 0;
        double a;
        double b;
        double c;
        double d;
        int nx;
        int ny;
        int N;
        int len;
        double (*f) (double, double);
        double* A = nullptr;
        int* I = nullptr;
        double* bb = nullptr;
        double* u = nullptr;
        double* v = nullptr;
        double* x = nullptr;
        double* r = nullptr;
        int task = task_status::WAIT;
        int curr_calcul = calculation::PROCESS;
        
        pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
        pthread_cond_t wait_condition = PTHREAD_COND_INITIALIZER;
        void deleting_mem();
        int initialize();
        theArgs() = default;

        
        
        ~theArgs () {
            deleting_mem();
        }

};


class Args {
	public:
	    theArgs * theAr = nullptr;
	    char* arguments;
		double eps = 1e-15;
		int maxit;
		int maxstep = 10;
		int p;
		int len;
		int N;
		int n_thread;
		double t1 = 0;
		double t2 = 0;
		double (*f) (double,double);
		int it = 0;
		io_status error_type=io_status::undef;
		double error_flag=0;
        io_status task;
        pthread_t tid;
        double* r1 = nullptr;
  	    double* r2 = nullptr;
	    double* r3 = nullptr;
	    double* r4 = nullptr;
		
		Args() = default;
		Args (const Args& x)=delete;
		Args (Args && x) = default;

};

void* thread_func(void* arg);
void* continue_thread (void* arg);



class Window : public QWidget { 
Q_OBJECT
    public:
        Args* ar;
        theArgs* theAr = nullptr;
        int func_id;
        const char* f_name;
        double a;
        double b;
        double c;
        double d;
        int nx;
        int ny;
        int mx = 10;
        int my = 10;
        int k = 0;
        double eps = 1e-15;
        int maxit;
        int maxstep = 10;
        int p = 2;
        int len;
        int N;
	    double* x = nullptr;
	    int n_thread;
	    double t1 = 0;
	    double t2 = 0;
	    double (*f) (double,double);
	    int it = 0;
        
        int state = 0;
        int zoom = 0;
        int perturbation = 0;
        double limit = 0;
        double limit1 = 0;
        double limit2 = 0;
        double l1 = 0;
        double l2 = 0;
        double* r1 = nullptr;
        double* r2 = nullptr;
        double* r3 = nullptr;
        double* r4 = nullptr;
        int (*action) (Window*) = nothing;
        QTimer* timer = nullptr;
        int return_p() {
            return p;
        }
        void thread_clean();
        void delete_mem();
        ~Window () {
            thread_clean();
            delete_mem();
            if (r1) delete [] r1;
            if (r2) delete [] r2;
            if (r3) delete [] r3;
            if (r4) delete [] r4;
        }
        void add_param(Args *ar);
        

    public:
      Window (QWidget *parent);

      QSize minimumSizeHint() const;
      QSize sizeHint() const;

      int commandLine(int argc, char **argv);
      QPointF l2g(double x_loc, double y_loc);
      
      void Polynomial();
      double count_value(double x);
      double count_spline_value(double x);
      double Lagrange_derivative(double x0, double x1, double x2, double x3,
                                          double y1, double y2, double y3);
      void Gauss_three_lines(int n, double* up, double* diag, double* down, double* d, double* b);
      void fill_diff(double* diff, double dist);
      void Spline(double* d, double* diff, double rev_dist, double rev_sq_dist);
      void init_r();
      void init_mem();
      void init_func ();
      void restart ();
      void LimitationsF();
      void LimitationsPf();
      void LimitationsError();
      void COLORS (double value, int &color1, int &color2, int &color3);
      void minicolor(QPainter* painter, QPen pen, QBrush brush, double value, QPointF* vert);
      void illustrateF (QPainter* painter);
      void illustratePf (QPainter* painter);
      void illustrateError (QPainter* painter);
      void connection_continuous();
      void closeEvent(QCloseEvent*);
      
    public slots:
        
        void change_func ();
        void change_state();
        void up_scale();
        void down_scale();
        void up_n();
        void down_n();
        void up_perturbation();
        void down_perturbation();
        void up_visualization();
        void down_visualization();
        void continuous();
        int close_tablo();
        
    protected:
      void paintEvent (QPaintEvent *event);
};


#endif
